// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Common
 * \brief Define helper functions for integration
 */
#ifndef DUMUX_COMMON_INTEGRATE_HH
#define DUMUX_COMMON_INTEGRATE_HH

#include <cmath>
#include <type_traits>

#include <dune/geometry/quadraturerules.hh>
#include <dune/common/concept.hh>

#if HAVE_DUNE_FUNCTIONS
#include <dune/functions/gridfunctions/gridfunction.hh>
#endif

#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {

// implementation details of the integrate functions
#ifndef DOXYGEN
namespace Impl {

struct HasLocalFunction
{
    template<class F>
    auto require(F&& f) -> decltype(
      localFunction(f),
      localFunction(f).unbind()
    );
};

template<class F>
static constexpr bool hasLocalFunction()
{ return Dune::models<HasLocalFunction, F>(); }

template<class Error,
         typename std::enable_if_t<IsIndexable<Error>::value, int> = 0>
Error sqrtNorm(const Error& error)
{
    using std::sqrt;
    auto e = error;
    for (int i = 0; i < error.size(); ++i)
        e[i] = sqrt(error[i]);
    return e;
}

template<class Error,
         typename std::enable_if_t<!IsIndexable<Error>::value, int> = 0>
Error sqrtNorm(const Error& error)
{
    using std::sqrt;
    return sqrt(error);
}

template<class T, typename = int>
struct FieldTypeImpl
{
    using type = T;
};

template<class T>
struct FieldTypeImpl<T, typename std::enable_if<(sizeof(std::declval<T>()[0]) > 0), int>::type>
{
    using type = typename FieldTypeImpl<std::decay_t<decltype(std::declval<T>()[0])>>::type;
};

template<class T>
using FieldType = typename FieldTypeImpl<T>::type;

} // end namespace Impl
#endif

/*!
 * \brief Integrate a grid function over a grid view
 * \param gg the grid geometry
 * \param f the grid function
 * \param order the order of the quadrature rule
 */
template<class GridGeometry, class SolutionVector,
         typename std::enable_if_t<!Impl::hasLocalFunction<SolutionVector>(), int> = 0>
auto integrateGridFunction(const GridGeometry& gg,
                           const SolutionVector& sol,
                           std::size_t order)
{
    using GridView = typename GridGeometry::GridView;
    using Scalar = typename Impl::FieldType< std::decay_t<decltype(sol[0])> >;

    Scalar integral(0.0);
    for (const auto& element : elements(gg.gridView()))
    {
        const auto elemSol = elementSolution(element, sol, gg);
        const auto geometry = element.geometry();
        const auto& quad = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(geometry.type(), order);
        for (auto&& qp : quad)
        {
            auto value = evalSolution(element, geometry, gg, elemSol, geometry.global(qp.position()));
            value *= qp.weight()*geometry.integrationElement(qp.position());
            integral += value;
        }
    }
    return integral;
}

/*!
 * \brief Integrate a function over a grid view
 * \param gv the grid view
 * \param f the first function
 * \param g the second function
 * \param order the order of the quadrature rule
 * \note dune functions currently doesn't support composing two functions
 */
template<class GridGeometry, class Sol1, class Sol2,
         typename std::enable_if_t<!Impl::hasLocalFunction<Sol1>(), int> = 0>
auto integrateL2Error(const GridGeometry& gg,
                      const Sol1& sol1,
                      const Sol2& sol2,
                      std::size_t order)
{
    using GridView = typename GridGeometry::GridView;
    using Scalar = typename Impl::FieldType< std::decay_t<decltype(sol1[0])> >;

    Scalar l2norm(0.0);
    for (const auto& element : elements(gg.gridView()))
    {
        const auto elemSol1 = elementSolution(element, sol1, gg);
        const auto elemSol2 = elementSolution(element, sol2, gg);

        const auto geometry = element.geometry();
        const auto& quad = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(geometry.type(), order);
        for (auto&& qp : quad)
        {
            const auto& globalPos = geometry.global(qp.position());
            const auto value1 = evalSolution(element, geometry, gg, elemSol1, globalPos);
            const auto value2 = evalSolution(element, geometry, gg, elemSol2, globalPos);
            const auto error = (value1 - value2);
            l2norm += (error*error)*qp.weight()*geometry.integrationElement(qp.position());
        }
    }

    using std::sqrt;
    return sqrt(l2norm);
}

#if HAVE_DUNE_FUNCTIONS

/*!
 * \brief Integrate a grid function over a grid view
 * \param gv the grid view
 * \param f the grid function
 * \param order the order of the quadrature rule
 * \note overload for a Dune::Funtions::GridFunction
 */
template<class GridView, class F,
         typename std::enable_if_t<Impl::hasLocalFunction<F>(), int> = 0>
auto integrateGridFunction(const GridView& gv,
                           const F& f,
                           std::size_t order)
{
    auto fLocal = localFunction(f);

    using Element = typename GridView::template Codim<0>::Entity;
    using LocalPosition = typename Element::Geometry::LocalCoordinate;
    using Scalar = typename Impl::FieldType< std::decay_t<decltype(fLocal(std::declval<LocalPosition>()))> >;

    Scalar integral(0.0);
    for (const auto& element : elements(gv))
    {
        fLocal.bind(element);

        const auto geometry = element.geometry();
        const auto& quad = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(geometry.type(), order);
        for (auto&& qp : quad)
        {
            auto value = fLocal(qp.position());
            value *= qp.weight()*geometry.integrationElement(qp.position());
            integral += value;
        }

        fLocal.unbind();
    }
    return integral;
}

/*!
 * \brief Integrate a function over a grid view
 * \param gv the grid view
 * \param f the first function
 * \param g the second function
 * \param order the order of the quadrature rule
 * \note overload for a Dune::Funtions::GridFunction
 * \note dune functions currently doesn't support composing two functions
 */
template<class GridView, class F, class G,
         typename std::enable_if_t<Impl::hasLocalFunction<F>(), int> = 0>
auto integrateL2Error(const GridView& gv,
                      const F& f,
                      const G& g,
                      std::size_t order)
{
    auto fLocal = localFunction(f);
    auto gLocal = localFunction(g);

    using Element = typename GridView::template Codim<0>::Entity;
    using LocalPosition = typename Element::Geometry::LocalCoordinate;
    using Scalar = typename Impl::FieldType< std::decay_t<decltype(fLocal(std::declval<LocalPosition>()))> >;

    Scalar l2norm(0.0);
    for (const auto& element : elements(gv))
    {
        fLocal.bind(element);
        gLocal.bind(element);

        const auto geometry = element.geometry();
        const auto& quad = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(geometry.type(), order);
        for (auto&& qp : quad)
        {
            const auto error = fLocal(qp.position()) - gLocal(qp.position());
            l2norm += (error*error)*qp.weight()*geometry.integrationElement(qp.position());
        }

        gLocal.unbind();
        fLocal.unbind();
    }

    using std::sqrt;
    return sqrt(l2norm);
}

/*!
 * \brief Integrate a function over a grid view
 * \param gv the grid view
 * \param f the first function
 * \param g the second function
 * \param order the order of the quadrature rule
 * \note overload for a Dune::Funtions::GridFunction
 * \note dune functions currently doesn't support composing two functions
 */
template<class GridView, class F, class G,
         typename std::enable_if_t<Impl::hasLocalFunction<F>(), int> = 0>
auto integrateFaceFluxError(const GridView& gv,
                            const F& f,
                            const G& g,
                            std::size_t order)
{
    auto fLocal = localFunction(f);
    auto gLocal = localFunction(g);

    using Element = typename GridView::template Codim<0>::Entity;
    using LocalPosition = typename Element::Geometry::LocalCoordinate;
    using Scalar = typename Impl::FieldType< std::decay_t<decltype(fLocal(std::declval<LocalPosition>()))> >;

    Scalar norm(0.0);
    Scalar ifNorm(0.0);
    Scalar stripNorm(0.0);
    Scalar stripIFNorm(0.0);
    for (const auto& element : elements(gv))
    {
        fLocal.bind(element);
        gLocal.bind(element);

        bool isInterfaceElement = false;
        if (element.hasBoundaryIntersections())
        {
            const auto y = element.geometry().center()[1];
            for (const auto& is : intersections(gv, element))
                if (is.boundary())
                {
                    if (y < 1.0 && is.centerUnitOuterNormal()[1] > 0.9)
                        isInterfaceElement = true;
                    else if (y > 1.0 && is.centerUnitOuterNormal()[1] < -0.9)
                        isInterfaceElement = true;
                }
        }

        Scalar localNorm = 0.0;
        Scalar localIFNorm = 0.0;
        for (const auto& is : intersections(gv, element))
        {
            const auto isGeometry = is.geometry();
            const auto normal = is.centerUnitOuterNormal();
            const auto gCenter = gLocal(isGeometry.center())*normal;

            Scalar errorSq = 0.0;
            const auto& quad = Dune::QuadratureRules<Scalar, GridView::dimension-1>::rule(isGeometry.type(), order);
            for (auto&& qp : quad)
            {
                const auto localPos = is.geometryInInside().global(qp.position());
                const auto error = fLocal(localPos)*normal - gCenter;
                errorSq += (error*error)*qp.weight()*isGeometry.integrationElement(qp.position());
            }

            localNorm += errorSq/isGeometry.volume();

            if (is.boundary())
            {
                const auto y = element.geometry().center()[1];
                if (y < 1.0 && normal[1] > 0.9)
                {
                    ifNorm += errorSq;
                    localIFNorm += errorSq/isGeometry.volume();
                }
                else if (y > 1.0 && normal[1] < -0.9)
                {
                    ifNorm += errorSq;
                    localIFNorm += errorSq/isGeometry.volume();
                }
            }
        }

        if (isInterfaceElement)
        {
            stripNorm += localNorm*element.geometry().volume();
            stripIFNorm += localIFNorm*element.geometry().volume();
        }
        norm += localNorm*element.geometry().volume();

        gLocal.unbind();
        fLocal.unbind();
    }

    std::cout << "Norm - strip norm - strip IF norm - if norm: " << std::sqrt(norm) << " - " << std::sqrt(stripNorm) << " - " << std::sqrt(stripIFNorm) << " - "<< std::sqrt(ifNorm) << std::endl;

    using std::sqrt;
    return sqrt(norm);
}

#endif

} // end namespace Dumux

#endif
