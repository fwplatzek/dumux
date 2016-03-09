// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Boundary types gathered on an element
 */
#ifndef DUMUX_BOX_ELEMENT_BOUNDARY_TYPES_HH
#define DUMUX_BOX_ELEMENT_BOUNDARY_TYPES_HH

#include "properties.hh"

#include <dumux/common/valgrind.hh>

namespace Dumux
{

/*!
 * \ingroup BoxModel
 * \ingroup ImplicitBoundaryTypes
 * \brief This class stores an array of BoundaryTypes objects
 */
template<class TypeTag>
class BoxElementBoundaryTypes : public std::vector<typename GET_PROP_TYPE(TypeTag, BoundaryTypes) >
{
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef std::vector<BoundaryTypes> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    enum { dim = GridView::dimension };
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    /*!
     * \brief Copy constructor.
     *
     * Copying a the boundary types of an element should be explicitly
     * requested
     */
    explicit BoxElementBoundaryTypes(const BoxElementBoundaryTypes &v)
        : ParentType(v)
    {}

    /*!
     * \brief Default constructor.
     */
    BoxElementBoundaryTypes()
    {
        hasDirichlet_ = false;
        hasNeumann_ = false;
        hasOutflow_ = false;
    }

    /*!
     * \brief Update the boundary types for all vertices of an element.
     *
     * \param problem The problem object which needs to be simulated
     * \param element The DUNE Codim<0> entity for which the boundary
     *                types should be collected
     * \param fvGeometry The element's finite volume geometry
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry)
    {
        int numVertices = element.subEntities(dim);

        this->resize(numVertices);

        hasDirichlet_ = false;
        hasNeumann_ = false;
        hasOutflow_ = false;

        for (auto&& scv : fvGeometry.scvs())
        {
            int scvIdxLocal = scv.indexInElement();
            (*this)[scvIdxLocal].reset();

            if (problem.model().onBoundary(scv))
            {
                (*this)[scvIdxLocal] = problem.boundaryTypes(element, scv);

                hasDirichlet_ = hasDirichlet_ || (*this)[scvIdxLocal].hasDirichlet();
                hasNeumann_ = hasNeumann_ || (*this)[scvIdxLocal].hasNeumann();
                hasOutflow_ = hasOutflow_ || (*this)[scvIdxLocal].hasOutflow();
            }
        }
    }

    /*!
     * \brief Update the boundary types for all vertices of an element.
     *
     * \param problem The problem object which needs to be simulated
     * \param element The DUNE Codim<0> entity for which the boundary
     *                types should be collected
     */
    void update(const Problem &problem,
                const Element &element)
    {
        const auto& fvGeometry = problem.model().fvGeometries(element);
        update(problem, element, fvGeometry);
    }

    /*!
     * \brief Returns whether the element has a vertex which contains
     *        a Dirichlet value.
     */
    bool hasDirichlet() const
    { return hasDirichlet_; }

    /*!
     * \brief Returns whether the element potentially features a
     *        Neumann boundary segment.
     */
    bool hasNeumann() const
    { return hasNeumann_; }

    /*!
     * \brief Returns whether the element potentially features an
     *        outflow boundary segment.
     */
    bool hasOutflow() const
    { return hasOutflow_; }

protected:
    bool hasDirichlet_;
    bool hasNeumann_;
    bool hasOutflow_;
};

} // namespace Dumux

#endif
