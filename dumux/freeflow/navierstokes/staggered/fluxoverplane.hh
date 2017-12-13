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
 *
 * \brief Class used calculate fluxes over planes. This only works for the staggered grid discretization.
 */
#ifndef DUMUX_FLUX_OVER_PLANE_STAGGERED_HH
#define DUMUX_FLUX_OVER_PLANE_STAGGERED_HH

#include <dune/common/fvector.hh>
#include <dumux/common/properties.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dumux/common/boundingboxtree.hh>



namespace Dumux
{


/*!
 * \brief  Class used to calculate fluxes over planes. This only works for the staggered grid discretization.
 */
template <class TypeTag>
class FluxOverPlane
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using Element = typename GridView::template Codim<0>::Entity;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dim>;

    static constexpr auto planeDim = dim - 1;
    using PlaneGeometryType = Dune::AffineGeometry< Scalar, planeDim, dim >;


    /*!
     * \brief Auxiliary class that holds the plane-specific data.
     */
    template<int mydim, int coordDim, class Scalar=double>
    class PlaneData
    {
        using Geo = Dune::AffineGeometry< Scalar, mydim, coordDim >;
    public:

        PlaneData() {}

        PlaneData(Geo&& geo)
        {
            values_.resize(1);
            geometries_.push_back(std::move(geo));
        }

        PlaneData(std::vector<Geo>&& geo)
        {
            values_.resize(geo.size());
            geometries_ = std::forward<decltype(geo)>(geo);
        }

        void addValue(int planeIdx, Scalar value)
        {
            values_[planeIdx] += value;
        }

        void addSubPlane(Geo&& geo)
        {
            values_.push_back(0.0);
            geometries_.push_back(std::move(geo));
        }

        auto& subPlanes() const
        {
            return geometries_;
        }

        Scalar value(int planeIdx) const
        {
            return values_[planeIdx];
        }

        auto& values() const
        {
            return values_;
        }

        void printPlaneBoundaries(int planeIdx) const
        {
            const auto& geometry = geometries_[planeIdx];
            for(int i = 0; i < geometry.corners(); ++i)
                std::cout << geometry.corner(i) << "  ";
        }

        void resetValues()
        {
            std::fill(values_.begin(), values_.end(), 0.0);
        }

    private:

        std::vector<Geo> geometries_;
        std::vector<Scalar> values_;
    };


public:

    using PlaneList = std::vector<PlaneGeometryType>;

    /*!
     * \brief The constructor
     *
     * \param assembler The assembler
     * \param sol The solution vector
     */
    template<class Assembler>
    FluxOverPlane(const Assembler& assembler,
                  const SolutionVector& sol)
    : problem_(assembler.problem())
    , gridVariables_(assembler.gridVariables())
    , sol_(sol)
    {
        verbose_  = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "FluxOverPlane.Verbose", false);
    }

    /*!
     * \brief Add a collection of sub planes under a given name
     *
     * \param name The name of the plane
     * \param sol The list of sub planes
     */
    void addPlane(const std::string& name, PlaneList&& planes )
    {
        planes_[name] = PlaneData<planeDim,dim>(std::forward<decltype(planes)>(planes));
    }

    /*!
     * \brief Add a plane under a given name, specifing the plane's corner points.
     *        This is a specialization for 2D, therefore the plane is actually a line.
     *
     * \param name The name of the plane
     * \param p0 The first corner
     * \param p1 The second corner
     */
    void addPlane(const std::string& name, const GlobalPosition& p0, const GlobalPosition& p1)
    {
        planes_[name].addSubPlane(makePlane(p0, p1));
    }

    /*!
     * \brief Add a plane under a given name, specifing the plane's corner points.
     *        This is a specialization for 3D.
     *
     * \param name The name of the plane
     * \param p0 The first corner
     * \param p1 The second corner
     * \param p2 The third corner
     * \param p3 The fourth corner
     */
    void addPlane(const std::string& name,
                  const GlobalPosition& p0,
                  const GlobalPosition& p1,
                  const GlobalPosition& p2,
                  const GlobalPosition& p3)
    {
        planes_[name].addSubPlane(makePlane(p0, p1, p2, p3));
    }

    /*!
     * \brief Creates a geometrical plane object.
     *        This is a specialization for 2D, therefore the plane is actually a line.
     *
     * \param p0 The first corner
     * \param p1 The second corner
     */
    static PlaneGeometryType makePlane(const GlobalPosition& p0, const GlobalPosition& p1)
    {
        const std::vector< Dune::FieldVector< Scalar, dim > > corners = {p0, p1};
        return PlaneGeometryType(Dune::GeometryTypes::line, corners);
    }

    /*!
     * \brief Creates a geometrical plane object.
     *        This is a specialization for 3D.
     *
     * \param p0 The first corner
     * \param p1 The second corner
     * \param p2 The third corner
     * \param p3 The fourth corner
     */
    static PlaneGeometryType makePlane(const GlobalPosition& p0,
                                       const GlobalPosition& p1,
                                       const GlobalPosition& p2,
                                       const GlobalPosition& p3)
    {
        const std::vector< Dune::FieldVector< Scalar, dim > > corners = {p0, p1, p2, p3};
        return PlaneGeometryType(Dune::GeometryTypes::quadrilateral, corners);
    }

    /*!
     * \brief Calculate the fluxes over all planes for a given flux type.
     *
     * \param fluxType The flux type. This can be a lambda of the following form:
     *                 [](const auto& element,
     *                    const auto& fvGeometry,
     *                    const auto& scvf,
     *                    const Scalar velocity) { return velocity * ... ;}
     */
    template<class FluxType>
    void calculateFluxes(const FluxType& fluxType)
    {
        // make sure to reset all the values of the planes, in case this method has been called already before
        for(auto&& plane : planes_)
            plane.second.resetValues();

        // make sure not to iterate over the same dofs twice
        std::vector<bool> dofVisited(problem_.fvGridGeometry().numFaceDofs(), false);

        for(auto&& element : elements(problem_.fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(problem_.fvGridGeometry());
            fvGeometry.bindElement(element);
            for(auto && scvf : scvfs(fvGeometry))
            {
                const auto dofIdx = scvf.dofIndex();
                // do nothing of the dof was already visited
                if(dofVisited[dofIdx])
                    continue;

                dofVisited[dofIdx] = true;

                // iterate through all planes and check if the flux at the given position
                // should be accounted for in the respective plane
                for(auto&& plane : planes_)
                {
                    const auto& subPlanes = plane.second.subPlanes();

                    for(int planeIdx = 0; planeIdx < subPlanes.size(); ++planeIdx)
                    {
                        if(BoundingBoxTreeHelper<dim>::pointInGeometry(subPlanes[planeIdx], scvf.center()))
                        {
                            const Scalar velocity = sol_[faceIdx][dofIdx][0];
                            const Scalar result = fluxType(element, fvGeometry, scvf, velocity);
                            plane.second.addValue(planeIdx, result);

                            if(verbose_)
                            {
                                std::cout << "Flux at face "  << scvf.center() << " (" << plane.first << "): " << result;
                                std::cout << " (directionIndex: " << scvf.directionIndex() << "; plane boundaries: ";
                                plane.second.printPlaneBoundaries(planeIdx);
                                std::cout << ", planeIdx " << planeIdx << ")" << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }

    /*!
     * \brief Return the fluxes of the individual sub planes of a given name.
     *
     * \param name The name of the plane
     */
    auto& values(const std::string& name) const
    {
        return planes_.at(name).values();
    }

    /*!
     * \brief Return the cumulative net fluxes of a plane of a given name.
     *
     * \param name The name of the plane
     */
    Scalar netFlux(const std::string& name) const
    {
        const auto& planeResults = values(name);
        return std::accumulate(planeResults.begin(), planeResults.end(), 0.0);
    }

private:

    std::map<std::string, PlaneData<planeDim ,dim> > planes_;
    const Problem& problem_;
    const GridVariables& gridVariables_;
    const SolutionVector& sol_;
    bool verbose_;
};

} //end namespace

#endif
