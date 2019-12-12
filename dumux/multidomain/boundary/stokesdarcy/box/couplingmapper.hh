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
 * \ingroup StokesDarcyCoupling
 * \copydoc Dumux::StokesDarcyCouplingMapper
 */

#ifndef DUMUX_STOKES_DARCY_BOX_COUPLINGMAPPER_HH
#define DUMUX_STOKES_DARCY_BOX_COUPLINGMAPPER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <memory>
#include <unordered_map>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/geometry/geometryintersection.hh>
#include <dumux/common/geometry/intersectingentities.hh>

#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling mapper for Stokes and Darcy domains with equal dimension.
 */
template<class MDTraits>
class StokesDarcyCouplingMapperBox
{
    using Scalar = typename MDTraits::Scalar;

public:
    static constexpr auto stokesCellCenterIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto stokesFaceIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto cellCenterIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto faceIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto stokesIdx = stokesCellCenterIdx;
    static constexpr auto darcyIdx = typename MDTraits::template SubDomain<2>::Index();


private:
    // obtain the type tags of the sub problems
    using StokesTypeTag = typename MDTraits::template SubDomain<0>::TypeTag;
    using DarcyTypeTag = typename MDTraits::template SubDomain<2>::TypeTag;

    // sub domain grid geometries & scvf geometries
    using StokesGG = GetPropType<StokesTypeTag, Properties::GridGeometry>;
    using DarcyGG = GetPropType<DarcyTypeTag, Properties::GridGeometry>;
    using StokesScvfGeometry = typename StokesGG::SubControlVolumeFace::Traits::Geometry;
    using DarcyScvfGeometry = typename DarcyGG::SubControlVolumeFace::Traits::Geometry;

    struct ElementMapInfo
    {
        std::size_t eIdx;
        std::size_t scvfIdx;
        std::size_t flipScvfIdx;
    };

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    using CouplingManager = GetPropType<StokesTypeTag, Properties::CouplingManager>;

    static_assert(StokesGG::discMethod == DiscretizationMethod::staggered, "The free flow domain must use the staggered discretization");
    static_assert(DarcyGG::discMethod == DiscretizationMethod::box, "The Darcy domain must use the Box discretization");
public:

    /*!
     * \brief Constructor
     */
    StokesDarcyCouplingMapperBox(const CouplingManager& couplingManager)
    : couplingManager_(couplingManager)
    {}

    /*!
     * \brief Main update routine
     */
    template<class Stencils, class StencilsB>
    void computeCouplingMapsAndStencils(Stencils& darcyToStokesCellCenterStencils,
                                        StencilsB& darcyToStokesFaceStencils,
                                        Stencils& stokesCellCenterToDarcyStencils,
                                        Stencils& stokesFaceToDarcyStencils)
    {
        computeCouplingMaps();

        const auto& stokesProblem = couplingManager_.problem(stokesIdx);
        const auto& darcyProblem = couplingManager_.problem(darcyIdx);

        const auto& stokesFvGridGeometry = stokesProblem.gridGeometry();
        const auto& darcyFvGridGeometry = darcyProblem.gridGeometry();

        auto darcyFvGeometry = localView(darcyFvGridGeometry);

        for(const auto& dataHandle : stokesElementToDarcyElementMap_)
        {
            if(dataHandle.second.size() > 1)
                DUNE_THROW(Dune::InvalidStateException, "Stokes face dof should only intersect with one Darcy element");

            const auto& data = dataHandle.second[0];
            const auto stokesElementIdx = dataHandle.first;
            const auto darcyEIdx = data.eIdx;
            const auto stokesScvfIdx = data.flipScvfIdx;
            const auto& stokesScvf = stokesFvGridGeometry.scvf(stokesScvfIdx);

            const auto& darcyElement = darcyFvGridGeometry.element(darcyEIdx);
            darcyFvGeometry.bind(darcyElement);

            darcyToStokesCellCenterStencils[darcyEIdx].push_back(stokesElementIdx);
            darcyToStokesFaceStencils[darcyEIdx].first.push_back(stokesScvf.dofIndex());
            darcyToStokesFaceStencils[darcyEIdx].second.push_back(stokesScvf.index());

            for (auto&& scv : scvs(darcyFvGeometry))
            {
                stokesCellCenterToDarcyStencils[stokesElementIdx].push_back(scv.dofIndex());
                stokesFaceToDarcyStencils[stokesScvf.dofIndex()].push_back(scv.dofIndex());
            }
        }
    }

    void computeCouplingMaps()
    {
        const auto& stokesProblem = couplingManager_.problem(stokesIdx);
        const auto& darcyProblem = couplingManager_.problem(darcyIdx);

        const auto& stokesFvGridGeometry = stokesProblem.gridGeometry();
        const auto& darcyFvGridGeometry = darcyProblem.gridGeometry();

        // find all darcy faces coupling to stokes
        isCoupledDarcyScvf_.resize(darcyFvGridGeometry.gridView().size(0));
        for (const auto& darcyElement : elements(darcyFvGridGeometry.gridView()))
        {
            const auto darcyEIdx = darcyFvGridGeometry.elementMapper().index(darcyElement);
            auto darcyFvGeometry = localView(darcyFvGridGeometry);
            darcyFvGeometry.bindElement(darcyElement);

            for (const auto& darcyScvf : scvfs(darcyFvGeometry))
            {
                if (!darcyScvf.boundary())
                    continue;

                // find all stokes elements that intersect with the face
                const auto& darcyScvfGeometry = darcyScvf.geometry();
                const auto rawIntersections = intersectingEntities(darcyScvfGeometry, stokesFvGridGeometry.boundingBoxTree());
                if (rawIntersections.empty())
                    continue;

                isCoupledDarcyScvf_[darcyEIdx].assign(darcyFvGeometry.numScvf(), false);
                for (const auto& rawIntersection : rawIntersections)
                {
                    const auto stokesEIdx = rawIntersection.second();
                    const auto stokesElement = stokesFvGridGeometry.element(stokesEIdx);
                    auto stokesFvGeometry = localView(stokesFvGridGeometry);
                    stokesFvGeometry.bindElement(stokesElement);

                    for(const auto& stokesScvf : scvfs(stokesFvGeometry))
                    {
                        if(!stokesScvf.boundary())
                            continue;

                        // intersect the geometries
                        using IntersectionAlgorithm = GeometryIntersection<DarcyScvfGeometry, StokesScvfGeometry>;
                        typename IntersectionAlgorithm::Intersection is;
                        if(IntersectionAlgorithm::intersection(darcyScvfGeometry, stokesScvf.geometry(), is))
                        {
                            isCoupledDarcyScvf_[darcyEIdx][darcyScvf.index()] = true;
                            darcyElementToStokesElementMap_[darcyEIdx].push_back({stokesEIdx, stokesScvf.index(), darcyScvf.index()});
                            stokesElementToDarcyElementMap_[stokesEIdx].push_back({darcyEIdx, darcyScvf.index(), stokesScvf.index()});
                        }
                    }
                }
            }
        }
    }

    /*!
     * \brief Returns whether a Darcy scvf is coupled to the other domain
     */
    bool isCoupledDarcyScvf(std::size_t eIdx, std::size_t scvfLocalIdx) const
    {
        if(isCoupledDarcyScvf_[eIdx].size() > 0)
            return isCoupledDarcyScvf_[eIdx][scvfLocalIdx];

        return false;
    }


    /*!
     * \brief A map that returns all Stokes elements coupled to a Darcy element
     */
    const auto& darcyElementToStokesElementMap() const
    {
        return darcyElementToStokesElementMap_;
    }

    /*!
     * \brief A map that returns all Darcy elements coupled to a Stokes element
     */
    const auto& stokesElementToDarcyElementMap() const
    {
        return stokesElementToDarcyElementMap_;
    }

private:
    std::unordered_map<std::size_t, std::vector<ElementMapInfo>> darcyElementToStokesElementMap_;
    std::unordered_map<std::size_t, std::vector<ElementMapInfo>> stokesElementToDarcyElementMap_;

    std::vector<std::vector<bool>> isCoupledDarcyScvf_;

    const CouplingManager& couplingManager_;
};

} // end namespace Dumux

#endif
