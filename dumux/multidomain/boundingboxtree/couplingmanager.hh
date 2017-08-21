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
 * \ingroup BoundaryCoupling
 * \brief Coupling manager for two domains with the same dimension.
 *        Intersection computation relies on boundingBoxTree and therefore
          does not require grid-glue.
 */
#ifndef DUMUX_COUPLINGMANAGER_STOKES_DARCY_HH
#define DUMUX_COUPLINGMANAGER_STOKES_DARCY_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/boundingboxtree.hh>
#include "couplingmapper.hh"
#include <dumux/common/exceptions.hh>

namespace Dumux
{

namespace Properties
{
// Property forward declarations
NEW_PROP_TAG(StokesProblemTypeTag);
NEW_PROP_TAG(DarcyProblemTypeTag);
NEW_PROP_TAG(GridView);
} // namespace Properties

/*!
 * \brief Manages the coupling between Stokes and Darcy elements
 * \ingroup BoundaryCoupling
 */
template<class TypeTag>
class CouplingManagerStokesDarcy
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);

    // obtain the type tags of the sub problems
    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using DarcyProblemTypeTag = typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);

    using StokesGridView = typename GET_PROP_TYPE(StokesProblemTypeTag, GridView);
    using DarcyGridView = typename GET_PROP_TYPE(DarcyProblemTypeTag, GridView);

    using StokesGrid = typename GET_PROP_TYPE(StokesProblemTypeTag, Grid);
    using DarcyGrid = typename GET_PROP_TYPE(DarcyProblemTypeTag, Grid);

    using StokesProblem = typename GET_PROP_TYPE(StokesProblemTypeTag, Problem);
    using DarcyProblem = typename GET_PROP_TYPE(DarcyProblemTypeTag, Problem);

    using DarcySubControlVolume = typename GET_PROP_TYPE(DarcyProblemTypeTag, SubControlVolume);
    using DarcySubControlVolumeFace = typename GET_PROP_TYPE(DarcyProblemTypeTag, SubControlVolumeFace);
    using StokesSubControlVolume = typename GET_PROP_TYPE(StokesProblemTypeTag, SubControlVolume);
    using StokesSubControlVolumeFace = typename GET_PROP_TYPE(StokesProblemTypeTag, SubControlVolumeFace);

    using DarcyElementVolumeVariables = typename GET_PROP_TYPE(DarcyProblemTypeTag, ElementVolumeVariables);
    using StokesElementVolumeVariables = typename GET_PROP_TYPE(StokesProblemTypeTag, ElementVolumeVariables);

    using DarcyFVElementGeometry = typename GET_PROP_TYPE(DarcyProblemTypeTag, FVElementGeometry);
    using StokesFVElementGeometry = typename GET_PROP_TYPE(StokesProblemTypeTag, FVElementGeometry);

    using StokesGlobalFaceVars = typename GET_PROP_TYPE(StokesProblemTypeTag, GlobalFaceVars);

    using FacePrimaryVariables = typename GET_PROP_TYPE(StokesProblemTypeTag, FacePrimaryVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(StokesProblemTypeTag, CellCenterPrimaryVariables);
    using DarcyPrimaryVariables = typename GET_PROP_TYPE(DarcyProblemTypeTag, PrimaryVariables);

    enum {
        dim = StokesGridView::dimension,
        dimWorld = StokesGridView::dimensionworld
    };

    using StokesElement = typename StokesGridView::template Codim<0>::Entity;
    using DarcyElement = typename DarcyGridView::template Codim<0>::Entity;

    using DarcyVertex = typename DarcyGrid::template Codim<dim>::Entity;

    using CouplingMapper = Dumux::CouplingMapperStokesDarcy<TypeTag>;

    using DofTypeIndices = typename GET_PROP(StokesProblemTypeTag, DofTypeIndices);
    using cellCenterIdx = typename DofTypeIndices::CellCenterIdx;
    using faceIdx = typename DofTypeIndices::FaceIdx;

public:

    /*!
     * \brief Constructor
     */
    CouplingManagerStokesDarcy(StokesProblem& stokesProblem, DarcyProblem& darcyProblem)
    : stokesProblem_(stokesProblem),
      darcyProblem_(darcyProblem),
      couplingMapper_(stokesProblem, darcyProblem, asImp_())
    {
        fluxStokesToDarcy_ = 0.0;

        // TODO find opposite Darcy elements for pressureInDarcyElement (all elements at once?)
    }


    /*!
     * \brief Return a reference to the stokes problem
     */
    const StokesProblem& stokesProblem() const
    {
        return stokesProblem_;
    }

    /*!
     * \brief Return a reference to the low dimensional problem
     */
    const DarcyProblem& darcyProblem() const
    {
        return darcyProblem_;
    }

    /*!
     * \brief Return a reference to the stokes gridview
     */
    const StokesGridView& stokesGridView() const
    {
        return stokesProblem().gridView();
    }

    /*!
     * \brief Return a reference to the low dimensional gridview
     */
    const DarcyGridView& darcyGridView() const
    {
        return darcyProblem().gridView();
    }

    void preInit()
    {}


    /*!
     * \brief Compute the coupling maps and stencil after the sub-problem have been initialized
     */
    void postInit()
    {
        couplingMapper_.computeCouplingMaps();
        computeStencils();
    }

    /*!
     * \brief Returns whether a stokes element sub control volume is considered for coupling
     *        with the darcy domain. This function is used for the boundaryType() method
     *        of the stokes problem which handles vertices.
     *
     * \param element The element
     */
    bool isStokesCouplingEntity(const StokesElement &element) const
    {
        return !(couplingStencil(element).empty());
    }


    /*!
     * \brief Returns whether a stokes element sub control volume is considered for coupling
     *        with the darcy domain. This function is used for the boundaryType() method
     *        of the stokes problem which handles vertices.
     *
     * \param vertex The vertex
     */
    bool isStokesCouplingEntity(const StokesElement &element, const StokesSubControlVolumeFace& scvf) const
    {
        return couplingMapper_.stokesFaceToDarcyMap().count(scvf.dofIndex());
    }

    /*!
     * \brief Returns whether a low dim vertex is considered for coupling
     *        with the stokes domain. This function is used for the boundaryType() method
     *        of the low dim problem which handles vertices.
     *
     * \param vertex The vertex
     */
    bool isDarcyCouplingEntity(const DarcyVertex &vertex) const
    {
        const auto darcyDofIdxGlobal = darcyGridView().indexSet().index(vertex);
        return couplingMapper_.darcyToStokesMap().count(darcyDofIdxGlobal);
    }

    /*!
     * \brief Returns whether a low dim vertex is considered for coupling
     *        with the stokes domain. This function is used for the boundaryType() method
     *        of the low dim problem which handles vertices.
     *
     * \param element The element
     */
    bool isDarcyCouplingEntity(const DarcyElement &element) const
    {
        return !(couplingStencil(element).empty());
    }

    /*!
     * \brief Returns whether a low dim vertex is considered for coupling
     *        with the stokes domain. This function is used for the boundaryType() method
     *        of the low dim problem which handles vertices.
     *
     * \param element The element
     */
    bool isDarcyCouplingEntity(const DarcySubControlVolume &scv) const
    {
        return couplingMapper_.darcyToStokesMap().count(scv.dofIndex());
    }

    /*!
     * \brief Returns whether a low dim vertex is considered for coupling
     *        with the stokes domain. This function is used for the boundaryType() method
     *        of the low dim problem which handles vertices.
     *
     * \param element The element
     */
    bool isDarcyCouplingEntity(const DarcySubControlVolumeFace &scvf) const
    {
        return couplingMapper_.darcyToStokesMap().count(scvf.insideScvIdx()); //dofIndex()); TODO
    }

    /*!
     * \brief Returns a reference to the coupling mapper.
     */
    const CouplingMapper &couplingMapper() const
    {
        return couplingMapper_;
    }

    /*!
     * \brief Returns a reference to the coupling mapper.
     */
    const auto &darcyToStokesMap() const
    {
        return couplingMapper_.darcyToStokesMap();
    }

    /*!
     * \brief Returns a reference to the coupling mapper.
     */
    const auto &stokesFaceToDarcyMap() const
    {
        return couplingMapper_.stokesFaceToDarcyMap();
    }

    void computeStencils()
    {
        const auto &stokesTree = stokesProblem_.boundingBoxTree();
        const auto &darcyTree = darcyProblem_.boundingBoxTree();

        // compute Stokes cell-center coupling stencil using the coupling map
        for (const auto &entry : couplingMapper_.stokesCCToDarcyMap())
        {
            const auto stokesElementIdx = entry.first;

            // get the darcy DOFs associated with the Stokes element
            auto &darcyInfo = couplingMapper_.stokesCCToDarcyMap().at(stokesElementIdx);

            Dune::dverb << "Stokes element " << stokesElementIdx <<  " is coupled to Darcy element:";
            std::vector<unsigned int> darcyDofs;

            const auto& stokesElement = stokesTree.entity(stokesElementIdx);
            auto stokesFvGeometry = localView(stokesProblem_.model().globalFvGeometry());
            stokesFvGeometry.bind(stokesElement);

            for (const auto &i :darcyInfo) // TODO here only one element
            {
                auto darcyElemIdx = i.darcyElementIdx;
                Dune::dverb << " " << darcyElemIdx ;
                stokesCCCouplingStencils_[stokesElementIdx].push_back(darcyElemIdx);


                for(auto&& stokesScvf : scvfs(stokesFvGeometry))
                    stokesFaceCouplingStencils_[stokesScvf.dofIndex()].push_back(darcyElemIdx);
            }
            Dune::dverb << std::endl;
        }

        Dune::dverb << std::endl;

        // compute Darcy coupling stencil using the coupling map
        for (const auto &entry : couplingMapper_.darcyToStokesMap())
        {
            const auto darcyElementIdx = entry.first;

            Dune::dverb << "Darcy element " << darcyElementIdx <<  " is coupled to Stokes element:";

            const auto& darcyElement = darcyTree.entity(darcyElementIdx);
            auto darcyFvGeometry = localView(darcyProblem_.model().globalFvGeometry());
            darcyFvGeometry.bind(darcyElement);

            const auto &stokesElements = couplingMapper_.darcyToStokesMap().at(darcyElementIdx).stokesElementIndices();

            for(const auto stokesElemIdx: stokesElements) // TODO here: only one element
            {
                const auto stokesElement = stokesTree.entity(stokesElemIdx);
                StokesFVElementGeometry fvGeometry = localView(stokesProblem_.model().globalFvGeometry());
                fvGeometry.bind(stokesElement);

                Dune::dverb << " Stokes element " << stokesElemIdx ;

                darcyToCCCouplingStencils_[darcyElementIdx].push_back(stokesElemIdx);

                for(auto&& stokesScvf : scvfs(fvGeometry))
                    darcyToFaceCouplingStencils_[darcyElementIdx].push_back(stokesScvf.dofIndex());
            }
            Dune::dverb << std::endl;
        }

        auto sortAndMakeUnique = [](auto&& stencils)
        {
            for(auto&& stencil : stencils)
            {
                std::sort(stencil.second.begin(), stencil.second.end());
                stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
            }
        };

        // sort and make unique
        sortAndMakeUnique(stokesCCCouplingStencils_);
        sortAndMakeUnique(stokesFaceCouplingStencils_);
        sortAndMakeUnique(darcyToCCCouplingStencils_);
        sortAndMakeUnique(darcyToFaceCouplingStencils_);
    }

    /*!
     * \brief Returns a coupling stencil for a Stokes element
     */
    const std::vector<unsigned int>& couplingStencil(const StokesElement& element) const
    {
        const unsigned int eIdx = stokesProblem().elementMapper().index(element);
        if (stokesCCCouplingStencils_.count(eIdx))
            return stokesCCCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Returns a coupling stencil for a Stokes element
     */
    const std::vector<unsigned int>& couplingStencil(const StokesSubControlVolumeFace& scvf) const
    {
        const unsigned int dofIdx = scvf.dofIndex();
        if (stokesFaceCouplingStencils_.count(dofIdx))
            return stokesFaceCouplingStencils_.at(dofIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Returns a coupling stencil for a Darcy element
     */
    const std::vector<unsigned int>& couplingStencil(const DarcyElement& element, cellCenterIdx) const
    {
        const unsigned int eIdx = darcyProblem().elementMapper().index(element);
        if (darcyToCCCouplingStencils_.count(eIdx))
            return darcyToCCCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Returns a coupling stencil for a Darcy element
     */
    const std::vector<unsigned int>& couplingStencil(const DarcyElement& element, faceIdx) const
    {
        const unsigned int eIdx = darcyProblem().elementMapper().index(element);
        if (darcyToFaceCouplingStencils_.count(eIdx))
            return darcyToFaceCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    //! evaluate coupling residual for the derivative Stokes DOF with respect to Darcy DOF
    auto evalStokesCCCouplingResidual(const StokesElement& element,
                              const StokesFVElementGeometry& fvGeometry,
                              const StokesElementVolumeVariables& elemVolVars,
                              const StokesGlobalFaceVars& globalFaceVars)
    {
        Scalar couplingScvfIdx = couplingSubControlVolumeFace(fvGeometry);
        Scalar outerNormalScalar = fvGeometry.scvf(couplingScvfIdx).outerNormalScalar();
        Scalar interfaceArea = fvGeometry.scvf(couplingScvfIdx).area();
        // TODO volVars = elemVolVars[scv], only 1 scv --> [0]?
        Scalar densityStokes = elemVolVars[0].density(); // TODO upwind?
        const Scalar velocity = globalFaceVars.faceVars(couplingScvfIdx).velocity();

        // mass coupling condition
        CellCenterPrimaryVariables stokesCCCouplingResidual(0.0);
        stokesCCCouplingResidual = densityStokes * velocity * outerNormalScalar * interfaceArea; // TODO nc --> rho_g^ff

        // TODO temperature

        // TODO transport

        fluxStokesToDarcy_ = stokesCCCouplingResidual;
        return stokesCCCouplingResidual; // fluxStokesCCToDarcy
    }

    //! evaluate coupling residual for the derivative stokes DOF with respect to low dim DOF
    auto evalStokesFaceCouplingResidual(const StokesElement& element,
                              const StokesFVElementGeometry& fvGeometry,
                              const StokesGlobalFaceVars& globalFaceVars,
                              const StokesSubControlVolumeFace& scvf)
    {
        FacePrimaryVariables stokesFaceCouplingResidual(0.0);
        if(stokesProblem_.onCouplingInterface(scvf.center()))
        {
            Scalar couplingScvfIdx = scvf.index();
            Scalar outerNormalScalar = fvGeometry.scvf(couplingScvfIdx).outerNormalScalar();
            Scalar interfaceArea = fvGeometry.scvf(couplingScvfIdx).area();
            Scalar pressureDarcy = pressureInDarcyElement(scvf);

            Scalar distanceCenters = (fvGeometry.scv(scvf.insideScvIdx()).center() - scvf.center())[1]; // TODO direction? 0,1,2?
            Scalar permeability = darcyProblem_.spatialParams().permeabilityAtPos(scvf.center());
            Scalar alphaBeaversJoseph = darcyProblem_.spatialParams().beaversJosephCoeffAtPos(scvf.center());
            Scalar beta = -1.0 * std::sqrt(permeability) / alphaBeaversJoseph / distanceCenters * outerNormalScalar;

            auto stokesElemVolVars = localView(stokesProblem_.model().curGlobalVolVars());
            stokesElemVolVars.bind(element, fvGeometry, stokesProblem_.model().curSol());
            const auto stokesVolVars = stokesElemVolVars[scvf.insideScvIdx()];

            Scalar effDynViscosity = stokesVolVars.viscosity();
//            Scalar effDynViscosity = BaseFluid::dynamicViscosity(stokesVolVars.pressure(), stokesVolVars.temperature(), stokesVolVars.massFrac());
//                                                                + (eddyKinematicViscosity_(elementInsideIdx) * stokesVolVars.density()); TODO
            const Scalar tangXVelocity = globalFaceVars.faceVars(scvf.dofIndex()).velocity();
            Scalar beaversJosephXVelocity = beta * tangXVelocity / (1.0 + beta);
//            Scalar tangZVelocity = faceSolutionVector[velocityZIdx];
//            Scalar beaversJosephZVelocity = beta * tangZVelocity / (1.0 + beta);

            // normal momentum coupling condition
            stokesFaceCouplingResidual[stokesProblem_.velocityYIdx] = pressureDarcy * outerNormalScalar * interfaceArea;

            // TODO only for 2D, add loop over tangDims
            // Fetzer2017 "classical BJ condition"
            // tangential momentum coupling condition --> b.c. for stokesProblem?
//            stokesFaceCouplingResidual[stokesProblem_.velocityXIdx] = 0.0; // TODO
            stokesFaceCouplingResidual[stokesProblem_.velocityXIdx] = -0.5 * effDynViscosity * (tangXVelocity - beaversJosephXVelocity)
                                                                           / distanceCenters * outerNormalScalar * interfaceArea;
            // TODO
//            stokesFaceCouplingResidual[stokesProblem_.velocityZIdx] = -0.5 * effDynViscosity * (tangZVelocity - beaversJosephZVelocity)
//                                                                           / distanceCenters * outerNormalScalar * interfaceArea;
        }
        return stokesFaceCouplingResidual;
    }

    //! evaluate coupling residual for the derivative Darcy DOF with respect to Stokes DOF
    auto evalDarcyCouplingResidual(const DarcyElement& element)
    {
        return (-1.0*fluxStokesToDarcy_); // fluxDarcyToStokesCC
    }

protected:
    // ! Returns the global index of the subcontrolvolumeface at the coupling interface
    const Scalar couplingSubControlVolumeFace(const StokesFVElementGeometry& fvGeometry)
    {
        Scalar couplingScvfIdx;
        // assumption: only one interface scvf (no corners in coupling interface)
        for(const auto &scvf : scvfs(fvGeometry))
        {
            if(stokesProblem_.onCouplingInterface(scvf.center()))
                couplingScvfIdx = scvf.index();
        }
        return couplingScvfIdx;
    }

    // ! Returns the pressure value of the adjacent Darcy element
    const Scalar pressureInDarcyElement(const StokesSubControlVolumeFace& scvf)
    {
        const auto& darcyTree = darcyProblem_.boundingBoxTree();

        // create a vector containing all Darcy elements coupled to the Stokes scvf
        const auto darcyCouplingInfo = stokesFaceToDarcyMap().at(scvf.dofIndex());

        // TODO assumption: same number of elements in each subdomain, only one coupling element on opposite side of interface
        const auto& darcyCouplingElement = darcyTree.entity(darcyCouplingInfo[0].darcyElementIdx);
        DarcyFVElementGeometry darcyFvGeometry = localView(darcyProblem_.model().globalFvGeometry());
        darcyFvGeometry.bind(darcyCouplingElement);

        const auto darcyDofIdx = darcyCouplingInfo[0].darcyDofIdx;

        auto darcyElemVolVars = localView(darcyProblem_.model().curGlobalVolVars());
        darcyElemVolVars.bind(darcyCouplingElement, darcyFvGeometry, darcyProblem_.model().curSol());
        const auto darcyVolVars = darcyElemVolVars[darcyDofIdx];

        return darcyVolVars.pressure();
    }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

private:
    StokesProblem& stokesProblem_;
    DarcyProblem& darcyProblem_;
    CouplingMapper couplingMapper_;

    DarcyPrimaryVariables fluxStokesToDarcy_;

    std::unordered_map<unsigned int, std::vector<unsigned int> > stokesCCCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > stokesFaceCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > darcyToCCCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > darcyToFaceCouplingStencils_;
    std::vector<unsigned int> emptyStencil_;
};

} // namespace Dumux

#endif
