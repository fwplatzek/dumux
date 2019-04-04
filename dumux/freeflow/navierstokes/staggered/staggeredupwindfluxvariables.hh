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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesUpwindVariables
 */
#ifndef DUMUX_NAVIERSTOKES_STAGGERED_UPWINDVARIABLES_HH
#define DUMUX_NAVIERSTOKES_STAGGERED_UPWINDVARIABLES_HH

#include <array>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/method.hh>
#include <dumux/freeflow/staggeredupwindmethods.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The upwinding variables class for the Navier-Stokes model using the staggered grid discretization.
 */
template<class TypeTag, int upwindSchemeOrder>
class StaggeredUpwindFluxVariables
{
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;

    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;
    using FaceVariables = typename GridFaceVariables::FaceVariables;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FacePrimaryVariables = GetPropType<TypeTag, Properties::FacePrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;

    static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();
    static constexpr auto faceIdx = FVGridGeometry::faceIdx();

public:
    /*!
     * \brief Returns the momentum in the frontal directon.
     *
     *        Checks if the model has higher order methods enabled and if the scvf in
     *        question is far enough from the boundary such that higher order methods can be employed.
     *        Then the corresponding set of momenta are collected and the prescribed
     *        upwinding method is used to calculate the momentum.
     */
    static FacePrimaryVariables computeUpwindedFrontalMomentum(const SubControlVolumeFace& scvf,
                                                               const ElementFaceVariables& elemFaceVars,
                                                               const ElementVolumeVariables& elemVolVars,
                                                               const GridFluxVariablesCache& gridFluxVarsCache,
                                                               const Scalar transportingVelocity)
    {
        Scalar momentum(0.0);
        const bool canHigherOrder = canFrontalSecondOrder_(scvf, transportingVelocity, std::integral_constant<bool, useHigherOrder>{});

        if (canHigherOrder)
        {
            const auto upwindingMomenta = getFrontalUpwindingMomenta_(scvf, elemFaceVars, elemVolVars[scvf.insideScvIdx()].density(),
                                                                      transportingVelocity, std::integral_constant<bool, useHigherOrder>{});
            momentum = doFrontalMomentumUpwinding_(scvf, upwindingMomenta, transportingVelocity, gridFluxVarsCache);
        }
        else
        {
            const auto upwindingMomenta = getFrontalUpwindingMomenta_(scvf, elemFaceVars, elemVolVars[scvf.insideScvIdx()].density(),
                                                                      transportingVelocity, std::integral_constant<bool, false>{});
            momentum = doFrontalMomentumUpwinding_(scvf, upwindingMomenta, transportingVelocity, gridFluxVarsCache);
        }

        return momentum;
    }

    /*!
     * \brief Returns the momentum in the lateral directon.
     *
     *        Evaluates which face is upstream.
     *        Checks if the model has higher order methods enabled and if the scvf in
     *        question is far enough from the boundary such that higher order methods can be employed.
     *        Then the corresponding set of momenta are collected and the prescribed
     *        upwinding method is used to calculate the momentum.
     */
    static FacePrimaryVariables computeUpwindedLateralMomentum(const Problem& problem,
                                                               const FVElementGeometry& fvGeometry,
                                                               const Element& element,
                                                               const SubControlVolumeFace& scvf,
                                                               const SubControlVolumeFace& normalFace,
                                                               const ElementVolumeVariables& elemVolVars,
                                                               const FaceVariables& faceVars,
                                                               const GridFluxVariablesCache& gridFluxVarsCache,
                                                               const int localSubFaceIdx,
                                                               const bool lateralFaceHasDirichletPressure,
                                                               const bool lateralFaceHasBJS)
    {
        // Get the transporting velocity, located at the scvf perpendicular to the current scvf where the dof
        // of interest is located.
        const Scalar transportingVelocity = faceVars.velocityNormalInside(localSubFaceIdx);

        // Check whether the own or the neighboring element is upstream.
        const bool selfIsUpstream = ( normalFace.directionSign() == sign(transportingVelocity) );

        Scalar momentum(0.0);
        const bool canHigherOrder = canLateralSecondOrder_(scvf, selfIsUpstream, localSubFaceIdx, std::integral_constant<bool, useHigherOrder>{});

        if (canHigherOrder)
        {
            const auto parallelUpwindingMomenta = getLateralUpwindingMomenta_(problem, fvGeometry, element, scvf, elemVolVars, faceVars,
                                                                              transportingVelocity, localSubFaceIdx,
                                                                              lateralFaceHasDirichletPressure, lateralFaceHasBJS,
                                                                              std::integral_constant<bool, useHigherOrder>{});
            momentum = doLateralMomentumUpwinding_(scvf, normalFace, parallelUpwindingMomenta, transportingVelocity, localSubFaceIdx, gridFluxVarsCache);
        }
        else
        {
            const auto parallelUpwindingMomenta = getLateralUpwindingMomenta_(problem, fvGeometry, element, scvf, elemVolVars, faceVars,
                                                                              transportingVelocity, localSubFaceIdx,
                                                                              lateralFaceHasDirichletPressure, lateralFaceHasBJS,
                                                                              std::integral_constant<bool, false>{});
            momentum = doLateralMomentumUpwinding_(scvf, normalFace, parallelUpwindingMomenta, transportingVelocity, localSubFaceIdx, gridFluxVarsCache);
        }

        return momentum;
    }

private:

  /*!
     * \brief Returns false as higher order methods are not enabled.
     *
     *        If higher order methods are not prescribed, this function will be called and false is returned..
     *
     * \param ownScvf the sub control volume face
     * \param transportingVelocity The average of the self and opposite velocities.
     * \param false_type False if higher order methods are not enabled.
     */
    static bool canFrontalSecondOrder_(const SubControlVolumeFace& ownScvf,
                                       const Scalar transportingVelocity,
                                       std::false_type)
    {
        return false;
    }

    /*!
     * \brief Returns whether or not the face in question is far enough from the wall to handle higher order methods.
     *
     *        Evaluates which face is upstream.
     *        If the face is upstream, and the scvf has a forward neighbor, higher order methods are possible.
     *        If the face is downstream, and the scvf has a backwards neighbor, higher order methods are possible.
     *        Otherwise, higher order methods are not possible.
     *        If higher order methods are not prescribed, this function will not be called.
     *
     * \param ownScvf the sub control volume face
     * \param transportingVelocity The average of the self and opposite velocities.
     * \param true_type True if higher order methods are enabled.
     */
    static bool canFrontalSecondOrder_(const SubControlVolumeFace& ownScvf,
                                       const Scalar transportingVelocity,
                                       std::true_type)
    {
        // Depending on selfIsUpstream I have to check if I have a forward or a backward neighbor to retrieve
        const bool selfIsUpstream = ownScvf.directionSign() != sign(transportingVelocity);

        if (selfIsUpstream)
        {
            if (ownScvf.hasForwardNeighbor(0))
                return true;
            else
                return false;
        }
        else
        {
            if (ownScvf.hasBackwardNeighbor(0))
                return true;
            else
                return false;
        }
    }

    /*!
     * \brief Returns an array of the two momenta needed for basic upwinding methods.
     *
     * \param scvf The sub control volume face
     * \param elemFaceVars The element face variables
     * \param density The given density \f$\mathrm{[kg/m^3]}\f$
     * \param transportingVelocity The average of the self and opposite velocities.
     * \param false_type False if higher order methods are not enabled.
     */
    static std::array<Scalar, 2> getFrontalUpwindingMomenta_(const SubControlVolumeFace& scvf,
                                                             const ElementFaceVariables& elemFaceVars,
                                                             const Scalar density,
                                                             const Scalar transportingVelocity,
                                                             std::false_type)
    {
        const Scalar momentumSelf = elemFaceVars[scvf].velocitySelf() * density;
        const Scalar momentumOpposite = elemFaceVars[scvf].velocityOpposite() * density;
        const bool selfIsUpstream = scvf.directionSign() != sign(transportingVelocity);

        return selfIsUpstream ? std::array<Scalar, 2>{momentumOpposite, momentumSelf}
                              : std::array<Scalar, 2>{momentumSelf, momentumOpposite};
    }

    /*!
     * \brief Returns an array of the three momenta needed for higher order upwinding methods.
     *
     * \param scvf The sub control volume face
     * \param elemFaceVars The element face variables
     * \param density The given density \f$\mathrm{[kg/m^3]}\f$
     * \param transportingVelocity The average of the self and opposite velocities.
     * \param true_type True if higher order methods are enabled.
     */
    static std::array<Scalar, 3> getFrontalUpwindingMomenta_(const SubControlVolumeFace& scvf,
                                                             const ElementFaceVariables& elemFaceVars,
                                                             const Scalar density,
                                                             const Scalar transportingVelocity,
                                                             std::true_type)
    {
        const bool selfIsUpstream = scvf.directionSign() != sign(transportingVelocity);
        std::array<Scalar, 3> momenta{0.0, 0.0, 0.0};

        if (selfIsUpstream)
        {
            momenta[0] = elemFaceVars[scvf].velocityOpposite() * density;
            momenta[1] = elemFaceVars[scvf].velocitySelf() * density;
            momenta[2] = elemFaceVars[scvf].velocityForward(0) * density;
        }
        else
        {
            momenta[0] = elemFaceVars[scvf].velocitySelf() * density;
            momenta[1] = elemFaceVars[scvf].velocityOpposite() * density;
            momenta[2] = elemFaceVars[scvf].velocityBackward(0) * density;
        }

        return momenta;
    }

    /*!
     * \brief Returns the upwinded momentum for basic upwind schemes
     *
     * \param scvf The sub control volume face
     * \param momenta The momenta to be upwinded
     * \param transportingVelocity The average of the self and opposite velocities.
     * \param gridFluxVarsCache The grid flux variables cache
     */
    static Scalar doFrontalMomentumUpwinding_(const SubControlVolumeFace& scvf,
                                              const std::array<Scalar, 2>& momenta,
                                              const Scalar transportingVelocity,
                                              const GridFluxVariablesCache& gridFluxVarsCache)
    {
        const auto& upwindScheme = gridFluxVarsCache.staggeredUpwindMethods();
        return upwindScheme.upwind(momenta[0], momenta[1]);
    }

    /*!
     * \brief Returns the upwinded momentum for higher order upwind schemes
     *
     * \param scvf The sub control volume face
     * \param momenta The momenta to be upwinded
     * \param transportingVelocity The average of the self and opposite velocities.
     * \param gridFluxVarsCache The grid flux variables cache
     */
    template<bool enable = useHigherOrder, std::enable_if_t<enable, int> = 0>
    static Scalar doFrontalMomentumUpwinding_(const SubControlVolumeFace& scvf,
                                              const std::array<Scalar, 3>& momenta,
                                              const Scalar transportingVelocity,
                                              const GridFluxVariablesCache& gridFluxVarsCache)
    {
        const auto& upwindScheme = gridFluxVarsCache.staggeredUpwindMethods();
        const bool selfIsUpstream = scvf.directionSign() != sign(transportingVelocity);
        std::array<Scalar,3> distances = getFrontalDistances_(scvf, selfIsUpstream);
        return upwindScheme.tvd(momenta, distances, selfIsUpstream, upwindScheme.tvdApproach());
    }

    /*!
     * \brief Returns an array of distances needed for non-uniform higher order upwind schemes
     *
     * \param ownScvf The sub control volume face
     * \param selfIsUpstream bool describing upstream face.
     */
    template<bool enable = useHigherOrder, std::enable_if_t<enable, int> = 0>
    static std::array<Scalar, 3> getFrontalDistances_(const SubControlVolumeFace& ownScvf,
                                                      const bool selfIsUpstream)
    {
        // Depending on selfIsUpstream the downstream and the (up)upstream distances are saved.
        // distances {upstream to downstream distance, up-upstream to upstream distance, downstream staggered cell size}
        std::array<Scalar, 3> distances{0.0, 0.0, 0.0};

        if (selfIsUpstream)
        {
            distances[0] = ownScvf.selfToOppositeDistance();
            distances[1] = ownScvf.axisData().inAxisForwardDistances[0];
            distances[2] = 0.5 * (ownScvf.axisData().selfToOppositeDistance + ownScvf.axisData().inAxisBackwardDistances[0]);
        }
        else
        {
            distances[0] = ownScvf.selfToOppositeDistance();
            distances[1] = ownScvf.axisData().inAxisBackwardDistances[0];
            distances[2] = 0.5 * (ownScvf.axisData().selfToOppositeDistance + ownScvf.axisData().inAxisForwardDistances[0]);
        }
        return distances;
    }

    /*!
     * \brief Check if a second order approximation for the lateral part of the advective term can be used
     *
     * This helper function checks if the scvf of interest is not too near to the
     * boundary so that a dof upstream with respect to the upstream dof is available.
     *
     * \param ownScvf The SubControlVolumeFace we are considering
     * \param selfIsUpstream @c true if the velocity at ownScvf is upstream wrt the transporting velocity
     * \param localSubFaceIdx The local subface index
     * \param true_type True when second order is enabled.
     */
    static bool canLateralSecondOrder_(const SubControlVolumeFace& ownScvf,
                                       const bool selfIsUpstream,
                                       const int localSubFaceIdx,
                                       std::true_type)
    {
        if (selfIsUpstream)
        {
            // The self velocity is upstream. The downstream velocity can be assigned or retrieved
            // from the boundary, even if there is no parallel neighbor.
            return true;
        }
        else
        {
            // The self velocity is downstream. If there is no parallel neighbor I cannot use a second order approximation.
            if (!ownScvf.hasParallelNeighbor(localSubFaceIdx, 0))
                return false;
            else
                return true;
        }
    }

    /*!
     * \brief Check if a second order approximation for the lateral part of the advective term can be used
     *        If higher order methods are not prescribed, this is called and false is returned.
     */
    static bool canLateralSecondOrder_(const SubControlVolumeFace& ownScvf,
                                       const bool selfIsUpstream,
                                       const int localSubFaceIdx,
                                       std::false_type)
    {   return false; }


    /*!
     * \brief Returns an array of the three momenta needed for higher order upwinding methods.
     *
     *  Only called if higher order methods are enabled and the scvf can use higher order methods.
     */
    static std::array<Scalar, 3> getLateralUpwindingMomenta_(const Problem& problem,
                                                             const FVElementGeometry& fvGeometry,
                                                             const Element& element,
                                                             const SubControlVolumeFace& ownScvf,
                                                             const ElementVolumeVariables& elemVolVars,
                                                             const FaceVariables& faceVars,
                                                             const Scalar transportingVelocity,
                                                             const int localSubFaceIdx,
                                                             bool lateralFaceHasDirichletPressure,
                                                             bool lateralFaceHasBJS,
                                                             std::true_type)
    {
        std::array<Scalar, 3> momenta{0.0, 0.0, 0.0};
        const SubControlVolumeFace& lateralFace = fvGeometry.scvf(ownScvf.insideScvIdx(), ownScvf.pairData(localSubFaceIdx).localNormalFaceIdx);

        // Check whether the own or the neighboring element is upstream.
        const bool selfIsUpstream = lateralFace.directionSign() == sign(transportingVelocity);

        // Get the volume variables of the own and the neighboring element
        const auto& insideVolVars = elemVolVars[lateralFace.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[lateralFace.outsideScvIdx()];

        // The local index of the faces that is opposite to localSubFaceIdx
        const int oppositeSubFaceIdx = localSubFaceIdx % 2 ? localSubFaceIdx - 1 : localSubFaceIdx + 1;

        if (selfIsUpstream)
        {
            momenta[1] = faceVars.velocitySelf() * insideVolVars.density();

            if(ownScvf.hasParallelNeighbor(localSubFaceIdx, 0))
                momenta[0] = faceVars.velocityParallel(localSubFaceIdx, 0) * insideVolVars.density();
            else
                momenta[0] = getParallelVelocityFromBoundary_(problem, ownScvf, lateralFace,
                                                              faceVars.velocitySelf(), localSubFaceIdx,
                                                              element, lateralFaceHasDirichletPressure,
                                                              lateralFaceHasBJS) * insideVolVars.density();

            // The "upstream-upstream" velocity is retrieved from the other parallel neighbor or from the boundary.
            if (ownScvf.hasParallelNeighbor(oppositeSubFaceIdx, 0))
                momenta[2] = faceVars.velocityParallel(oppositeSubFaceIdx, 0) * insideVolVars.density();
            else
                momenta[2] = getParallelVelocityFromOtherBoundary_(problem, fvGeometry, ownScvf,
                                                                   oppositeSubFaceIdx, element,
                                                                   (momenta[1]/insideVolVars.density())) * insideVolVars.density();
        }
        else
        {
            momenta[0] = faceVars.velocitySelf() * outsideVolVars.density();
            momenta[1] = faceVars.velocityParallel(localSubFaceIdx, 0) * outsideVolVars.density();

            // If there is another parallel neighbor I can assign the "upstream-upstream" velocity, otherwise I retrieve it from the boundary.
            if (ownScvf.hasParallelNeighbor(localSubFaceIdx, 1))
                momenta[2] = faceVars.velocityParallel(localSubFaceIdx, 1) * outsideVolVars.density();
            else
            {
                const Element& elementParallel = fvGeometry.fvGridGeometry().element(fvGeometry.scv(lateralFace.outsideScvIdx()));
                const SubControlVolumeFace& firstParallelScvf = fvGeometry.scvf(lateralFace.outsideScvIdx(), ownScvf.localFaceIdx());
                momenta[2] = getParallelVelocityFromOtherBoundary_(problem, fvGeometry, firstParallelScvf,
                                                                   localSubFaceIdx, elementParallel,
                                                                   (momenta[1]/outsideVolVars.density())) * outsideVolVars.density();
            }
        }
        return momenta;

    }

    /*!
     * \brief Returns an array of the two momenta needed for basic upwinding methods.
     *
     *  Called if higher order methods are not enabled of if the scvf can not use higher order methods.
     */
    static std::array<Scalar, 2> getLateralUpwindingMomenta_(const Problem& problem,
                                                             const FVElementGeometry& fvGeometry,
                                                             const Element& element,
                                                             const SubControlVolumeFace& ownScvf,
                                                             const ElementVolumeVariables& elemVolVars,
                                                             const FaceVariables& faceVars,
                                                             const Scalar transportingVelocity,
                                                             const int localSubFaceIdx,
                                                             bool lateralFaceHasDirichletPressure,
                                                             bool lateralFaceHasBJS,
                                                             std::false_type)
    {
         // Check whether the own or the neighboring element is upstream.
        const SubControlVolumeFace& lateralFace = fvGeometry.scvf(ownScvf.insideScvIdx(), ownScvf.pairData(localSubFaceIdx).localNormalFaceIdx);
        const bool selfIsUpstream = lateralFace.directionSign() == sign(transportingVelocity);

        // Get the volume variables of the own and the neighboring element
        const auto& insideVolVars = elemVolVars[lateralFace.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[lateralFace.outsideScvIdx()];

        const Scalar momentumSelf = faceVars.velocitySelf() * insideVolVars.density();

        const Scalar momentumParallel = ownScvf.hasParallelNeighbor(localSubFaceIdx, 0)
                                      ? faceVars.velocityParallel(localSubFaceIdx, 0) * outsideVolVars.density()
                                      : (getParallelVelocityFromBoundary_(problem, ownScvf, lateralFace,
                                                                         faceVars.velocitySelf(), localSubFaceIdx, element,
                                                                         lateralFaceHasDirichletPressure, lateralFaceHasBJS)
                                      * insideVolVars.density());

        return selfIsUpstream ? std::array<Scalar, 2>{momentumParallel, momentumSelf}
                              : std::array<Scalar, 2>{momentumSelf, momentumParallel};
    }

    /*!
     * \brief Returns the upwinded momentum for basic upwind schemes
     *
     *        Fowards to Frontal Momentum Upwinding method
     */
    static Scalar doLateralMomentumUpwinding_(const SubControlVolumeFace& scvf,
                                              const SubControlVolumeFace& normalScvf,
                                              const std::array<Scalar, 2>& momenta,
                                              const Scalar transportingVelocity,
                                              const int localSubFaceIdx,
                                              const GridFluxVariablesCache& gridFluxVarsCache)
    {
        return doFrontalMomentumUpwinding_(scvf, momenta, transportingVelocity, gridFluxVarsCache);
    }

    /*!
     * \brief Returns the upwinded momentum for higher order upwind schemes
     *
     * \param scvf The sub control volume face
     * \param normalScvf The normal sub control volume face
     * \param momenta The momenta to be upwinded
     * \param transportingVelocity The average of the self and opposite velocities.
     * \param localSubFaceIdx  The local index of the subface
     * \param gridFluxVarsCache The grid flux variables cache
     */
    static Scalar doLateralMomentumUpwinding_(const SubControlVolumeFace& scvf,
                                              const SubControlVolumeFace& normalScvf,
                                              const std::array<Scalar, 3>& momenta,
                                              const Scalar transportingVelocity,
                                              const int localSubFaceIdx,
                                              const GridFluxVariablesCache& gridFluxVarsCache)
    {
        const bool selfIsUpstream = ( normalScvf.directionSign() == sign(transportingVelocity) );
        const auto& upwindScheme = gridFluxVarsCache.staggeredUpwindMethods();
        std::array<Scalar,3> distances = getLateralDistances_(scvf, localSubFaceIdx, selfIsUpstream);
        return upwindScheme.tvd(momenta, distances, selfIsUpstream, upwindScheme.tvdApproach());
    }

    /*!
     * \brief Returns an array of distances needed for non-uniform higher order upwind schemes
     *
     * \param ownScvf The sub control volume face
     * \param localSubFaceIdx  The local index of the subface
     * \param selfIsUpstream bool describing upstream face.
     */
    template<bool enable = useHigherOrder, std::enable_if_t<enable, int> = 0>
    static std::array<Scalar, 3> getLateralDistances_(const SubControlVolumeFace& ownScvf,
                                                      const int localSubFaceIdx,
                                                      const bool selfIsUpstream)
    {
        // distances {upstream to downstream distance, up-upstream to upstream distance, downstream staggered cell size}
        std::array<Scalar, 3> distances{0.0, 0.0, 0.0};

        // The local index of the faces that is opposite to localSubFaceIdx
        const int oppositeSubFaceIdx = localSubFaceIdx % 2 ? localSubFaceIdx - 1 : localSubFaceIdx + 1;

        if (selfIsUpstream)
        {
            distances[0] = ownScvf.cellCenteredParallelDistance(localSubFaceIdx, 0);
            distances[1] = ownScvf.cellCenteredParallelDistance(oppositeSubFaceIdx, 0);
            if(ownScvf.hasParallelNeighbor(localSubFaceIdx, 0))
                distances[2] = ownScvf.pairData(localSubFaceIdx).parallelCellWidths[0];
            else
                distances[2] = ownScvf.area() / 2.0;
        }
        else
        {
            distances[0] = ownScvf.cellCenteredParallelDistance(localSubFaceIdx, 0);
            distances[1] = ownScvf.cellCenteredParallelDistance(localSubFaceIdx, 1);
            distances[2] = ownScvf.pairData(localSubFaceIdx).parallelCellWidths[0];
        }

        return distances;
    }

    /*!
     * \brief Return the outer parallel velocity for normal faces that are on the boundary and therefore have no neighbor.
     *
     * Calls the problem to retrieve a fixed value set on the boundary.
     *
     * \param problem The problem
     * \param scvf The SubControlVolumeFace that is normal to the boundary
     * \param normalFace The face at the boundary
     * \param velocitySelf the velocity at scvf
     * \param localSubFaceIdx The local index of the face that is on the boundary
     * \param element The element that is on the boundary
     * \param lateralFaceHasDirichletPressure @c true if there is a dirichlet condition for the pressure on the boundary
     * \param lateralFaceHasBJS @c true if there is a BJS condition fot the velocity on the boundary
     */
    static Scalar getParallelVelocityFromBoundary_(const Problem& problem,
                                                   const SubControlVolumeFace& scvf,
                                                   const SubControlVolumeFace& normalFace,
                                                   const Scalar velocitySelf,
                                                   const int localSubFaceIdx,
                                                   const Element& element,
                                                   const bool lateralFaceHasDirichletPressure,
                                                   const bool lateralFaceHasBJS)
    {
        // If there is a Dirichlet condition for the pressure we assume zero gradient for the velocity,
        // so the velocity at the boundary equal to that on the scvf.
        if (lateralFaceHasDirichletPressure)
            return velocitySelf;

        const auto ghostFace = makeParallelGhostFace_(scvf, localSubFaceIdx);
        if (lateralFaceHasBJS)
            return problem.bjsVelocity(element, scvf, normalFace, localSubFaceIdx, velocitySelf);
        return problem.dirichlet(element, ghostFace)[Indices::velocity(scvf.directionIndex())];
    }

    /*!
     * \brief Return a velocity value from a boundary for which the boundary conditions have to be checked.
     *
     * \param problem The problem
     * \param scvf The SubControlVolumeFace that is normal to the boundary
     * \param localIdx The local index of the face that is on the boundary
     * \param boundaryElement The element that is on the boundary
     * \param parallelVelocity The velocity over scvf
     */
    static Scalar getParallelVelocityFromOtherBoundary_(const Problem& problem,
                                                        const FVElementGeometry& fvGeometry,
                                                        const SubControlVolumeFace& scvf,
                                                        const int localIdx,
                                                        const Element& boundaryElement,
                                                        const Scalar parallelVelocity)
    {
        // A ghost subface at the boundary is created, featuring the location of the sub face's center
        const SubControlVolumeFace& boundaryNormalFace = fvGeometry.scvf(scvf.insideScvIdx(), scvf.pairData(localIdx).localNormalFaceIdx);
        GlobalPosition boundarySubFaceCenter = scvf.pairData(localIdx).virtualFirstParallelFaceDofPos + boundaryNormalFace.center();
        boundarySubFaceCenter *= 0.5;
        const SubControlVolumeFace boundarySubFace = makeGhostFace_(boundaryNormalFace, boundarySubFaceCenter);

        // The boundary condition is checked, in case of symmetry or Dirichlet for the pressure
        // a gradient of zero is assumed in the direction normal to the bounadry, while if there is
        // Dirichlet of BJS for the velocity the related values are exploited.
        const auto bcTypes = problem.boundaryTypes(boundaryElement, boundarySubFace);

        if (bcTypes.isDirichlet(Indices::velocity(scvf.directionIndex())))
        {
            const SubControlVolumeFace ghostFace = makeParallelGhostFace_(scvf, localIdx);
            return problem.dirichlet(boundaryElement, ghostFace)[Indices::velocity(scvf.directionIndex())];
        }
        else if (bcTypes.isSymmetry() || bcTypes.isDirichlet(Indices::pressureIdx))
            return parallelVelocity;
        else if (bcTypes.isBJS(Indices::velocity(scvf.directionIndex())))
        {
            const SubControlVolumeFace ghostFace = makeParallelGhostFace_(scvf, localIdx);
            return problem.bjsVelocity(boundaryElement, scvf, boundaryNormalFace, localIdx, parallelVelocity);
        }
        else
        {
            // Neumann conditions are not well implemented
            DUNE_THROW(Dune::InvalidStateException, "Something went wrong with the boundary conditions for the momentum equations at global position " << boundarySubFaceCenter);
        }
    }


    //! helper function to conveniently create a ghost face used to retrieve boundary values from the problem
    static SubControlVolumeFace makeGhostFace_(const SubControlVolumeFace& ownScvf, const GlobalPosition& pos)
    {
        return SubControlVolumeFace(pos, std::vector<unsigned int>{ownScvf.insideScvIdx(), ownScvf.outsideScvIdx()},
                                    ownScvf.directionIndex(), ownScvf.axisData().selfDof, ownScvf.index());
    };

    //! helper function to conveniently create a ghost face which is outside the domain, parallel to the scvf of interest
    static SubControlVolumeFace makeParallelGhostFace_(const SubControlVolumeFace& ownScvf, const int localSubFaceIdx)
    {
        return makeGhostFace_(ownScvf, ownScvf.pairData(localSubFaceIdx).virtualFirstParallelFaceDofPos);
    };
};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_UPWINDVARIABLES_HH