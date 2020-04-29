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

#ifndef DUMUX_NAVIERSTOKES_FLUXHELPER_HH
#define DUMUX_NAVIERSTOKES_FLUXHELPER_HH

#include <dumux/common/math.hh>

namespace Dumux {
namespace NavierStokes {

template<class Traits, bool enableEneryBalance, bool isCompositional>
struct NavierStokesAdvectiveFluxHelper
{
    template<class NumEqVector, class VolumeVariables, class SubControlVolumeFace, class Scalar>
    static void advectiveFlux(NumEqVector& flux,
                              const VolumeVariables& insideVolVars,
                              const VolumeVariables& outsideVolVars,
                              const SubControlVolumeFace& scvf,
                              const Scalar velocity,
                              const Scalar upwindWeight)
    {
        advectiveEnergyFlux(flux, insideVolVars, outsideVolVars, velocity, scvf, upwindWeight);
        advectiveMassFlux(flux, insideVolVars, outsideVolVars, velocity, scvf, upwindWeight);
    }

    template <bool enable = enableEneryBalance, typename std::enable_if_t<!enable, int> = 0, class ...Args>
    static void advectiveEnergyFlux(Args&&... args) {}

    template <bool enable = enableEneryBalance, typename std::enable_if_t<enable, int> = 0,
              class NumEqVector, class VolumeVariables, class SubControlVolumeFace, class Scalar>
    static void advectiveEnergyFlux(NumEqVector& flux,
                                    const VolumeVariables& insideVolVars,
                                    const VolumeVariables& outsideVolVars,
                                    const SubControlVolumeFace& scvf,
                                    const Scalar velocity,
                                    const Scalar upwindWeight)
    {
        static constexpr auto localEnergyBalanceIdx = NumEqVector::dimension - 1;
        auto upwindTerm = [](const auto& volVars) { return volVars.density() * volVars.enthalpy(); };
        flux[localEnergyBalanceIdx] += advectiveUpwindFlux(insideVolVars,
                                                           outsideVolVars,
                                                           scvf,
                                                           velocity,
                                                           upwindWeight,
                                                           upwindTerm);
    }

    template <bool enable = isCompositional, typename std::enable_if_t<!enable, int> = 0,
              class NumEqVector, class VolumeVariables, class SubControlVolumeFace, class Scalar>
    static void advectiveMassFlux(NumEqVector& flux,
                                const VolumeVariables& insideVolVars,
                                const VolumeVariables& outsideVolVars,
                                const SubControlVolumeFace& scvf,
                                const Scalar velocity,
                                const Scalar upwindWeight)
    {
        auto upwindTerm = [](const auto& volVars) { return volVars.density(); };
        flux[Traits::Indices::conti0EqIdx - Traits::dim()] = advectiveUpwindFlux(insideVolVars,
                                                                                 outsideVolVars,
                                                                                 scvf,
                                                                                 velocity,
                                                                                 upwindWeight,
                                                                                 upwindTerm);
    }

    template <bool enable = isCompositional, typename std::enable_if_t<enable, int> = 0,
              class NumEqVector, class VolumeVariables, class SubControlVolumeFace, class Scalar>
    static void advectiveMassFlux(NumEqVector& flux,
                                  const VolumeVariables& insideVolVars,
                                  const VolumeVariables& outsideVolVars,
                                  const SubControlVolumeFace& scvf,
                                  const Scalar velocity,
                                  const Scalar upwindWeight)
    {
        static constexpr auto numComponents = Traits::numFluidComponents();
        static constexpr bool useMoles = Traits::useMoles();

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            auto upwindTerm = [compIdx](const auto& volVars)
            {
                const auto density = useMoles ? volVars.molarDensity() : volVars.density();
                const auto fraction =  useMoles ? volVars.moleFraction(compIdx) : volVars.massFraction(compIdx);
                return density * fraction;
            };

            flux[compIdx] = advectiveUpwindFlux(insideVolVars,
                                                outsideVolVars,
                                                scvf,
                                                velocity,
                                                upwindWeight,
                                                upwindTerm);
        }

        // in case one balance is substituted by the total mass balance
        if (Traits::replaceCompEqIdx() < numComponents)
        {
            flux[Traits::replaceCompEqIdx()] = std::accumulate(flux.begin(), flux.end(), 0.0);
        }
    }

    template<class VolumeVariables, class SubControlVolumeFace, class Scalar, class UpwindTerm>
    static Scalar advectiveUpwindFlux(const VolumeVariables& insideVolVars,
                                      const VolumeVariables& outsideVolVars,
                                      const SubControlVolumeFace& scvf,
                                      const Scalar velocity,
                                      const Scalar upwindWeight,
                                      UpwindTerm upwindTerm)
    {
        const bool insideIsUpstream = scvf.directionSign() == sign(velocity);

        const auto& upstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;
        const auto& downstreamVolVars = insideIsUpstream ? outsideVolVars : insideVolVars;

        const Scalar flux = (upwindWeight * upwindTerm(upstreamVolVars) +
                            (1.0 - upwindWeight) * upwindTerm(downstreamVolVars))
                            * velocity * scvf.directionSign();
        return flux;
    }

};

} // end namespace NavierStokes
} // end namespace Dumux

#endif
