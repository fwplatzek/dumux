// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Holger Class                                 *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief   This file contains the data which is required to calculate
 *          all fluxes of components over a face of a finite volume for
 *          the two-phase, two-component model.
 */
/*!
 * \ingroup ThreePThreeCModel
 */
#ifndef DUMUX_3P3C_FLUX_VARIABLES_HH
#define DUMUX_3P3C_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/spline.hh>

#include "3p3cproperties.hh"

namespace Dumux
{

/*!
 * \brief This template class contains the data which is required to
 *        calculate all fluxes of components over a face of a finite
 *        volume for the two-phase, two-component model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */
template <class TypeTag>
class ThreePThreeCFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef Dune::FieldVector<Scalar, dim> 	DimVector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        comp0Idx = Indices::comp0Idx,
        comp1Idx = Indices::comp1Idx,
        comp2Idx = Indices::comp2Idx
    };

public:
    /*
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param faceIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    ThreePThreeCFluxVariables(const Problem &problem,
                              const Element &element,
                              const FVElementGeometry &fvGeometry,
                              const int faceIdx,
                              const ElementVolumeVariables &elemVolVars,
                              const bool onBoundary = false)
        : fvGeometry_(fvGeometry),scvfIdx_(faceIdx), onBoundary_(onBoundary)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            density_[phaseIdx] = Scalar(0);
            molarDensity_[phaseIdx] = Scalar(0);
            potentialGrad_[phaseIdx] = Scalar(0);
            massFractionComp0Grad_[phaseIdx] = Scalar(0);
            massFractionComp1Grad_[phaseIdx] = Scalar(0);
            massFractionComp2Grad_[phaseIdx] = Scalar(0);
            moleFractionComp0Grad_[phaseIdx] = Scalar(0);
            moleFractionComp1Grad_[phaseIdx] = Scalar(0);
            moleFractionComp2Grad_[phaseIdx] = Scalar(0);
        }

        calculateGradients_(problem, element, elemVolVars);
        calculateVelocities_(problem, element, elemVolVars);
        calculateporousDiffCoeff_(problem, element, elemVolVars);
    };

private:
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        // calculate gradients
    	DimVector tmp(0.0);
        for (int idx = 0;
             idx < fvGeometry_.numFAP;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const DimVector &feGrad = face().grad[idx];

            // index for the element volume variables 
            int volVarsIdx = face().fapIndices[idx];

            // compute sum of pressure gradients for each phase
            for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
            {
                // the pressure gradient
                tmp = feGrad;
                tmp *= elemVolVars[volVarsIdx].pressure(phaseIdx);
                potentialGrad_[phaseIdx] += tmp;
            }

            // the concentration gradient of the components
            // component in the phases
            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().massFraction(wPhaseIdx, comp0Idx);
            massFractionComp0Grad_[wPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().massFraction(nPhaseIdx, comp0Idx);
            massFractionComp0Grad_[nPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().massFraction(gPhaseIdx, comp0Idx);
            massFractionComp0Grad_[gPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().massFraction(wPhaseIdx, comp1Idx);
            massFractionComp1Grad_[wPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().massFraction(nPhaseIdx, comp1Idx);
            massFractionComp1Grad_[nPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().massFraction(gPhaseIdx, comp1Idx);
            massFractionComp1Grad_[gPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().massFraction(wPhaseIdx, comp2Idx);
            massFractionComp2Grad_[wPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().massFraction(nPhaseIdx, comp2Idx);
            massFractionComp2Grad_[nPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().massFraction(gPhaseIdx, comp2Idx);
            massFractionComp2Grad_[gPhaseIdx] += tmp;

            // the molar concentration gradients of the components
            // in the phases
            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().moleFraction(wPhaseIdx, comp0Idx);
            moleFractionComp0Grad_[wPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().moleFraction(nPhaseIdx, comp0Idx);
            moleFractionComp0Grad_[nPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().moleFraction(gPhaseIdx, comp0Idx);
            moleFractionComp0Grad_[gPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().moleFraction(wPhaseIdx, comp1Idx);
            moleFractionComp1Grad_[wPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().moleFraction(nPhaseIdx, comp1Idx);
            moleFractionComp1Grad_[nPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().moleFraction(gPhaseIdx, comp1Idx);
            moleFractionComp1Grad_[gPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().moleFraction(wPhaseIdx, comp2Idx);
            moleFractionComp2Grad_[wPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().moleFraction(nPhaseIdx, comp2Idx);
            moleFractionComp2Grad_[nPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().moleFraction(gPhaseIdx, comp2Idx);
            moleFractionComp2Grad_[gPhaseIdx] += tmp;
        }

        // correct the pressure gradients by the hydrostatic
        // pressure due to gravity
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            int i = face().i;
            int j = face().j;
            Scalar fI = rhoFactor_(phaseIdx, i, elemVolVars);
            Scalar fJ = rhoFactor_(phaseIdx, j, elemVolVars);
            if (fI + fJ <= 0)
                fI = fJ = 0.5; // doesn't matter because no phase is
            // present in both cells!
            density_[phaseIdx] =
                (fI*elemVolVars[i].density(phaseIdx) +
                 fJ*elemVolVars[j].density(phaseIdx))
                /
                (fI + fJ);
            // phase density
            molarDensity_[phaseIdx]
                =
                (fI*elemVolVars[i].molarDensity(phaseIdx) +
                 fJ*elemVolVars[j].molarDensity(phaseIdx))
                /
                (fI + fJ);

            tmp = problem.gravity();
            tmp *= density_[phaseIdx];

            potentialGrad_[phaseIdx] -= tmp;
        }
    }

    Scalar rhoFactor_(int phaseIdx, int scvIdx, const ElementVolumeVariables &elemVolVars)
    {
        static const Scalar eps = 1e-2;
        const Scalar sat = elemVolVars[scvIdx].density(phaseIdx);
        if (sat > eps)
            return 0.5;
        if (sat <= 0)
            return 0;

        static const Dumux::Spline<Scalar> sp(0, eps, // x0, x1
                                              0, 0.5, // y0, y1
                                              0, 0); // m0, m1
        return sp.eval(sat);
    }

    void calculateVelocities_(const Problem &problem,
                              const Element &element,
                              const ElementVolumeVariables &elemVolVars)
    {
        const SpatialParameters &spatialParams = problem.spatialParams();
        // multiply the pressure potential with the intrinsic
        // permeability
        DimMatrix K;
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            spatialParams.meanK(K,
                                spatialParams.intrinsicPermeability(element,
                                                                    fvGeometry_,
                                                                    face().i),
                                spatialParams.intrinsicPermeability(element,
                                                                    fvGeometry_,
                                                                    face().j));
            K.mv(potentialGrad_[phaseIdx], Kmvp_[phaseIdx]);
            KmvpNormal_[phaseIdx] = - (Kmvp_[phaseIdx] * face().normal);
        }

        // set the upstream and downstream vertices
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            upstreamIdx_[phaseIdx] = face().i;
            downstreamIdx_[phaseIdx] = face().j;

            if (KmvpNormal_[phaseIdx] < 0) {
                std::swap(upstreamIdx_[phaseIdx],
                          downstreamIdx_[phaseIdx]);
            }
        }
    }

    void calculateporousDiffCoeff_(const Problem &problem,
                               const Element &element,
                               const ElementVolumeVariables &elemVolVars)
    {

        const VolumeVariables &volVarsI = elemVolVars[face().i];
        const VolumeVariables &volVarsJ = elemVolVars[face().j];

        Dune::FieldMatrix<Scalar, numPhases, numComponents> diffusionCoefficientMatrix_i = volVarsI.diffusionCoefficient();
        Dune::FieldMatrix<Scalar, numPhases, numComponents> diffusionCoefficientMatrix_j = volVarsJ.diffusionCoefficient();

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // make sure to calculate only diffusion coefficents
            // for phases which exist in both finite volumes
            /* \todo take care: This should be discussed once again
             * as long as a meaningful value can be found for the required mole fraction
             * diffusion should work even without this one here */
            if (volVarsI.saturation(phaseIdx) <= 0 ||
                volVarsJ.saturation(phaseIdx) <= 0)
            {
                porousDiffCoeff_[phaseIdx][comp0Idx] = 0.0;
                porousDiffCoeff_[phaseIdx][comp1Idx] = 0.0;
                porousDiffCoeff_[phaseIdx][comp2Idx] = 0.0;
                continue;
            }

            // calculate tortuosity at the nodes i and j needed
            // for porous media diffusion coefficient

            Scalar tauI =
                1.0/(volVarsI.porosity() * volVarsI.porosity()) *
                pow(volVarsI.porosity() * volVarsI.saturation(phaseIdx), 7.0/3);
            Scalar tauJ =
                1.0/(volVarsJ.porosity() * volVarsJ.porosity()) *
                pow(volVarsJ.porosity() * volVarsJ.saturation(phaseIdx), 7.0/3);
            // Diffusion coefficient in the porous medium

            // -> harmonic mean
            porousDiffCoeff_[phaseIdx][comp0Idx] = harmonicMean(volVarsI.porosity() * volVarsI.saturation(phaseIdx) * tauI * diffusionCoefficientMatrix_i[phaseIdx][comp0Idx],
                                                                volVarsJ.porosity() * volVarsJ.saturation(phaseIdx) * tauJ * diffusionCoefficientMatrix_j[phaseIdx][comp0Idx]);
            porousDiffCoeff_[phaseIdx][comp1Idx] = harmonicMean(volVarsI.porosity() * volVarsI.saturation(phaseIdx) * tauI * diffusionCoefficientMatrix_i[phaseIdx][comp1Idx],
                                                                volVarsJ.porosity() * volVarsJ.saturation(phaseIdx) * tauJ * diffusionCoefficientMatrix_j[phaseIdx][comp1Idx]);
            porousDiffCoeff_[phaseIdx][comp2Idx] = harmonicMean(volVarsI.porosity() * volVarsI.saturation(phaseIdx) * tauI * diffusionCoefficientMatrix_i[phaseIdx][comp2Idx],
                                                                volVarsJ.porosity() * volVarsJ.saturation(phaseIdx) * tauJ * diffusionCoefficientMatrix_j[phaseIdx][comp2Idx]);

        }
    }

public:
    /*!
     * \brief Return the pressure potential multiplied with the
     *        intrinsic permeability which goes from vertex i to
     *        vertex j.
     *
     * Note that the length of the face's normal is the area of the
     * phase, so this is not the actual velocity by the integral of
     * the velocity over the face's area. Also note that the phase
     * mobility is not yet included here since this would require a
     * decision on the upwinding approach (which is done in the
     * actual model).
     */
    Scalar KmvpNormal(int phaseIdx) const
    { return KmvpNormal_[phaseIdx]; }

    /*!
     * \brief Return the pressure potential multiplied with the
     *        intrinsic permeability as DimVector (for velocity output)
     */
    DimVector Kmvp(int phaseIdx) const
    { return Kmvp_[phaseIdx]; }

    /*!
     * \brief Return the local index of the upstream control volume
     *        for a given phase.
     */
    int upstreamIdx(int phaseIdx) const
    { return upstreamIdx_[phaseIdx]; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase.
     */
    int downstreamIdx(int phaseIdx) const
    { return downstreamIdx_[phaseIdx]; }

    /*!
     * \brief The diffusivity matrix
     */
    Dune::FieldMatrix<Scalar, numPhases, numComponents> porousDiffCoeff() const
    { return porousDiffCoeff_; };

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ of a phase at the integration
     *        point.
     */
    DUMUX_DEPRECATED_MSG("use density instead")
    Scalar densityAtIP(int phaseIdx) const
    {
      density(phaseIdx);
    }

    /*!
	 * \brief Return density \f$\mathrm{[kg/m^3]}\f$ of a phase.
	 */
	Scalar density(int phaseIdx) const
	{ return density_[phaseIdx]; }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ of a phase at the integration
     *        point.
     */
	DUMUX_DEPRECATED_MSG("use molarDensity instead")
    Scalar molarDensityAtIP(int phaseIdx) const
    {
      molarDensity(phaseIdx);
    }
    /*!
	 * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ of a phase.
	 */
	Scalar molarDensity(int phaseIdx) const
	{ return molarDensity_[phaseIdx]; }

    /*!
     * \brief The mass fraction gradients of the components in a phase.
     */
	DUMUX_DEPRECATED_MSG("use massFractionComp0Grad instead")
    const DimVector &wConcentrationGrad(int phaseIdx) const
    { massFractionComp0Grad(phaseIdx); };

	const DimVector &massFractionComp0Grad(int phaseIdx) const
	{return massFractionComp0Grad_[phaseIdx];}

	DUMUX_DEPRECATED_MSG("use massFractionComp1Grad instead")
    const DimVector &cConcentrationGrad(int phaseIdx) const
    { massFractionComp1Grad(phaseIdx); };

	const DimVector &massFractionComp1Grad(int phaseIdx) const
	{ return massFractionComp1Grad_[phaseIdx]; };

	DUMUX_DEPRECATED_MSG("use massFractionComp2Grad instead")
    const DimVector &aConcentrationGrad(int phaseIdx) const
    {  massFractionComp2Grad(phaseIdx); };

	const DimVector &massFractionComp2Grad(int phaseIdx) const
	{ return massFractionComp2Grad_[phaseIdx]; };


    /*!
     * \brief The molar concentration gradients of the components in a phase.
     */
	DUMUX_DEPRECATED_MSG("use moleFractionComp0Grad instead")
    const DimVector &molarWConcGrad(int phaseIdx) const
    {  moleFractionComp0Grad(phaseIdx); };

	const DimVector &moleFractionComp0Grad(int phaseIdx) const
	{ return moleFractionComp0Grad_[phaseIdx]; };

	DUMUX_DEPRECATED_MSG("use moleFractionComp1Grad instead")
    const DimVector &molarCConcGrad(int phaseIdx) const
    { moleFractionComp1Grad(phaseIdx); };

	const DimVector &moleFractionComp1Grad(int phaseIdx) const
	{ return moleFractionComp1Grad_[phaseIdx]; };

	DUMUX_DEPRECATED_MSG("use moleFractionComp2Grad instead")
    const DimVector &molarAConcGrad(int phaseIdx) const
    { moleFractionComp2Grad(phaseIdx); };

	const DimVector &moleFractionComp2Grad(int phaseIdx) const
    { return moleFractionComp2Grad_[phaseIdx]; };

    const SCVFace &face() const
    {
        if (onBoundary_)
            return fvGeometry_.boundaryFace[scvfIdx_];
        else
            return fvGeometry_.subContVolFace[scvfIdx_];
    }

protected:
    const FVElementGeometry &fvGeometry_;
    int scvfIdx_;
    const bool onBoundary_;

    // gradients
    DimVector potentialGrad_[numPhases];
    DimVector massFractionComp0Grad_[numPhases];
    DimVector massFractionComp1Grad_[numPhases];
    DimVector massFractionComp2Grad_[numPhases];
    DimVector moleFractionComp0Grad_[numPhases];
    DimVector moleFractionComp1Grad_[numPhases];
    DimVector moleFractionComp2Grad_[numPhases];

    // density of each face at the integration point
    Scalar density_[numPhases], molarDensity_[numPhases];

    // intrinsic permeability times pressure potential gradient
    DimVector Kmvp_[numPhases];
    // projected on the face normal
    Scalar KmvpNormal_[numPhases];

    // local index of the upwind vertex for each phase
    int upstreamIdx_[numPhases];
    // local index of the downwind vertex for each phase
    int downstreamIdx_[numPhases];

    // the diffusivity matrix for the porous medium
    Dune::FieldMatrix<Scalar, numPhases, numComponents> porousDiffCoeff_;
};

} // end namepace

#endif
