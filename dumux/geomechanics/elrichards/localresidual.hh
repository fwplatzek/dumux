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
 * \brief Element-wise calculation of the residual for the linear elastic,
 * two-phase model in the fully implicit scheme.
 */
#ifndef DUMUX_ELRICHARDS_LOCAL_RESIDUAL_HH
#define DUMUX_ELRICHARDS_LOCAL_RESIDUAL_HH

#include <dumux/implicit/box/localresidual.hh>
#include "properties.hh"


namespace Dumux {
/*!
 * \ingroup ElRichardsModel
 * \ingroup ImplicitLocalResidual
 *
 * \brief Element-wise calculation of the Jacobian matrix for problems
 * using the two-phase linear-elasticity fully implicit model.
 */
template<class TypeTag>
class ElRichardsLocalResidual: public BoxLocalResidual<TypeTag>
//class ElRichardsLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;


    enum {
        dim = GridView::dimension
    };

    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

     enum {
        //contiWEqIdx = Indices::contiWEqIdx,
        contiEqIdx = Indices::contiEqIdx,
        //contiNEqIdx = Indices::contiNEqIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        //nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, EffectivePermeabilityModel) EffectivePermeabilityModel;

    enum {
        numFluidPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

   //typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    static const bool useHead = GET_PROP_VALUE(TypeTag, UseHead);


public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    ElRichardsLocalResidual() {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = 0.5;
        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit,
                MassUpwindWeight);
    }
    ;

    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a finite sub-control volume.
     *
     *  \param storage The phase mass within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &storage, int scvIdx,
            bool usePrevSol) const {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit Euler method.
      //  const ElementVolumeVariables &elemVolVars =
         //       usePrevSol ? this->prevVolVars_() : this->curVolVars_();

/*          const VolumeVariables &volVars =
            usePrevSol ?
            this->prevVolVars_(scvIdx) :
            this->curVolVars_(scvIdx);


        // partial time derivative of the wetting phase mass
        //pressure head formulation
        storage[contiEqIdx] =
                volVars.saturation(wPhaseIdx)
                * volVars.porosity();
        // pressure formulation
        if (!useHead)
            storage[contiEqIdx] *= volVars.density(wPhaseIdx);
    }  */
        const ElementVolumeVariables &elemVolVars =
                usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        storage = Scalar(0);

        // wetting phase mass
        storage[contiEqIdx] = volVars.density(wPhaseIdx)
                * volVars.saturation(wPhaseIdx) * volVars.effPorosity;// effPorosity is effected by geometry and coupling. this part has been coupled to use just porosity insted of effporosity in line 269 elementvolvariables.
//         // non-wetting phase mass
//         storage[contiNEqIdx] = volVars.density(nPhaseIdx)
//                 * volVars.saturation(nPhaseIdx) * volVars.effPorosity;
    }

    /*!
     * \brief Evaluates the mass flux over a face of a sub-control
     *        volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face for each phase
     * \param fIdx The index of the SCV face
     * \param onBoundary A boolean variable to specify whether the flux variables
     *        are calculated for interior SCV faces or boundary faces, default=false
     */
    void computeFlux(PrimaryVariables &flux, int fIdx, const bool onBoundary = false) const
    {
        FluxVariables fluxVars;
        fluxVars.update(this->problem_(),
                        this->element_(),
                        this->fvGeometry_(),
                        fIdx,
                        this->curVolVars_(),
                        onBoundary);

        flux = 0;
//        this->computeAdvectiveFlux(flux, fluxVars);
        asImp_()->computeAdvectiveFlux(flux, fluxVars);
        asImp_()->computeDiffusiveFlux(flux, fluxVars);
    }



    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each phase
     * \param fluxVars The flux variables at the current SCV
     *
     * This method is called by compute flux and is mainly there for
     * derived models to ease adding equations selectively.
     */
    void computeAdvectiveFlux(PrimaryVariables &flux,
            const FluxVariables &fluxVars) const {

         // calculate effective permeability based on effective porosity
        // according to the relation given in Rutqvist and Tsang (2002)
        // this evaluation should be moved to another location
        DimVector tmpVec;

        DimMatrix Keff, Keff_i, Keff_j;
        Keff_i = EffectivePermeabilityModel::effectivePermeability(this->curVolVars_()[fluxVars.face().i],
                              this->problem_().spatialParams(),
                              this->element_(),
                              this->fvGeometry_(),
                              fluxVars.face().i);
        Keff_j = EffectivePermeabilityModel::effectivePermeability(this->curVolVars_()[fluxVars.face().j],
                              this->problem_().spatialParams(),
                              this->element_(),
                              this->fvGeometry_(),
                              fluxVars.face().j);

        this->problem_().spatialParams().meanK(Keff, Keff_i, Keff_j);
        // loop over all phases
        for (int phaseIdx = 0; phaseIdx < numFluidPhases; ++phaseIdx) {
            // data attached to upstream and the downstream vertices
            // of the current phase
            // calculate the flux in the normal direction of the
            // current sub control volume face

            // if geomechanical feedback on flow is taken into account the effective permeability is
            // applied for the flux calculations
            if (this->problem_().coupled() == true) {
                Keff.mv(fluxVars.potentialGrad(phaseIdx), tmpVec);
            } else {
                fluxVars.intrinsicPermeability().mv(
                        fluxVars.potentialGrad(phaseIdx), tmpVec);
            }
            Scalar normalFlux = -(tmpVec * fluxVars.face().normal);

            // data attached to upstream and the downstream vertices
            // of the current phase
            const VolumeVariables &up = this->curVolVars_(
                    fluxVars.upstreamIdx(phaseIdx));
            const VolumeVariables &dn = this->curVolVars_(
                    fluxVars.downstreamIdx(phaseIdx));

            // add advective flux of current phase
//             int eqIdx = (phaseIdx == wPhaseIdx) ? contiEqIdx : contiNEqIdx; //in my case, equation index(eqIdx) is always equal to contiEqIdx
//             flux[eqIdx] += normalFlux
//                     * ((massUpwindWeight_) * up.density(phaseIdx)
//                             * up.mobility(phaseIdx)
//                             + (massUpwindWeight_) * dn.density(phaseIdx)
//                                     * dn.mobility(phaseIdx));

            //from richards

        //pressure head formulation
        flux[contiEqIdx] =
            fluxVars.volumeFlux(wPhaseIdx);
        //pressure formulation
        if(!useHead)
            flux[contiEqIdx] *=((    massUpwindWeight_)*up.density(wPhaseIdx)
                                +
                                (1 - massUpwindWeight_)*dn.density(wPhaseIdx));

//             if (this->problem_().coupled() == true)
//             {
//                 if (std::abs(flux[contiEqIdx]) > 1.e-15)
//                 {
//                     std::cout << "massUpwindWeight_" << massUpwindWeight_ << std::endl;
//                     std::cout << "up.density(phaseIdx)" << up.density(phaseIdx) << std::endl;
//                     std::cout << "volumeFlux" << fluxVars.volumeFlux(wPhaseIdx) << std::endl;
//                 }
//             }

            // if geomechanical feedback on flow is taken into account add the flux contribution
            // of the displacement velocity

//             if (this->problem_().coupled() == true) { // here means when we couple it with geomechanics (after the initialization, instead of porosity, we should consider the effective porosity)
//                 // use upwind displacement velocity to calculate phase transport (?)
//                 flux[contiEqIdx] += up.effPorosity * up.saturation(phaseIdx)
//                         * up.density(phaseIdx) * fluxVars.timeDerivUNormal();
//             }

        }
    }


    /*!
     * \brief Adds the diffusive flux to the flux vector over
     *        the face of a sub-control volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each phase
     * \param fluxVars The flux variables at the current SCV
     *
     * This function doesn't do anything but may be used by the
     * non-isothermal three-phase models to calculate diffusive heat
     * fluxes
     */
    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        // diffusive fluxes
        flux += 0.0;
    }


    /*!
     * \brief Calculate the source term of the equation
     *
     * \param q The source/sink in the SCV for each phase
     * \param scvIdx The index of the SCV
     *
     */
    void computeSource(PrimaryVariables &q, int scvIdx) {
        // retrieve the source term intrinsic to the problem
        this->problem_().source(q, this->element_(), this->fvGeometry_(),
                scvIdx);
    }

protected:
    Implementation *asImp_() {
        return static_cast<Implementation *>(this);
    }
    const Implementation *asImp_() const {
        return static_cast<const Implementation *>(this);
    }

private:
    Scalar massUpwindWeight_;
};
}
#endif