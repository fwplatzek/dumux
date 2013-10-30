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
 * \brief Contains the quantities which are constant within a
 *        finite volume in the non-isothermal three-phase three-component
 *        model.
 */
#ifndef DUMUX_3P3CNI_VOLUME_VARIABLES_HH
#define DUMUX_3P3CNI_VOLUME_VARIABLES_HH

#include <dumux/implicit/3p3c/3p3cvolumevariables.hh>


namespace Dumux
{

/*!
 * \ingroup ThreePThreeCNIModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the non-isothermal three-phase three-component
 *        model.
 */
template <class TypeTag>
class ThreePThreeCNIVolumeVariables : public ThreePThreeCVolumeVariables<TypeTag>
{
    //! \cond false
    typedef ThreePThreeCVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { temperatureIdx = Indices::temperatureIdx };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    //! \endcond

public:
    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int scvIdx,
                bool isOldSol)
    {
        // vertex update data for the mass balance
        ParentType::update(priVars,
                           problem,
                           element,
                           fvGeometry,
                           scvIdx,
                           isOldSol);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(this->fluidState());

        // the internal energies and the enthalpies
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar h = FluidSystem::enthalpy(this->fluidState_, paramCache, phaseIdx);
            this->fluidState_.setEnthalpy(phaseIdx, h);
        }
    };

    /*!
     * \brief Returns the total internal energy of a phase in the
     *        sub-control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar internalEnergy(int phaseIdx) const
    { return this->fluidState_.internalEnergy(phaseIdx); };

    /*!
     * \brief Returns the total enthalpy of a phase in the sub-control
     *        volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar enthalpy(int phaseIdx) const
    { return this->fluidState_.enthalpy(phaseIdx); };

    /*!
     * \brief Returns the total heat capacity \f$\mathrm{[J/(K*m^3]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar heatCapacity() const
    { return heatCapacity_; };

protected:
    // this method gets called by the parent class. since this method
    // is protected, we are friends with our parent..
    friend class ThreePThreeCVolumeVariables<TypeTag>;

    static Scalar temperature_(const PrimaryVariables &priVars,
                               const Problem& problem,
                               const Element &element,
                               const FVElementGeometry &fvGeometry,
                               const int scvIdx)
    {
        return priVars[temperatureIdx];
    }

    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param priVars The solution primary variables
     * \param problem The problem
     * \param element The element
     * \param fvGeometry Evaluate function with solution of current or previous time step
     * \param scvIdx The local index of the SCV (sub-control volume)
     * \param isOldSol Evaluate function with solution of current or previous time step
     */
    void updateEnergy_(const PrimaryVariables &priVars,
                       const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &fvGeometry,
                       const int scvIdx,
                       bool isOldSol)
    {
        // copmute and set the heat capacity of the solid phase
        heatCapacity_ = problem.spatialParams().heatCapacity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(heatCapacity_);
    };

    Scalar heatCapacity_;
};

} // end namespace

#endif
