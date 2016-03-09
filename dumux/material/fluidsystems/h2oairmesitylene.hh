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
 * \brief @copybrief Dumux::FluidSystems::H2OAirMesitylene
 */
#ifndef DUMUX_H2O_AIR_MESITYLENE_FLUID_SYSTEM_HH
#define DUMUX_H2O_AIR_MESITYLENE_FLUID_SYSTEM_HH

#include <dumux/material/idealgas.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/mesitylene.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

#include <dumux/material/binarycoefficients/h2o_air.hh>
#include <dumux/material/binarycoefficients/h2o_mesitylene.hh>
#include <dumux/material/binarycoefficients/air_mesitylene.hh>

#include <dumux/material/fluidsystems/base.hh>

namespace Dumux
{
namespace FluidSystems
{

/*!
 * \ingroup Fluidsystems
 * \brief A three-phase fluid system featuring gas, NAPL and water as phases and
 * distilled water \f$(\mathrm{H_2O})\f$ and air (Pseudo component composed of
 * \f$\mathrm{79\%\;N_2}\f$, \f$\mathrm{20\%\;O_2}\f$ and Mesitylene \f$(\mathrm{C_6H_3(CH_3)_3})\f$ as components. It assumes all phases to be ideal mixtures.
 */
template <class Scalar,
          class H2OType = Dumux::TabulatedComponent<Scalar, Dumux::H2O<Scalar> > >
class H2OAirMesitylene
    : public BaseFluidSystem<Scalar, H2OAirMesitylene<Scalar, H2OType> >
{
    typedef H2OAirMesitylene<Scalar, H2OType> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

public:
    typedef Dumux::Mesitylene<Scalar> NAPL;
    typedef Dumux::Air<Scalar> Air;
    typedef H2OType H2O;


    static const int numPhases = 3;
    static const int numComponents = 3;

    static const int wPhaseIdx = 0; // index of the water phase
    static const int nPhaseIdx = 1; // index of the NAPL phase
    static const int gPhaseIdx = 2; // index of the gas phase

    static const int H2OIdx = 0;
    static const int NAPLIdx = 1;
    static const int airIdx = 2;

    // export component indices to indicate the main component
    // of the corresponding phase at atmospheric pressure 1 bar
    // and room temperature 20°C:
    static const int wCompIdx = H2OIdx;
    static const int nCompIdx = NAPLIdx;
    static const int gCompIdx = airIdx;

    /*!
     * \brief Initialize the fluid system's static parameters generically
     *
     * If a tabulated H2O component is used, we do our best to create
     * tables that always work.
     */
    static void init()
    {
        init(/*tempMin=*/273.15,
             /*tempMax=*/623.15,
             /*numTemp=*/100,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/200);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water \f$\mathrm{[K]}\f$
     * \param tempMax The maximum temperature used for tabulation of water \f$\mathrm{[K]}\f$
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
     * \param pressMax The maximum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        if (H2O::isTabulated) {
            std::cout << "Initializing tables for the H2O fluid properties ("
                      << nTemp*nPress
                      << " entries).\n";

            H2O::init(tempMin, tempMax, nTemp,
                      pressMin, pressMax, nPress);
        }
    }


    /*!
     * \brief Return whether a phase is liquid
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gPhaseIdx;
    }

    static bool isIdealGas(int phaseIdx)
    { return phaseIdx == gPhaseIdx && H2O::gasIsIdeal() && Air::gasIsIdeal() && NAPL::gasIsIdeal(); }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are indepent on the fluid composition. This assumtion is true
     * if Henry's law and Rault's law apply. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Rault's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
        return true;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be compressible.
     *
     * Compressible means that the partial derivative of the density
     * to the fluid pressure is always larger than zero.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // gases are always compressible
        if (phaseIdx == gPhaseIdx)
            return true;
        else if (phaseIdx == wPhaseIdx)
            // the water component decides for the water phase...
            return H2O::liquidIsCompressible();

        // the NAPL component decides for the napl phase...
        return NAPL::liquidIsCompressible();
    }

    /*!
     * \brief Return the human readable name of a phase (used in indices)
     */
    static const char *phaseName(int phaseIdx)
    {
        switch (phaseIdx) {
        case wPhaseIdx: return "w";
        case nPhaseIdx: return "n";
        case gPhaseIdx: return "g";
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the human readable name of a component (used in indices)
     */
    static const char *componentName(int compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::name();
        case airIdx: return Air::name();
        case NAPLIdx: return NAPL::name();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     * \param compIdx The index of the component
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::molarMass();
        case airIdx: return Air::molarMass();
        case NAPLIdx: return NAPL::molarMass();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Given all mole fractions in a phase, return the phase
     *        density \f$\mathrm{[kg/m^3]}\f$.
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     */
    using Base::density;
    template <class FluidState>
    static Scalar density(const FluidState &fluidState, int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx) {
            // See: Ochs 2008
            // \todo: proper citation
            Scalar rholH2O = H2O::liquidDensity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
            Scalar clH2O = rholH2O/H2O::molarMass();

            // this assumes each dissolved molecule displaces exactly one
            // water molecule in the liquid
            return
                clH2O*(H2O::molarMass()*fluidState.moleFraction(wPhaseIdx, H2OIdx)
                       +
                       Air::molarMass()*fluidState.moleFraction(wPhaseIdx, airIdx)
                       +
                       NAPL::molarMass()*fluidState.moleFraction(wPhaseIdx, NAPLIdx));
        }
        else if (phaseIdx == nPhaseIdx) {
            // assume pure NAPL for the NAPL phase
            Scalar pressure = NAPL::liquidIsCompressible()?fluidState.pressure(phaseIdx):1e100;
            return NAPL::liquidDensity(fluidState.temperature(phaseIdx), pressure);
        }

        assert (phaseIdx == gPhaseIdx);
        Scalar pH2O =
            fluidState.moleFraction(gPhaseIdx, H2OIdx)  *
            fluidState.pressure(gPhaseIdx);
        Scalar pAir =
            fluidState.moleFraction(gPhaseIdx, airIdx)  *
            fluidState.pressure(gPhaseIdx);
        Scalar pNAPL =
            fluidState.moleFraction(gPhaseIdx, NAPLIdx)  *
            fluidState.pressure(gPhaseIdx);
        return
            H2O::gasDensity(fluidState.temperature(phaseIdx), pH2O) +
            Air::gasDensity(fluidState.temperature(phaseIdx), pAir) +
            NAPL::gasDensity(fluidState.temperature(phaseIdx), pNAPL);
    }

    /*!
     * \brief Return the viscosity of a phase \f$\mathrm{[Pa s]}\f$.
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     */
    using Base::viscosity;
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx) {
            // assume pure water viscosity
            return H2O::liquidViscosity(fluidState.temperature(phaseIdx),
                                        fluidState.pressure(phaseIdx));
        }
        else if (phaseIdx == nPhaseIdx) {
            // assume pure NAPL viscosity
            return NAPL::liquidViscosity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        }

        assert (phaseIdx == gPhaseIdx);

        /* Wilke method. See:
         *
         * See: R. Reid, et al.: The Properties of Gases and Liquids,
         * 4th edition, McGraw-Hill, 1987, 407-410
         * 5th edition, McGraw-Hill, 2001, p. 9.21/22
         *
         * in this case, we use a simplified version in order to avoid
         * computationally costly evaluation of sqrt and pow functions and
         * divisions
         * -- compare e.g. with Promo Class p. 32/33
         */
        Scalar muResult;
        const Scalar mu[numComponents] = {
            H2O::gasViscosity(fluidState.temperature(phaseIdx), H2O::vaporPressure(fluidState.temperature(phaseIdx))),
            Air::simpleGasViscosity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx)),
            NAPL::gasViscosity(fluidState.temperature(phaseIdx), NAPL::vaporPressure(fluidState.temperature(phaseIdx)))
        };
        // molar masses
        const Scalar M[numComponents] = {
            H2O::molarMass(),
            Air::molarMass(),
            NAPL::molarMass()
        };

        Scalar muAW = mu[airIdx]*fluidState.moleFraction(gPhaseIdx, airIdx)
            + mu[H2OIdx]*fluidState.moleFraction(gPhaseIdx, H2OIdx)
            / (fluidState.moleFraction(gPhaseIdx, airIdx)
               + fluidState.moleFraction(gPhaseIdx, H2OIdx));
        Scalar xAW = fluidState.moleFraction(gPhaseIdx, airIdx)
            + fluidState.moleFraction(gPhaseIdx, H2OIdx);

        Scalar MAW = (fluidState.moleFraction(gPhaseIdx, airIdx)*Air::molarMass()
                      + fluidState.moleFraction(gPhaseIdx, H2OIdx)*H2O::molarMass())
            / xAW;

        Scalar phiCAW = 0.3; // simplification for this particular system
        /* actually like this
         * Scalar phiCAW = std::pow(1.+std::sqrt(mu[NAPLIdx]/muAW)*std::pow(MAW/M[NAPLIdx],0.25),2)
         *                 / std::sqrt(8.*(1.+M[NAPLIdx]/MAW));
         */
        Scalar phiAWC = phiCAW * muAW*M[NAPLIdx]/(mu[NAPLIdx]*MAW);

        muResult = (xAW*muAW)/(xAW+fluidState.moleFraction(gPhaseIdx, NAPLIdx)*phiAWC)
            + (fluidState.moleFraction(gPhaseIdx, NAPLIdx) * mu[NAPLIdx])
            / (fluidState.moleFraction(gPhaseIdx, NAPLIdx) + xAW*phiCAW);
        return muResult;
    }


    /*!
     * \brief Given all mole fractions, return the diffusion
     *        coefficent in \f$\mathrm{[m^2/s]}\f$ of a component in a phase.
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     * \param compIdx The index of the component
     */
    using Base::diffusionCoefficient;
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        switch (phaseIdx)
        {
            case gPhaseIdx:
            {
                switch (compIdx)
                {
                    case NAPLIdx:
                    {
                        Scalar diffWC = Dumux::BinaryCoeff::H2O_Mesitylene::gasDiffCoeff(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
                        Scalar diffAW = Dumux::BinaryCoeff::H2O_Air::gasDiffCoeff(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
                        const Scalar xga = fluidState.moleFraction(gPhaseIdx, airIdx);
                        const Scalar xgw = fluidState.moleFraction(gPhaseIdx, H2OIdx);
                        const Scalar xgc = fluidState.moleFraction(gPhaseIdx, NAPLIdx);
                        return (1.- xgw)/(xga/diffAW + xgc/diffWC);
                    }
                    case H2OIdx:
                    {
                        Scalar diffAC = Dumux::BinaryCoeff::Air_Mesitylene::gasDiffCoeff(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
                        Scalar diffWC = Dumux::BinaryCoeff::H2O_Mesitylene::gasDiffCoeff(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
                        const Scalar xga = fluidState.moleFraction(gPhaseIdx, airIdx);
                        const Scalar xgw = fluidState.moleFraction(gPhaseIdx, H2OIdx);
                        const Scalar xgc = fluidState.moleFraction(gPhaseIdx, NAPLIdx);
                        return (1.- xgc)/(xgw/diffWC + xga/diffAC);
                    }
                    case airIdx:
                        DUNE_THROW(Dune::InvalidStateException, "Diffusivity of air in the gas phase is constraint by sum of diffusive fluxes = 0 !");
                }
            }
            case wPhaseIdx:
            {
                Scalar diffACl = 1.e-9; // BinaryCoeff::Air_Mesitylene::liquidDiffCoeff(temperature, pressure);
                Scalar diffWCl = 1.e-9; // BinaryCoeff::H2O_Mesitylene::liquidDiffCoeff(temperature, pressure);
                Scalar diffAWl = 1.e-9; // BinaryCoeff::H2O_Air::liquidDiffCoeff(temperature, pressure);

                Scalar xwa = fluidState.moleFraction(wPhaseIdx, airIdx);
                Scalar xww = fluidState.moleFraction(wPhaseIdx, H2OIdx);
                Scalar xwc = fluidState.moleFraction(wPhaseIdx, NAPLIdx);

                switch (compIdx)
                {
                    case NAPLIdx:
                        return (1.- xww)/(xwa/diffAWl + xwc/diffWCl);
                    case airIdx:
                        return (1.- xwc)/(xww/diffWCl + xwa/diffACl);
                    case H2OIdx:
                        DUNE_THROW(Dune::InvalidStateException,
                                   "Diffusivity of water in the water phase "
                                   "is constraint by sum of diffusive fluxes = 0 !\n");
                }
            }
            case nPhaseIdx:
                DUNE_THROW(Dune::InvalidStateException, "Diffusion coefficients of substances in non-wetting liquid phase are undefined!\n");
        }
        return 0;
    }

    using Base::binaryDiffusionCoefficient;
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::H2OAirMesitylene::binaryDiffusionCoefficient()");
    }

    /*!
     * \brief Returns the fugacity coefficient \f$\mathrm{[-]}\f$ of a component in a
     *        phase.
     *
     * In this case, things are actually pretty simple. We have an ideal
     * solution. Thus, the fugacity coefficient is 1 in the gas phase
     * (fugacity equals the partial pressure of the component in the gas phase
     * respectively in the liquid phases it is the inverse of the
     * Henry coefficients scaled by pressure
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     * \param compIdx The index of the component
     */
    using Base::fugacityCoefficient;
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == wPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::vaporPressure(T)/p;
            else if (compIdx == airIdx)
                return Dumux::BinaryCoeff::H2O_Air::henry(T)/p;
            else if (compIdx == NAPLIdx)
                return Dumux::BinaryCoeff::H2O_Mesitylene::henry(T)/p;
        }

        // for the NAPL phase, we assume currently that nothing is
        // dissolved. this means that the affinity of the NAPL
        // component to the NAPL phase is much higher than for the
        // other components, i.e. the fugacity cofficient is much
        // smaller.
        if (phaseIdx == nPhaseIdx) {
            Scalar phiNapl = NAPL::vaporPressure(T)/p;
            if (compIdx == NAPLIdx)
                return phiNapl;
            else if (compIdx == airIdx)
                return 1e6*phiNapl;
            else if (compIdx == H2OIdx)
                return 1e6*phiNapl;
        }

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        assert(phaseIdx == gPhaseIdx);
        return 1.0;
    }


    /*!
     * \brief Given all mole fractions in a phase, return the specific
     *        phase enthalpy \f$\mathrm{[J/kg]}\f$.
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     */
    /*!
     *  \todo This system neglects the contribution of gas-molecules in the liquid phase.
     *        This contribution is probably not big. Somebody would have to find out the enthalpy of solution for this system. ...
     */
    using Base::enthalpy;
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx) {
            return H2O::liquidEnthalpy(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        }
        else if (phaseIdx == nPhaseIdx) {
            return NAPL::liquidEnthalpy(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        }
        else if (phaseIdx == gPhaseIdx) {  // gas phase enthalpy depends strongly on composition
            Scalar hgc = NAPL::gasEnthalpy(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));
            Scalar hgw = H2O::gasEnthalpy(fluidState.temperature(phaseIdx),
                                          fluidState.pressure(phaseIdx));
            // pressure is only a dummy here (not dependent on pressure, just temperature)
            Scalar hga = Air::gasEnthalpy(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));

            Scalar result = 0;
            result += hgw * fluidState.massFraction(gPhaseIdx, H2OIdx);
            result += hga * fluidState.massFraction(gPhaseIdx, airIdx);
            result += hgc * fluidState.massFraction(gPhaseIdx, NAPLIdx);

            return result;
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }
 /*!
     * \brief Return the heat capacity in \f$\mathrm{[J/(kg K)]}\f$.
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     */
    using Base::heatCapacity;
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::H2OAirMesitylene::heatCapacity()");
    }

 /*!
     * \brief Return the thermal conductivity \f$\mathrm{[W/(m K)]}\f$.
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     */
    using Base::thermalConductivity;
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx){
            const Scalar temperature  = fluidState.temperature(phaseIdx) ;
            const Scalar pressure = fluidState.pressure(phaseIdx);
            return H2O::liquidThermalConductivity(temperature, pressure);
        }
        else if (phaseIdx == nPhaseIdx)
        {
            //Thermal Conductivity of p-Xylene
            //taken from the Dortmund Data Bank
            //see http://www.ddbst.de/en/EED/PCP/TCN_C176.php
            const Scalar lambdaPureXylene = 0.13;
            return lambdaPureXylene; // conductivity of p-xylene [W/(m K)]
        }
        else if (phaseIdx == gPhaseIdx)
        {
            //gPhaseIdx i.e. air
            // Isobaric Properties for Nitrogen in: NIST Standard
            // see http://webbook.nist.gov/chemistry/fluid/
            // evaluated at p=.1 MPa, T=20°C
            // Nitrogen: 0.025398
            // Oxygen: 0.026105
            // lambda_air is approximately 0.78*lambda_N2+0.22*lambda_O2
            const Scalar lambdaPureAir = 0.0255535;

            return lambdaPureAir; // conductivity of pure air [W/(m K)]
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

private:
    static Scalar waterPhaseDensity_(Scalar T, Scalar pw, Scalar xww, Scalar xwa, Scalar xwc)
    {
        Scalar rholH2O = H2O::liquidDensity(T, pw);
        Scalar clH2O = rholH2O/H2O::molarMass();

        // this assumes each dissolved molecule displaces exactly one
        // water molecule in the liquid
        return
            clH2O*(xww*H2O::molarMass() + xwa*Air::molarMass() + xwc*NAPL::molarMass());
    }

    static Scalar gasPhaseDensity_(Scalar T, Scalar pg, Scalar xgw, Scalar xga, Scalar xgc)
    {
        return H2O::gasDensity(T, pg*xgw) + Air::gasDensity(T, pg*xga) + NAPL::gasDensity(T, pg*xgc);
    }

    static Scalar NAPLPhaseDensity_(Scalar T, Scalar pn)
    {
        return
            NAPL::liquidDensity(T, pn);
    }

};
} // end namespace FluidSystems

#ifdef DUMUX_PROPERTIES_HH
// forward defintions of the property tags
namespace Properties {
    NEW_PROP_TAG(Scalar);
    NEW_PROP_TAG(Components);
}

/*!
 * \brief A threephase fluid system with water, air and mesitylene as components.
 *
 * This is an adapter to use Dumux::H2OAirMesityleneFluidSystem<TypeTag>, as is
 * done with most other classes in Dumux.
 *  This fluidsystem is applied by default with the tabulated version of
 *  water of the IAPWS-formulation.
 *
 *  To change the component formulation (ie to use nontabulated or
 *  incompressible water), or to switch on verbosity of tabulation,
 *  use the property system and the property "Components":
 *
 *        // Select desired version of the component
 *        SET_PROP(myApplicationProperty, Components) : public GET_PROP(TypeTag, DefaultComponents)
 *        {
 *            typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 *
 *        // Do not use the defaults !
 *        //    typedef Dumux::TabulatedComponent<Scalar, Dumux::H2O<Scalar> > H2O;
 *
 *        // Apply e.g. untabulated water:
 *        typedef Dumux::H2O<Scalar> H2O;
 *        };
 *
 *   Also remember to initialize tabulated components (FluidSystem::init()), while this
 *   is not necessary for non-tabularized ones.
 */
template<class TypeTag>
class H2OAirMesityleneFluidSystem
: public FluidSystems::H2OAirMesitylene<typename GET_PROP_TYPE(TypeTag, Scalar),
                                        typename GET_PROP(TypeTag, Components)::H2O>
{};
#endif

} // end namespace Dumux

#endif
