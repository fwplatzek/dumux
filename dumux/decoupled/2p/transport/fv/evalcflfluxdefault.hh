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
#ifndef DUMUX_EVALCFLFLUX_DEFAULT_HH
#define DUMUX_EVALCFLFLUX_DEFAULT_HH

/**
 * @file
 * @brief  Fluxes to evaluate a CFL-Condition
 */

#include <dumux/decoupled/common/impetproperties.hh> 
#include "evalcflflux.hh"

namespace Dumux
{
/*!\ingroup IMPES
 * @brief  Default implementation of cfl-fluxes to evaluate a CFL-Condition
 *
 * Compares the maximum of inflow and outflow to the element volume weighted by relative permeability and viscosity ratios.
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class EvalCflFluxDefault: public EvalCflFlux<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

      typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    enum
    {
        vt = Indices::velocityTotal
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

public:

    /*! \brief adds a flux to the cfl-criterion evaluation
     *
     * \copydetails EvalCflFlux::addFlux(Scalar&,Scalar&,Scalar&,Scalar&,Scalar,const Element&,int)
     */
    void addFlux(Scalar& lambdaW, Scalar& lambdaNw, Scalar& viscosityW, Scalar& viscosityNw, Scalar flux, const Element& element, int phaseIdx = -1)
    {
        addFlux(lambdaW, lambdaNw, viscosityW, viscosityNw, flux, phaseIdx);
    }

    /*! \brief adds a flux to the cfl-criterion evaluation
     *
     * \copydetails EvalCflFlux::addFlux(Scalar&,Scalar&,Scalar&,Scalar&,Scalar,const Intersection&,int)
     */
    void addFlux(Scalar& lambdaW, Scalar& lambdaNw, Scalar& viscosityW, Scalar& viscosityNw, Scalar flux, const Intersection& intersection, int phaseIdx = -1)
    {
        addFlux(lambdaW, lambdaNw, viscosityW, viscosityNw, flux, phaseIdx);
    }

    /*! \brief Returns the CFL flux-function
     *
     * \copydetails EvalCflFlux::getCflFluxFunction(const Element&)
     */
    Scalar getCflFluxFunction(const Element& element);

    DUNE_DEPRECATED_MSG("use getCflFluxFunction() (uncapitalized 'fl') instead")
    Scalar getCFLFluxFunction(const Element& element)
    { return getCflFluxFunction(element); }

    /*! \brief  Returns the CFL time-step
     *
     * \copydetails EvalCflFlux::getDt(const Element&)
     */
    Scalar getDt(const Element& element)
    {
        Scalar porosity = std::max(problem_.spatialParams().porosity(element), porosityThreshold_);

        return (getCflFluxFunction(element) * porosity * element.geometry().volume());
    }

    //! resets the accumulated CFL-fluxes to zero
    void reset()
    {
        fluxWettingOut_ = 0;
        fluxNonwettingOut_ = 0;
        fluxIn_ = 0;
        fluxOut_ = 0;
    }

    /*! \brief Constructs an EvalCflFluxDefault object
     *
     * \param problem A problem type object
     */
    EvalCflFluxDefault (Problem& problem)
    : problem_(problem)
    {
        porosityThreshold_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Impet, PorosityThreshold);
        reset();
    }

private:

    void addFlux(Scalar& lambdaW, Scalar& lambdaNw, Scalar& viscosityW, Scalar& viscosityNw, Scalar flux, int phaseIdx = -1)
    {
        Scalar krSum = lambdaW * viscosityW + lambdaNw * viscosityNw;
        Scalar viscosityRatio = 1 - fabs(0.5 - viscosityNw / (viscosityW + viscosityNw));//1 - fabs(viscosityWI-viscosityNwI)/(viscosityWI+viscosityNwI);

        switch (phaseIdx)
         {
         case wPhaseIdx:
         {
             //for time step criterion
             if (flux >= 0)
             {
                 fluxWettingOut_ += flux / (krSum * viscosityRatio);
             }
             if (flux < 0)
             {
                 fluxIn_ -= flux / (krSum * viscosityRatio);
             }

             break;
         }

             //for time step criterion if the non-wetting phase velocity is used
         case nPhaseIdx:
         {
             if (flux >= 0)
             {
                 fluxNonwettingOut_ += flux / (krSum * viscosityRatio);
             }
             if (flux < 0)
             {
                 fluxIn_ -= flux / (krSum * viscosityRatio);
             }

             break;
         }
         default:
         {
             if (flux >= 0)
             {
                 fluxOut_ += flux / (krSum * viscosityRatio);
             }
             if (flux < 0)
             {
                 fluxIn_ -= flux / (krSum * viscosityRatio);
             }

             break;
         }
         }
    }

    Scalar getCFLFluxIn(int phaseIdx = 0)
    {
        if (std::isnan(fluxIn_) || std::isinf(fluxIn_))
        {
            fluxIn_ = 1e-100;
        }

            return fluxIn_;
    }

    Scalar getCFLFluxOut(int phaseIdx = 0)
    {
        if (std::isnan(fluxWettingOut_) || std::isinf(fluxWettingOut_))
        {
            fluxWettingOut_ = 1e-100;
        }

        if (std::isnan(fluxNonwettingOut_) || std::isinf(fluxNonwettingOut_))
        {
            fluxNonwettingOut_ = 1e-100;
        }
        if (std::isnan(fluxOut_) || std::isinf(fluxOut_))
        {
            fluxOut_ = 1e-100;
        }

        if (phaseIdx == wPhaseIdx)
            return fluxWettingOut_;
        else if (phaseIdx == nPhaseIdx)
            return fluxNonwettingOut_;
        else
            return fluxOut_;
    }

    Problem& problem_;//problem data
    Scalar fluxWettingOut_;
    Scalar fluxNonwettingOut_;
    Scalar fluxOut_;
    Scalar fluxIn_;
    Scalar porosityThreshold_;
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, VelocityFormulation);
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation);
};

// Returns the CFL flux-function
template<class TypeTag>
typename EvalCflFluxDefault<TypeTag>::Scalar EvalCflFluxDefault<TypeTag>::getCflFluxFunction(const Element& element)
{
    Scalar residualSatW = problem_.spatialParams().materialLawParams(element).swr();
    Scalar residualSatNw = problem_.spatialParams().materialLawParams(element).snr();

    // compute dt restriction
    Scalar volumeCorrectionFactor = 1 - residualSatW - residualSatNw;
    Scalar volumeCorrectionFactorOutW = 0;
    Scalar volumeCorrectionFactorOutNw = 0;

        Scalar satW = problem_.variables().cellData(problem_.variables().index(element)).saturation(wPhaseIdx);
        volumeCorrectionFactorOutW = std::max((satW - residualSatW), 1e-2);
        volumeCorrectionFactorOutNw = std::max((1 - satW - residualSatNw), 1e-2);

    //make sure correction is in the right range. If not: force dt to be not min-dt!
    if (volumeCorrectionFactorOutW <= 0)
    {
        volumeCorrectionFactorOutW = 1e100;
    }
    if (volumeCorrectionFactorOutNw <= 0)
    {
        volumeCorrectionFactorOutNw = 1e100;
    }

    //correct volume
    Scalar cFLFluxIn = volumeCorrectionFactor / getCFLFluxIn();
    Scalar cFLFluxOut = 0;

    if (velocityType_ == vt)
    {
        cFLFluxOut = volumeCorrectionFactor / getCFLFluxOut(-1);
    }
    else
    {
        cFLFluxOut = std::min(volumeCorrectionFactorOutW / getCFLFluxOut(wPhaseIdx), volumeCorrectionFactorOutNw / getCFLFluxOut(nPhaseIdx));
    }

    //determine timestep
    Scalar cFLFluxFunction = std::min(cFLFluxIn, cFLFluxOut);

    return cFLFluxFunction;
}

}

#endif
