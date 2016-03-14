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
 * \brief Specification of the material params for the van Genuchten
 *        capillary pressure model.
 *
 * In comparison to the 2p version, this parameter container also includes
 * the residual saturations, as their inclusion is very model-specific.
 */
#ifndef PARKERVANGEN_PARAMS_3P_HH
#define PARKERVANGEN_PARAMS_3P_HH

#include <dune/common/fvector.hh>
#include <memory>

namespace Dumux
{
/*!
 * \brief Reference implementation of a van Genuchten params
 */
template<class ScalarT>
class ParkerVanGen3PParams
{
public:
    typedef ScalarT Scalar;
    typedef std::array<Scalar,500> SimpleTable;

    ParkerVanGen3PParams()
    {betaGw_ = betaNw_ = betaGn_ = 1.0; tableKrwInitialized_ = false;}

    ParkerVanGen3PParams(Scalar vgAlpha, Scalar vgn, Scalar KdNAPL, Scalar rhoBulk,
                         Dune::FieldVector<Scalar, 4> residualSaturation, Scalar betaNw = 1.0,
                         Scalar betaGn = 1.0, Scalar betaGw = 1.0, bool regardSnr=false)
    {
        setVgAlpha(vgAlpha);
        setVgn(vgn);
        setSwr(residualSaturation[0]);
        setSnr(residualSaturation[1]);
        setSgr(residualSaturation[2]);
        setKrRegardsSnr(regardSnr);
        setKdNAPL(KdNAPL);
        setBetaNw(betaNw);
        setBetaGn(betaGn);
        setBetaGw(betaGw);
        setRhoBulk(rhoBulk);
        tableKrwInitialized_ = false;
    }

    /*!
     * \brief Return the \f$\mathrm{\alpha}\f$ shape parameter of van Genuchten's
     *        curve.
     */
    Scalar vgAlpha() const
    { return vgAlpha_; }

    /*!
     * \brief Set the \f$\mathrm{\alpha}\f$ shape parameter of van Genuchten's
     *        curve.
     * \param v Set shape parameter
     */
    void setVgAlpha(Scalar v)
    { vgAlpha_ = v; }

    /*!
     * \brief Return the \f$\mathrm{m}\f$ shape parameter of van Genuchten's
     *        curve.
     */
    Scalar vgm() const
    { return vgm_; }

    /*!
     * \brief Set the \f$\mathrm{m}\f$ shape parameter of van Genuchten's
     *        curve.
     * \param m Set shape parameter
     *
     * The \f$\mathrm{n}\f$ shape parameter is set to \f$n = \frac{1}{1 - m}\f$
     */
    void setVgm(Scalar m)
    { vgm_ = m; vgn_ = 1/(1 - vgm_); }

    /*!
     * \brief Return the \f$\mathrm{n}\f$ shape parameter of van Genuchten's
     *        curve.
     */
    Scalar vgn() const
    { return vgn_; }

    /*!
     * \brief Set the \f$\mathrm{n}\f$ shape parameter of van Genuchten's
     *        curve.
     * \param n Set shape parameter
     *
     * The \f$\mathrm{n}\f$ shape parameter is set to \f$\mathrm{m = 1 - \frac{1}{n}}\f$
     */
    void setVgn(Scalar n)
    { vgn_ = n; vgm_ = 1 - 1/vgn_; }

    /*!
     * \brief Return the residual saturation.
     * \param phaseIdx Indicator, The saturation of phases
     */
    Scalar satResidual(int phaseIdx) const
    {
        switch (phaseIdx)
        {
        case 0:
            return swr_;
        case 1:
            return snr_;
        case 2:
            DUNE_THROW(Dune::NotImplemented, "sgr for three phases not required and therefore not implemented");
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Set all residual saturations.
     */
    void setResiduals(Dune::FieldVector<Scalar, 3> residualSaturation)
    {
        setSwr(residualSaturation[0]);
        setSnr(residualSaturation[1]);
        setSgr(residualSaturation[2]);
    }


    /*!
     * \brief Return the residual wetting saturation.
     */
    Scalar swr() const
    { return swr_; }

    /*!
     * \brief Set the residual wetting saturation.
     * \param input Set residual wetting saturation
     */
    void setSwr(Scalar input)
    { swr_ = input; }

    /*!
     * \brief Return the residual non-wetting saturation.
     */
    Scalar snr() const
    { return snr_; }

    /*!
     * \brief Set the residual non-wetting saturation.
     * \param input Set the resiudal non-wetting saturation
     */
    void setSnr(Scalar input)
    { snr_ = input; }

     /*!
     * \brief Return the residual gas saturation.
     */
    Scalar sgr() const
    {
        return sgr_;
    }

     /*!
     * \brief Set the residual gas saturation.
     * \param input Set the resiudal gas saturation
     */
    void setSgr(Scalar input)
    {
         sgr_ = input;
    }

     /*!
     * \brief Set the residual total liquid saturation.
     */
    Scalar swrx() const
    {
         std::cerr << "swrx for three phases not implemented anymore. Equals swr" << std::endl;
         return swr_;
    }

     /*!
     * \brief Set the residual total liquid saturation.
     * \param v Set the resiudal gas saturation
     */
    void setSwrx(Scalar v)
    {
         std::cerr << "swrx for three phases not implemented anymore. Equals swr" << std::endl;
    }
    /*!
     * \brief defines the scaling parameters of capillary pressure between the phases (=1 for Gas-Water)
     */
    void setBetaNw(Scalar input)
    { betaNw_ = input; }

    void setBetaGn(Scalar input)
    { betaGn_ = input; }

    void setBetaGw(Scalar input)
    { betaGw_ = input; }

    /*!
     * \brief Return the values for the beta scaling parameters of capillary pressure between the phases
     */
    Scalar betaNw() const
    { return betaNw_; }

    Scalar betaGn() const
    { return betaGn_; }

    Scalar betaGw() const
    { return betaGw_; }

    /*!
     * \brief defines if residual n-phase saturation should be regarded in its relative permeability.
     * \param input Regard residual n-phase saturation
     */
    void setKrRegardsSnr(bool input)
    { krRegardsSnr_ = input; }

    /*!
     * \brief Calls if residual n-phase saturation should be regarded in its relative permeability.
     */
    bool krRegardsSnr() const
    { return krRegardsSnr_; }


    /*!
     * \brief Return the bulk density of the porous medium in \f$\mathrm{[kg/m^3]}\f$
     */
    Scalar rhoBulk() const
    { return rhoBulk_; }

    /*!
     * \brief Set the bulk density of the porous medium
     * \param input Density of the porous medium in \f$\mathrm{[kg/m^3]}\f$
     */
    void setRhoBulk(Scalar input)
    { rhoBulk_ = input; }

    /*!
     * \brief Return the adsorption coefficient
     */
    Scalar KdNAPL() const
    { return KdNAPL_; }

    /*!
     * \brief Set the adsorption coefficient
     * \param input Set the adsorption coefficient
     */
    void setKdNAPL(Scalar input)
    { KdNAPL_ = input; }

    void initTableKrw(SimpleTable &input)
    {
        tableKrw_ = input;
    }

    void initTableKrg(SimpleTable &input)
    {
        tableKrg_ = input;
    }

    const SimpleTable &tableKrw() const
    {
        return tableKrw_;
    }

    const SimpleTable &tableKrg() const
    {
        return tableKrg_;
    }

    bool tableKrwInitialized() const
    {
        return tableKrwInitialized_;
    }

    void initTableKrw(std::function<Scalar(Scalar)> f)
    {
        tableKrw_ = std::unique_ptr<SimpleTable>(new SimpleTable());
        Scalar swLow = -1.0;
        Scalar swHigh = 2.0;
        Scalar increment = (swHigh - swLow)/499;
        for(auto &i : *tableKrw_)
        {
            i = f(swLow);
            swLow += increment;
        }
        tableKrwInitialized_ = true;
    }

    void printTableKrw()
    {
        Scalar swLow = -1.0;
        Scalar swHigh = 2.0;
        Scalar increment = (swHigh - swLow)/499;
        for(auto &i : *tableKrw_)
        {
          std::cout << "sw " << swLow << ", krw " << i << std::endl;
          swLow += increment;
        }
    }

    Scalar getKrw(const Scalar sw) const
    {
        Scalar swLow = -1.0;
        Scalar swHigh = 2.0;
//         Scalar m = 1000/(swHigh-swLow);
//         Scalar c = -m;
//         Scalar y = m*sw +c;
        Scalar y = 499/(swHigh-swLow) * (sw-swLow);
        int up = y;
        int dn = up;
        if(up >= y)
            --dn;
        else
        {
            ++dn;
            std::swap(up,dn);
        }
        const Scalar tableUp = (*tableKrw_)[up];
        const Scalar tableDn = (*tableKrw_)[dn];


//         if(up < dn)
//             DUNE_THROW(Dune::InvalidStateException, "up < dn");
//         if(up == dn)
//             DUNE_THROW(Dune::InvalidStateException, "up == dn");
//         if(up>1000)
//             DUNE_THROW(Dune::InvalidStateException, "up>1000");
//         if(dn<0)
//             DUNE_THROW(Dune::InvalidStateException, "dn<0");

//         Scalar interpolated = ((*tableKrw_)[up]- (*tableKrw_)[dn])/(up-dn) * (y-dn) + (*tableKrw_)[dn];
        return (tableUp- tableDn)/(up-dn) * (y-dn) + tableDn;
//         std::cout << interpolated << std::endl;
//         if(interpolated >= (*tableKrw_)[0] &&  interpolated <= (*tableKrw_)[499])
//             return interpolated;
//         else
//         {
//             std::cout << "up: " << up << " dn: " << dn << " result: " << interpolated << std::endl;
//             std::cout << "(*tableKrw_)[0] " << (*tableKrw_)[0] << " (*tableKrw_)[1000] " << (*tableKrw_)[999]  << std::endl;
//             DUNE_THROW(Dune::InvalidStateException, "bla");
//         }

    }

private:
    Scalar vgAlpha_;
    Scalar vgm_;
    Scalar vgn_;
    Scalar swr_;
    Scalar snr_;
    Scalar sgr_;

    Scalar KdNAPL_;
    Scalar rhoBulk_;

    Scalar betaNw_;
    Scalar betaGn_;
    Scalar betaGw_;

    std::unique_ptr<SimpleTable> tableKrw_;
    bool tableKrwInitialized_;
    SimpleTable tableKrg_;

    bool krRegardsSnr_ ;
};
} // namespace Dumux

#endif
