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
 * \brief This material law takes a material law defined for effective
 *        saturations and converts it to a material law defined on
 *        absolute saturations.
 *
 */
#ifndef DUMUX_EFF_TO_ABS_LAW_HH
#define DUMUX_EFF_TO_ABS_LAW_HH

#include "efftoabslawparams.hh"

namespace Dumux
{
/*!
 * \ingroup fluidmatrixinteractionslaws
 *
 * \brief This material law takes a material law defined for effective
 *        saturations and converts it to a material law defined on absolute
 *        saturations.
 *
 *        The idea: "material laws" (like VanGenuchten or BrooksCorey) are defined for effective saturations.
 *        The numeric calculations however are performed with absolute saturations. The EffToAbsLaw class gets
 *        the "material laws" actually used as well as the corresponding parameter container as template arguments.
 *
 *        Subsequently, the desired function (pc, sw... ) of the actually used "material laws" are called but with the
 *        saturations already converted from absolute to effective.
 *
 *        This approach makes sure that in the "material laws" only effective saturations are considered, which makes sense,
 *        as these laws only deal with effective saturations. This also allows for changing the calculation of the effective
 *        saturations easily, as this is subject of discussion / may be problem specific.
 *
 *        Additionally, handing over effective saturations to the "material laws" in stead of them calculating effective
 *        saturations prevents accidently "converting twice".
 *
 *        This boils down to:
 *        - the actual material laws (linear, VanGenuchten...) do not need to deal with any kind of conversion
 *        - the definition of the material law in the spatial parameters is not really intuitive, but using it is:
 *          Hand in values, get back values, do not deal with conversion.
 */
template <class EffLawT, class AbsParamsT = EffToAbsLawParams<typename EffLawT::Params> >
class EffToAbsLaw
{
    typedef EffLawT EffLaw;

public:
    typedef AbsParamsT Params;
    typedef typename EffLaw::Scalar Scalar;

    /*!
     * \brief The capillary pressure-saturation curve.
     *
     *
     * \param sw        Absolute saturation of the wetting phase \f$\overline{S}_w\f$. It is converted to effective saturation
     *                  and then handed over to the material law actually used for calculation.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Capillary pressure calculated by specific constitutive relation (EffLaw e.g. Brooks & Corey, van Genuchten, linear...)
     *
     */
    static Scalar pc(const Params &params, Scalar sw)
    {
        return EffLaw::pc(params, swToSwe(params, sw));
    }

    DUNE_DEPRECATED_MSG("use pc() (uncapitalized 'c') instead")
    static Scalar pC(const Params &params, Scalar sw)
    {
        return pc(params, sw);
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * \param pc        Capillary pressure \f$p_C\f$:
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *\return           Absolute wetting phase saturation calculated as inverse of (EffLaw e.g. Brooks & Corey, van Genuchten, linear...) constitutive relation.
     *
     * \return The absolute saturation of the wetting phase \f$S_w\f$
     */
    static Scalar sw(const Params &params, Scalar pc)
    {
        return sweToSw_(params, EffLaw::sw(params, pc));
    }

    DUNE_DEPRECATED_MSG("use sw() (uncapitalized 's') instead")
    static Scalar Sw(const Params &params, Scalar pc)
    {
        return sw(params, pc);
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure w.r.t the absolute saturation.
     *
     *        In this case the chain rule needs to be applied:
     \f[
             p_c = p_c( \overline S_w (S_w))
             \rightarrow p_c ^\prime = \frac{\partial  p_c}{\partial \overline S_w} \frac{\partial \overline S_w}{\partial S_w}
     \f]
     * \param sw        Absolute saturation of the wetting phase \f$\overline{S}_w\f$.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Partial derivative of \f$p_c\f$ w.r.t. effective saturation according to EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
    */
    static Scalar dpc_dsw(const Params &params, Scalar sw)
    {
        return EffLaw::dpc_dsw(params, swToSwe(params, sw) )*dswe_dsw_(params);
    }

    DUNE_DEPRECATED_MSG("use dpc_dsw() (uncapitalized 'c', 's') instead")
    static Scalar dpC_dSw(const Params &params, Scalar sw)
    {
        return dpc_dsw(params, sw);
    }

    /*!
     * \brief Returns the partial derivative of the absolute
     *        saturation w.r.t. the capillary pressure.
     *
     * In this case the chain rule needs to be applied:
     \f[
            S_w = S_w(\overline{S}_w (p_c) )
            \rightarrow S_w^\prime = \frac{\partial S_w}{\partial \overline S_w} \frac{\partial \overline S_w}{\partial p_c}
     \f]
     *
     *
     * \param pc        Capillary pressure \f$p_C\f$:
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Partial derivative of effective saturation w.r.t. \f$p_c\f$ according to EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     */
    static Scalar dsw_dpc(const Params &params, Scalar pc)
    {
        return EffLaw::dsw_dpc(params, pc)*dsw_dswe_(params);
    }

    DUNE_DEPRECATED_MSG("use dsw_dpc() (uncapitalized 's', 'c') instead")
    static Scalar dSw_dpC(const Params &params, Scalar pc)
    {
        return dsw_dpc(params, pc);
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param sw        Absolute saturation of the wetting phase \f$\overline{S}_w\f$. It is converted to effective saturation
     *                  and then handed over to the material law actually used for calculation.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Relative permeability of the wetting phase calculated as implied by EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     *
     */
    static Scalar krw(const Params &params, Scalar sw)
    {
        return EffLaw::krw(params, swToSwe(params, sw));
    };

    /*!
     * \brief The relative permeability for the non-wetting phase.
     *
     * \param sw        Absolute saturation of the wetting phase \f${S}_w\f$. It is converted to effective saturation
     *                  and then handed over to the material law actually used for calculation.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Relative permeability of the non-wetting phase calculated as implied by EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     */
    static Scalar krn(const Params &params, Scalar sw)
    {
        return EffLaw::krn(params, swToSwe(params, sw));
    }

    /*!
     * \brief Convert an absolute wetting saturation to an effective one.
     *
     * \param sw        Absolute saturation of the wetting phase \f${S}_w\f$.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Effective saturation of the wetting phase.
     */
    static Scalar swToSwe(const Params &params, Scalar sw)
    {
        return (sw - params.swr())/(1 - params.swr() - params.snr());
    }

    DUNE_DEPRECATED_MSG("use swToSwe() (uncapitalized 's') instead")
    static Scalar SwToSwe(const Params &params, Scalar sw)
    {
        return swToSwe(params, sw);
    }

    /*!
     * \brief Convert an absolute non-wetting saturation to an effective one.
     *
     * \param sn        Absolute saturation of the non-wetting phase \f${S}_n\f$.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Effective saturation of the non-wetting phase.
     */
    static Scalar snToSne(const Params &params, Scalar sn)
    {
        return (sn - params.snr())/(1 - params.swr() - params.snr());
    }

    DUNE_DEPRECATED_MSG("use snToSne() (uncapitalized 's') instead")
    static Scalar SnToSne(const Params &params, Scalar sn)
    {
        return snToSne(params, sn);
    }

private:
    /*!
     * \brief Convert an effective wetting saturation to an absolute one.
     *
     * \param swe       Effective saturation of the non-wetting phase \f$\overline{S}_n\f$.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Absolute saturation of the non-wetting phase.
     */
    static Scalar sweToSw_(const Params &params, Scalar swe)
    {
        return swe*(1 - params.swr() - params.snr()) + params.swr();
    }

    /*!
     * \brief           Derivative of the effective saturation w.r.t. the absolute saturation.
     *
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Derivative of the effective saturation w.r.t. the absolute saturation.
     */
    static Scalar dswe_dsw_(const Params &params)
    { return 1.0/(1 - params.swr() - params.snr()); }

    /*!
     * \brief           Derivative of the absolute saturation w.r.t. the effective saturation.
     *
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Derivative of the absolute saturation w.r.t. the effective saturation.
     */
    static Scalar dsw_dswe_(const Params &params)
    { return 1 - params.swr() - params.snr(); }
};
}

#endif
