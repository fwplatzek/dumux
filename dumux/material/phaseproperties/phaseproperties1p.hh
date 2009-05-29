// $Id$

#ifndef DUNE_PHASEPROPERTIES1P_HH
#define DUNE_PHASEPROPERTIES1P_HH

#include <dumux/material/property_baseclasses.hh>
#include <dune/common/exceptions.hh>

namespace Dune
{
/** \ingroup properties
 * @brief Fluid properties of the Interstitial Fluid
 */
class InterstitialFluid : public Liquid_GL
{
public:
    InterstitialFluid()
    {}

    double viscosity (double T, double p, double X = 0.) const
    {
        return 0.00069152;
    }


    double density (double T, double p, double X = 0.) const
    {
        return 0;
    }


    double enthalpy (double T, double p, double Xa = 0.) const
    {
        return 0;
    }


    double intEnergy(double T, double p, double Xa = 0.) const
    {
        return 0;
    }


    double diffCoeff(double T=283.15, double p=1e5) const
    {
        return 3.7378e-6;
    }


    double Xa_Max(double T, double p) const
    {
        return 0;
    }


    double p_vap(double T) const
    {
        return 0;
    }


    double henry(double T) const
    {
        return 0;
    }


    double T_vap(double p) const
    {
        return 0;
    }
};

//fluid properties of blood

class Blood: public Liquid_GL
{
public:
    Blood()
    {}

    double viscosity (double T, double p, double X = 0.) const
    {
        return 0.0069152;
    }


    double density (double T, double p, double X = 0.) const
    {
        return 0;
    }


    double enthalpy (double T, double p, double Xa = 0.) const
    {
        return 0;
    }


    double intEnergy(double T, double p, double Xa = 0.) const
    {
        return 0;
    }


    double diffCoeff(double T=283.15, double p=1e5) const
    {
        return 2.7378e-6;
    }


    double Xa_Max(double T, double p) const
    {
        return 0;
    }


    double p_vap(double T) const
    {
        return 0;
    }


    double henry(double T) const
    {
        return 0;
    }


    double T_vap(double p) const
    {
        return 0;
    }
};

}
#endif

