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
 * \ingroup ThreePWaterOilModel
 * \brief Adds vtk output fields specific to the twop model
 */
#ifndef DUMUX_3P2CNI_VTK_OUTPUT_FIELDS_HH
#define DUMUX_3P2CNI_VTK_OUTPUT_FIELDS_HH

#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ThreePWaterOilModel
 * \brief Adds vtk output fields specific to the three-phase three-component model
 */
template<class FluidSystem, class Indices>
class ThreePWaterOilVtkOutputFields
{

public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        // register standardized vtk output fields
        vtk.addVolumeVariable( [](const auto& v){ return v.saturation(Indices::wPhaseIdx); }, "sw");
        vtk.addVolumeVariable( [](const auto& v){ return v.saturation(Indices::nPhaseIdx); },"sn");
        vtk.addVolumeVariable( [](const auto& v){ return v.saturation(Indices::gPhaseIdx); },"sg");
        vtk.addVolumeVariable( [](const auto& v){ return v.pressure(Indices::wPhaseIdx); },"pw");
        vtk.addVolumeVariable( [](const auto& v){ return v.pressure(Indices::nPhaseIdx); },"pn");
        vtk.addVolumeVariable( [](const auto& v){ return v.pressure(Indices::gPhaseIdx); },"pg");
        vtk.addVolumeVariable( [](const auto& v){ return v.density(Indices::wPhaseIdx); },"rhow");
        vtk.addVolumeVariable( [](const auto& v){ return v.density(Indices::nPhaseIdx); },"rhon");
        vtk.addVolumeVariable( [](const auto& v){ return v.density(Indices::gPhaseIdx); },"rhog");
        vtk.addVolumeVariable( [](const auto& v){ return
        v.mobility(Indices::wPhaseIdx); },"MobW");
        vtk.addVolumeVariable( [](const auto& v){ return v.mobility(Indices::nPhaseIdx); },"MobN");
        vtk.addVolumeVariable( [](const auto& v){ return v.mobility(Indices::gPhaseIdx); },"MobG");
        vtk.addVolumeVariable( [](const auto& v){ return
        v.viscosity(Indices::wPhaseIdx); },"ViscosW");
        vtk.addVolumeVariable( [](const auto& v){ return v.viscosity(Indices::nPhaseIdx); },"ViscosN");
        vtk.addVolumeVariable( [](const auto& v){ return v.viscosity(Indices::gPhaseIdx); },"ViscosG");

        for (int i = 0; i < FluidSystem::numPhases; ++i)
            for (int j = 0; j < FluidSystem::numComponents; ++j)
                vtk.addVolumeVariable([i,j](const auto& v){ return v.moleFraction(i,j); },
                                      "x^" + FluidSystem::phaseName(i) + "_" +  FluidSystem::componentName(j));

        vtk.addVolumeVariable( [](const auto& v){ return v.porosity(); },"porosity");
        vtk.addVolumeVariable( [](const auto& v){ return v.priVars().state(); },"phase presence");
        vtk.addVolumeVariable( [](const auto& v){ return v.permeability(); },"permeability");
    }
};

} // end namespace Dumux

#endif