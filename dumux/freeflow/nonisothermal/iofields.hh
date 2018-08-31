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
 * \ingroup FreeflowNIModel
 * \copydoc Dumux::FreeflowNonIsothermalIOFields
 */
#ifndef DUMUX_FREEFLOW_NI_IO_FIELDS_HH
#define DUMUX_FREEFLOW_NI_IO_FIELDS_HH

#include <dune/common/deprecated.hh>

#include <dumux/io/name.hh>

namespace Dumux
{

/*!
 * \ingroup FreeflowNIModel
 * \brief Adds I/O fields specific to non-isothermal free-flow models
 */
template<class IsothermalIOFields, class ModelTraits>
class FreeflowNonIsothermalIOFields
{

public:

    //! Initialize the non-isothermal specific output fields.
    template <class OutputModule>
    DUNE_DEPRECATED_MSG("use initOutputModule instead")
    static void init(OutputModule& out)
    {
        initOutputModule(out);
    }

    //! Add the non-isothermal specific output fields.
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        IsothermalIOFields::initOutputModule(out);

        out.addVolumeVariable([](const auto& v){ return v.temperature(); }, IOName::temperature());
        out.addVolumeVariable([](const auto& v){ return v.thermalConductivity(); }, "lambda");
        if (ModelTraits::usesTurbulenceModel())
            out.addVolumeVariable([](const auto& v){ return v.effectiveThermalConductivity() - v.thermalConductivity(); }, "lambda_t");
    }

    //! return the names of the primary variables
    template <class FluidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        if (pvIdx < ModelTraits::numEq() - 1)
            return IsothermalIOFields::template primaryVariableName<FluidSystem>(pvIdx, state);
        else
            return IOName::temperature();
    }
};

} // end namespace Dumux

#endif
