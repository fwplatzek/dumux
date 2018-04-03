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
 * \ingroup ZeroEqModel
 * \copydoc Dumux::EddyViscosityModels
 */
#ifndef DUMUX_EDDYVISCOSITY_MODELS_HH
#define DUMUX_EDDYVISCOSITY_MODELS_HH

namespace Dumux {

/*!
 * \ingroup ZeroEqModel
 * \brief The available eddy viscosity models
 *
 * The following models are available:
 *  -# Prandtl's mixing length, e.g. \cite Oertel2012a
 *  -# Van-Driest modification, \cite vanDriest1956a and \cite Hanna1981a
 *  -# Baldwin-Lomax, \cite Baldwin1978a
 */
class EddyViscosityModels
{
public:
    static constexpr int none = 0;
    static constexpr int prandtl = 1;
    static constexpr int modifiedVanDriest = 2;
    static constexpr int baldwinLomax = 3;
};

} // end namespace Dumux

#endif