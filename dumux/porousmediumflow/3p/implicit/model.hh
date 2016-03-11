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
 * \brief Adaption of the fully implicit scheme to the three-phase flow model.
 *
 * The model is designed for simulating three fluid phases with water, gas, and
 * a liquid contaminant (NAPL - non-aqueous phase liquid)
 */
#ifndef DUMUX_3P_MODEL_HH
#define DUMUX_3P_MODEL_HH

#include <dumux/porousmediumflow/implicit/velocityoutput.hh>
#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup ThreePModel
 * \brief Adaption of the fully implicit scheme to the three-phase flow model.
 *
 * This model implements three-phase flow of three fluid phases
 * \f$\alpha \in \{ water, gas, NAPL \}\f$
 * The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum.
 *
 * By inserting this into the equations for the conservation of the
 * components, the well-known multiphase flow equation is obtained.
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * The model uses commonly applied auxiliary conditions like
 * \f$S_w + S_n + S_g = 1\f$ for the saturations.
 * Furthermore, the phase pressures are related to each other via
 * capillary pressures between the fluid phases, which are functions of
 * the saturation, e.g. according to the approach of Parker et al.
 *
 * The used primary variables are gas phase pressure \f$p_g\f$,
 * water saturation \f$S_w\f$ and NAPL saturation \f$S_n\f$.
 */
template<class TypeTag>
class ThreePModel: public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,
    };


    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     *
     * \param sol The solution vector
     * \param writer The writer for multi-file VTK datasets
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dimWorld> > VectorField;

        // get the number of degrees of freedom
        unsigned numDofs = this->numDofs();

        // create the required scalar fields
        ScalarField *saturation[numPhases];
        ScalarField *pressure[numPhases];
        ScalarField *density[numPhases];

        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            saturation[phaseIdx] = writer.allocateManagedBuffer(numDofs);
            pressure[phaseIdx] = writer.allocateManagedBuffer(numDofs);
            density[phaseIdx] = writer.allocateManagedBuffer(numDofs);
        }

        ScalarField *temperature = writer.allocateManagedBuffer (numDofs);
        ScalarField *poro = writer.allocateManagedBuffer(numDofs);
        ScalarField *perm = writer.allocateManagedBuffer(numDofs);
        // VectorField *velocityN = writer.template allocateManagedBuffer<double, dimWorld>(numDofs);
        // VectorField *velocityW = writer.template allocateManagedBuffer<double, dimWorld>(numDofs);
        // VectorField *velocityG = writer.template allocateManagedBuffer<double, dimWorld>(numDofs);
        // ImplicitVelocityOutput<TypeTag> velocityOutput(this->problem_());

        // if (velocityOutput.enableOutput()) // check if velocity output is demanded
        // {
        //     // initialize velocity fields
        //     for (unsigned int i = 0; i < numDofs; ++i)
        //     {
        //         (*velocityN)[i] = Scalar(0);
        //         (*velocityW)[i] = Scalar(0);
        //         (*velocityG)[i] = Scalar(0);
        //     }
        // }

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank = writer.allocateManagedBuffer (numElements);


        for (const auto& element : elements(this->gridView_()))
        {
            if(element.partitionType() == Dune::InteriorEntity)
            {
                const auto& fvGeometry = this->fvGeometries(element);

                for (auto&& scv : fvGeometry.scvs())
                {
                    auto eIdx = scv.elementIndex();
                    auto dofIdxGlobal = scv.dofIndex();
                    (*rank)[eIdx] = this->gridView_().comm().rank();

                    const auto& volVars = this->curVolVars(scv);

                    for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
                    {
                        (*saturation[phaseIdx])[dofIdxGlobal] = volVars.saturation(phaseIdx);
                        (*pressure[phaseIdx])[dofIdxGlobal] = volVars.pressure(phaseIdx);
                        (*density[phaseIdx])[dofIdxGlobal] = volVars.density(phaseIdx);
                    }

                    (*poro)[dofIdxGlobal] = volVars.porosity();
                    (*perm)[dofIdxGlobal] = volVars.permeability();
                    (*temperature)[dofIdxGlobal] = volVars.temperature();

                    // velocity output
                    // velocityOutput.calculateVelocity(*velocityW, elemVolVars, fvGeometry, element, wPhaseIdx);
                    // velocityOutput.calculateVelocity(*velocityN, elemVolVars, fvGeometry, element, nPhaseIdx);
                    // velocityOutput.calculateVelocity(*velocityN, elemVolVars, fvGeometry, element, gPhaseIdx);
                }
            }
        }

        writer.attachDofData(*saturation[wPhaseIdx], "sw", isBox);
        writer.attachDofData(*saturation[nPhaseIdx], "sn", isBox);
        writer.attachDofData(*saturation[gPhaseIdx], "sg", isBox);
        writer.attachDofData(*pressure[wPhaseIdx], "pw", isBox);
        writer.attachDofData(*pressure[nPhaseIdx], "pn", isBox);
        writer.attachDofData(*pressure[gPhaseIdx], "pg", isBox);
        writer.attachDofData(*density[wPhaseIdx], "rhow", isBox);
        writer.attachDofData(*density[nPhaseIdx], "rhon", isBox);
        writer.attachDofData(*density[gPhaseIdx], "rhog", isBox);

        writer.attachDofData(*poro, "porosity", isBox);
        writer.attachDofData(*perm, "permeability", isBox);
        writer.attachDofData(*temperature, "temperature", isBox);

        // if (velocityOutput.enableOutput()) // check if velocity output is demanded
        // {
        //     writer.attachDofData(*velocityW,  "velocityW", isBox, dim);
        //     writer.attachDofData(*velocityN,  "velocityN", isBox, dim);
        //     writer.attachDofData(*velocityG,  "velocityG", isBox, dim);
        // }

        writer.attachCellData(*rank, "process rank");
    }

};

}

#include "propertydefaults.hh"

#endif
