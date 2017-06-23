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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/

/*!
* \file
*
* \brief Adaption of the fully implicit scheme to the two-phase flow model.
*/

#ifndef DUMUX_TWOPMINC_MODEL_HH
#define DUMUX_TWOPMINC_MODEL_HH

#include <dune/common/exceptions.hh>

#include <dumux/porousmediumflow/implicit/velocityoutput.hh>
#include <dumux/porousmediumflow/2p/implicit/model.hh>

#include "properties.hh"
#include "localresidual.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPMincModel
 */
template<class TypeTag >
class TwoPMincModel : public TwoPModel<TypeTag>
{
    typedef TwoPMincModel<TypeTag> ThisType;
    typedef TwoPModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        nPhaseIdx = Indices::nPhaseIdx,
        wPhaseIdx = Indices::wPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum{ numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum{ numContinua = GET_PROP_VALUE(TypeTag, NumContinua) };
    enum{ dim = GridView::dimension };
    enum{ dimWorld = GridView::dimensionworld };

    enum {
        // Options for choosing the nested volume elements
        constantVolFraction_=0,
        equidistantNestedVolElements_ = 1,
        DFM_volFraction_ = 2
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<int, dim> GridResolution;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    /*!
     * \brief The constructor.
     */
    TwoPMincModel()
        :ParentType(),
         interactingContinuaType_(GET_PARAM_FROM_GROUP(TypeTag, int, Problem, InteractingContinuaType))
    {
        static_assert(dim==2, "The TwoPMincModel only works for two dimensional grids");
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     *
     * \param sol The global solution vector
     * \param writer The writer for multi-file VTK datasets
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;

        // get the number of degrees of freedom
        unsigned numDofs = this->numDofs();

        ScalarField *pw[numContinua];
        ScalarField *pn[numContinua];
        ScalarField *pc[numContinua];
        ScalarField *sw[numContinua];
        ScalarField *sn[numContinua];
        ScalarField *rhoW[numContinua];
        ScalarField *rhoN[numContinua];
        ScalarField *mobW[numContinua];
        ScalarField *mobN[numContinua];
        ScalarField *poro[numContinua];
        ScalarField *perm[numContinua];
        ScalarField *Te;

        // create the required scalar fields
        for (int nC=0; nC<numContinua; nC++)
        {
            pw[nC] = writer.allocateManagedBuffer(numDofs);
            pn[nC] = writer.allocateManagedBuffer(numDofs);
            pc[nC] = writer.allocateManagedBuffer(numDofs);
            sw[nC] = writer.allocateManagedBuffer(numDofs);
            sn[nC] = writer.allocateManagedBuffer(numDofs);
            rhoW[nC] = writer.allocateManagedBuffer(numDofs);
            rhoN[nC] = writer.allocateManagedBuffer(numDofs);
            mobW[nC] = writer.allocateManagedBuffer(numDofs);
            mobN[nC] = writer.allocateManagedBuffer(numDofs);
            poro[nC] = writer.allocateManagedBuffer(numDofs);
            perm[nC] = writer.allocateManagedBuffer(numDofs);

        }
        Te = writer.allocateManagedBuffer(numDofs);
        VectorField *velocityN = writer.template allocateManagedBuffer<double, dim>(numDofs);
        VectorField *velocityW = writer.template allocateManagedBuffer<double, dim>(numDofs);
        ImplicitVelocityOutput<TypeTag> velocityOutput(this->problem_());

        if (velocityOutput.enableOutput()) // check if velocity output is demanded
        {
            // initialize velocity fields
            for (unsigned int i = 0; i < numDofs; ++i)
            {
                (*velocityN)[i] = Scalar(0);
                (*velocityW)[i] = Scalar(0);
            }
        }

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank = writer.allocateManagedBuffer(numElements);

        for (const auto& element : elements(this->gridView_()))
        {
            int eIdx = this->elementMapper().index(element);
            (*rank)[eIdx] = this->gridView_().comm().rank();

            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView_(), element);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                               element,
                               fvGeometry,
                               false /* TODO oldSol? */);

            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                int globalIdx = this->dofMapper().subIndex(element, scvIdx, dofCodim);
                for (int nC=0; nC < numContinua; nC++) {
                    (*pw[nC])[globalIdx] = elemVolVars[scvIdx].pressure(wPhaseIdx, nC);
                    (*pn[nC])[globalIdx] = elemVolVars[scvIdx].pressure(nPhaseIdx, nC);
                    (*pc[nC])[globalIdx] = elemVolVars[scvIdx].capillaryPressure(nC);
                    (*sw[nC])[globalIdx] = elemVolVars[scvIdx].saturation(wPhaseIdx, nC);
                    (*sn[nC])[globalIdx] = elemVolVars[scvIdx].saturation(nPhaseIdx, nC);
                    (*rhoW[nC])[globalIdx] = elemVolVars[scvIdx].density(wPhaseIdx, nC);
                    (*rhoN[nC])[globalIdx] = elemVolVars[scvIdx].density(nPhaseIdx, nC);
                    (*mobW[nC])[globalIdx] = elemVolVars[scvIdx].mobility(wPhaseIdx, nC);
                    (*mobN[nC])[globalIdx] = elemVolVars[scvIdx].mobility(nPhaseIdx, nC);
                    (*poro[nC])[globalIdx] = elemVolVars[scvIdx].porosity(nC);
                    (*perm[nC])[globalIdx] = elemVolVars[scvIdx].intrinsicPermeability(nC);
                }
                (*Te)[globalIdx] = elemVolVars[scvIdx].temperature();
            }

            // velocity output
            velocityOutput.calculateVelocity(*velocityW, elemVolVars, fvGeometry, element, wPhaseIdx);
            velocityOutput.calculateVelocity(*velocityN, elemVolVars, fvGeometry, element, nPhaseIdx);
        }
        for (int nC=0; nC<numContinua; nC++)
        {
            const std::string num = std::to_string(nC);
            writer.attachDofData(*pw[nC], "pW_" + num, isBox);
            writer.attachDofData(*pn[nC], "pN_" + num, isBox);
            writer.attachDofData(*pc[nC], "pC_" + num, isBox);
            writer.attachDofData(*sw[nC], "Sw_" + num, isBox);
            writer.attachDofData(*sn[nC], "Sn_" + num, isBox);
            writer.attachDofData(*rhoW[nC], "rhoW_" + num, isBox);
            writer.attachDofData(*rhoN[nC], "rhoN_" + num, isBox);
            writer.attachDofData(*mobW[nC], "mobW_" + num, isBox);
            writer.attachDofData(*mobN[nC], "mobM_" + num, isBox);
            writer.attachDofData(*poro[nC], "por_" + num, isBox);
            writer.attachDofData(*perm[nC], "K_" + num, isBox);
        }
        writer.attachDofData(*Te, "temperature", isBox);

        if (velocityOutput.enableOutput()) // check if velocity output is demanded
        {
            writer.attachDofData(*velocityW,  "velocityW", isBox, dim);
            writer.attachDofData(*velocityN,  "velocityN", isBox, dim);
        }
        writer.attachCellData(*rank, "process rank");
    }


    /*
     * \brief Geometric implementation of classic MINC (Pruess and Narasimhan 1982, 1985)
     */
    void calculateMincGeometricParameters(const GridResolution& res,
                                          const GlobalPosition& bboxMin,
                                          const GlobalPosition& bboxMax,
                                          const std::string& filename="empty.csv")
    {
        const Scalar length_x = bboxMax[0]- bboxMin[0];
        const Scalar length_y = bboxMax[1]- bboxMin[1];

        cellSize_[0] = length_x/res[0];
        cellSize_[1] = length_y/res[1];

        /*
         * Intialize values
         */
        for (int nC = 0; nC<numContinua; nC++)
        {
            volFraction_[nC] = 0;
            interfaceArea_[nC] = 0.0;
            distNestedContinua_[nC] =0.0;
        }

        /*
         * Equally distanced nested continua
         * Physically it is not correct to have different discretization
         * distances between x, y or z directions. The classic implementation
         * is to have equally distanced nested volumes from the nearest fracture
         */
        switch(interactingContinuaType_)
        {
        case equidistantNestedVolElements_:
            {
            /*
             * The current implentation of the MINC model
             * assumes that there is only one single distance between
             * adjacent continua faces, i.e., the equidistant
             * case only works for quadratic elements.
             */
            assert( Dune::FloatCmp::eq(cellSize_[0],cellSize_[1]) );

            distNestedContinua_ = cellSize_[0]/(2.0*numContinua-1.0);

            Dune::FieldVector <Scalar, dim> delta;
            delta[0] = cellSize_[0]/(2.0*numContinua-1.0); //the increments of the equidistanced nested elements
            delta[1] = cellSize_[1]/(2.0*numContinua-1.0);
            Dune::FieldVector <Scalar, numContinua> volContinuum;

            //LENGTH of the nested volume elements
            for (int nC=0; nC<numContinua; nC++)
            {
                for (int d=0; d<dim; d++)
                {
                    length_[d][nC] = cellSize_[d] - 2*delta[d]*nC;
                }

                volContinuum[nC] = length_[0][nC];
                for (int d=1; d<dim; d++)
                {
                    volContinuum[nC] *= length_[d][nC];
                }
            }

            //VOLUME FRACTION
            for (int nC= 0; nC<numContinua-1; nC++)
            {
                volFraction_[nC] = volContinuum[nC] - volContinuum[nC+1];
                volFraction_[nC] /= volContinuum[0];
            }
            volFraction_[numContinua-1] = volContinuum[numContinua-1]/volContinuum[0];

            }
            break;
        // EquidistantNestedVolumeElements
        /*
         * Constant volume fractions - and the interface area is determined
         */
        case constantVolFraction_:{
            Scalar x1 = 0;
            Scalar x2 = 0;

            Dune::FieldVector <Scalar, numContinua> x1_test(0.0);
            Dune::FieldVector <Scalar, numContinua> x2_test(0.0);

            length_[0][0] = cellSize_[0];
            length_[1][0] = cellSize_[1];
            Scalar sumDistNestedContinua_ = 0;
            volFraction_ = 1.0/numContinua;

            assert(numContinua>1);
            //assuming that nC=0 is the outermost frame that surrounds the matrix continua
            for (int nC=1; nC<numContinua; nC++)
            {
                /*
                 * In the following the distance between the interacting
                 * continua faces is calculated (x1,x2).
                 *
                 * For 2 continua in a quadrativ element
                 * the relation can be expressed via
                 * the volume equality as
                 *
                 * (length - 2x)^2 = length^2 - (length-2x)^2
                 *
                 * For rectangular cells the distance between the different
                 * continua faces is the same in all spatial directions.
                 *
                 * Quadratic equation aX^2 + bX + c = 0;
                 */
                const Scalar a = 4;
                const Scalar b = (-2)*(length_[0][nC-1]+length_[1][nC-1]);
                const Scalar c = (volFraction_[nC] * length_[0][0] * length_[1][0]);
                x1 = (-b-sqrt(pow(b,2)-4*a*c))/(2*a);
                x2 = (-b+sqrt(pow(b,2)-4*a*c))/(2*a);
                x1_test[nC] = x1;
                x2_test[nC] = x2;
                /*
                 * Eliminate the edges which are negative
                 * the ones which have the double bigger than the edge
                 * the ones which are bigger than the previous edge
                 */
                using std::min;
                if ((    (min<Scalar>(cellSize_[0],cellSize_[1])- 2* x1) < 0.0
                        || (x1 < 0.0)
                        || (x1_test[nC]<x1_test[nC-1])
                )
                        && ((min<Scalar>(cellSize_[0],cellSize_[1])- 2* x2) > 0.0
                                && x2 > 0.0
                                && (x2_test[nC]>x2_test[nC-1])
                        ))
                {
                    distNestedContinua_[nC] = x2;
                }
                else if (( (min<Scalar>(cellSize_[0],cellSize_[1])- 2* x2) < 0
                        || (x2 < 0)
                        || (x2_test[nC]<x2_test[nC-1])
                )
                        && ((min<Scalar>(cellSize_[0],cellSize_[1])- 2* x1) > 0
                                && (x1 > 0)
                                && (x1_test[nC]>x1_test[nC-1])
                        ))
                {
                    distNestedContinua_[nC] = x1;
                }
                else
                {
                    distNestedContinua_[nC] = min<Scalar>(x1, x2);
                    DUNE_THROW(Dune::InvalidStateException, "Check the solution of the nested continua geometric parameters");
                }

                length_[0][nC] = length_[0][nC-1] - 2*distNestedContinua_[nC];
                length_[1][nC] = length_[1][nC-1] - 2*distNestedContinua_[nC];
                sumDistNestedContinua_ += distNestedContinua_[nC];
            }

            }
            break; //ConstantVolumeFractions

        default:
            DUNE_THROW(Dune::InvalidStateException, "InteractingContinuaType not set!");
        }

        //INTEFACE AREA
        /*
         * Calculate the inteface areas between the interacting continua
         * in the same coarse block
         */
        for(int nC=0; nC<numContinua; nC++)
        {

            if (dim == 2)
            {
                interfaceArea_[nC] = 2*(length_[0][nC]+length_[1][nC]);
            }

            else if(dim == 3)
            {
                interfaceArea_[nC] = 2*(length_[0][nC]*length_[1][nC]
                                      +length_[0][nC]*length_[2][nC]
                                      +length_[1][nC]*length_[2][nC]);
            }
        }

        Valgrind::CheckDefined(interfaceArea_);
        Valgrind::CheckDefined(volFraction_);
        Valgrind::CheckDefined(cellSize_);
        Valgrind::CheckDefined(distNestedContinua_);
    }

    const Dune::FieldVector<Scalar, numContinua>& getDistNestedContinua() const
    { return distNestedContinua_; }

    const Dune::FieldVector<Scalar, numContinua>& getInterfaceArea() const
    { return interfaceArea_; }

    const Dune::FieldVector<Scalar, numContinua>& getVolFraction() const
    { return volFraction_; }


private:
     Dune::FieldVector <Scalar, dim> cellSize_;
     Dune::FieldMatrix <Scalar, dim, numContinua> length_;//dimensions of the nested elements length width
     Dune::FieldVector <Scalar, numContinua> volFraction_;
     Dune::FieldVector <Scalar, numContinua> interfaceArea_;
     Dune::FieldVector <Scalar, numContinua> distNestedContinua_; //distance between two nested continua
     int interactingContinuaType_;
};

}

#include "propertydefaults.hh"

#endif