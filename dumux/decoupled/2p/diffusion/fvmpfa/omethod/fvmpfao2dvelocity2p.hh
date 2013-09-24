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
#ifndef DUMUX_FVMPFAO2DVELOCITY2P_HH
#define DUMUX_FVMPFAO2DVELOCITY2P_HH


#include <dune/grid/common/gridenums.hh>
#include <dumux/decoupled/2p/diffusion/diffusionproperties2p.hh>
#include <dumux/decoupled/common/fv/mpfa/fvmpfaproperties.hh>
#include <dumux/decoupled/common/fv/mpfa/mpfaointeractionvolume.hh>

/**
 * @file
 * @brief  Velocity calculation using a 2-d MPFA O-method
 */

namespace Dumux
{
template<class TypeTag> class FvMpfaO2dVelocity2P
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum
        {
            dim = GridView::dimension, dimWorld = GridView::dimensionworld
        };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;
#else
    typedef typename Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<Scalar, dim> ReferenceElement;
#endif

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename SpatialParams::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;

    typedef typename Element::Geometry Geometry;
    #if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
        typedef typename Geometry::JacobianTransposed JacobianTransposed;
    #else
        typedef typename Geometry::Jacobian JacobianTransposed;
    #endif

    typedef typename GET_PROP_TYPE(TypeTag, GridTypeIndices) GridTypeIndices;

    typedef typename Dumux::FVMPFAOInteractionVolume<TypeTag> InteractionVolume;
    typedef std::vector<Dune::FieldVector<bool, 2 * dim> > InnerBoundaryVolumeFaces;

    enum
        {
            pw = Indices::pressureW,
            pn = Indices::pressureNw,
            pglobal = Indices::pressureGlobal,
            sw = Indices::saturationW,
            sn = Indices::saturationNw,
            vw = Indices::velocityW,
            vn = Indices::velocityNw,
            vt = Indices::velocityTotal
        };
    enum
        {
            wPhaseIdx = Indices::wPhaseIdx,
            nPhaseIdx = Indices::nPhaseIdx,
            pressureIdx = Indices::pressureIdx,
            saturationIdx = Indices::saturationIdx,
            pressureEqIdx = Indices::pressureEqIdx,
            satEqIdx = Indices::satEqIdx,
            numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
        };

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:
    FvMpfaO2dVelocity2P(Problem& problem) :
        problem_(problem), gravity_(problem.gravity())
    {
        density_[wPhaseIdx] = 0.;
        density_[nPhaseIdx] = 0.;
        viscosity_[wPhaseIdx] = 0.;
        viscosity_[nPhaseIdx] = 0.;

        vtkOutputLevel_ = GET_PARAM_FROM_GROUP(TypeTag, int, Vtk, OutputLevel);
    }

    void calculateInnerInteractionVolumeVelocity(InteractionVolume& interactionVolume, CellData& cellData1, CellData& cellData2, CellData& cellData3, CellData& cellData4, InnerBoundaryVolumeFaces& innerBoundaryVolumeFaces);
    void calculateBoundaryInteractionVolumeVelocity(InteractionVolume& interactionVolume, CellData& cellData, int elemIdx);

    void initialize(bool solveTwice = true)
    {
        ElementIterator element = problem_.gridView().template begin<0>();
        FluidState fluidState;
        fluidState.setPressure(wPhaseIdx, problem_.referencePressure(*element));
        fluidState.setPressure(nPhaseIdx, problem_.referencePressure(*element));
        fluidState.setTemperature(problem_.temperature(*element));
        fluidState.setSaturation(wPhaseIdx, 1.);
        fluidState.setSaturation(nPhaseIdx, 0.);
        density_[wPhaseIdx] = FluidSystem::density(fluidState, wPhaseIdx);
        density_[nPhaseIdx] = FluidSystem::density(fluidState, nPhaseIdx);
        viscosity_[wPhaseIdx] = FluidSystem::viscosity(fluidState, wPhaseIdx);
        viscosity_[nPhaseIdx] = FluidSystem::viscosity(fluidState, nPhaseIdx);

        return;
    }

    //! \brief Write data files
    /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        if (vtkOutputLevel_ > 0)
        {
            Dune::BlockVector < DimVector > &velocityWetting = *(writer.template allocateManagedBuffer<Scalar, dim>(
                                                                                                                    problem_.gridView().size(0)));
            Dune::BlockVector < DimVector > &velocityNonwetting = *(writer.template allocateManagedBuffer<Scalar, dim>(
                                                                                                                       problem_.gridView().size(0)));

            // compute update vector
            ElementIterator eItEnd = problem_.gridView().template end<0>();
            for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
            {
                // cell index
                int globalIdx = problem_.variables().index(*eIt);

                CellData & cellData = problem_.variables().cellData(globalIdx);

                Dune::FieldVector < Scalar, 2 * dim > fluxW(0);
                Dune::FieldVector < Scalar, 2 * dim > fluxNw(0);

                // run through all intersections with neighbors and boundary
                IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
                for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
                {
                    int isIndex = isIt->indexInInside();

                    fluxW[isIndex] += isIt->geometry().volume()
                        * (isIt->centerUnitOuterNormal() * cellData.fluxData().velocity(wPhaseIdx, isIndex));
                    fluxNw[isIndex] += isIt->geometry().volume()
                        * (isIt->centerUnitOuterNormal() * cellData.fluxData().velocity(nPhaseIdx, isIndex));
                }

                DimVector refVelocity(0);
                refVelocity[0] = 0.5 * (fluxW[1] - fluxW[0]);
                refVelocity[1] = 0.5 * (fluxW[3] - fluxW[2]);

                const DimVector& localPos = ReferenceElements::general(eIt->geometry().type()).position(0, 0);

                // get the transposed Jacobian of the element mapping
                const JacobianTransposed jacobianT = eIt->geometry().jacobianTransposed(localPos);

                // calculate the element velocity by the Piola transformation
                DimVector elementVelocity(0);
                jacobianT.umtv(refVelocity, elementVelocity);
                elementVelocity /= eIt->geometry().integrationElement(localPos);

                velocityWetting[globalIdx] = elementVelocity;

                refVelocity = 0;
                refVelocity[0] = 0.5 * (fluxNw[1] - fluxNw[0]);
                refVelocity[1] = 0.5 * (fluxNw[3] - fluxNw[2]);

                // calculate the element velocity by the Piola transformation
                elementVelocity = 0;
                jacobianT.umtv(refVelocity, elementVelocity);
                elementVelocity /= eIt->geometry().integrationElement(localPos);

                velocityNonwetting[globalIdx] = elementVelocity;
            }

            writer.attachCellData(velocityWetting, "wetting-velocity", dim);
            writer.attachCellData(velocityNonwetting, "non-wetting-velocity", dim);
        }

        return;
    }

private:
    Problem& problem_;
    const GlobalPosition& gravity_; //!< vector including the gravity constant

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    int vtkOutputLevel_;

    static constexpr Scalar threshold_ = 1e-15;
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, VelocityFormulation);//!< gives kind of velocity used (\f$ 0 = v_w\f$, \f$ 1 = v_n\f$, \f$ 2 = v_t\f$)
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation);//!< gives kind of pressure used (\f$ 0 = p_w\f$, \f$ 1 = p_n\f$, \f$ 2 = p_{global}\f$)
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation);//!< gives kind of saturation used (\f$ 0 = S_w\f$, \f$ 1 = S_n\f$)
};
// end of template

template<class TypeTag>
void FvMpfaO2dVelocity2P<TypeTag>::calculateInnerInteractionVolumeVelocity(InteractionVolume& interactionVolume, CellData& cellData1, CellData& cellData2, CellData& cellData3, CellData& cellData4, InnerBoundaryVolumeFaces& innerBoundaryVolumeFaces)
{
    ElementPointer & elementPointer1 = interactionVolume.getSubVolumeElement(0);
    ElementPointer & elementPointer2 = interactionVolume.getSubVolumeElement(1);
    ElementPointer & elementPointer3 = interactionVolume.getSubVolumeElement(2);
    ElementPointer & elementPointer4 = interactionVolume.getSubVolumeElement(3);

    // cell index
    int globalIdx1 = problem_.variables().index(*elementPointer1);
    int globalIdx2 = problem_.variables().index(*elementPointer2);
    int globalIdx3 = problem_.variables().index(*elementPointer3);
    int globalIdx4 = problem_.variables().index(*elementPointer4);

    // get pressure values
    Dune::FieldVector < Scalar, 2 * dim > potW(0);
    Dune::FieldVector < Scalar, 2 * dim > potNw(0);

    potW[0] = cellData1.potential(wPhaseIdx);
    potW[1] = cellData2.potential(wPhaseIdx);
    potW[2] = cellData3.potential(wPhaseIdx);
    potW[3] = cellData4.potential(wPhaseIdx);

    potNw[0] = cellData1.potential(nPhaseIdx);
    potNw[1] = cellData2.potential(nPhaseIdx);
    potNw[2] = cellData3.potential(nPhaseIdx);
    potNw[3] = cellData4.potential(nPhaseIdx);

    //get mobilities of the phases
    Dune::FieldVector < Scalar, numPhases > lambda1(cellData1.mobility(wPhaseIdx));
    lambda1[nPhaseIdx] = cellData1.mobility(nPhaseIdx);

    //compute total mobility of cell 1
    Scalar lambdaTotal1 = lambda1[wPhaseIdx] + lambda1[nPhaseIdx];

    //get mobilities of the phases
    Dune::FieldVector < Scalar, numPhases > lambda2(cellData2.mobility(wPhaseIdx));
    lambda2[nPhaseIdx] = cellData2.mobility(nPhaseIdx);

    //compute total mobility of cell 1
    Scalar lambdaTotal2 = lambda2[wPhaseIdx] + lambda2[nPhaseIdx];

    //get mobilities of the phases
    Dune::FieldVector < Scalar, numPhases > lambda3(cellData3.mobility(wPhaseIdx));
    lambda3[nPhaseIdx] = cellData3.mobility(nPhaseIdx);

    //compute total mobility of cell 1
    Scalar lambdaTotal3 = lambda3[wPhaseIdx] + lambda3[nPhaseIdx];

    //get mobilities of the phases
    Dune::FieldVector < Scalar, numPhases > lambda4(cellData4.mobility(wPhaseIdx));
    lambda4[nPhaseIdx] = cellData4.mobility(nPhaseIdx);

    //compute total mobility of cell 1
    Scalar lambdaTotal4 = lambda4[wPhaseIdx] + lambda4[nPhaseIdx];

    Scalar gn12nu14 = interactionVolume.getNtkrkNu_df(lambdaTotal1, 0, 0, 1);
    Scalar gn12nu12 = interactionVolume.getNtkrkNu_df(lambdaTotal1, 0, 0, 0);
    Scalar gn14nu14 = interactionVolume.getNtkrkNu_df(lambdaTotal1, 0, 1, 1);
    Scalar gn14nu12 = interactionVolume.getNtkrkNu_df(lambdaTotal1, 0, 1, 0);
    Scalar gn12nu23 = interactionVolume.getNtkrkNu_df(lambdaTotal2, 1, 1, 0);
    Scalar gn12nu21 = interactionVolume.getNtkrkNu_df(lambdaTotal2, 1, 1, 1);
    Scalar gn23nu23 = interactionVolume.getNtkrkNu_df(lambdaTotal2, 1, 0, 0);
    Scalar gn23nu21 = interactionVolume.getNtkrkNu_df(lambdaTotal2, 1, 0, 1);
    Scalar gn43nu32 = interactionVolume.getNtkrkNu_df(lambdaTotal3, 2, 0, 1);
    Scalar gn43nu34 = interactionVolume.getNtkrkNu_df(lambdaTotal3, 2, 0, 0);
    Scalar gn23nu32 = interactionVolume.getNtkrkNu_df(lambdaTotal3, 2, 1, 1);
    Scalar gn23nu34 = interactionVolume.getNtkrkNu_df(lambdaTotal3, 2, 1, 0);
    Scalar gn43nu41 = interactionVolume.getNtkrkNu_df(lambdaTotal4, 3, 1, 0);
    Scalar gn43nu43 = interactionVolume.getNtkrkNu_df(lambdaTotal4, 3, 1, 1);
    Scalar gn14nu41 = interactionVolume.getNtkrkNu_df(lambdaTotal4, 3, 0, 0);
    Scalar gn14nu43 = interactionVolume.getNtkrkNu_df(lambdaTotal4, 3, 0, 1);

    // compute transmissibility matrix T = CA^{-1}B+F
    Dune::FieldMatrix < Scalar, 2 * dim, 2 * dim > C(0), F(0), A(0), B(0);

    // evaluate matrix C, F, A, B
    C[0][0] = -gn12nu12;
    C[0][3] = -gn12nu14;
    C[1][0] = gn23nu21;
    C[1][1] = -gn23nu23;
    C[2][1] = gn43nu32;
    C[2][2] = gn43nu34;
    C[3][2] = -gn14nu43;
    C[3][3] = gn14nu41;

    F[0][0] = gn12nu12 + gn12nu14;
    F[1][1] = -gn23nu21 + gn23nu23;
    F[2][2] = -gn43nu34 - gn43nu32;
    F[3][3] = gn14nu43 - gn14nu41;

    A[0][0] = gn12nu12 + gn12nu21;
    A[0][1] = -gn12nu23;
    A[0][3] = gn12nu14;
    A[1][0] = -gn23nu21;
    A[1][1] = gn23nu23 + gn23nu32;
    A[1][2] = gn23nu34;
    A[2][1] = -gn43nu32;
    A[2][2] = -gn43nu34 - gn43nu43;
    A[2][3] = gn43nu41;
    A[3][0] = -gn14nu12;
    A[3][2] = gn14nu43;
    A[3][3] = -gn14nu41 - gn14nu14;

    //                        std::cout << A << "\n";

    B[0][0] = gn12nu12 + gn12nu14;
    B[0][1] = gn12nu21 - gn12nu23;
    B[1][1] = -gn23nu21 + gn23nu23;
    B[1][2] = gn23nu34 + gn23nu32;
    B[2][2] = -gn43nu34 - gn43nu32;
    B[2][3] = -gn43nu43 + gn43nu41;
    B[3][0] = -gn14nu12 - gn14nu14;
    B[3][3] = gn14nu43 - gn14nu41;

    //flux vector
    Dune::FieldVector < Scalar, 2 * dim > fluxW(0);
    Dune::FieldVector < Scalar, 2 * dim > fluxNw(0);

    // compute T
    A.invert();
    F += C.rightmultiply(B.leftmultiply(A));
    Dune::FieldMatrix < Scalar, 2 * dim, 2 * dim > T(F);

    T.mv(potW, fluxW);
    T.mv(potNw, fluxNw);

    Scalar potentialDiffW12 = fluxW[0];
    Scalar potentialDiffW14 = fluxW[3];
    Scalar potentialDiffW32 = -fluxW[1];
    Scalar potentialDiffW34 = -fluxW[2];

    Scalar potentialDiffNw12 = fluxNw[0];
    Scalar potentialDiffNw14 = fluxNw[3];
    Scalar potentialDiffNw32 = -fluxNw[1];
    Scalar potentialDiffNw34 = -fluxNw[2];

    //store potentials for further calculations (saturation, ...)
    cellData1.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(0, 0), fluxW[0]);
    cellData1.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(0, 0), fluxNw[0]);
    cellData1.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(0, 1), fluxW[3]);
    cellData1.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(0, 1), fluxNw[3]);
    cellData2.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(1, 0), fluxW[1]);
    cellData2.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(1, 0), fluxNw[1]);
    cellData2.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(1, 1), -fluxW[0]);
    cellData2.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(1, 1), -fluxNw[0]);
    cellData3.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(2, 0), -fluxW[2]);
    cellData3.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(2, 0), -fluxNw[2]);
    cellData3.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(2, 1), -fluxW[1]);
    cellData3.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(2, 1), -fluxNw[1]);
    cellData4.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(3, 0), -fluxW[3]);
    cellData4.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(3, 0), -fluxNw[3]);
    cellData4.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(3, 1), fluxW[2]);
    cellData4.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(3, 1), fluxNw[2]);

    //compute mobilities of face 1
    Dune::FieldVector < Scalar, numPhases > lambda12Upw(0.0);
    lambda12Upw[wPhaseIdx] = (potentialDiffW12 >= 0) ? lambda1[wPhaseIdx] : lambda2[wPhaseIdx];
    lambda12Upw[nPhaseIdx] = (potentialDiffNw12 >= 0) ? lambda1[nPhaseIdx] : lambda2[nPhaseIdx];

    //compute mobilities of face 4
    Dune::FieldVector < Scalar, numPhases > lambda14Upw(0.0);
    lambda14Upw[wPhaseIdx] = (potentialDiffW14 >= 0) ? lambda1[wPhaseIdx] : lambda4[wPhaseIdx];
    lambda14Upw[nPhaseIdx] = (potentialDiffNw14 >= 0) ? lambda1[nPhaseIdx] : lambda4[nPhaseIdx];

    //compute mobilities of face 2
    Dune::FieldVector < Scalar, numPhases > lambda32Upw(0.0);
    lambda32Upw[wPhaseIdx] = (potentialDiffW32 >= 0) ? lambda3[wPhaseIdx] : lambda2[wPhaseIdx];
    lambda32Upw[nPhaseIdx] = (potentialDiffNw32 >= 0) ? lambda3[nPhaseIdx] : lambda2[nPhaseIdx];

    //compute mobilities of face 3
    Dune::FieldVector < Scalar, numPhases > lambda34Upw(0.0);
    lambda34Upw[wPhaseIdx] = (potentialDiffW34 >= 0) ? lambda3[wPhaseIdx] : lambda4[wPhaseIdx];
    lambda34Upw[nPhaseIdx] = (potentialDiffNw34 >= 0) ? lambda3[nPhaseIdx] : lambda4[nPhaseIdx];

    for (int i = 0; i < numPhases; i++)
    {
        // evaluate parts of velocity
        DimVector vel12 = interactionVolume.getNormal(0, 0);
        DimVector vel14 = interactionVolume.getNormal(0, 1);
        DimVector vel23 = interactionVolume.getNormal(1, 0);
        DimVector vel21 = interactionVolume.getNormal(1, 1);
        DimVector vel34 = interactionVolume.getNormal(2, 0);
        DimVector vel32 = interactionVolume.getNormal(2, 1);
        DimVector vel41 = interactionVolume.getNormal(3, 0);
        DimVector vel43 = interactionVolume.getNormal(3, 1);

        Dune::FieldVector < Scalar, 2 * dim > flux(0);
        switch (i)
        {
        case wPhaseIdx:
        {
            flux = fluxW;
            break;
        }
        case nPhaseIdx:
        {
            flux = fluxNw;
            break;
        }
        }

        vel12 *= flux[0] / (2 * interactionVolume.getFaceArea(0, 0)); //divide by 2 because the flux is related to the half face!
        vel14 *= flux[3] / (2 * interactionVolume.getFaceArea(0, 1));
        vel23 *= flux[1] / (2 * interactionVolume.getFaceArea(1, 0));
        vel21 *= flux[0] / (2 * interactionVolume.getFaceArea(1, 1));
        vel34 *= flux[2] / (2 * interactionVolume.getFaceArea(2, 0));
        vel32 *= flux[1] / (2 * interactionVolume.getFaceArea(2, 1));
        vel41 *= flux[3] / (2 * interactionVolume.getFaceArea(3, 0));
        vel43 *= flux[2] / (2 * interactionVolume.getFaceArea(3, 1));

        Scalar lambdaT12 = lambda12Upw[wPhaseIdx] + lambda12Upw[nPhaseIdx];
        Scalar lambdaT14 = lambda14Upw[wPhaseIdx] + lambda14Upw[nPhaseIdx];
        Scalar lambdaT32 = lambda32Upw[wPhaseIdx] + lambda32Upw[nPhaseIdx];
        Scalar lambdaT34 = lambda34Upw[wPhaseIdx] + lambda34Upw[nPhaseIdx];
        Scalar fracFlow12 = (lambdaT12 > threshold_) ? lambda12Upw[i] / (lambdaT12) : 0.0;
        Scalar fracFlow14 = (lambdaT14 > threshold_) ? lambda14Upw[i] / (lambdaT14) : 0.0;
        Scalar fracFlow32 = (lambdaT32 > threshold_) ? lambda32Upw[i] / (lambdaT32) : 0.0;
        Scalar fracFlow34 = (lambdaT34 > threshold_) ? lambda34Upw[i] / (lambdaT34) : 0.0;

        vel12 *= fracFlow12;
        vel14 *= fracFlow14;
        vel23 *= fracFlow32;
        vel21 *= fracFlow12;
        vel34 *= fracFlow34;
        vel32 *= fracFlow32;
        vel41 *= fracFlow14;
        vel43 *= fracFlow34;

        if (innerBoundaryVolumeFaces[globalIdx1][interactionVolume.getIndexOnElement(0, 0)])
        {
            vel12 *= 2;
        }
        if (innerBoundaryVolumeFaces[globalIdx1][interactionVolume.getIndexOnElement(0, 1)])
        {
            vel14 *= 2;
        }
        if (innerBoundaryVolumeFaces[globalIdx2][interactionVolume.getIndexOnElement(1, 0)])
        {
            vel23 *= 2;
        }
        if (innerBoundaryVolumeFaces[globalIdx2][interactionVolume.getIndexOnElement(1, 1)])
        {
            vel21 *= 2;
        }
        if (innerBoundaryVolumeFaces[globalIdx3][interactionVolume.getIndexOnElement(2, 0)])
        {
            vel34 *= 2;
        }
        if (innerBoundaryVolumeFaces[globalIdx3][interactionVolume.getIndexOnElement(2, 1)])
        {
            vel32 *= 2;
        }
        if (innerBoundaryVolumeFaces[globalIdx4][interactionVolume.getIndexOnElement(3, 0)])
        {
            vel41 *= 2;
        }
        if (innerBoundaryVolumeFaces[globalIdx4][interactionVolume.getIndexOnElement(3, 1)])
        {
            vel43 *= 2;
        }

        //store velocities
        cellData1.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(0, 0), vel12);
        cellData1.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(0, 1), vel14);
        cellData2.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(1, 0), vel23);
        cellData2.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(1, 1), vel21);
        cellData3.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(2, 0), vel34);
        cellData3.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(2, 1), vel32);
        cellData4.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(3, 0), vel41);
        cellData4.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(3, 1), vel43);
    }
    //set velocity marker
    cellData1.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(0, 0));
    cellData1.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(0, 1));
    cellData2.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(1, 0));
    cellData2.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(1, 1));
    cellData3.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(2, 0));
    cellData3.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(2, 1));
    cellData4.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(3, 0));
    cellData4.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(3, 1));
}

template<class TypeTag>
void FvMpfaO2dVelocity2P<TypeTag>::calculateBoundaryInteractionVolumeVelocity(InteractionVolume& interactionVolume, CellData& cellData, int elemIdx)
{
        ElementPointer & elementPointer = interactionVolume.getSubVolumeElement(elemIdx);

        // get global coordinate of cell centers
        const GlobalPosition& globalPos = elementPointer->geometry().center();

        //permeability vector at boundary
        DimMatrix permeability(problem_.spatialParams().intrinsicPermeability(*elementPointer));

        //get mobilities of the phases
        Dune::FieldVector < Scalar, numPhases > lambda(cellData.mobility(wPhaseIdx));
        lambda[nPhaseIdx] = cellData.mobility(nPhaseIdx);

        for (int faceIdx = 0; faceIdx < dim; faceIdx++)
        {
            int intVolFaceIdx = interactionVolume.getFaceIndexFromSubVolume(elemIdx, faceIdx);

            if (interactionVolume.isBoundaryFace(intVolFaceIdx))
            {
                if (interactionVolume.getBoundaryType(intVolFaceIdx).isDirichlet(pressureEqIdx))
                {
                    int boundaryFaceIdx = interactionVolume.getIndexOnElement(elemIdx, faceIdx);

                    const ReferenceElement& referenceElement = ReferenceElements::general(
                            elementPointer->geometry().type());

                    const LocalPosition& localPos = referenceElement.position(boundaryFaceIdx, 1);

                    const GlobalPosition& globalPosFace = elementPointer->geometry().global(localPos);

                    DimVector distVec(globalPosFace - globalPos);
                    Scalar dist = distVec.two_norm();
                    DimVector unitDistVec(distVec);
                    unitDistVec /= dist;

                    // get pc and lambda at the boundary
                    Scalar satWBound = cellData.saturation(wPhaseIdx);
                    //check boundary sat at face 1
                    if (interactionVolume.getBoundaryType(intVolFaceIdx).isDirichlet(satEqIdx))
                    {
                        Scalar satBound = interactionVolume.getDirichletValues(intVolFaceIdx)[saturationIdx];
                        switch (saturationType_)
                        {
                        case sw:
                        {
                            satWBound = satBound;
                              break;
                        }
                        case sn:
                        {
                            satWBound = 1 - satBound;
                            break;
                        }
                        }

                    }

                    Scalar pcBound = MaterialLaw::pc(
                            problem_.spatialParams().materialLawParams(*elementPointer), satWBound);

                    Scalar gravityDiffBound = (problem_.bBoxMax() - globalPosFace) * gravity_
                            * (density_[nPhaseIdx] - density_[wPhaseIdx]);

                    pcBound += gravityDiffBound;

                    Dune::FieldVector < Scalar, numPhases
                            > lambdaBound(
                                    MaterialLaw::krw(
                                            problem_.spatialParams().materialLawParams(*elementPointer),
                                            satWBound));
                    lambdaBound[nPhaseIdx] = MaterialLaw::krn(
                            problem_.spatialParams().materialLawParams(*elementPointer), satWBound);
                    lambdaBound[wPhaseIdx] /= viscosity_[wPhaseIdx];
                    lambdaBound[nPhaseIdx] /= viscosity_[nPhaseIdx];

                    Scalar gdeltaZ = (problem_.bBoxMax()-globalPosFace) * gravity_;
                    Scalar potentialBoundW = interactionVolume.getDirichletValues(intVolFaceIdx)[pressureIdx] + density_[wPhaseIdx]*gdeltaZ;
                    Scalar potentialBoundNw = potentialBoundW;

                    //calculate potential gradients
                    switch (pressureType_)
                    {
                    case pw:
                    {
                        potentialBoundNw += pcBound;
                        break;
                    }
                    case pn:
                    {
                        //calculate potential gradients
                        potentialBoundW -= pcBound;
                        break;
                    }
                    }

                    Scalar potentialDiffW = (cellData.potential(wPhaseIdx) - potentialBoundW) / dist;
                    Scalar  potentialDiffNw = (cellData.potential(nPhaseIdx) - potentialBoundNw) / dist;

                    //store potentials for further calculations (saturation, ...)
                    cellData.fluxData().addUpwindPotential(wPhaseIdx, boundaryFaceIdx, potentialDiffW);
                    cellData.fluxData().addUpwindPotential(nPhaseIdx, boundaryFaceIdx, potentialDiffNw);

                    //calculated phase velocities from advective velocities -> capillary pressure velocity already added in pressure part!
                    DimVector velocityW(0);
                    DimVector velocityNw(0);

                    // calculate capillary pressure gradient
                    DimVector pressGradient = unitDistVec;
                    pressGradient *= (cellData.potential(wPhaseIdx) - potentialBoundW) / dist;
                    permeability.mv(pressGradient, velocityW);

                    pressGradient = unitDistVec;
                    pressGradient *= (cellData.potential(nPhaseIdx) - potentialBoundNw) / dist;
                    permeability.mv(pressGradient, velocityNw);

                    velocityW *= (potentialDiffW >= 0.) ? lambda[wPhaseIdx] : lambdaBound[wPhaseIdx];
                    velocityNw *= (potentialDiffNw >= 0.) ? lambda[nPhaseIdx] : lambdaBound[nPhaseIdx];

                    //velocity is calculated from two vertices of one intersection!
                    velocityW *= 0.5;
                    velocityNw *= 0.5;

                    //store velocities
                        velocityW += cellData.fluxData().velocity(wPhaseIdx, boundaryFaceIdx);
                        velocityNw += cellData.fluxData().velocity(nPhaseIdx, boundaryFaceIdx);
                        cellData.fluxData().setVelocity(wPhaseIdx, boundaryFaceIdx, velocityW);
                        cellData.fluxData().setVelocity(nPhaseIdx, boundaryFaceIdx, velocityNw);
                        cellData.fluxData().setVelocityMarker(boundaryFaceIdx);
                }
                else if (interactionVolume.getBoundaryType(intVolFaceIdx).isNeumann(pressureEqIdx))
                {
                    int boundaryFaceIdx = interactionVolume.getIndexOnElement(elemIdx, faceIdx);

                    const ReferenceElement& referenceElement = ReferenceElements::general(
                            elementPointer->geometry().type());

                    const LocalPosition& localPos = referenceElement.position(boundaryFaceIdx, 1);

                    const GlobalPosition& globalPosFace = elementPointer->geometry().global(localPos);

                    DimVector distVec(globalPosFace - globalPos);
                    Scalar dist = distVec.two_norm();
                    DimVector unitDistVec(distVec);
                    unitDistVec /= dist;

                    // get neumann boundary value
                    PrimaryVariables boundValues(interactionVolume.getNeumannValues(intVolFaceIdx));

                    boundValues[wPhaseIdx] /= density_[wPhaseIdx];
                    boundValues[nPhaseIdx] /= density_[nPhaseIdx];

                    DimVector velocityW(unitDistVec);
                    DimVector velocityNw(unitDistVec);

                    velocityW *= boundValues[wPhaseIdx] / (2 * interactionVolume.getFaceArea(elemIdx, faceIdx));
                    velocityNw *= boundValues[nPhaseIdx]
                            / (2 * interactionVolume.getFaceArea(elemIdx, faceIdx));

                    //store potentials for further calculations (saturation, ...)
                    cellData.fluxData().addUpwindPotential(wPhaseIdx, boundaryFaceIdx, boundValues[wPhaseIdx]);
                    cellData.fluxData().addUpwindPotential(nPhaseIdx, boundaryFaceIdx, boundValues[nPhaseIdx]);

                    //store velocities
                    velocityW += cellData.fluxData().velocity(wPhaseIdx, boundaryFaceIdx);
                    velocityNw += cellData.fluxData().velocity(nPhaseIdx, boundaryFaceIdx);
                    cellData.fluxData().setVelocity(wPhaseIdx, boundaryFaceIdx, velocityW);
                    cellData.fluxData().setVelocity(nPhaseIdx, boundaryFaceIdx, velocityNw);
                    cellData.fluxData().setVelocityMarker(boundaryFaceIdx);
                }
                else
                {
                    DUNE_THROW(Dune::NotImplemented,
                            "No valid boundary condition type defined for pressure equation!");
                }
        }
    }
}

}
// end of Dune namespace
#endif
