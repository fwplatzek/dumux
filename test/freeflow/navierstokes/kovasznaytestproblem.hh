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
 * \ingroup NavierStokesTests
 * \brief Test for the staggered grid Navier-Stokes model with analytical solution (Kovasznay 1947)
 */
#ifndef DUMUX_KOVASZNAY_TEST_PROBLEM_HH
#define DUMUX_KOVASZNAY_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include "l2error.hh"

namespace Dumux
{
template <class TypeTag>
class KovasznayTestProblem;

namespace Properties
{
NEW_TYPE_TAG(KovasznayTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokes));

// the fluid system
SET_PROP(KovasznayTestTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
SET_TYPE_PROP(KovasznayTestTypeTag, Grid, Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the problem property
SET_TYPE_PROP(KovasznayTestTypeTag, Problem, Dumux::KovasznayTestProblem<TypeTag> );

SET_BOOL_PROP(KovasznayTestTypeTag, EnableFVGridGeometryCache, true);

SET_BOOL_PROP(KovasznayTestTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(KovasznayTestTypeTag, EnableGridVolumeVariablesCache, true);

SET_BOOL_PROP(KovasznayTestTypeTag, EnableInertiaTerms, true);
}

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the staggered grid (Kovasznay 1947)
 * \todo doc me!
 */
template <class TypeTag>
class KovasznayTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

    static constexpr auto dimWorld = GET_PROP_TYPE(TypeTag, GridView)::dimensionworld;
    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

public:
    KovasznayTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        printL2Error_ = getParam<bool>("Problem.PrintL2Error");

        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
        Scalar reynoldsNumber = 1.0 / kinematicViscosity_;
        lambda_ = 0.5 * reynoldsNumber
                        - std::sqrt(reynoldsNumber * reynoldsNumber * 0.25 + 4.0 * M_PI * M_PI);

        using CellArray = std::array<unsigned int, dimWorld>;
        const auto numCells = getParam<CellArray>("Grid.Cells");

        cellSizeX_ = (this->fvGridGeometry().bBoxMax()[0] - this->fvGridGeometry().bBoxMin()[0]) / numCells[0];
        cellSizeY_ = (this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1]) / numCells[1];

        createAnalyticalSolution_();
    }

   /*!
     * \name Problem parameters
     */
    // \{

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    void postTimeStep(const SolutionVector& curSol) const
    {
        if(printL2Error_)
        {
            using L2Error = NavierStokesTestL2Error<Scalar, ModelTraits, PrimaryVariables>;
            const auto l2error = L2Error::calculateL2Error(*this, curSol);
            const int numCellCenterDofs = this->fvGridGeometry().numCellCenterDofs();
            const int numFaceDofs = this->fvGridGeometry().numFaceDofs();
            std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
                    << std::setw(6) << numCellCenterDofs << " cc dofs and " << numFaceDofs << " face dofs (total: " << numCellCenterDofs + numFaceDofs << "): "
                    << std::scientific
                    << "L2(p) = " << l2error.first[Indices::pressureIdx] << " / " << l2error.second[Indices::pressureIdx]
                    << " , L2(vx) = " << l2error.first[Indices::velocityXIdx] << " / " << l2error.second[Indices::velocityXIdx]
                    << " , L2(vy) = " << l2error.first[Indices::velocityYIdx] << " / " << l2error.second[Indices::velocityYIdx]
                    << std::endl;
        }
    }

   /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 298.0; }


   /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

    // \}
   /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // set Dirichlet values for the velocity everywhere
        values.setDirichlet(Indices::velocityXIdx);
        values.setDirichlet(Indices::velocityYIdx);

        // set a fixed pressure in one cell
        if (isLowerLeftCell_(globalPos))
            values.setDirichletCell(Indices::pressureIdx);

        return values;
    }

   /*!
     * \brief Return dirichlet boundary values at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition & globalPos) const
    {
        // use the values of the analytical solution
        return analyticalSolution(globalPos);
    }

   /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        PrimaryVariables values;
        values[Indices::pressureIdx] = 0.5 * (1.0 - std::exp(2.0 * lambda_ * x));
        values[Indices::velocityXIdx] = 1.0 - std::exp(lambda_ * x) * std::cos(2.0 * M_PI * y);
        values[Indices::velocityYIdx] = 0.5 * lambda_ / M_PI * std::exp(lambda_ * x) * std::sin(2.0 * M_PI * y);

        return values;
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values[Indices::pressureIdx] = 0.0;
        values[Indices::velocityXIdx] = 0.0;
        values[Indices::velocityYIdx] = 0.0;

        if (globalPos[0] < (this->fvGridGeometry().bBoxMin()[0] + eps_) ||
            globalPos[0] > (this->fvGridGeometry().bBoxMax()[0] - eps_) ||
            globalPos[1] < (this->fvGridGeometry().bBoxMin()[1] + eps_) ||
            globalPos[1] > (this->fvGridGeometry().bBoxMax()[1] - eps_))
        {
            values[Indices::velocityXIdx] = dirichletAtPos(globalPos)[Indices::velocityXIdx];
            values[Indices::velocityYIdx] = dirichletAtPos(globalPos)[Indices::velocityYIdx];
        }

        if (isLowerLeftCell_(globalPos))
            values[Indices::pressureIdx] = dirichletAtPos(globalPos)[Indices::pressureIdx];

        return values;
    }

   /*!
     * \brief Returns the analytical solution for the pressure
     */
    auto& getAnalyticalPressureSolution() const
    {
        return analyticalPressure_;
    }

   /*!
     * \brief Returns the analytical solution for the velocity
     */
    auto& getAnalyticalVelocitySolution() const
    {
        return analyticalVelocity_;
    }

   /*!
     * \brief Returns the analytical solution for the velocity at the faces
     */
    auto& getAnalyticalVelocitySolutionOnFace() const
    {
        return analyticalVelocityOnFace_;
    }

private:

   /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    void createAnalyticalSolution_()
    {
        analyticalPressure_.resize(this->fvGridGeometry().numCellCenterDofs());
        analyticalVelocity_.resize(this->fvGridGeometry().numCellCenterDofs());
        analyticalVelocityOnFace_.resize(this->fvGridGeometry().numFaceDofs());

        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto ccDofIdx = scv.dofIndex();
                auto ccDofPosition = scv.dofPosition();
                auto analyticalSolutionAtCc = analyticalSolution(ccDofPosition);

                // velocities on faces
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    const auto faceDofIdx = scvf.dofIndex();
                    const auto faceDofPosition = scvf.center();
                    const auto dirIdx = scvf.directionIndex();
                    const auto analyticalSolutionAtFace = analyticalSolution(faceDofPosition);
                    analyticalVelocityOnFace_[faceDofIdx][dirIdx] = analyticalSolutionAtFace[Indices::velocity(dirIdx)];
                }

                analyticalPressure_[ccDofIdx] = analyticalSolutionAtCc[Indices::pressureIdx];

                for(int dirIdx = 0; dirIdx < ModelTraits::dim(); ++dirIdx)
                    analyticalVelocity_[ccDofIdx][dirIdx] = analyticalSolutionAtCc[Indices::velocity(dirIdx)];
            }
        }
    }

    bool isLowerLeftCell_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < (this->fvGridGeometry().bBoxMin()[0] + 0.5*cellSizeX_ + eps_);
    }

    Scalar eps_;
    Scalar cellSizeX_;
    Scalar cellSizeY_;

    Scalar kinematicViscosity_;
    Scalar lambda_;
    bool printL2Error_;
    std::vector<Scalar> analyticalPressure_;
    std::vector<VelocityVector> analyticalVelocity_;
    std::vector<VelocityVector> analyticalVelocityOnFace_;
};
} //end namespace

#endif
