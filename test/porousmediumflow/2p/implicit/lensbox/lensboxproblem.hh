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
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_SatConditionProblem_HH
#define DUMUX_SatConditionProblem_HH

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>
#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/implicit/cellcentered/propertydefaults.hh>
#include <dumux/porousmediumflow/2p/implicit/gridadaptindicator.hh>
#include <dumux/implicit/adaptive/gridadaptinitializationindicator.hh>
#include <dumux/implicit/box/vertidxtoscvneighbormapper.hh>
#include <dumux/porousmediumflow/2p/implicit/vertextominpcelemmapper.hh>
#include "lensboxspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class SatConditionProblem;

//////////
// Specify the properties for the satcondition problem
//////////
namespace Properties
{
NEW_TYPE_TAG(SatConditionProblem, INHERITS_FROM(TwoP, SatConditionSpatialParams));
NEW_TYPE_TAG(SatConditionBoxProblem, INHERITS_FROM(BoxModel, SatConditionProblem));
NEW_TYPE_TAG(SatConditionBoxAdaptiveProblem, INHERITS_FROM(BoxModel, SatConditionProblem));
NEW_TYPE_TAG(SatConditionCCProblem, INHERITS_FROM(CCModel, SatConditionProblem));
NEW_TYPE_TAG(SatConditionCCAdaptiveProblem, INHERITS_FROM(CCModel, SatConditionProblem));

#if HAVE_UG
SET_TYPE_PROP(SatConditionCCProblem, Grid, Dune::UGGrid<2>);
SET_TYPE_PROP(SatConditionBoxProblem, Grid, Dune::UGGrid<2>);
#else
SET_TYPE_PROP(SatConditionCCProblem, Grid, Dune::YaspGrid<2>);
SET_TYPE_PROP(SatConditionBoxProblem, Grid, Dune::YaspGrid<2>);
#endif

#if HAVE_DUNE_ALUGRID
SET_TYPE_PROP(SatConditionBoxAdaptiveProblem, Grid, Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>);
SET_TYPE_PROP(SatConditionCCAdaptiveProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);
#endif

// Set the problem property
SET_TYPE_PROP(SatConditionProblem, Problem, Dumux::SatConditionProblem<TypeTag>);

// Set the wetting phase
SET_PROP(SatConditionProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::FluidSystems::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(SatConditionProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::FluidSystems::LiquidPhase<Scalar, Dumux::DNAPL<Scalar> > type;
};

// Linear solver settings
SET_TYPE_PROP(SatConditionCCProblem, LinearSolver, Dumux::ILU0BiCGSTABBackend<TypeTag> );
SET_TYPE_PROP(SatConditionBoxProblem, LinearSolver, Dumux::ILU0BiCGSTABBackend<TypeTag> );
#if HAVE_DUNE_ALUGRID
SET_TYPE_PROP(SatConditionCCAdaptiveProblem, LinearSolver, Dumux::ILU0BiCGSTABBackend<TypeTag> );
SET_TYPE_PROP(SatConditionBoxAdaptiveProblem, LinearSolver, Dumux::ILU0BiCGSTABBackend<TypeTag> );

SET_BOOL_PROP(SatConditionCCAdaptiveProblem, AdaptiveGrid, true);
SET_TYPE_PROP(SatConditionCCAdaptiveProblem, AdaptionIndicator, TwoPImplicitGridAdaptIndicator<TypeTag>);
SET_TYPE_PROP(SatConditionCCAdaptiveProblem,  AdaptionInitializationIndicator, ImplicitGridAdaptInitializationIndicator<TypeTag>);

SET_BOOL_PROP(SatConditionBoxAdaptiveProblem, AdaptiveGrid, true);
SET_TYPE_PROP(SatConditionBoxAdaptiveProblem, AdaptionIndicator, TwoPImplicitGridAdaptIndicator<TypeTag>);
SET_TYPE_PROP(SatConditionBoxAdaptiveProblem,  AdaptionInitializationIndicator, ImplicitGridAdaptInitializationIndicator<TypeTag>);
#endif

NEW_PROP_TAG(BaseProblem);
SET_TYPE_PROP(SatConditionBoxProblem, BaseProblem, ImplicitPorousMediaProblem<TypeTag>);
SET_TYPE_PROP(SatConditionCCProblem, BaseProblem, ImplicitPorousMediaProblem<TypeTag>);
#if HAVE_DUNE_ALUGRID
SET_TYPE_PROP(SatConditionCCAdaptiveProblem, BaseProblem, ImplicitPorousMediaProblem<TypeTag>);
SET_TYPE_PROP(SatConditionBoxAdaptiveProblem, BaseProblem, ImplicitPorousMediaProblem<TypeTag>);
#endif
}

/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is a vertical column sized 0.1m times 0.5m and features a
 * low permeablility zone which spans from a height of 0.155 m to 0.355m
 * and is bounded above and below by a medium with higher permability.
 * The problem is simulated quasi 1D so there is only one cell in
 * x-direction.
 *
 * On the top, left and right of the domain neumann boundary
 * conditions are used, while dirichlet conditions apply on the bottom.
 *
 * DNAPL is injected at the top boundary at a rate of
 * 0.05 kg/(s m^2), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 */
template <class TypeTag >
class SatConditionProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;
    typedef typename GridView::Grid Grid;
    enum {dim=Grid::dimension};
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;

    enum {

        // primary variable indices
        pwIdx = Indices::pwIdx,
        snIdx = Indices::snIdx,

        // equation indices
        contiNEqIdx = Indices::contiNEqIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,


        // world dimension
        dimWorld = GridView::dimensionworld
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dumux::VertIdxToScvNeighborMapper<GridView> VertIdxToScvNeighborMapper;
    typedef typename Dumux::VertIdxToMinPcMapper<TypeTag> VertIdxToMinPcMapper;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    SatConditionProblem(TimeManager &timeManager,
                const GridView &gridView)
    : ParentType(timeManager, gridView),
    vertIdxToScvNeighborMapper_(gridView),
    vertIdxToMinPcMapper_(gridView, this->spatialParams())
    {
        eps_ = 3e-6;
        temperature_ = 273.15 + 20; // -> 20°C

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);
        vertIdxToMinPcMapper_.update();
        vertIdxToScvNeighborMapper_.update();

    }

    /*!
     * \name Problem parameters
     */
    // \{

        /*!
     * \brief Returns the vertIdxToScvNeighborMapper
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const VertIdxToScvNeighborMapper &vertIdxToScvNeighborMapper() const
    {
        return vertIdxToScvNeighborMapper_;
    }

    const VertIdxToMinPcMapper &vertIdxToMinPcMapper() const
    {
        return vertIdxToMinPcMapper_;
    }


    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    {
        return name_.c_str();
    }

    /*!
     * \brief User defined output after the time integration
     *
     * Will be called diretly after the time integration.
     */
    void postTimeStep()
    {
        // Calculate storage terms
        PrimaryVariables storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout<<"Storage: " << storage << std::endl;
        }
    }

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * This problem assumes a uniform temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; }

    /*!
     * \brief Returns the source term
     *
     * \param values Stores the source values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param globalPos The global position
     */
    void sourceAtPos(PrimaryVariables &values,
                const GlobalPosition &globalPos) const
    {
        values = 0.0;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param values Stores the value of the boundary type
     * \param globalPos The global position
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
            const GlobalPosition &globalPos) const
    {
        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)) {
            values.setAllDirichlet();
        }
        else {
            values.setAllNeumann();
        }
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;
        fluidState.setTemperature(temperature_);
        fluidState.setPressure(FluidSystem::wPhaseIdx, /*pressure=*/1.0e5);
        fluidState.setPressure(FluidSystem::nPhaseIdx, /*pressure=*/1.0e5);

        Scalar densityW = FluidSystem::density(fluidState, FluidSystem::wPhaseIdx);

        Scalar height = this->bBoxMax()[1] - this->bBoxMin()[1];
        Scalar depth = this->bBoxMax()[1] - globalPos[1];
        Scalar alpha = 1 + 1.5/height;
        Scalar width = this->bBoxMax()[0] - this->bBoxMin()[0];
        Scalar factor = (width*alpha + (1.0 - alpha)*globalPos[0])/width;

        // hydrostatic pressure scaled by alpha
        values[pwIdx] = 1.0e5 - factor*densityW*this->gravity()[1]*depth;
        values[snIdx] = 0.0;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumannAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        values = 0.0;
        if (onInlet_(globalPos)) {
            values[contiNEqIdx] = -0.04; // kg / (m * s)
        }
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{


    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;
        fluidState.setTemperature(temperature_);
        fluidState.setPressure(FluidSystem::wPhaseIdx, /*pressure=*/1.0e5);
        fluidState.setPressure(FluidSystem::nPhaseIdx, /*pressure=*/1.0e5);

        Scalar densityW = FluidSystem::density(fluidState, FluidSystem::wPhaseIdx);

        Scalar depth = this->bBoxMax()[1] - globalPos[1];

        // hydrostatic pressure
        values[pwIdx] = 1.0e5 - densityW*this->gravity()[1]*depth;
        values[snIdx] = 0.0;
    }
    // \}

private:

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->bBoxMax()[1] - eps_;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->bBoxMax()[0] - this->bBoxMin()[0];
        Scalar lambda = (this->bBoxMax()[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda && lambda < 2.0/3.0;
    }

    Scalar temperature_;
    Scalar eps_;
    std::string name_;
    VertIdxToScvNeighborMapper vertIdxToScvNeighborMapper_;
    VertIdxToMinPcMapper vertIdxToMinPcMapper_;
};
} //end namespace

#endif
