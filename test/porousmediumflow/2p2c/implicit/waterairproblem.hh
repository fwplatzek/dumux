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
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium.
 */
#ifndef DUMUX_WATER_AIR_PROBLEM_HH
#define DUMUX_WATER_AIR_PROBLEM_HH

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/material/fluidsystems/h2on2.hh>

#include <dumux/porousmediumflow/2p2c/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include "waterairspatialparams.hh"

#define ISOTHERMAL 0

namespace Dumux
{
template <class TypeTag>
class WaterAirProblem;

namespace Properties
{
NEW_TYPE_TAG(WaterAirProblem, INHERITS_FROM(TwoPTwoCNI, WaterAirSpatialParams));
NEW_TYPE_TAG(WaterAirBoxProblem, INHERITS_FROM(BoxModel, WaterAirProblem));
NEW_TYPE_TAG(WaterAirCCProblem, INHERITS_FROM(CCTpfaModel, WaterAirProblem));

// Set the grid type
SET_TYPE_PROP(WaterAirProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(WaterAirProblem, Problem, WaterAirProblem<TypeTag>);

// Set the wetting phase
SET_TYPE_PROP(WaterAirProblem, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(WaterAirProblem, UseMoles, true);
}


/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium. During
 *        buoyancy driven upward migration the gas passes a high
 *        temperature area.
 *
 * The domain is sized 40 m times 40 m in a depth of 1000 m. The rectangular area
 * with the increased temperature (380 K) starts at (20 m, 1 m) and ends at
 * (30 m, 30 m).
 *
 * For the mass conservation equation neumann boundary conditions are used on
 * the top and on the bottom of the domain, while dirichlet conditions
 * apply on the left and the right boundary.
 * For the energy conservation equation dirichlet boundary conditions are applied
 * on all boundaries.
 *
 * Gas is injected at the bottom boundary from 15 m to 25 m at a rate of
 * 0.001 kg/(s m), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * At the dirichlet boundaries a hydrostatic pressure, a gas saturation of zero and
 * a geothermal temperature gradient of 0.03 K/m are applied.
 *
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. The default setting for useMoles is true.
 *
 * This problem uses the \ref TwoPTwoCModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2p2cni</tt> or
 * <tt>./test_cc2p2cni</tt>
 *  */
template <class TypeTag >
class WaterAirProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Grid = typename GridView::Grid;
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

#if !ISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif

        // Phase State
        wPhaseOnly = Indices::wPhaseOnly,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        conti0EqIdx = Indices::conti0EqIdx,
        contiNEqIdx = conti0EqIdx + Indices::nCompIdx
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Intersection = typename GridView::Intersection;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \brief The constructor.
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    WaterAirProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        maxDepth_ = 1000.0; // [m]

        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM(TypeTag, std::string, Problem.Name);

        //stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout << "The problem uses mole-fractions" << std::endl;
        else
            std::cout << "The problem uses mass-fractions" << std::endl;

        this->spatialParams().plotMaterialLaw();
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

#if ISOTHERMAL
    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param element The element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param scvIdx The local vertex index (SCV index)
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    {
        return 273.15 + 10; // -> 10°C
    }
#endif

    //! \copydoc Dumux::ImplicitProblem::source()
    PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        return PrimaryVariables(0.0);
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if(globalPos[0] > 40 - eps_ || globalPos[0] < eps_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();

#if !ISOTHERMAL
        bcTypes.setDirichlet(temperatureIdx);
#endif
        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local vertex index
     * \param boundaryFaceIdx The index of the boundary face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^2*s) or kg/(m^2*s))
     */
    PrimaryVariables neumann(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);

        const auto& globalPos = scvf.ipGlobal();

        // negative values for injection
        if (globalPos[0] > 15 - eps_ && globalPos[0] < 25 + eps_ && globalPos[1] < eps_)
        {
            values[contiNEqIdx] = -1e-3/FluidSystem::molarMass(nCompIdx); //(kg/(m^2*s) or mole/(m^2*s) )
        }
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
     * \param values The initial values for the primary variables
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        auto priVars = initial_(globalPos);

#if !ISOTHERMAL
        if (globalPos[0] > 20 - eps_ && globalPos[0] < 30 + eps_ && globalPos[1] < 30 + eps_)
            priVars[temperatureIdx] = 380;
#endif
        return priVars;
    }

    /*!
     * \brief Return the initial phase state inside a control volume.
     *
     * \param vertex The vertex
     * \param vIdxGlobal The global index of the vertex
     * \param globalPos The global position
     */
    int initialPhasePresence(const SubControlVolume& scv) const
    {
        return wPhaseOnly;
    }

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        Scalar densityW = 1000.0;
        priVars[pressureIdx] = 1e5 + (maxDepth_ - globalPos[1])*densityW*9.81;
        priVars[switchIdx] = 0.0;
#if !ISOTHERMAL
        priVars[temperatureIdx] = 283.0 + (maxDepth_ - globalPos[1])*0.03;
#endif
        return priVars;
    }

    Scalar maxDepth_;
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
};

} //end namespace Dumux

#endif
