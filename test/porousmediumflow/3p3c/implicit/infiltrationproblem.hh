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
 * \brief Isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 */
#ifndef DUMUX_INFILTRATIONPROBLEM_HH
#define DUMUX_INFILTRATIONPROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidsystems/h2oairmesitylene.hh>

#include <dumux/porousmediumflow/3p3c/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include "infiltrationspatialparameters.hh"

namespace Dumux
{
template <class TypeTag>
class InfiltrationProblem;

namespace Properties
{
NEW_TYPE_TAG(InfiltrationProblem, INHERITS_FROM(ThreePThreeC, InfiltrationSpatialParams));
NEW_TYPE_TAG(InfiltrationBoxProblem, INHERITS_FROM(BoxModel, InfiltrationProblem));
NEW_TYPE_TAG(InfiltrationCCProblem, INHERITS_FROM(CCTpfaModel, InfiltrationProblem));

// Set the grid type
SET_TYPE_PROP(InfiltrationProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(InfiltrationProblem, Problem, Dumux::InfiltrationProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(InfiltrationProblem,
              FluidSystem,
              Dumux::FluidSystems::H2OAirMesitylene<typename GET_PROP_TYPE(TypeTag, Scalar)>);


// Maximum tolerated relative error in the Newton method
SET_SCALAR_PROP(InfiltrationProblem, NewtonMaxRelativeShift, 1e-8);

// -1 backward differences, 0: central differences, +1: forward differences
SET_INT_PROP(InfiltrationProblem, ImplicitNumericDifferenceMethod, 0);
}

/*!
 * \ingroup ThreePThreeCModel
 * \ingroup ImplicitTestProblems
 * \brief Isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 *
 * The 2D domain of this test problem is 500 m long and 10 m deep, where
 * the lower part represents a slightly inclined groundwater table, and the
 * upper part is the vadose zone.
 * A LNAPL (Non-Aqueous Phase Liquid which is lighter than water) infiltrates
 * (modelled with a Neumann boundary condition) into the vadose zone. Upon
 * reaching the water table, it spreads (since lighter than water) and migrates
 * on top of the water table in the direction of the slope.
 * On its way through the vadose zone, it leaves a trace of residually trapped
 * immobile NAPL, which can in the following dissolve and evaporate slowly,
 * and eventually be transported by advection and diffusion.
 *
 * Left and right boundaries are constant head boundaries (Dirichlet),
 * Top and bottom are Neumann boundaries, all no-flow except for the small
 * infiltration zone in the upper left part.
 *
 * This problem uses the \ref ThreePThreeCModel.
 *
 * This problem should typically be simulated for 30 days.
 * A good choice for the initial time step size is 60 s.
 * To adjust the simulation time it is necessary to edit the file test_box3p3c.input
 * or test_cc3p3c.input.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box3p3c test_box3p3c.input</tt> or
 * <tt>./test_cc3p3c test_cc3p3c.input</tt>
 */
template <class TypeTag >
class InfiltrationProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,

        // Phase State
        wgPhaseOnly = Indices::wgPhaseOnly,

        contiWEqIdx = Indices::conti0EqIdx, //!< Index of the mass conservation equation for the water component
        contiNEqIdx = Indices::conti1EqIdx,//!< Index of the mass conservation equation for the contaminant component
        contiAEqIdx = Indices::conti2EqIdx,//!< Index of the mass conservation equation for the gas component

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    InfiltrationProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView), eps_(1e-6)
    {
        temperature_ = 273.15 + 10.0; // -> 10 degrees Celsius
        FluidSystem::init(/*tempMin=*/temperature_ - 1,
                          /*tempMax=*/temperature_ + 1,
                          /*nTemp=*/3,
                          /*pressMin=*/0.8*1e5,
                          /*pressMax=*/3*1e5,
                          /*nPress=*/200);

        name_ = GET_RUNTIME_PARAM(TypeTag, std::string, Problem.Name);
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
    const std::string name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param globalPos The position
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return temperature_;
    }

    /*!
     * \brief Returns the source term at specific position in the domain.
     *
     * \param values The source values for the primary variables
     * \param globalPos The position
     */
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
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
        BoundaryTypes values;
        if(globalPos[0] > 500. - eps_)
            values.setAllDirichlet();
        else if(globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
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
        PrimaryVariables values;

        Scalar y = globalPos[1];
        Scalar x = globalPos[0];
        Scalar sw, swr=0.12, sgr=0.03;

        if(y >(-1.E-3*x+5) )
        {
            Scalar pc = 9.81 * 1000.0 * (y - (-5E-4*x+5));
            if (pc < 0.0) pc = 0.0;

            sw = invertPcgw_(pc, this->spatialParams().materialLawParamsAtPos(globalPos));
            if (sw < swr) sw = swr;
            if (sw > 1.-sgr) sw = 1.-sgr;

            values[pressureIdx] = 1e5 ;
            values[switch1Idx] = sw;
            values[switch2Idx] = 1.e-6;
        }else {
            values[pressureIdx] = 1e5 + 9.81 * 1000.0 * ((-5E-4*x+5) - y);
            values[switch1Idx] = 1.-sgr;
            values[switch2Idx] = 1.e-6;
        }

        return values;
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
     */
    PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        // negative values for injection
        if ((globalPos[0] <= 75.0+eps_) && (globalPos[0] >= 50.0+eps_) && (globalPos[1] >= 10.0-eps_))
        {
            values[contiWEqIdx] = -0.0;
            //mole flow conversion to mass flow with molar mass M(Mesit.)=0,120 kg/mol --> 1.2e-4 kg/(sm)
            //the 3p3c model uses mole fractions
            values[contiNEqIdx] = -0.001;
            values[contiAEqIdx] = -0.0;
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
        PrimaryVariables values;
        initial_(values, globalPos);
        return values;
    }

    /*!
     * \brief Return the initial phase state inside a control volume.
     *
     * \param globalPos The global position
     */
    int initialPhasePresence(const SubControlVolume& scv) const
    {
        return wgPhaseOnly;
    }

        /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer. Adjust this in case of anisotropic permeabilities.
     */
    void addOutputVtkFields()
    {
        // get the number of degrees of freedom
        unsigned numDofs = this->model().numDofs();

        // create the scalar field required for the permeabilities
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        ScalarField *Kxx = this->resultWriter().allocateManagedBuffer(numDofs);

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = this->model().fvGeometries(element);
            for (auto&& scv : fvGeometry.scvs())
            {
                auto dofIdxGlobal = scv.dofIndex();
                (*Kxx)[dofIdxGlobal] = this->spatialParams().intrinsicPermeability(scv);
            }
        }

        this->resultWriter().attachDofData(*Kxx, "permeability", isBox);
    }

private:
    // internal method for the initial condition
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar y = globalPos[1];
        Scalar x = globalPos[0];
        Scalar sw, swr=0.12, sgr=0.03;

        if(y >(-1.E-3*x+5) )
        {
            Scalar pc = 9.81 * 1000.0 * (y - (-5E-4*x+5));
            if (pc < 0.0) pc = 0.0;

            sw = invertPcgw_(pc, this->spatialParams().materialLawParamsAtPos(globalPos));
            if (sw < swr) sw = swr;
            if (sw > 1.-sgr) sw = 1.-sgr;

            values[pressureIdx] = 1e5 ;
            values[switch1Idx] = sw;
            values[switch2Idx] = 1.e-6;
        }else {
            values[pressureIdx] = 1e5 + 9.81 * 1000.0 * ((-5E-4*x+5) - y);
            values[switch1Idx] = 1.-sgr;
            values[switch2Idx] = 1.e-6;
        }
    }

    static Scalar invertPcgw_(Scalar pcIn, const MaterialLawParams &pcParams)
    {
        Scalar lower,upper;
        int k;
        int maxIt = 50;
        Scalar bisLimit = 1.;
        Scalar sw, pcgw;
        lower=0.0; upper=1.0;
        for (k=1; k<=25; k++)
        {
            sw = 0.5*(upper+lower);
            pcgw = MaterialLaw::pcgw(pcParams, sw);
            Scalar delta = pcgw-pcIn;
            if (delta<0.) delta*=-1.;
            if (delta<bisLimit)
            {
                return(sw);
            }
            if (k==maxIt) {
                return(sw);
            }
            if (pcgw>pcIn) lower=sw;
            else upper=sw;
        }
        return(sw);
    }

    Scalar temperature_;
    Scalar eps_;
    std::string name_;
};
} //end namespace

#endif
