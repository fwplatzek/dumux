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
/**
 * \file
 * \ingroup OnePTests
 * \brief Definition of a problem, for the 1p2c problem:
 * Component transport of oxygen in interstitial fluid.
 */
#ifndef DUMUX_TISSUE_PROBLEM_HH
#define DUMUX_TISSUE_PROBLEM_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include "soilspatialparams.hh"

namespace Dumux {

template <class TypeTag>
class SoilProblem;

namespace Properties {

NEW_TYPE_TAG(SoilTypeTag, INHERITS_FROM(CCTpfaModel, Richards));

// Set the grid type
SET_TYPE_PROP(SoilTypeTag, Grid, Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 3> >);

SET_BOOL_PROP(SoilTypeTag, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(SoilTypeTag, EnableGridVolumeVariablesCache, true);
SET_BOOL_PROP(SoilTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(SoilTypeTag, SolutionDependentAdvection, false);
SET_BOOL_PROP(SoilTypeTag, SolutionDependentMolecularDiffusion, false);
SET_BOOL_PROP(SoilTypeTag, SolutionDependentHeatConduction, false);

// Set the problem property
SET_TYPE_PROP(SoilTypeTag, Problem, SoilProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(SoilTypeTag, SpatialParams, SoilSpatialParams<TypeTag>);

// Set the material law
SET_PROP(SoilTypeTag, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<RegularizedVanGenuchten<Scalar>>;
};

} // end namespace Properties


/*!
 * \ingroup OnePTests
 */
template <class TypeTag>
class SoilProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {
        // world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

public:
    SoilProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry)
    , couplingManager_(couplingManager)
    {
        //read parameters from input file
        name_ = getParam<std::string>("Problem.Name") + "_3d";
        initPressure_ = getParam<Scalar>("BoundaryConditions.InitialSoilPressure");
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

    /*!
     * \brief Returns the temperature within the domain [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // in [K]

    /*
      * \brief Returns the reference pressure [Pa] of the non-wetting
     *        fluid phase within a finite volume
     *
     * This problem assumes a constant reference pressure of 1 bar.
     */
    Scalar nonWettingReferencePressure() const
    { return 1.0e5; }


    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initialAtPos(globalPos); }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of Dumux::PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = this->couplingManager().bulkPointSources(); }

    /*!
     * \brief Evaluate the point sources (added by addPointSources)
     *        for all phases within a given sub-control-volume.
     *
     * This is the method for the case where the point source is
     * solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param pointSource A single point source
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub-control volume within the element
     *
     * For this method, the \a values() method of the point sources returns
     * the absolute rate mass generated or annihilate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void pointSource(PointSource& source,
                     const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume &scv) const
    {
        // compute source at every integration point
        const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[Indices::pressureIdx];
        const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[Indices::pressureIdx];

        const auto& spatialParams = this->couplingManager().problem(Dune::index_constant<1>{}).spatialParams();
        const auto lowDimElementIdx = this->couplingManager().pointSourceData(source.id()).lowDimElementIdx();
        const Scalar Kr = spatialParams.Kr(lowDimElementIdx);
        const Scalar rootRadius = spatialParams.radius(lowDimElementIdx);

        // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
        const auto density = 1000;
        const Scalar sourceValue = 2* M_PI *rootRadius * Kr *(pressure1D - pressure3D)*density;
        source = sourceValue*source.quadratureWeight()*source.integrationElement();

    }

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
        PrimaryVariables priVars({initPressure_});
        priVars.setState(Indices::bothPhases);
        return priVars;
    }

    //! Called after every time step
    //! Output the total global exchange term
    void computeSourceIntegral(const SolutionVector& sol, const GridVariables& gridVars)
    {
        PrimaryVariables source(0.0);
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, sol);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto pointSources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
                pointSources *= scv.volume()*elemVolVars[scv].extrusionFactor();
                source += pointSources;
            }
        }

        std::cout << "Global integrated source (soil): " << source << " (kg/s) / "
                  <<                           source*3600*24*1000 << " (g/day)" << '\n';

    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    Scalar initPressure_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace Dumux

#endif
