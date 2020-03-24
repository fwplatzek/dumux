This tutorial was copied from dumux/test/porousmediumflow/tracer/1ptracer.

# One-phase flow with random permeability distribution and a tracer model

## Problem set-up
This example contains a contaminant transported by a base groundwater flow in a randomly distributed permeability field. The figure below shows the simulation set-up. The permeability values range between 6.12e-15 and 1.5 e-7 $`m^2`$. A pressure gradient between the top and the bottom boundary leads to a groundwater flux from the bottom to the top. Neumann no-flow boundaries are assigned to the left and right boundary. Initially, there is a contaminant concentration at the bottom of the domain.

![](./img/setup.png)

## Model description
Two different models are applied to simulate the system: In a first step, the groundwater velocity is evaluated under stationary conditions using the single phase model.
In a second step, the contaminant is transported with the groundwater velocity field. It is assumed, that the dissolved contaminant does not affect density and viscosity of the groundwater, and thus, it is handled as a tracer by the tracer model. The tracer model is then solved instationarily.

### 1p Model
The single phase model uses Darcy's law as the equation for the momentum conservation:

```math
\textbf v = - \frac{\textbf K}{\mu} \left(\textbf{grad}\, p - \varrho {\textbf g} \right),
```

with the darcy velocity $` \textbf v `$, the permeability $` \textbf K`$, the dynamic viscosity $` \mu`$, the pressure $`p`$, the density $`\rho`$ and the gravity $`\textbf g`$.

Darcy's law is inserted into the mass balance equation:

```math
\phi \frac{\partial \varrho}{\partial t} + \text{div} \textbf v = 0,
```

where $`\phi`$ is the porosity.

The equation is discretized using cell-centered finite volumes with two-point flux approximation as spatial discretization scheme for the pressure as primary variable. For details on the discretization schemes available in DuMuX, have a look at the [handbook](https://dumux.org/handbook).

### Tracer Model
The transport of the contaminant component $`\kappa`$ occurs with the velocity field $`\textbf v`$
that is computed with the __1p model__ (see above):

```math
\phi \frac{ \partial \varrho X^\kappa}{\partial t} - \text{div} \left\lbrace \varrho X^\kappa {\textbf v} + \varrho D^\kappa_\text{pm} \textbf{grad} X^\kappa \right\rbrace = 0,
```

where $`X^\kappa`$ is the mass fraction of the contaminant component $`\kappa`$ and $` D^\kappa_\text{pm} `$ is the effective diffusivity.
The effective diffusivity is a function of the diffusion coefficient of the component $`D^\kappa`$ and the porosity and tortuosity $`\tau`$ of the porous medium (see [dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh)):

```math
D^\kappa_\text{pm}= \phi \tau D^\kappa.
```

The primary variable of this model is the mass fraction $`X^\kappa`$. We apply the same spatial discretization as in the single phase model and use the implicit Euler method for time discretization. For more information, have a look at the dumux [handbook](https://dumux.org/handbook).

In the following, we take a close look at the files containing the set-up: The boundary conditions and spatially distributed parameters for the single phase model are set in `problem_1p.hh` and `spatialparams_1p.hh`.
For the tracer model, this is done in the files `problem_tracer.hh` and `spatialparams_tracer.hh`, respectively. Afterwards, we show the different steps for solving the model in the source file `main.cc`. Finally, some simulation results are shown.


## The file `spatialparams_1p.hh`


This file contains the __spatial parameters class__ which defines the
distributions for the porous medium parameters permeability and porosity
over the computational grid

In this example, we use a randomly generated and element-wise distributed
permeability field. For this, we use the random number generation facilitie
provided by the C++ standard library.
```cpp
#include <random>
```
We include the spatial parameters class for single-phase models discretized
by finite volume schemes, from which the spatial parameters defined for this
example will inherit.
```cpp
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {
```
In the `OnePTestSpatialParams` class, we define all functions needed to describe
the porous medium, e.g. porosity and permeability, for the 1p_problem.
```cpp
template<class GridGeometry, class Scalar>
class OnePTestSpatialParams
: public FVSpatialParamsOneP<GridGeometry, Scalar,
                             OnePTestSpatialParams<GridGeometry, Scalar>>
{
```
The following convenience aliases will be used throughout this class:
```cpp
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar,
                                           OnePTestSpatialParams<GridGeometry, Scalar>>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

public:
```
The spatial parameters must export the type used to define permeabilities.
Here, we are using scalar permeabilities, but tensors are also supported.
```cpp
    using PermeabilityType = Scalar;
```
### Generation of the random permeability field
We generate the random permeability field upon construction of the spatial parameters class
```cpp
    OnePTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), K_(gridGeometry->gridView().size(0), 0.0)
    {
```
The permeability of the domain and the lens are obtained from the `params.input` file.
```cpp
        permeability_ = getParam<Scalar>("SpatialParams.Permeability");
        permeabilityLens_ = getParam<Scalar>("SpatialParams.PermeabilityLens");
```
Furthermore, the position of the lens, which is defined by the position of the lower left and the upper right corners, are obtained from the input file.
```cpp
        lensLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensLowerLeft");
        lensUpperRight_ =getParam<GlobalPosition>("SpatialParams.LensUpperRight");
```
We generate random fields for the permeability using lognormal distributions, with `permeability_` as mean value and 10 % of it as standard deviation.
A separate distribution is used for the lens using `permeabilityLens_`. A permeability value is created for each element of the grid and is stored in the vector `K_`.
```cpp
        std::mt19937 rand(0);
        std::lognormal_distribution<Scalar> K(std::log(permeability_), std::log(permeability_)*0.1);
        std::lognormal_distribution<Scalar> KLens(std::log(permeabilityLens_), std::log(permeabilityLens_)*0.1);
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            const auto eIdx = gridGeometry->elementMapper().index(element);
            const auto globalPos = element.geometry().center();
            K_[eIdx] = isInLens_(globalPos) ? KLens(rand) : K(rand);
        }
    }
```
### Properties of the porous matrix
This function returns the permeability $`[m^2]`$ to be used within a sub-control volume (`scv`) inside the element `element`.
One can define the permeability as function of the primary variables on the element, which are given in the provided `ElementSolution`.
Here, we use element-wise distributed permeabilities that were randomly generated in the constructor (see above).
```cpp
    template<class ElementSolution>
    const PermeabilityType& permeability(const Element& element,
                                         const SubControlVolume& scv,
                                         const ElementSolution& elemSol) const
    {
        return K_[scv.elementIndex()];
    }

```
We set the porosity $`[-]`$ for the whole domain to a value of $`20 \%`$.
```cpp
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.2; }
```
We reference to the permeability field. This is used in the main function to write an output for the permeability field.
```cpp
    const std::vector<Scalar>& getKField() const
    { return K_; }

private:
```
The following function returns true if a given position is inside the lens.
We use an epsilon of 1.5e-7 here for floating point comparisons.
```cpp
    bool isInLens_(const GlobalPosition& globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i)
        {
            if (globalPos[i] < lensLowerLeft_[i] + 1.5e-7
                || globalPos[i] > lensUpperRight_[i] - 1.5e-7)
                return false;
        }

        return true;
    }

    GlobalPosition lensLowerLeft_, lensUpperRight_;
    Scalar permeability_, permeabilityLens_;
    std::vector<Scalar> K_;
};

} // end namespace Dumux

```



## The file `problem_1p.hh`


This file contains the __problem class__ which defines the initial and boundary
conditions for the single-phase flow simulation.

### Include files
This header contains the porous medium problem class that this class is derived from:
```cpp
#include <dumux/porousmediumflow/problem.hh>
```
This header contains the class that specifies all spatially variable parameters
related to this problem.
```cpp
#include "spatialparams_1p.hh"
```
### The problem class
We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
As this is a porous medium flow problem, we inherit from the base class `PorousMediumFlowProblem`.
```cpp
namespace Dumux {

template<class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
```
We use convenient declarations that we derive from the property system.
```cpp
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
```
This is the constructor of our problem class:
```cpp
    OnePTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry) {}
```
First, we define the type of boundary conditions depending on the location. Two types of boundary conditions
can be specified: Dirichlet or Neumann boundary condition. On a Dirichlet boundary, the values of the
primary variables need to be fixed. On a Neumann boundary condition, values for derivatives need to be fixed.
Mixed boundary conditions (different types for different equations on the same boundary) are not accepted for
cell-centered finite volume schemes.
```cpp
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;
```
we retrieve the global position, i.e. the vector with the global coordinates,
of the integration point on the boundary sub-control volume face `scvf`
```cpp
        const auto globalPos = scvf.ipGlobal();
```
we define a small epsilon value
```cpp
        Scalar eps = 1.0e-6;
```
We specify Dirichlet boundaries on the top and bottom of our domain:
```cpp
        if (globalPos[dimWorld-1] < eps || globalPos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps)
            values.setAllDirichlet();
```
The top and bottom of our domain are Neumann boundaries:
```cpp
        else
            values.setAllNeumann();

        return values;
    }
```
Second, we specify the values for the Dirichlet boundaries. We need to fix values of our primary variable
```cpp
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
```
we retreive again the global position
```cpp
        const auto& pos = scvf.ipGlobal();
        PrimaryVariables values(0);
```
and assign pressure values in [Pa] according to a pressure gradient to 1e5 Pa at the top and 1.1e5 Pa at the bottom.
```cpp
        values[0] = 1.0e+5*(1.1 - pos[dimWorld-1]*0.1);
        return values;
    }
```
We need to specify a constant temperature for our isothermal problem.
Fluid properties that depend on temperature will be calculated with this value.
```cpp
    Scalar temperature() const
    {
        return 283.15; // 10°C
    }
```
This is everything the one phase problem class contains.
```cpp
};
```
We leave the namespace Dumux.
```cpp
} // end namespace Dumux
```



## The file `properties_1p.hh`


This file defines the `TypeTag` used for the single-phase simulation, for
which we then define the necessary properties.

### Include files
The `TypeTag` defined for this simulation will inherit all properties from the
`OneP` type tag, a convenience type tag that predefines most of the required
properties for single-phase flow simulations in DuMuX. The properties that are
defined in this file are those that depend on user choices and no meaningful
default can be set.
```cpp
#include <dumux/porousmediumflow/1p/model.hh>
```
We want to use `YaspGrid`, an implementation of the dune grid interface for structured grids:
```cpp
#include <dune/grid/yaspgrid.hh>
```
In this example, we want to discretize the equations with the cell centered finite volume
scheme using two-point-flux approximation:
```cpp
#include <dumux/discretization/cctpfa.hh>
```
The fluid properties are specified in the following headers (we use liquid water as the fluid phase):
```cpp
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
```
The local residual for incompressible flow is included.
The one-phase flow model (included above) uses a default implementation of the
local residual for single-phase flow. However, in this example we are using an
incompressible fluid phase. Therefore, we are including the specialized local
residual which contains functionality to analytically compute the entries of
the Jacobian matrix. We will use this in the main file.
```cpp
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>
```
We include the problem and spatial parameters headers used for this simulation.
```cpp
#include "problem_1p.hh"
#include "spatialparams_1p.hh"
```
### Basic property definitions for the 1p problem
We enter the namespace Dumux::Properties in order to import the entire Dumux namespace for general use:
```cpp
namespace Dumux:: Properties {
```
A `TypeTag` for our simulation is created which inherits from the one-phase flow model
and the cell centered finite volume scheme with two-point-flux discretization scheme:
```cpp
namespace TTag {
struct IncompressibleTest { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
}
```
We use a structured 2D grid:
```cpp
template<class TypeTag>
struct Grid<TypeTag, TTag::IncompressibleTest> { using type = Dune::YaspGrid<2>; };
```
The problem class specifies initial and boundary conditions:
```cpp
template<class TypeTag>
struct Problem<TypeTag, TTag::IncompressibleTest> { using type = OnePTestProblem<TypeTag>; };
```
We define the spatial parameters for our simulation:
```cpp
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::IncompressibleTest>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = OnePTestSpatialParams<GridGeometry, Scalar>;
};
```
We use the local residual that contains analytic derivative methods for incompressible flow:
```cpp
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::IncompressibleTest> { using type = OnePIncompressibleLocalResidual<TypeTag>; };
```
In the following we define the fluid system to be used:
```cpp
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::IncompressibleTest>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};
```
This enables grid-wide caching of the volume variables.
```cpp
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
```
This enables grid wide caching for the flux variables.
```cpp
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
```
This enables grid-wide caching for the finite volume grid geometry
```cpp
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
```
The caches store values that were already calculated for later usage. This increases the memory demand but makes the simulation faster.
```cpp
```
We leave the namespace Dumux::Properties.
```cpp
} // end namespace Dumux::Properties
```



## The file `spatialparams_tracer.hh`


In this file, we define spatial properties of the porous medium such as the permeability and the porosity in various functions for the tracer problem.
Furthermore, spatial dependent properties of the tracer fluid system are defined and in the end two functions handle the calculated volume fluxes from the solution of the 1p problem.
We use the properties for porous medium flow models, declared in the file `properties.hh`.
```cpp
#include <dumux/porousmediumflow/properties.hh>
```
As in the 1p spatialparams, we inherit from the spatial parameters for single-phase models using finite volumes, which we include here.
```cpp
#include <dumux/material/spatialparams/fv1p.hh>
```
We enter the namespace Dumux
```cpp
namespace Dumux {
```
In the `TracerTestSpatialParams` class, we define all functions needed to describe spatially dependent parameters for the `tracer_problem`.
```cpp
template<class GridGeometry, class Scalar>
class TracerTestSpatialParams
: public FVSpatialParamsOneP<GridGeometry, Scalar,
                             TracerTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar,
                                           TracerTestSpatialParams<GridGeometry, Scalar>>;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:

    TracerTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry) {}
```
### Properties of the porous matrix
We define the same porosity for the whole domain as in the 1p spatialparams.
```cpp
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.2; }
```
We do not consider dispersivity for the tracer transport. So we set the dispersivity coefficient to zero.
```cpp
    template<class ElementSolution>
    Scalar dispersivity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    { return 0; }
```
### Properties of the fluid system
In the following, we define fluid properties that are spatial parameters in the tracer model.
They can possible vary in space but are usually constants.
Furthermore, spatially constant values of the fluid system are defined in the `TracerFluidSystem` class in `problem.hh`.
We define the fluid density to a constant value of 1000 $`\frac{kg}{m^3}`$.
```cpp
    Scalar fluidDensity(const Element &element,
                        const SubControlVolume& scv) const
    { return 1000; }
```
This interface defines the fluid molar mass within the sub-control volume `scv`.
```cpp
    Scalar fluidMolarMass(const Element &element,
                          const SubControlVolume& scv) const
    { return fluidMolarMassAtPos(scv.dofPosition()); }
```
This interface defines the fluid molar mass depending on the position in the domain.
```cpp
    Scalar fluidMolarMassAtPos(const GlobalPosition &globalPos) const
    { return 18.0; }
```
### The volume fluxes
We define a function which returns the volume flux across the given sub-control volume face `scvf`.
This flux is obtained from the vector `volumeFlux_` that contains the fluxes across al sub-control volume faces of the discretization.
This vector can be set using the `setVolumeFlux` function.
```cpp
    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    {
        return volumeFlux_[scvf.index()];
    }
```
We define a function that allows setting the volume fluxes for all sub-control volume faces of the discretization.
This is used in the main function after these fluxes have been based on the pressure solution obtained with the single-phase model.
```cpp
    void setVolumeFlux(const std::vector<Scalar>& f)
    { volumeFlux_ = f; }

private:
    std::vector<Scalar> volumeFlux_;
};

} // end namespace Dumux

```



## The file `problem_tracer.hh`


This file contains the __problem class__ which defines the initial and boundary
conditions for the tracer transport simulation.
```cpp
```
### Include files
This header contains the porous medium problem class that this class is derived from:
```cpp
#include <dumux/porousmediumflow/problem.hh>
```
This header contains the class that specifies all spatially variable parameters
related to this problem.
```cpp
#include "spatialparams_tracer.hh"
```
### The problem class
We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
As this is a porous medium flow problem, we inherit from the base class `PorousMediumFlowProblem`.
```cpp
namespace Dumux {

template <class TypeTag>
class TracerTestProblem : public PorousMediumFlowProblem<TypeTag>
{
```
We use convenient declarations that we derive from the property system.
```cpp
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
```
We create a convenience bool stating whether mole or mass fractions are used
```cpp
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
```
We create additional convenience integers to make dimWorld and numComponents available in the problem
```cpp
    static constexpr int dimWorld = GridView::dimensionworld;
    static const int numComponents = FluidSystem::numComponents;

public:
```
This is the constructor of our problem class:
```cpp
    TracerTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
```
We print to the terminal whether mole or mass fractions are used
```cpp
        if(useMoles)
            std::cout<<"problem uses mole fractions" << '\n';
        else
            std::cout<<"problem uses mass fractions" << '\n';
    }
```
We define the type of boundary conditions depending on the location.
All boundaries are set to a neumann-type flow boundary condition.
```cpp
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }
```
We specify the initial conditions for the primary variable (tracer concentration) depending on the location.
```cpp
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables initialValues(0.0);
```
The initial contamination is located at the bottom of the domain:
```cpp
        if (globalPos[1] < 0.1 + eps_)
        {
```
We chose a mole fraction of $`1e-9`$, but in case the mass fractions
are used by the model, we have to convert this value:
```cpp
            if (useMoles)
                initialValues = 1e-9;
            else
                initialValues = 1e-9*FluidSystem::molarMass(0)
                                    /this->spatialParams().fluidMolarMassAtPos(globalPos);
        }
        return initialValues;
    }
```
We implement an outflow boundary on the top of the domain and prescribe zero-flux Neumann boundary conditions on all other boundaries.
```cpp
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];
        const auto& globalPos = scvf.center();
```
This is the outflow boundary, where tracer is transported by advection with the given flux field.
```cpp
        if (globalPos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps_)
        {
            values = this->spatialParams().volumeFlux(element, fvGeometry, elemVolVars, scvf)
                     * volVars.massFraction(0, 0) * volVars.density(0)
                     / scvf.area();
            assert(values>=0.0 && "Volume flux at outflow boundary is expected to have a positive sign");
        }
```
Prescribe zero-flux Neumann boundary conditions elsewhere
```cpp
        else
            values = 0.0;

        return values;
    }

private:
```
We assign a private global variable for the epsilon:
```cpp
    static constexpr Scalar eps_ = 1e-6;
```
This is everything the tracer problem class contains.
```cpp
};
```
We leave the namespace Dumux here.
```cpp
}

```



## The file `properties_tracer.hh`


This file defines the `TypeTag` used for the tracer transport simulation, for
which we then define the necessary properties.

### Include files
As for the single-phase problem, a`TypeTag` is defined for this simulation.
Here, we inherit all properties from the `Tracer` type tag, a convenience type tag
that predefines most of the required properties for tracer transport flow simulations in DuMuX.
```cpp
#include <dumux/porousmediumflow/tracer/model.hh>
```
Again, we use YaspGrid, an implementation of the dune grid interface for structured grids:
```cpp
#include <dune/grid/yaspgrid.hh>
```
and the cell centered, two-point-flux discretization.
```cpp
#include <dumux/discretization/cctpfa.hh>
```
This includes the base class for fluid systems. We will define a custom fluid
system that inherits from that class.
```cpp
#include <dumux/material/fluidsystems/base.hh>
```
We include the problem and spatial parameters headers used for this simulation.
```cpp
#include "problem_tracer.hh"
#include "spatialparams_tracer.hh"
```
### Basic property definitions for the tracer transport problem
We enter the namespace Dumux
```cpp
namespace Dumux {
```
In the following, we create a new tracer fluid system and derive from the base fluid system.
```cpp
template<class TypeTag>
class TracerFluidSystem : public FluidSystems::Base<GetPropType<TypeTag, Properties::Scalar>,
                                                               TracerFluidSystem<TypeTag>>
{
```
We define some convenience aliases to be used inside this class.
```cpp
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
```
We specify that the fluid system only contains tracer components,
```cpp
    static constexpr bool isTracerFluidSystem()
    { return true; }
```
and that no component is the main component
```cpp
    static constexpr int getMainComponent(int phaseIdx)
    { return -1; }
```
We define the number of components of this fluid system (one single tracer component)
```cpp
    static constexpr int numComponents = 1;
```
This interface is designed to define the names of the components of the fluid system.
Here, we only have a single component, so `compIdx` should always be 0.
The component name is used for the vtk output.
```cpp
    static std::string componentName(int compIdx = 0)
    { return "tracer_" + std::to_string(compIdx); }
```
We set the phase name for the phase index (`phaseIdx`) for velocity vtk output:
Here, we only have a single phase, so `phaseIdx` should always be zero.
```cpp
    static std::string phaseName(int phaseIdx = 0)
    { return "Groundwater"; }
```
We set the molar mass of the tracer component with index `compIdx` (should again always be zero here).
```cpp
    static Scalar molarMass(unsigned int compIdx = 0)
    { return 0.300; }
```
We set the value for the binary diffusion coefficient. This
might depend on spatial parameters like pressure / temperature.
But, in this case we neglect diffusion and return 0.0:
```cpp
    static Scalar binaryDiffusionCoefficient(unsigned int compIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    { return 0.0; }
};
```
We enter the namespace Properties
```cpp
namespace Properties {
```
A `TypeTag` for our simulation is created which inherits from the tracer model and the
cell centered discretization scheme using two-point flux approximation.
```cpp
namespace TTag {
struct TracerTest { using InheritsFrom = std::tuple<Tracer>; };
struct TracerTestCC { using InheritsFrom = std::tuple<TracerTest, CCTpfaModel>; };
}
```
We enable caching for the grid volume variables, the flux variables and the FV grid geometry.
```cpp
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };
```
We use the same grid as in the stationary one-phase model, a structured 2D grid:
```cpp
template<class TypeTag>
struct Grid<TypeTag, TTag::TracerTest> { using type = Dune::YaspGrid<2>; };
```
The problem class that specifies initial and boundary conditions:
```cpp
template<class TypeTag>
struct Problem<TypeTag, TTag::TracerTest> { using type = TracerTestProblem<TypeTag>; };
```
We define the spatial parameters for our tracer simulation:
```cpp
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TracerTest>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = TracerTestSpatialParams<GridGeometry, Scalar>;
};
```
One can choose between a formulation in terms of mass or mole fractions.
Here, we are using mass fractions.
```cpp
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TracerTest> { static constexpr bool value = false; };
```
We use solution-independent molecular diffusion coefficients. Per default, solution-dependent
diffusion coefficients are assumed during the computation of the jacobian matrix entries. Specifying
solution-independent diffusion coefficients can speed up computations:
```cpp
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::TracerTestCC>
{ static constexpr bool value = false; };
```
We set the above created tracer fluid system:
```cpp
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TracerTest> { using type = TracerFluidSystem<TypeTag>; };
```
We leave the namespace Properties and Dumux.
```cpp
} // end namespace Properties
} // end namespace Dumux

```



## The file `main.cc`


We look now at the main file for the tracer problem. We set up two problems in this file and solve them sequentially, first the 1p problem and afterwards the tracer problem. The result of the 1p problem is the pressure distribution in the problem domain. We use it to calculate the volume fluxes, which act as an input for the tracer problem. Based on this volume fluxes, we calculate the transport of a tracer in the following tracer problem.
### Includes
```cpp
#include <config.h>
```
This includes the `TypeTags` and properties to be used for the single-phase
and the tracer simulations.
```cpp
#include "properties_1p.hh"
#include "properties_tracer.hh"
```
Further, we include a standard header file for C++, to get time and date information
```cpp
#include <ctime>
```
and another one for in- and output.
```cpp
#include <iostream>
```
Dumux is based on DUNE, the Distributed and Unified Numerics Environment, which provides several grid managers and linear solvers.
Here, we include classes related to parallel computations, time measurements and file I/O.
```cpp
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
```
In Dumux, the property system is used to specify classes and compile-time options to be used by the model.
For this, different properties are defined containing type definitions, values and methods.
All properties are declared in the file `properties.hh`.
```cpp
#include <dumux/common/properties.hh>
```
The following file contains the parameter class, which manages the definition and retrieval of input
parameters by a default value, the inputfile or the command line.
```cpp
#include <dumux/common/parameters.hh>
```
The file `dumuxmessage.hh` contains the class defining the start and end message of the simulation.
```cpp
#include <dumux/common/dumuxmessage.hh>
```
The following file contains the class, which defines the sequential linear solver backends.
```cpp
#include <dumux/linear/seqsolverbackend.hh>
```
Further we include the assembler, which assembles the linear systems for finite volume schemes (box-scheme, tpfa-approximation, mpfa-approximation).
```cpp
#include <dumux/assembly/fvassembler.hh>
```
The containing class in the following file defines the different differentiation methods used to compute the derivatives of the residual.
```cpp
#include <dumux/assembly/diffmethod.hh>
```
We need the following class to simplify the writing of dumux simulation data to VTK format.
```cpp
#include <dumux/io/vtkoutputmodule.hh>
```
The gridmanager constructs a grid from the information in the input or grid file. There is a specification for the different supported grid managers.
```cpp
#include <dumux/io/grid/gridmanager.hh>
```
### Beginning of the main function
```cpp
int main(int argc, char** argv) try
{
    using namespace Dumux;
```
Convenience aliases for the type tags of the two problems, which are defined in the individual problem files.
```cpp
    using OnePTypeTag = Properties::TTag::IncompressibleTest;
    using TracerTypeTag = Properties::TTag::TracerTestCC;
```
We initialize MPI. Finalization is done automatically on exit.
```cpp
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);
```
We print the dumux start message.
```cpp
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);
```
We parse the command line arguments.
```cpp
    Parameters::init(argc, argv);
```
### Create the grid
The `GridManager` class creates the grid from information given in the input file.
This can either be a grid file, or in the case of structured grids, by specifying the coordinates
of the corners of the grid and the number of cells to be used to discretize each spatial direction.
Here, we solve both the single-phase and the tracer problem on the same grid.
Hence, the grid is only created once using the grid type defined by the type tag of the 1p problem.
```cpp
    GridManager<GetPropType<OnePTypeTag, Properties::Grid>> gridManager;
    gridManager.init();
```
We compute on the leaf grid view.
```cpp
    const auto& leafGridView = gridManager.grid().leafGridView();
```
### Set-up and solving of the 1p problem
In the following section, we set up and solve the 1p problem. As the result of this problem, we obtain the pressure distribution in the domain.
#### Set-up
We create and initialize the finite volume grid geometry, the problem, the linear system, including the jacobian matrix, the residual and the solution vector and the gridvariables.
We need the finite volume geometry to build up the subcontrolvolumes (scv) and subcontrolvolume faces (scvf) for each element of the grid partition.
```cpp
    using GridGeometry = GetPropType<OnePTypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();
```
In the problem, we define the boundary and initial conditions.
```cpp
    using OnePProblem = GetPropType<OnePTypeTag, Properties::Problem>;
    auto problemOneP = std::make_shared<OnePProblem>(gridGeometry);
```
The jacobian matrix (`A`), the solution vector (`p`) and the residual (`r`) are parts of the linear system.
```cpp
    using JacobianMatrix = GetPropType<OnePTypeTag, Properties::JacobianMatrix>;
    using SolutionVector = GetPropType<OnePTypeTag, Properties::SolutionVector>;
    SolutionVector p(leafGridView.size(0));

    auto A = std::make_shared<JacobianMatrix>();
    auto r = std::make_shared<SolutionVector>();
```
The grid variables store variables (primary and secondary variables) on sub-control volumes and faces (volume and flux variables).
```cpp
    using OnePGridVariables = GetPropType<OnePTypeTag, Properties::GridVariables>;
    auto onePGridVariables = std::make_shared<OnePGridVariables>(problemOneP, gridGeometry);
    onePGridVariables->init(p);
```
#### Assembling the linear system
We create and inizialize the assembler.
```cpp
    using OnePAssembler = FVAssembler<OnePTypeTag, DiffMethod::analytic>;
    auto assemblerOneP = std::make_shared<OnePAssembler>(problemOneP, gridGeometry, onePGridVariables);
    assemblerOneP->setLinearSystem(A, r);
```
We assemble the local jacobian and the residual and stop the time needed, which is displayed in the terminal output, using the `assemblyTimer`. Further, we start the timer to evaluate the total time of the assembly, solving and updating.
```cpp
    Dune::Timer timer;
    Dune::Timer assemblyTimer; std::cout << "Assembling linear system ..." << std::flush;
    assemblerOneP->assembleJacobianAndResidual(p);
    assemblyTimer.stop(); std::cout << " took " << assemblyTimer.elapsed() << " seconds." << std::endl;
```
We want to solve `Ax = -r`.
```cpp
    (*r) *= -1.0;
```
#### Solution
We set the linear solver "UMFPack" as the linear solver. Afterwards we solve the linear system. The time needed to solve the system is recorded by the `solverTimer` and displayed in the terminal output.
```cpp
    using LinearSolver = UMFPackBackend;
    Dune::Timer solverTimer; std::cout << "Solving linear system ..." << std::flush;
    auto linearSolver = std::make_shared<LinearSolver>();
    linearSolver->solve(*A, p, *r);
    solverTimer.stop(); std::cout << " took " << solverTimer.elapsed() << " seconds." << std::endl;
```
#### Update and output
We update the grid variables with the new solution.
```cpp
    Dune::Timer updateTimer; std::cout << "Updating variables ..." << std::flush;
    onePGridVariables->update(p);
    updateTimer.elapsed(); std::cout << " took " << updateTimer.elapsed() << std::endl;

```
We initialize the vtkoutput. Each model has a predefined model specific output with relevant parameters for that model. We add the pressure data from the solution vector (`p`) and the permeability field as output data.
```cpp
    using GridView = GetPropType<OnePTypeTag, Properties::GridView>;
    Dune::VTKWriter<GridView> onepWriter(leafGridView);
    onepWriter.addCellData(p, "p");
    const auto& k = problemOneP->spatialParams().getKField();
    onepWriter.addCellData(k, "permeability");
    onepWriter.write("1p");
```
We stop the timer and display the total time of the simulation as well as the cumulative CPU time.
```cpp
    timer.stop();

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

```
### Computation of the volume fluxes
We use the results of the 1p problem to calculate the volume fluxes in the model domain.
```cpp
    using Scalar =  GetPropType<OnePTypeTag, Properties::Scalar>;
    std::vector<Scalar> volumeFlux(gridGeometry->numScvf(), 0.0);

    using FluxVariables =  GetPropType<OnePTypeTag, Properties::FluxVariables>;
    auto upwindTerm = [](const auto& volVars) { return volVars.mobility(0); };
```
We iterate over all elements
```cpp
    for (const auto& element : elements(leafGridView))
    {
```
Compute the element-local views on geometry, primary and secondary variables
as well as variables needed for flux computations
```cpp
        auto fvGeometry = localView(*gridGeometry);
        fvGeometry.bind(element);

        auto elemVolVars = localView(onePGridVariables->curGridVolVars());
        elemVolVars.bind(element, fvGeometry, p);

        auto elemFluxVars = localView(onePGridVariables->gridFluxVarsCache());
        elemFluxVars.bind(element, fvGeometry, elemVolVars);
```
We calculate the volume fluxes for all sub-control volume faces except for Neumann boundary faces
```cpp
        for (const auto& scvf : scvfs(fvGeometry))
        {
```
skip Neumann boundary faces
```cpp
            if (scvf.boundary() && problemOneP->boundaryTypes(element, scvf).hasNeumann())
                continue;
```
let the `FluxVariables` class do the flux computation.
```cpp
            FluxVariables fluxVars;
            fluxVars.init(*problemOneP, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
            volumeFlux[scvf.index()] = fluxVars.advectiveFlux(0, upwindTerm);
        }
    }

```
### Set-up and solving of the tracer problem
#### Set-up
Similar to the 1p problem, we first create and initialize the problem.
```cpp
    using TracerProblem = GetPropType<TracerTypeTag, Properties::Problem>;
    auto tracerProblem = std::make_shared<TracerProblem>(gridGeometry);
```
We use the volume fluxes calculated in the previous section as input for the tracer model.
```cpp
    tracerProblem->spatialParams().setVolumeFlux(volumeFlux);
```
We create and initialize the solution vector. As the tracer problem is transient, the initial solution defined in the problem is applied to the solution vector.
```cpp
    SolutionVector x(leafGridView.size(0));
    tracerProblem->applyInitialSolution(x);
    auto xOld = x;
```
We create and initialize the grid variables.
```cpp
    using GridVariables = GetPropType<TracerTypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(tracerProblem, gridGeometry);
    gridVariables->init(x);
```
We read in some time loop parameters from the input file. The parameter `tEnd` defines the duration of the simulation, dt the initial time step size and `maxDt` the maximal time step size.
```cpp
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
```
We instantiate the time loop.
```cpp
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
```
We create and inizialize the assembler with time loop for the instationary problem.
```cpp
    using TracerAssembler = FVAssembler<TracerTypeTag, DiffMethod::analytic, /*implicit=*/false>;
    auto assembler = std::make_shared<TracerAssembler>(tracerProblem, gridGeometry, gridVariables, timeLoop);
    assembler->setLinearSystem(A, r);
```
We initialize the vtk output module and add a velocity output.
```cpp
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, tracerProblem->name());
    using IOFields = GetPropType<TracerTypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);
    using VelocityOutput = GetPropType<TracerTypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    vtkWriter.write(0.0);

```
We define 10 check points in the time loop at which we will write the solution to vtk files.
```cpp
    timeLoop->setPeriodicCheckPoint(tEnd/10.0);
```
#### The time loop
We start the time loop and calculate a new time step as long as `tEnd` is not reached. In every single time step, the problem is assembled and solved.
```cpp
    timeLoop->start(); do
    {
```
First we define the old solution as the solution of the previous time step for storage evaluations.
```cpp
        assembler->setPreviousSolution(xOld);
```
Then the linear system is assembled.
```cpp
        Dune::Timer assembleTimer;
        assembler->assembleJacobianAndResidual(x);
        assembleTimer.stop();
```
We solve the linear system `A(xOld-xNew) = r`.
```cpp
        Dune::Timer solveTimer;
        SolutionVector xDelta(x);
        linearSolver->solve(*A, xDelta, *r);
        solveTimer.stop();
```
We calculate the actual solution and update it in the grid variables.
```cpp
        updateTimer.reset();
        x -= xDelta;
        gridVariables->update(x);
        updateTimer.stop();
```
We display the statistics of the actual time step.
```cpp
        const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
        std::cout << "Assemble/solve/update time: "
                  <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                  <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                  <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
                  <<  std::endl;
```
The new solution is defined as the old solution.
```cpp
        xOld = x;
        gridVariables->advanceTimeStep();
```
We advance the time loop to the next time step.
```cpp
        timeLoop->advanceTimeStep();
```
We write the Vtk output on check points.
```cpp
        if (timeLoop->isCheckPoint())
            vtkWriter.write(timeLoop->time());
```
We report the statistics of this time step.
```cpp
        timeLoop->reportTimeStep();
```
We set the time step size dt of the next time step.
```cpp
        timeLoop->setTimeStepSize(dt);

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

```
### Final Output
```cpp
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;

}
```
### Exception handling
In this part of the main file we catch and print possible exceptions that could
occur during the simulation.
```cpp
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
```

## Results

The 1p-model calculated a stationary pressure distribution. It is shown in the following figure:
![](./img/pressure.png)


The random permeability distribution generates the velocity profile shown in the left plot of the next figure. The image in the middle illustrates the tracer distribution after 2500s and the image on the right after 5000s.

| ![](img/velocityprofile.png)| ![](img/tracer_2500.png) | ![](img/tracer_5000.png)|
|:---:|:---:|:---:|
| velocity profile| tracer concentration after 2500s | tracer concentration after 5000s |
