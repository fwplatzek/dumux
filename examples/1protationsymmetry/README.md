<!-- Important: This file has been automatically generated by generate_example_docs.py. Do not edit this file directly! -->

# Rotation-symmetric pressure distribution

__In this example, you will learn how to__

* solve a rotation-symmetric problem one-dimensionally
* perform a convergence test against an analytical solution
* apply the `Rotational Extrusion` filters in [ParaView](https://www.paraview.org/) for a two-dimensional visualization of the one-dimensional results


__Result__. With the `Rotational Extrusion` and the `Warp By Scalar` filter in [ParaView](https://www.paraview.org/),
the pressure distribution of this example looks as shown in the following picture:

<figure>
    <center>
        <img src="img/result.png" alt="Rotation-symmetric pressure distribution" width="60%"/>
        <figcaption> <b> Fig.1 </b> - Result TODO.</figcaption>
    </center>
</figure>


__Table of contents__. This description is structured as follows:

[[_TOC_]]


## Problem setup

We consider a single-phase problem that leads to a rotation-symmetric pressure distribution.
The following figure illustrates the setup:

<figure>
    <center>
        <img src="img/setup.svg" alt="Rotation-symmetric setup" width="60%"/>
        <figcaption> <b> Fig.1 </b> - Setup for the rotation-symmetric problem. The pressure boundary conditions are shown by the colored lines and the simulation domain is depicted in grey.</figcaption>
    </center>
</figure>

This could, for example, represent a cross section of an injection/extraction well in a homogeneous
and isotropic porous medium, where the well with radius $`r_1`$ is cut out and the
injection/extraction pressure $`p_1`$ is prescribed as a Dirichlet boundary condition. At the outer
radius $`r_2`$, we set the pressure $`p_2`$. In the polar coordinates $`r`$ and $`\varphi`$, the
solution to this problem is independent of the angular coordinate $`\varphi`$ and can be reduced to
a one-dimensional problem in the radial coordinate $`r`$. Therefore, in this example, we want to
solve the problem on a one-dimensional computational domain as illustrated by the orange line in
the above figure.

## Mathematical model

In this example we are using the single-phase model of DuMuX, which considers Darcy's law to relate
the Darcy velocity $`\textbf v`$ to gradients of the pressure $`p`$. For an isotropic medium and
neglecting gravitational forces, this can be written as:

```math
\textbf v = - \frac{k}{\mu} \text{grad} p.
```

Here, $`k`$ is the permeability of the porous medium and $`\varrho`$ and $`\mu`$ are the density
and the dynamic viscosity of the fluid. In the model, the mass balance equation for the fluid
phase is solved:

```math
\phi \frac{\partial \varrho}{\partial t} + \text{div} \left( \varrho \textbf v \right) = 0,
```

where $`\phi`$ is the porosity of the porous medium. Let us now introduce the transformation
$`(x, y)^T = \Phi ( r, \varphi )`$ from polar into cartesian coordinates (see e.g.
[wikipedia.org](https://en.wikipedia.org/wiki/Polar_coordinate_system#Converting_between_polar_and_Cartesian_coordinates)),
and denote with

```math
\tilde{p} \left( r, \varphi \right)
    = \tilde{p} \left( r \right)
    = p \left( \Phi^{-1}(x, y) \right)
```

and

```math
\tilde{\mathbf{v}} \left( r, \varphi \right)
    = \tilde{v}_r \left( r \right)
    = \mathbf{v} \left( \Phi^{-1}(x, y) \right)
```

the pressure and velocity distributions expressed in polar coordinates. The first identity
in the two above equations originates from the rotational symmetry of the problem and the
resulting independence of pressure and velocity on $`\varphi`$ in polar coordinates. Thus, in
polar coordinates we can write the mass balance equation as:

```math
\phi \frac{\partial \varrho}{\partial t}
   - \frac{\partial}{\partial r} \left( \frac{k}{\mu} \frac{\partial \tilde{p}}{\partial r} \right)
   = 0.
```

## Discretization

We employ a finite-volume scheme to spatially discretize the mass balance equation shown above.
The discrete equation describing mass conservation inside a control volume $`K`$ is obtained
by integration and reads:

```math
    | K | \left( \phi \, \partial \varrho / \partial t \right)_K
    + \sum_{\sigma \in \mathcal{S}_K} | \sigma | \left( \varrho v_r \right)_\sigma
    = 0,
```

where $`\sigma`$ are the faces of the control volume such that
$`\bigcup_{\sigma \in \mathcal{S}_K} \sigma \equiv \partial K`$ and where the notation $`( \cdot )_K`$
and $`( \cdot )_\sigma`$ was used to denote quantities evaluated for the control volume $`K`$ or a
face $`\sigma`$, respectively. The volume of the control volume is denoted with $`| K |`$ and
$`| \sigma |`$ is the area of a face.

Integration over polar coordinates requires taking into account the Jacobian determinant of the
coordinate transformation from polar to cartesian coordinates (see e.g.
[wikipedia.org](https://en.wikipedia.org/wiki/Polar_coordinate_system#Generalization)).
Let us discretize the domain by the intervals
$`K_i = (i\Delta r, (i+1)\Delta r) \times (0, 2 \Pi)`$,
$`i \in \{1, \dots, N \}`$, as control volumes.
As a result, their volumes are

```math
| K_i | = \Pi \left( ((i+1)*\Delta r)^2 - (i*\Delta r)^2 \right)
```

and the area of a face $`\sigma \in \mathcal{S}_{K_i}`$ is

```math
| \sigma | = 2 \Pi r_\sigma,
```

where $`r_\sigma`$ is the radius at which the face is defined.

# Implementation & Post processing

## Part 1: Rotation-symmetric one-phase flow simulation setup

| [:arrow_right: Click to continue with part 1 of the documentation](doc/problem.md) |
|---:|


## Part 2: Main program flow

| [:arrow_right: Click to continue with part 2 of the documentation](doc/main.md) |
|---:|


## Part 3: Post-processing with ParaView

| [:arrow_right: Click to continue with part 3 of the documentation](doc/paraview.md) |
|---:|