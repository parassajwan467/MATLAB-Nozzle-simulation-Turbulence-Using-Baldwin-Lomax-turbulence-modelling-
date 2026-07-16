# MATLAB-Nozzle-simulation-Turbulence-Using-Baldwin-Lomax-turbulence-modelling-


A finite-volume MATLAB framework developed for the numerical simulation of compressible flow through axisymmetric convergent-divergent rocket nozzles.

The solver has been developed progressively from the compressible Euler equations to the viscous compressible Navier-Stokes equations with algebraic turbulence closure using the Baldwin-Lomax turbulence model. The numerical methodology follows classical finite-volume techniques and draws inspiration from the computational strategy presented by Mehta and Jayachandran for viscous rocket nozzle flows.

---

## Overview

The present framework is intended for the simulation of internal compressible nozzle flows encountered in rocket propulsion applications. The code architecture is modular and designed to facilitate future extension towards multiphase gas-particle flow simulations.

Current capabilities include:

- Axisymmetric compressible flow formulation
- Structured body-fitted finite volume discretization
- Inviscid and viscous flow simulations
- Algebraic turbulence modelling
- Shock-capturing artificial dissipation
- Local time stepping for accelerated convergence

---

## Governing Equations

The solver employs the conservative form of the axisymmetric compressible governing equations:

### Inviscid Formulation
- Continuity equation
- Axial momentum equation
- Radial momentum equation
- Total energy equation

### Viscous Formulation
The viscous extension includes:

- Molecular viscous stresses
- Fourier heat conduction
- Axisymmetric viscous source terms
- Effective viscosity treatment for turbulent flows

The governing equations are solved in integral finite-volume form over quadrilateral control volumes.

---

## Numerical Methodology

### Spatial Discretization
- Cell-centered Finite Volume Method (FVM)
- Structured body-fitted mesh
- Face-based flux evaluation

### Time Integration
- Explicit multi-stage Runge-Kutta time marching
- Local CFL-controlled pseudo time stepping

### Inviscid Flux Treatment
- Central differencing formulation
- Jameson second-order and fourth-order artificial dissipation
- Pressure-based shock sensor for adaptive stabilization

### Viscous Flux Treatment
- Gradient reconstruction at cell centres and faces
- Viscous stress tensor evaluation
- Heat flux computation using Fourier's law

---

## Physical Models

### Molecular Viscosity
Dynamic viscosity is evaluated using Sutherland's law:

- Reference viscosity:
  μ_ref = 1.716 × 10⁻⁵ kg/(m·s)
- Reference temperature:
  T_ref = 273.15 K
- Sutherland constant:
  S = 111 K

### Turbulence Closure
The solver employs the Baldwin-Lomax algebraic turbulence model:

- Inner layer eddy viscosity formulation
- Outer layer wake formulation
- Van Driest damping function
- Effective viscosity computation

Model constants include:

| Parameter | Value |
|----------|------:|
| κ (Von Kármán constant) | 0.40 |
| A⁺ | 26 |
| K | 0.0168 |
| Ccp | 1.60 |
| Cwk | 0.25 |

---

## Boundary Conditions

### Inlet
- Supersonic inflow initialization using analytical isentropic nozzle solution

### Exit
- Supersonic outflow with zero-gradient extrapolation

### Axis of Symmetry
- Symmetry boundary condition
- Zero radial velocity
- Zero normal gradients

### Wall Boundary
- No-slip wall condition
- Adiabatic wall assumption

---

## Repository Structure

```text
MATLAB-Nozzle-Simulation
│
├── Initialization
│   ├── InitializeFlow.m
│   ├── PrimitiveFromW.m
│   └── IsentropicInitialization.m
│
├── PreProcessing
│   ├── GeometryGeneration.m
│   ├── GridGeneration.m
│   ├── FaceMetrics.m
│   └── MeshUtilities.m
│
├── InviscidSolver
│   ├── BuildResidual.m
│   ├── ComputePressureSensor.m
│   ├── ComputeArtificialDissipation.m
│   ├── ApplyBoundaryConditions.m
│   └── TimeIntegration.m
│
├── TurbulenceModelling
│   ├── ComputeGradients.m
│   ├── ComputeViscousProperties.m
│   ├── ComputeViscousFlux.m
│   ├── BuildViscousResidual.m
│   └── ComputeBaldwinLomax.m
│
├── PostProcessing
│   ├── PlotMachContours.m
│   ├── PlotPressureContours.m
│   ├── PlotResidualHistory.m
│   ├── PlotWallShearStress.m
│   └── PlotTurbulentViscosity.m
│
├── main.m
│
└── README.md
```

---

## Verification and Validation

The solver has been verified and assessed using:

- One-dimensional isentropic nozzle theory
- Area-Mach number relation
- Conservation of mass flow rate
- Expected boundary layer development
- Wall shear stress behaviour
- Residual convergence characteristics

---

## Research Motivation

The long-term objective of this framework is the numerical simulation of turbulent gas-particle flows inside solid rocket motor nozzles using an Eulerian-Eulerian formulation.

Future developments include:

- Gas-particle momentum coupling
- Particle energy exchange
- Two-phase compressible flow modelling
- Particle loading effects
- Combustion source terms
- Advanced turbulence closures
- Higher-order spatial reconstruction schemes

---

## References

1. Mehta, R. C., and Jayachandran, R.

   *A Fast Algorithm to Solve Viscous Two-Phase Flow in an Axisymmetric Rocket Nozzle.*

2. Baldwin, B. S., and Lomax, H. (1978).

   *Thin Layer Approximation and Algebraic Model for Separated Turbulent Flows.*

3. Anderson, J. D.

   *Computational Fluid Dynamics: The Basics with Applications.*

4. Hirsch, C.

   *Numerical Computation of Internal and External Flows.*

---

## Author

**Paras**  
Undergraduate Student, Department of Aerospace Engineering  
Indian Institute of Space Science and Technology (IIST)

---

## Disclaimer

This code has been developed primarily for academic research and educational purposes. The implementation continues to evolve as additional physical models and numerical capabilities are incorporated.
