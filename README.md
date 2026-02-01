## DESCRIPTION 
- This FORTRAN code simulates the diffusion of Iron inside the INTRA CLUSTER MEDIUM; The workflow is the following:
1.  Obtain the gas density profile for the ICM in hydrostatic equilibrium;
2.  Use the profiles obtained to feed the diffusion equation and see the evolution in 5Gyrs

The following modules are implemented: 
* algoritmi.mod, containing the algorithms used to resolve the equations
  The algorithms used are *euler* and *ftcs*

* equazioni.mod, containing the equations to be solved
  The equations to be solved are:
    * Hydrostatic Equilibrium
      * *eq_idro1* perfect and isothermal gas;
      * *eq_idro2* gas with a radial temperature profile as the one from literature;
      * *eq_idro3* previous case + turbulence;
      * *diff* diffusion equation for iron in the ICM;

* funzioni.mod, cointains:
  * *NFW* the Navarro-Frank-White DM-density profile;
  * *massa* the computation of the mass up to a given radius of a density profile;
  * *t_r* the temperature profile from (name)

## GRID
- A staggered mesh is implemented.
  The width, spacing and deformation factor can be set;

## REQUIREMENTS
- Install the *gfortran* compiler

## INSTRUCTIONS
- run "make" to compile
- Run the simulation changing the simulation paramaters as you want
  Simulation parameters:
