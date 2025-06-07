
# 📌 compositionalIGFoam Repository Documentation

[![DOI](https://img.shields.io/badge/DOI-10.3390/en17143401-blue)](https://www.mdpi.com/1996-1073/17/14/3401)
[![DOI](https://img.shields.io/badge/DOI-InterPore2024-blue)](https://events.interpore.org/event/46/contributions/7161)

This code has been published in the `Journal of Energies` and presented in `InterPore 2024 Conference - 16th Annual Meeting & Conference Courses`.

# CompositionalIGFoam for OpenFOAM
About
-----

This code is developed on top of the `PorousMultiphaseFoam(PMF)` project for compositional flow modeling in porous media. It has been developed based on OpenFOAM-v10 and includes a new solver (compositionalIGFoam) for H2O / CO2 interactions by **[`Papi et al. (2024)`](https://www.mdpi.com/1996-1073/17/14/3401)**. The impesFoam solver of the PMF package has been modified to accommodate a compositional solver. The Henry's and Raoult's laws are implemented to account for the concentration of CO2 and H2O components in the liquid and gas phases by assuming an equilibrium condition. So no species transport equation is needed to be solved.
The developed model is the first of its kind in OpenFOAM to address CO2/water interactions in a multi-phase flow at the Darcy scale.

## Origin

This code is based on porousMultiphaseFoam(PMF), originally developed by Horgue et al. (2015).

PMF is an open-source toolbox dedicated to simulation of flow and transport processes in porous media.

It is based on the OpenFOAM environment and therefore benefits from its multiple feature such as pre-and post-processing tools or parallel efficiency for example.


# General Information

- This toolbox needs only a standard OpenFOAM installation and is compatible with the OpenFOAM-v10
  (see www.openfoam.org or www.openfoam.com).

- Please cite the related paper in the `doc` folder if you are using this
  toolbox.

- Read the COPYING_OPENFOAM file for information about OpenFOAM and this
  toolbox Copyrights.

# Installation instructions :


- First, source the OpenFOAM configuration file, i.e. (example for ubuntu
  version) :

> source /opt/openfoamv10/etc/bashrc

- then in the "compositionalIGFoam" directory, run :

> ./Allwmake -jX

  to install the package (with X the number of processors).

- Dynamic libraries are compiled in the standard OpenFOAM user directory :

> $FOAM_USER_LIBBIN

- The executable solver "compositionalIGFoam" is placed in the standard OpenFOAM user
  directory $FOAM_USER_APPBIN.

- Each tutorial directory contains "run" and "clean" files to test installation
  and validate the solver.

- A python script runTutorials.py can be used to test all components.

- To remove compilation and temporary files, run :

> ./Allwclean


---


## 🧪 Test Case Documentation

An individual test case can be found in the following folder:

- 📄 CO2 injection into a water-saturated core model
  ↳ [View Markdown](tutorials/compositionalIGFoam-tutorials/coreInjection/co2InjectionCore)

---

How to cite
---------------

If your are using the compositionalIGFoam, please cite as follows:

- [A New Straightforward Darcy-Scale Compositional Solver in OpenFOAM for CO2/Water Mutual Solubility in CO2 Storage Processes in Aquifers](https://www.mdpi.com/1996-1073/17/14/3401)
  A. Papi, A. Jahanbakhsh, M.M. Maroto-Valer;
  *MDPI Energies  17(14), 3401*, **2024**.
