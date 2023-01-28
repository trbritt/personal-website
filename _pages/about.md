---
layout: about
title: Tristan Britt
permalink: /
description:  PhD Student at the <a href="https://cpm.research.mcgill.ca">Centre for Physics of Materials</a>, McGill University.
profile:
  align: right
  image: headshot.jpg
  address: >
    <p>801 Rue Sherbrooke Ouest</p>
    <p>025 Otto Maass Chemistry Building</p>
    <p>Montreal, QC, Canada H3A 0B8</p>

education: true
news: true  # includes a list of news items
selected_papers: true # includes a list of papers marked as "selected={true}"
social: true  # includes social icons at the bottom of the page
---
### biography

Having worked extensively in the managerial teams of large scale organisations in the time before completing my undergraduate degree, I’ve developed team leadership skills from the beginning of my professional career. My industry research includes simulation and optimisation of custom radiofrequency (RF) waveguides for use in cryogenic axion dark matter experiments (ADMX) at the Korea Advanced Institute of Science and Technology (KAIST). In addition, I have also saw through the co-design and optimisation of the bending magnet currently in use in the LEReC beamline in the electron Relativistic Heavy Ion Collider (eRHIC) at Brookhaven National Lab (BNL). These projects have given me the experience and preparedness to handle high quality real-time research.

My current focus is on examining the unique properties of 2D materials as revealed through momentum and time resolved ultrafast electron diffraction spectroscopy at McGill University, under the advisorship of Professor Bradley J Siwick, including both experiment and numerical simulations using time dependent perturbative density functional theory (TDPDFT). Experiments under progress include the observation of chiral phonons in monolayer molybdenum disulphide (MoS$$_2$$), as well as the thermoelectric properties of tin selenide (SnSe) approaching the PNMA phase transition at 600K.

---

### education

###### **McGill University** - Montreal, QC, CA | *Experimental Condensed Matter PhD* • In Progress 
The thesis is titled 'Ultrafast phonon dynamics at the 2D limit: a comparison of ultrafast electron scattering to *ab-initio* calculations'. 

###### **McGill University** - Montreal, QC, CA | *M.Sc.* • June 2019 - August 2021
My master's degreed focused on instrument design and optimisation with high-performance computations. Firstly, simulations were done to verify the results of a time-resolved electron energy loss spectroscopy (trEELS) setup that is to be included in upcoming upgrades to the UEDS beamline at McGill. These simulations are written in `C`, and rely on OpenMP and MPI for parallel computation. The rest of the work focused on the design of a high performance TM$$_{010}$$ radiofrequency compression cavity, used to temporally focus the Coulomb expanded ultrafast electron bunch used in the scattering experiments. 
###### **Indiana University** - Bloomington, IN, USA | *B.Sc.* with Summa Cum Laude • May 2019
I received B.S. degrees in both Physics and Applied Mathematics. Research conducted while pursuing these degrees included neutron interferometry experiments, RF cavity design, cryogenic magnet optimisation, and experimental fiber bundle computations.
<div class="row">

  <div class="columns download">
     <p>
        <a href="assets/pdf/Tristan_Britt_Resume_Professional.pdf" class="button"><i class="fa fa-download"></i>  Download Resume</a>
     </p>
  </div>

</div>

---

### software

I have proficiency in the following software packages, suites, or applications.

##### `COMSOL`
A finite-element multiphysics solver. Among the many packages provided by this suite, I am proficient in primarily the radiofrequency (RF) electronics, Particle Tracing, Optimisation, Heat Transfer, and E&M solvers.

##### `OPERA`
A finite-element electromagnetic field solver. This was used in the [published](/publications/) work from BNL where a cryogenic bending magnet was designed. It is currently installed and in use at the electron cooling beamline of the electron Relativistic Heavy Ion Collider (eRHIC).

##### `ROXIE`
A finite-element beam magnet field solver. Designed by CERN, this `C` suite of programmes allows for highly efficient computation of various multipole magnetic fields, with primary applications being the use of such complicated fields for beam manipulation inside syncotrons and other light or electron sources. This was used to compliment the `OPERA` simulations discussed above.

##### `CERN ROOT`
A data analysis suite, used to examine large datasets efficiently. Originally created by CERN to process complex collision data, this suite has functionality for multivariate fitting, visualisation, and more, for 'big data' sets.

##### Autodesk's `AutoCAD` Suite + `Solidworks` + `CREO`
These CAD design softwares allow for 3D design optimisation in otherwise untractable physical geometries. I have used Solidworks and `AutoCAD` to design and optimise the RF cavity designed for the Siwick Lab, as well as the cavity used in the ADMX experiment in Korea. 

##### `ANSYS`
This multiphysics simulation engine in my use has been directly coupled to the CAD work previously described. The designs are made in CAD and then brought into `ANSYS` for processing. An example is the computation of all electromagnetic standing modes inside the RF cavity, to identify which geometry provides the largest electric field strength at the volume of interest.

##### `HDF` Data format
I am proficient in HDF5 data, as the experimental data of the instrument at McGill is stored as such. Owing to its compact memory architecture, I have used HDF5 as the default data format for all simulations I've performed over the last few years. The density functional theory calculations require extremely expensive computation on big data and as such, HDF5 is a well-presented solution to the corresponding storage of the data. I am proficient in using HDF5 with Python, and have made my own modules and interfaces for parallel I/O to HDF5 with `Fortran`.

---

### programming experience

My programmes have utilised the following techniques/approaches.
* Density Functional Theory
* Perturbative Density Functional Theory
* Object-oriented C++ computational electromagnetics simulations
* Finite Element Method
* Integral Equation Method
* Finite Differnce Time Domain (FDTD)
* High Frequency Methods
* RF Design and Analysis
* OpenMP
* CUDA-accelerated implementations

I am proficient in the following programming languages.

##### Fortran - `F77`, `F90`, modern Fortran (`F08`)
The expensive nature of the computations done in density functional theory require an arrayoptimised language that is easily parallelisable. As such, I have become a programmer whose primary language is optimised Fortran. The standard density functional theory simulations done by the open-source suite `Quantum Espresso` are all in modern fortran. I have expanded this suite for my own custom calculations, the details of which are given in [this project](/projects/2_project).

##### `Python`
Python is the primary lagnuage used to process the experimental scattering data collected for my work at McGill University. I have taught classes on Python at my time at McGill, ranging from basic language syntax to nonlinear curve fitting, to data visualisation, to basic machine learning. I am certified on Linkedin as proficient in `Python`.

##### `C` flavour languages
My simulations of the magnets in my time at BNL, my electron bunch simulations in my masters, and other academic projects have made my skills in `C` flavour languages workable. I am more familiar with `C++` than C, owing to its convenient structure-oriented programming approaches which have shown to be more convenient for the projects I've used `C` for. I am certified on Linkedin as proficient in these languages.

##### Markdown and $$\LaTeX$$
These two typesetting languages have been the default method through which I have written all my academic published work, my class assignments during my previous degrees, and this website.

##### Javascript
I have recently aimed to be more creative with my coding, and to get a perspective on a programming language new and different to that I've developed during my PhD so far. Javascript has been my entrance into generative coding, and I will start posting more projects about it soon! (see my projects page [here](/projects/))

##### Matlab
A great programming language to teach new students programming, I've used `MATLAB` a few times before in classes to introduce students to an interpreted language with rgeat flexibility and low barrier-of-entry.

---
