---
title: 'obscura: A modular C++ tool and library for the direct detection of (sub-GeV) dark matter via nuclear and electron recoils'
tags: 
  - c++
  - astroparticle physics
  - dark matter
  - direct detection
authors: 
  - name: Timon Emken
    orcid: 0000-0002-4251-2229
    affiliation: 1
affiliations: 
  - name: The Oskar Klein Centre, Department of Physics, Stockholm University, AlbaNova, SE-10691 Stockholm, Sweden 
    index: 1
date: 15 July 2021
bibliography: paper.bib
---

# Summary

The observation of a large number of gravitational anomalies on astrophysical and cosmological scales have convinced us that the majority of matter in the Universe is invisible [@Bertone:2004pz;@Bertone:2018krk].
This *dark matter* (DM) must be fundamentally different from the known matter we can describe using the Standard Model of Particle Physics (SM).
Its only established property is that it interacts gravitationally and indeed dominates the gravitational potential of galaxies and galaxy clusters.
One of the leading hypothesis is that DM is made up of one or more new particles and that galaxies such as our Milky Way are embedded in gigantic haloes of these as of yet undetected particles.
Our planet would at any moment be penetrated by a stream of these particles without much of an effect.
If these dark particles interact with nuclei and/or electrons via some new force besides gravity, they would on occasion collide with a terrestrial particle.
*Direct detection experiments* search for these kind of interactions and aim to observe DM events within a detector caused by an interaction with target nuclei [@Goodman:1984dc;@Drukier:1986tm;@Wasserman:1986hh] or electrons [@Kopp:2009et;@Essig:2011nj].
These experiments are typically placed deep underground to shield them from possible backgrounds e.g. due to cosmic rays.

In order to interpret the outcome of direct detection experiments, we need to make predictions for the expected events caused by the incoming DM particles.
In all cases, this requires making a number of assumptions about the possible particle attributes of DM (e.g. mass and interaction strength) and the properties of the galactic DM halo (e.g. the local DM density and their energy distribution) [@Lewin:1995rx;@DelNobile:2021icc].

`obscura` is a tool to make quantitative predictions for direct DM searches, analyse experimental data, and derive e.g. exclusion limits as seen in \autoref{fig:constraints}.
`obscura` can e.g. be used to compute the expected event rates in terrestrial detectors looking for rare interactions between the DM and nuclei or electrons.
There are many different experimental techniques and targets proposed and applied for direct detection experiments [@Griffin:2019mvc].
Additionally, due to our ignorance about the particle physics of DM there exists a plethora of viable assumptions and models.
The vast variety of viable assumptions is reflected by the modular, polymorphic structure of all modules of the `obscura` library which allows to easily extend `obscura`'s functionality to the users' new idea on the fundamental nature of DM particles, or on a new detection technology.
For example, the library can handle any kind of DM particles of any mass, provided that the scattering is well-described by non-relativistic dynamics, and that the differential (nucleus and/or electron) scattering cross sections depend only on the momentum transfer, the relative speed between DM and target, and at most one additional dynamic parameter such as the center-of-mass energy or the local temperature of the target.
Furthermore, a generic structure also allows applications of (a subset of) the `obscura` classes in a variety of DM research projects even beyond an the context of direct detection, e.g. to compute DM capture rates in the Sun [@Emken:2021lgc].

For more details on `obscura` and its implementation in C++, we refer to the [documentation](https://obscura.readthedocs.io)[^1].

[^1]: The latest version of the documentation can be found under [https://obscura.readthedocs.io](https://obscura.readthedocs.io).

![Excluded regions (90% confidence level) of the DM parameter space given by the $(m_\mathrm{DM},\sigma_i)$ plane, where $m_\mathrm{DM}$ is the assumed DM mass and $\sigma_i$ is the interaction cross section with target $i$. For comparison, the dashed lines denote the official results published by the experimental collaborations. Some of the `obscura` results are conservative due to a simplified analysis.\label{fig:constraints}](obscura_DD_Constraints.png){ width=85% } 



# The modular structure of direct detection computations

Making predictions and performing analyses for direct detection experiments involves methods and results from statistics, astrophysics, particle physics, nuclear and atomic physics, and condensed matter physics.
For each of these fields, we need to make choices and assumptions which will affect our interpretation of DM searches.

As an example, let us look at the energy spectrum of DM induced ionization events, as derived in [@Essig:2015cda].
\begin{equation}
 \frac{\mathrm{d} R_\mathrm{ion}}{\mathrm{d} E_e} = N_T \frac{\rho_\chi}{m_\mathrm{DM}}\sum_{n,\ell} \int \mathrm{d}q^2\int \mathrm{d}v\; v f_\chi(v) \frac{1}{4E_e}\frac{\mathrm{d}\sigma_e}{\mathrm{d}q^2} \left|f_\mathrm{ion}^{n\ell}(q,E_e)\right|^2\, .
\end{equation}
The DM mass $m_\mathrm{DM}$ and the differential DM-electron scattering cross section $\frac{\mathrm{d}\sigma_e}{\mathrm{d}q^2}$ are defined by the assumed particle physics of the hypothetical DM particle the experiment is probing.
The velocity distribution $f_\chi(v)$ and the local DM energy density $\rho_\chi$ are important inputs from astrophysics and cosmology.
Lastly, the ionization form factor $f_\mathrm{ion}^{n\ell}(q,E_e)$ encapsulates the atomic physics of the electronic bound states and describes the probability of an electron with quantum numbers $(n\ell)$ to get ionized by an incoming DM particle.
As we can see, the evaluation of this expression for the electron recoil spectrum is highly modular combining inputs from various fields of research.
This modularity should be reflected in the structure of corresponding research software.

It is our ambition for the `obscura` code that the basic functionality does not rely on specific choices and that no particular assumption is hard-coded.
Instead the basic code's setup is polymorphic and written in terms of generic base classes widely agnostic to specific assumptions.
Of course, a number of standard ideas and models are implemented as derived classes, which also illustrate the usage of the base classes.
The classes are described in more detail in the [documentation](https://obscura.readthedocs.io).

![The class structure of `obscura`.\label{fig:flowchart}](FlowChart.png){ width=70% } 

External research software can use `obscura` by implementing its classes following the dependencies indicated by the flow chart of \autoref{fig:flowchart} and computing standard quantities in the context of direct detection of dark matter.
In addition, these classes are meant to be general-purpose and can be applied in other contexts depending on the research project's main objectives.
It is also possible to exploit the polymorphic structure and extend its functionality by creating new derived classes based on the users' own ideas.
As a final benefit of the polymorphic structure, any research software that is formulated entirely in terms of the abstract base classes can later on be used with any derived classes and allows analyses and research for a broad range of alternative assumptions without changing the core of the scientific code. 
As an example, the `DaMaSCUS-SUN` code uses `obscura` in the context of Monte Carlo simulations [@Emken:2021lgc;@Emken2021].

# Statement of need

For the interpretation of past and future direct searches for DM particles, it is important to be able to provide accurate predictions for event rates and spectra under a variety of possible and viable assumptions in a computationally efficient way.
While there exists a few tools to compute DM induced nuclear recoil spectra, such as DDCalc [@GAMBITDarkMatterWorkgroup:2017fax;@GAMBIT:2018eea] or WimPyDD [@Jeong:2021bpl], `obscura` is not limited to nuclear targets.
Instead its main focus lies on sub-GeV DM searches probing electron recoils which typically requires methods from atomic and condensed matter physics, see e.g. [@Essig:2015cda;@Catena:2019gfa;@Catena:2021qsr].
In the context of sub-GeV DM searches, new ideas such as target materials or detection techniques are being proposed regularly, and the theoretical modelling of these are getting improved continuosly, see e.g [@Griffin:2021znd].
At the same time, currently running experiments continue to publish their results and analyses, setting increasingly strict bounds on the DM parameter space.
In such a dynamic field, `obscura` can be an invaluable tool due to its high level of adaptability and facilitate and accelerate the development of new, reliable research software for the preparation of a DM discovery in the hopefully near future.


# Acknowledgements
The author thanks Radovan Bast for  valuable  discussions and support regarding research software engineering.
The author was supported by the Knut & Alice Wallenberg Foundation (PI, Jan Conrad).

# References