---
title: 'obscura: A modular C++ tool and library for the direct detection of (sub-GeV) dark matter particles via nuclear and electron recoils'
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

<!-- Notes:
  - link to documentation
  - already used in damascus-sun
  - Your paper should include:
  - A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline.
  - Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it. 
  - Emphasize polymorphisms and modular structure
-->

# Summary

The observation of a large number of gravitational anomalies on astrophysical and cosmological scales have convinced us that the majority of matter in the Universe is invisible [Bertone:2004pz;@Bertone:2018krk].
This *dark matter* (DM) must be fundamentally different from the visible matter we can describe using the Standard Model of Particle Physics (SM).
Its only established property is that it interacts gravitationally and indeed dominates the gravitational potential of galaxies and galaxy clusters.
One of the leading hypothesis is that DM is made up of one or more new particles and that galaxies such as our Milky Way are embedded in gigantic haloes of these as of yet undetected particles.
Our planet would at any moment be penetrated by a stream of these particles without much of an effect.
If these dark particles interact with nuclei and/or electrons via some new force besides gravity, they would on occasion collide with a terrestrial atom.
*Direct detection experiments* search for these kind of interactions and aim to observe DM events within a detector caused by an interaction with target nuclei [@Goodman:1984dc;@Drukier:1986tm;@Wasserman:1986hh] or electrons [@Kopp:2009et;@Essig:2011nj].
These experiments are typically placed deep underground to shield them from possible backgrounds e.g. due to cosmic rays.

In order to interpret the outcome of direct detection experiments, we need to make predictions for the expected events caused by the incoming DM particles.
This cannot be done without making a number of assumptions about the possible particle attributes of DM (e.g. mass and interaction strength) and the properties of the galactic DM halo (e.g. the local DM density and their energy distribution) [@Lewin:1995rx;@DelNobile:2021icc].

`obscura` is a tool to make quantitative predictions for direct DM searches, analyse experimental data, and derive e.g. exclusion limits as seen in Figure \autoref{fig:constraints}.
It can be used to compute the expected event rates in terrestrial detectors looking for rare interactions between the DM and nuclei or electrons.
There are many different experimental techniques and targets proposed and applied for direct detection experiments [@Griffin:2019mvc].
Additionally, due to our ignorance about the particle physics of DM there exists a plethora of viable assumptions and models.
This is reflected by the modular, polymorphic structure of all modules of the `obscura` library which allows to easily extend `obscura`'s functionality to the users' new idea on what the DM particles could be like or behave, or to a new detection technology.
This generic structure also allows applications of (a subset of) the `obscura` classes in a variety of DM research projects even beyond direct detection.

For more details on `obscura` and its implementation in C++, we refer to the [documentation](https://obscura.readthedocs.io)[^1].

[^1]: The latest version of the documentation can be found under [https://obscura.readthedocs.io](https://obscura.readthedocs.io).

![Excluded regions (95% confidence level) of the DM parameter space given by the $(m_\mathrm{DM},\sigma_i)$ plane, where $m_\mathrm{DM}$ is the assumed DM mass and $\sigma_i$ is the interaction cross section with target $i$. \label{fig:constraints}](obscura_DD_Constraints.png){ width=80% } 

# Statement of need

While there exists a few tools to compute DM induced nuclear recoil spectra, such as DDCalc [@GAMBITDarkMatterWorkgroup:2017fax;@GAMBIT:2018eea] or WimPyDD [@Jeong:2021bpl], `obscura` is not limited to nuclear targets.
Instead the main focus lies on sub-GeV DM searches probing electron recoils which typically requires methods from atomic and condensed matter physics, see e.g. [@Essig:2015cda;@Catena:2019gfa;Catena:2021qsr].

 

# The modular structure of direct detection computations

Making predictions for direct detection experiments involves methods and results from statistics, astrophysics, particle physics, nuclear and atomic physics, and condensed matter physics.
For each of these fields, we need to make choices and assumptions which affect our interpretation of DM searches.
It is our ambition for the `obscura` code that the basic functionality does not rely on specific choices and that no particular assumption is hard-coded.
Instead the basic code's setup is polymorphic and written in terms of generic base classes widely agnostic to specific assumptions.
Of course a number of standard ideas are implemented as derived classes.
The user is free to include their own ideas as their own derived classes.
The classes are described in more detail in the [documentation](https://obscura.readthedocs.io).

As an example, let us look at the energy spectrum of DM induced ionization events, as derived in [@Essig:2015cda].
\begin{equation}
 \frac{\mathrm{d} R_\mathrm{ion}}{\mathrm{d} E_e} = N_T \frac{\rho_\chi}{m_\chi}\sum_{n,\ell} \int \mathrm{d}q^2\int \mathrm{d}v\; v f_\chi(v) \frac{1}{4E_e}\frac{\mathrm{d}\sigma_e}{\mathrm{d}q^2} \left|f_\mathrm{ion}^{n\ell}(q,E_e)\right|^2\, .
\end{equation}


![The class structure of `obscura`.\label{fig:flowchart}](FlowChart.png){ width=70% } 

Research software can implement `obscura`'s classes following the dependencies indicated by the flow chart of figure \autoref{fig:flowchart}.
Alternatively, they can also be used in all kinds of functions depending on the research project's aim.
As an example, the `DaMaSCUS-SUN` code uses `obscura` in the context of Monte Carlo simulations [@Emken:2021lgc;@Emken2021].

# Acknowledgements
The author thanks Radovan Bast for  valuable  discussions and support regarding research software engineering.
The author was supported by the Knut & Alice Wallenberg Foundation.

# References