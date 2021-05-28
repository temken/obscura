[![Build Status](https://github.com/temken/obscura/workflows/Build%20Status/badge.svg)](https://github.com/temken/obscura/actions)
[![codecov](https://codecov.io/gh/temken/obscura/branch/master/graph/badge.svg)](https://codecov.io/gh/temken/obscura)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

# obscura - Direct detection of dark matter with nucleus and electron recoil experiments

[![DOI](https://zenodo.org/badge/250298432.svg)](https://zenodo.org/badge/latestdoi/250298432)

A C++ tool and library for dark matter direct detection computations for both nuclear and electron recoil experiments.

<img src="paper/FlowChart.png" width="500">

## Dependencies

- [libphysica](https://github.com/temken/libphysica)

## Included experiments

<img src="paper/obscura_DD_Constraints.png" width="500">

The following nuclear and electron recoil direct detection experiments are implemented.

### Nuclear recoil experiments

#### CRESST-II

- *Results on light dark matter particles with a low-threshold CRESST-II detector*  
CRESST Collaboration (G. Angloher et al.)  
[![Eur.Phys.J. C76 (2016) no.1, 25](https://img.shields.io/badge/Eur.Phys.J.-C76(2016)no.1,25-255773.svg)](https://link.springer.com/article/10.1140/epjc/s10052-016-3877-3)
[![[arXiv:1509.01515]](https://img.shields.io/badge/arXiv-1509.01515-B31B1B.svg)](https://arxiv.org/abs/1509.01515)

- *Description of CRESST-II data*  
CRESST Collaboration (G. Angloher et al.)  
[![[arXiv:1701.08157]](https://img.shields.io/badge/arXiv-1701.08157-B31B1B.svg)](https://arxiv.org/abs/1701.08157)

#### CRESST-III

- *First results on low-mass dark matter from the CRESST-III experiment*  
CRESST Collaboration (F. Petricca et al.)  
[![J.Phys.Conf.Ser. 1342 (2020) no.1, 012076](https://img.shields.io/badge/J.Phys.Conf.Ser.-1342(2020)no.1,012076-255773.svg)](https://iopscience.iop.org/article/10.1088/1742-6596/1342/1/012076)
[![[arXiv:1711.07692]](https://img.shields.io/badge/arXiv-1711.07692-B31B1B.svg)](https://arxiv.org/abs/1711.07692)

- *Description of CRESST-III data*  
CRESST Collaboration (A.H. Abdelhameed et al.)  
[![[arXiv:1905.07335]](https://img.shields.io/badge/arXiv-1905.07335-B31B1B.svg)](https://arxiv.org/abs/1905.07335)

#### CRESST-surface

- *Results on MeV-scale dark matter from a gram-scale cryogenic calorimeter operated above ground*  
CRESST Collaboration (G. Angloher et al.)  
[![Eur.Phys.J. C77 (2017) no.9, 637](https://img.shields.io/badge/Eur.Phys.J.-C77(2017)no.9,637-255773.svg)](https://link.springer.com/article/10.1140%2Fepjc%2Fs10052-017-5223-9)
[![[arXiv:1707.06749]](https://img.shields.io/badge/arXiv-1707.06749-B31B1B.svg)](https://arxiv.org/abs/1707.06749)

#### DAMIC_N_2012

- *Direct Search for Low Mass Dark Matter Particles with CCDs*  
DAMIC Collaboration (J. Barreto et al.)  
[![Phys.Lett. B711 (2012) 264](https://img.shields.io/badge/Phys.Lett.B-711(2012)264-255773.svg)](https://www.sciencedirect.com/science/article/pii/S0370269312003887?via%3Dihub)
[![[arXiv:1105.5191]](https://img.shields.io/badge/arXiv-1105.5191-B31B1B.svg)](https://arxiv.org/abs/1105.5191)

#### XENON1T_N_2017

- *First Dark Matter Search Results from the XENON1T Experiment*  
XENON Collaboration (E. Aprile et al.)  
[![Phys.Rev.Lett. 119 (2017) no.18, 181301](https://img.shields.io/badge/Phys.Rev.Lett.-119(2017)no.18,181301-255773.svg)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.181301)
[![[arXiv:1705.06655]](https://img.shields.io/badge/arXiv-1705.06655-B31B1B.svg)](https://arxiv.org/abs/1705.06655)


### Electron recoil experiments

#### CDMS-HVeV_2018

- *First Dark Matter Constraints from a SuperCDMS Single-Charge Sensitive Detector*  
SuperCDMS Collaboration (R. Agnese et al.)  
[![Phys.Rev.Lett. 121 (2018) no.5, 051301](https://img.shields.io/badge/Phys.Rev.Lett.-121(2018)no.5,051301-255773.svg)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.051301)
[![[arXiv:1804.10697]](https://img.shields.io/badge/arXiv-1804.10697-B31B1B.svg)](https://arxiv.org/abs/1804.10697)


- *Constraints on low-mass, relic dark matter candidates from a surface-operated SuperCDMS single-charge sensitive detector*  
SuperCDMS Collaboration (D.W. Amaral et al.)    
<!-- [![Phys.Rev.Lett. 121 (2018) no.5, 051301](https://img.shields.io/badge/Phys.Rev.Lett.-121(2018)no.5,051301-255773.svg)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.051301) -->
[![[arXiv:2005.14067]](https://img.shields.io/badge/arXiv-2005.14067-B31B1B.svg)](https://arxiv.org/abs/2005.14067)

<!-- #### DAMIC-e -->
<!-- 1907.12628 -->

#### DarkSide-50_S2

- *Constraints on Sub-GeV Dark-Matter–Electron Scattering from the DarkSide-50 Experiment*  
DarkSide Collaboration (P. Agnes et al.)  
[![Phys.Rev.Lett. 121 (2018) no.11, 111303](https://img.shields.io/badge/Phys.Rev.Lett.-121(2018)no.11,111303-255773.svg)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.111303)
[![[arXiv:1802.06998]](https://img.shields.io/badge/arXiv-1802.06998-B31B1B.svg)](https://arxiv.org/abs/1802.06998)

#### protoSENSEI@surface

- *SENSEI: First Direct-Detection Constraints on sub-GeV Dark Matter from a Surface Run*  
SENSEI Collaboration (Michael Crisler et al.)   
[![Phys.Rev.Lett. 121 (2018) no.6, 061803](https://img.shields.io/badge/Phys.Rev.Lett.-121(2018)no.6-255773.svg)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.061803)
[![[arXiv:1804.00088]](https://img.shields.io/badge/arXiv-1804.00088-B31B1B.svg)](https://arxiv.org/abs/1804.00088)

#### protoSENSEI@MINOS

- *SENSEI: Direct-Detection Constraints on Sub-GeV Dark Matter from a Shallow Underground Run Using a Prototype Skipper-CCD*  
SENSEI Collaboration (Orr Abramoff et al.)   
[![Phys.Rev.Lett. 122 (2019) no.16, 161801](https://img.shields.io/badge/Phys.Rev.Lett.-122(2019)no.16,161801-255773.svg)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.161801)
[![[arXiv:1901.10478]](https://img.shields.io/badge/arXiv-1901.10478-B31B1B.svg)](https://arxiv.org/abs/1901.10478)

#### SENSEI@MINOS

- *SENSEI: Direct-Detection Results on sub-GeV Dark Matter from a New Skipper-CCD*  
SENSEI Collaboration (Liron Barak et al.) 
[![[arXiv:2004.11378]](https://img.shields.io/badge/arXiv-2004.11378-B31B1B.svg)](https://arxiv.org/abs/2004.11378)

#### XENON10_S2

- *A search for light dark matter in XENON10 data*  
XENON10 Collaboration (J. Angle et al.)  
[![Phys.Rev.Lett. 107 (2011) 051301](https://img.shields.io/badge/Phys.Rev.Lett.-107(2011)051301-255773.svg)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.107.051301)
[![[arXiv:1104.3088]](https://img.shields.io/badge/arXiv-1104.3088-B31B1B.svg)](https://arxiv.org/abs/1104.3088)

- *First Direct Detection Limits on sub-GeV Dark Matter from XENON10*  
Rouven Essig, Aaron Manalaysay, Jeremy Mardon, Peter Sorensen, Tomer Volansky.  
[![Phys.Rev.Lett. 109 (2012) 021301](https://img.shields.io/badge/Phys.Rev.Lett.-109(2012)021301-255773.svg)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.109.021301)
[![[arXiv:1206.2644]](https://img.shields.io/badge/arXiv-1206.2644-B31B1B.svg)](https://arxiv.org/abs/1206.2644)

- *New Constraints and Prospects for sub-GeV Dark Matter Scattering off Electrons in Xenon*  
Rouven Essig, Tomer Volansky, Tien-Tien Yu  
[![Phys.Rev. D96 (2017) no.4, 043017](https://img.shields.io/badge/Phys.Rev.D-96(2017)no.4-255773.svg)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.96.043017)
[![[arXiv:1703.00910]](https://img.shields.io/badge/arXiv-1703.00910-B31B1B.svg)](https://arxiv.org/abs/1703.00910)


#### XENON100_S2

- *Low-mass dark matter search using ionization signals in XENON100*  
XENON Collaboration (E. Aprile et al.)   
[![Phys.Rev. D94 (2016) no.9, 092001](https://img.shields.io/badge/Phys.Rev.D-94(2016)no.9-255773.svg)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.092001)
[![[arXiv:1605.06262]](https://img.shields.io/badge/arXiv-1605.06262-B31B1B.svg)](https://arxiv.org/abs/1605.06262)

- *New Constraints and Prospects for sub-GeV Dark Matter Scattering off Electrons in Xenon*  
Rouven Essig, Tomer Volansky, Tien-Tien Yu   
[![Phys.Rev. D96 (2017) no.4, 043017](https://img.shields.io/badge/Phys.Rev.D-96(2017)no.4-255773.svg)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.96.043017)
[![[arXiv:1703.00910]](https://img.shields.io/badge/arXiv-1703.00910-B31B1B.svg)](https://arxiv.org/abs/1703.00910)

#### XENON1T_S2

- *Light Dark Matter Search with Ionization Signals in XENON1T*  
XENON Collaboration (E. Aprile et al.)  
[![Phys.Rev.Lett. 123 (2019) no.25, 251801](https://img.shields.io/badge/Phys.Rev.Lett.-123(2019)no.25,251801-255773.svg)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.251801)
[![[arXiv:1907.11485]](https://img.shields.io/badge/arXiv-1907.11485-B31B1B.svg)](https://arxiv.org/abs/1907.11485)

## CITATION

```
@software{Emken2021-2,
  author = {Emken, Timon},
  title = {{obscura - A C++ library for dark matter detection computations [Code, v0.1.0]}},
  year         = {2021},
  publisher    = {Zenodo},
  version      = {v0.1.0},
  doi          = {DOI:10.5281/zenodo.4557188},
  url          = {https://doi.org/10.5281/zenodo.4557188},
  howpublished={The code can be found under \url{https://github.com/temken/obscura}. Version 0.1.0 is archived as \href{https://doi.org/10.5281/zenodo.4557188}{DOI:10.5281/zenodo.4557188}}
}
```

## VERSION HISTORY

- 23.02.2021: Release of version 0.1.0

## AUTHORS & CONTACT

The author of this tool/library is Timon Emken.

For questions, bug reports or other suggestions please contact [timon.emken@fysik.su.se](mailto:timon.emken@fysik.su.se).


## LICENSE

This project is licensed under the MIT License - see the LICENSE file.
