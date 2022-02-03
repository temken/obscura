.. image:: https://github.com/temken/obscura/actions/workflows/main.yml/badge.svg?branch=master
   :target: https://github.com/temken/obscura/actions/workflows/main.yml
   :alt: Build Status
.. image:: https://codecov.io/gh/temken/obscura/branch/master/graph/badge.svg?token=1Pe1QMcngr
   :target: https://codecov.io/gh/temken/obscura
   :alt: Code Coverage 
.. image:: https://readthedocs.org/projects/obscura/badge/?version=latest
   :target: https://obscura.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

========================================================================================
*obscura* - Direct detection of dark matter with nucleus and electron recoil experiments
========================================================================================

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4557187.svg
   :target: https://doi.org/10.5281/zenodo.4557187
   :alt: DOI
.. image:: https://joss.theoj.org/papers/fd8076268036956d3bf08193c4fc2db9/status.svg
   :target: https://joss.theoj.org/papers/fd8076268036956d3bf08193c4fc2db9
   :alt: JOSS paper

A modular C++ tool and library for dark matter direct detection computations for both nuclear and electron recoil experiments.

The purpose of this documentation or manual is to provide insight into the polymorphic class structure of *obscura* and how it can be applied in different contexts.
It should also serve as a guide and describe the usage of *obscura* via code examples.

The documentation does not contain a review of the physics implemented in the library.
For more physics details, we refer to e.g. chapter 3 of [Emken2019]_ or [Nobile2021]_.

If you want to contribute to `obscura`, please check out the `contribution guidelines <https://github.com/temken/obscura/blob/master/docs/CONTRIBUTING.md>`_.

.. image:: https://raw.githubusercontent.com/temken/obscura/master/paper/FlowChart.png
   :width: 500

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   01_Getting_Started
   02_Main_Modules
   03_Targets
   04_DM_Particle
   05_DM_Distribution
   06_DM_Detector
   07_Examples
   08_Experiments
   09_Citations
   10_Release_History
   11_License
   12_Contact
   References

For the interpretation of past and future direct searches for DM particles, it is important to be able to provide accurate predictions for event rates and spectra under a variety of possible and viable assumptions in a computationally efficient way.
While there exists a few tools to compute DM induced nuclear recoil spectra, such as `DDCalc <https://ddcalc.hepforge.org/>`_ or `WimPyDD <https://wimpydd.hepforge.org/>`_, `obscura` is not limited to nuclear targets.
Instead its main focus lies on sub-GeV DM searches probing electron recoils which typically requires methods from atomic and condensed matter physics, see e.g. [Essig2012]_ or [Catena2019]_.
In the context of sub-GeV DM searches, new ideas such as target materials or detection techniques are being proposed regularly, and the theoretical modelling of these are getting improved continuosly.
At the same time, currently running experiments continue to publish their results and analyses, setting increasingly strict bounds on the DM parameter space.
In such a dynamic field, `obscura` can be an invaluable tool due to its high level of adaptability and facilitate and accelerate the development of new, reliable research software for the preparation of a DM discovery in the hopefully near future.