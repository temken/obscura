==============================
6. The ``DM_Detector`` classes
==============================

The details of a direct detection experiment are summarized in the ``DM_Detector`` class declared in `/include/obscura/Direct_Detection.hpp <https://github.com/temken/obscura/blob/master/include/obscura/Direct_Detection.hpp>`_.
In particular, it responsible for:

1. The statistical methods to compute likelihoods and exclusion limits. Since these are independent of the type of experiment, this functionality is part of the base class ``DM_Detector``. As of now, *obscura* implements the following statistical analyses.
   1. Poisson statistics
   2. Binned Poisson statistics
   3. Maximum gap following [Yellin2002]_. 
2. The detector details, such as detection efficiencies, energy resolution, target particles, etc. These can be very specific and are implemented in classes derived from ``DM_Detector``, e.g. ``DM_Detector_Nucleus``.


We provide a number of examples of how to construct different instances of derived classes of ``DM_Detector``.

--------------------------
Nuclear recoil experiments
--------------------------

For experiments looking for DM induced nuclear recoils, *obscura* contains the ``DM_Detector_Nucleus`` class that is declared in `/include/obscura/Direct_Detection_Nucleus.hpp <https://github.com/temken/obscura/blob/master/include/obscura/Direct_Detection_Nucleus.hpp>`_.

For example, assume we have a nuclear recoil experiment with :math:`\mathrm{CaWO}_4` crystals, an energy threshold of 500 eV, and an exposure of 100 kg days.
This information suffices to define a toy experiment.

.. code-block:: c++

   #include "libphysica/Natural_Units.hpp"

   #include "obscura/Target_Nucleus.hpp"
   #include "obscura/DM_Detector_Nucleus.hpp"

   using namespace libphysica::natural_units;

   // ...

   double exposure = 100.0 * kg * day;
   std::vector<Nucleus> nuclear_targets = {obscura::Get_Nucleus(8), obscura::Get_Nucleus(20), obscura::Get_Nucleus(74)};
   std::vector<double> target_ratios = {4, 1, 1};
   double energy_threshold = 500 * eV;
   obscura::DM_Detector_Nucleus detector("Nuclear recoil experiment", exposure, nuclear_targets, target_ratios);

---------------------------
Electron recoil experiments
---------------------------

For electron recoil experiments with atomic targets, we have to use the ``DM_Detector_Ionization_ER`` class that can be found in `/include/obscura/Direct_Detection_ER.hpp <https://github.com/temken/obscura/blob/master/include/obscura/Direct_Detection_ER.hpp>`_

Here is an example of a xenon target experiment probing DM-electron interactions and DM induced ionizations. As exposure we choose 100 kg days, and we furthermore assume that only events with at least 4 ionized electrons can be detectred.

.. code-block:: c++

   #include "libphysica/Natural_Units.hpp"

   #include "obscura/DM_Detector_ER.hpp"

   using namespace libphysica::natural_units;

   // ...

   double exposure = 100.0 * kg * day;
   obscura::DM_Detector_Ionization_ER xenon_experiment("Electron recoil experiment", exposure, "Xe");
   argon_experiment.Use_Electron_Threshold(4);

Alternatively, many experiments looking for sub-GeV DM use semiconductor crystals as targets. In this case, there is another derived class, ``DM_Detector_Crystal``.

Again we construct an example toy experiment. This time, we assume a silicon crystal target, choose an exposure of 10 g year, and assume that only events with at least 2 electron-hole pairs can trigger the detector.

.. code-block:: c++

   #include "libphysica/Natural_Units.hpp"

   #include "obscura/DM_Detector_Crystal.hpp"

   using namespace libphysica::natural_units;

   // ...

   double exposure = 10.0 * gram * year;
   obscura::DM_Detector_Crystal silicon_experiment("Crystal target experiment", exposure, "Si");
   silicon_experiment.Use_Q_Threshold(2);

