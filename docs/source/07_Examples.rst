====================================
7. Examples: Putting it all together
====================================


-------------------------------------------------------------
Computing the recoil spectrum of SI & SD nuclear interactions
-------------------------------------------------------------

As an example for nuclear recoils that also illustrates nicely the modular structure of *obscura*, we compute the nuclear recoil spectrum :math:`\frac{ \mathrm{d}R}{\mathrm{d}E_R}` for a 10 GeV DM particle interacting with xenon nuclei via spin-independent and spin-dependent interactions.

For the definition and details of the nuclear recoil spectrum, see e.g. chapter 3.5 of [Emken2019]_.

1. First we define the DM particle objects that describe SI and SD interactions

.. code-block:: c++

   // 1. DM particle (SI and Sd)
   obscura::DM_Particle_SI dm_SI(10.0 * GeV);
   dm_SI.Set_Sigma_Proton(1.0e-40 * cm * cm);
   dm_SI.Print_Summary();

   obscura::DM_Particle_SD dm_SD(10.0 * GeV);
   dm_SD.Set_Sigma_Proton(1.0e-40 * cm * cm);
   dm_SD.Print_Summary();    

The ``Print_Summary()`` function is a member of many of the classes and provides a terminal output that summarizes the object.

2. For the DM distribution we use the standard halo model with default parameters.

.. code-block:: c++

   // 2. DM distribution
   obscura::Standard_Halo_Model shm;
   shm.Print_Summary();

3. As target nuclei, we choose xenon and import the nuclear data.

.. code-block:: c++

   // 3. Direct detection targets
   obscura::Nucleus xenon = obscura::Get_Nucleus("Xe");
   xenon.Print_Summary();

4. With these three objects, we can compute the differential nuclear recoil spectrum for a given recoil energy :math:`E_R`.

.. code-block:: c++

   double E_R		= 1.0 * keV;
   double dRdER_SI = obscura::dRdER_Nucleus(E_R, dm_SI, shm, xenon);
   double dRdER_SD = obscura::dRdER_Nucleus(E_R, dm_SD, shm, xenon);

5. The results are given in natural units in powers of GeV. To convert it to another unit, we can use the unit functionality of the *libphysica* library.

.. code-block:: c++

   std::cout << "SI-interactions: \tdR/dER (1 keV) = " << In_Units(dRdER_SI, 1.0 / kg / year / keV) << " events / kg / year / keV" << std::endl;
   std::cout << "SD-interactions: \tdR/dER (1 keV) = " << In_Units(dRdER_SD, 1.0 / kg / year / keV) << " events / kg / year / keV" << std::endl;

.. raw:: html

	<details>
	<summary><a>The full main.cpp</a></summary>
 
.. code-block:: c++

   #include <iostream>

   #include "libphysica/Natural_Units.hpp"

   #include "obscura/DM_Halo_Models.hpp"
   #include "obscura/DM_Particle_Standard.hpp"
   #include "obscura/Direct_Detection_Nucleus.hpp"
   #include "obscura/Target_Nucleus.hpp"

   using namespace libphysica::natural_units;

   int main()
   {

   	// 1. DM particle (SI and Sd)
   	obscura::DM_Particle_SI dm_SI(10.0 * GeV);
   	dm_SI.Set_Sigma_Proton(1.0e-40 * cm * cm);
   	dm_SI.Print_Summary();

   	obscura::DM_Particle_SD dm_SD(10.0 * GeV);
   	dm_SD.Set_Sigma_Proton(1.0e-40 * cm * cm);
   	dm_SD.Print_Summary();

   	// 2. DM distribution
   	obscura::Standard_Halo_Model shm;
   	shm.Print_Summary();

   	// 3. Direct detection targets
   	obscura::Nucleus xenon = obscura::Get_Nucleus("Xe");
   	xenon.Print_Summary();

   	// 4. Evalute the nuclear recoil spectrum
   	double E_R		= 1.0 * keV;
   	double dRdER_SI = obscura::dRdER_Nucleus(E_R, dm_SI, shm, xenon);
   	double dRdER_SD = obscura::dRdER_Nucleus(E_R, dm_SD, shm, xenon);

   	std::cout << "SI-interactions: \tdR/dER (1 keV) = " << In_Units(dRdER_SI, 1.0 / kg / year / keV) << " events / kg / year / keV" << std::endl;
   	std::cout << "SD-interactions: \tdR/dER (1 keV) = " << In_Units(dRdER_SD, 1.0 / kg / year / keV) << " events / kg / year / keV" << std::endl;

   	return 0;
   }

.. raw:: html

	</details>

.. raw:: html

	<details>
	<summary><a>The terminal output</a></summary>

.. code-block::

   ----------------------------------------
   DM particle summary:
           Mass:                   10 GeV
           Spin:                   0.5
           Low mass:               [ ]

           Interaction:            Spin-Independent (SI)

           Coupling ratio fixed:   [x]
           Isospin conservation:   [x]
           Coupling ratio:         fn/fp = 1

           Sigma_P[cm^2]:          1e-40
           Sigma_N[cm^2]:          1e-40
           Sigma_E[cm^2]:          1e-40

           Interaction type:       Contact
   ----------------------------------------

   ----------------------------------------
   DM particle summary:
           Mass:                   10 GeV
           Spin:                   0.5
           Low mass:               [ ]
   Interaction:            Spin-Dependent (SD)

           Coupling ratio fixed:   [x]
           Isospin conservation:   [x]
           Coupling ratio:         fn/fp = 1

           Sigma_P[cm^2]:          1e-40
           Sigma_N[cm^2]:          1e-40
           Sigma_E[cm^2]:          1e-40

   ----------------------------------------

   Dark matter distribution - Summary
           Standard halo model (SHM)

           Local DM density[GeV/cm^3]:     0.4
           Speed domain [km/sec]:          [0,777]
           Average DM velocity [km/sec]:   (-11.1 , -232 , -7.3)
           Average DM speed [km/sec]:      330

           Speed dispersion v_0[km/sec]:   220
           Gal. escape velocity [km/sec]:  544
           Observer's velocity [km/sec]:   (11.1 , 232 , 7.3)
           Observer's speed [km/sec]:      233


   Xe
   Isotope Z       A       Abund.[%]       Spin    <sp>    <sn>
   ------------------------------------------------------------
   Xe-124  54      124     0.095           0       0       0
   Xe-126  54      126     0.089           0       0       0
   Xe-128  54      128     1.91            0       0       0
   Xe-129  54      129     26.4            0.5     0.01    0.329
   Xe-130  54      130     4.07            0       0       0
   Xe-131  54      131     21.2            1.5     -0.009  -0.272
   Xe-132  54      132     26.9            0       0       0
   Xe-134  54      134     10.4            0       0       0
   Xe-136  54      136     8.86            0       0       0
   Total:          131     99.999

   SI-interactions:        dR/dER (1 keV) = 13621.8 events / kg / year / keV
   SD-interactions:        dR/dER (1 keV) = 0.132525 events / kg / year / keV

.. raw:: html

   </details>

--------------------------------------------------------------------------
Exclusion limits for a sub-GeV DM particle via electron recoil experiments
--------------------------------------------------------------------------

As a second example for an application of *obscura*, we will compute the 95% confidence level exclusion limit on the DM-electron cross section for a sub-GeV DM particle.

We assume a DM mass of 100 MeV, and two different direct detection experiments.

1. An argon based experiment with an exposure of 100 kg years and an observational threshold of at least 4 ionized electrons.
2. A semiconductor experiment with Si crystal targets, an exposure of 10 gram years, and an observational threshold of minimum 2 electron-hole pairs.

Let us set up the different objects to obtain the limits.

1. First we define the DM particle object with 100 MeV mass.

.. code-block:: c++

   // 1. DM particle
   obscura::DM_Particle_SI dm(100.0 * MeV);
   dm.Print_Summary();   

2. For the DM distribution we again use the standard halo model with default parameters.

.. code-block:: c++

   // 2. DM distribution
   obscura::Standard_Halo_Model shm;
   shm.Print_Summary();

3. For the first experiment, we create an instance of the ``DM_Detector_Ionization_ER`` class and specify the desired detector properties of the toy experiment.

.. code-block:: c++

   // 3. Argon target experiment
   obscura::DM_Detector_Ionization_ER argon_experiment("Argon toy experiment", 100.0 * kg * year, "Ar");
   argon_experiment.Use_Electron_Threshold(4);
   argon_experiment.Print_Summary();

4. The same for the semiconductor experiment:

.. code-block:: c++

   // 4. Si target experiment
   obscura::DM_Detector_Crystal silicon_experiment("Silicon toy experiment", 10.0 * gram * year, "Si");
   silicon_experiment.Use_Q_Threshold(2);
   silicon_experiment.Print_Summary();

4. With these three objects, we can compute the limit on the DM-electron cross section.

.. code-block:: c++

   // 5. Compute the 95% CL exclusion limits for m = 100.0 MeV
   double limit_Ar = argon_experiment.Upper_Limit(dm, shm, 0.95);
   double limit_Si = silicon_experiment.Upper_Limit(dm, shm, 0.95);

1. As in the previous example, the results are given in natural units in powers of GeV. We convert it to :math:`\mathrm{cm}^2`, and print the result on the terminal.

.. code-block:: c++

   std::cout << "Argon experiment: \tsigma_e < " << In_Units(limit_Ar, cm * cm) << " cm^2 (95%CL)" << std::endl;
   std::cout << "Silicon experiment: \tsigma_e < " << In_Units(limit_Si, cm * cm) << " cm^2 (95%CL)" << std::endl;

.. raw:: html

	<details>
	<summary><a>The full main.cpp</a></summary>
 
.. code-block:: c++

   #include <iostream>

   #include "libphysica/Natural_Units.hpp"

   #include "obscura/DM_Halo_Models.hpp"
   #include "obscura/DM_Particle_Standard.hpp"
   #include "obscura/Direct_Detection_Crystal.hpp"
   #include "obscura/Direct_Detection_ER.hpp"
   #include "obscura/Target_Atom.hpp"
   #include "obscura/Target_Crystal.hpp"

   using namespace libphysica::natural_units;

   int main()
   {

   	// 1. DM particle
   	obscura::DM_Particle_SI dm(100.0 * MeV);
   	dm.Print_Summary();

   	// 2. DM distribution
   	obscura::Standard_Halo_Model shm;
   	shm.Print_Summary();

   	// 3. Argon target experiment
   	obscura::DM_Detector_Ionization_ER argon_experiment("Argon toy experiment", 100.0 * kg * year, "Ar");
   	argon_experiment.Use_Electron_Threshold(4);
   	argon_experiment.Print_Summary();

   	// 4. Si target experiment
   	obscura::DM_Detector_Crystal silicon_experiment("Silicon toy experiment", 10.0 * gram * year, "Si");
   	silicon_experiment.Use_Q_Threshold(2);
   	silicon_experiment.Print_Summary();

   	// 5. Compute the 95% CL exclusion limits for m = 100.0 MeV
   	double limit_Ar = argon_experiment.Upper_Limit(dm, shm, 0.95);
   	double limit_Si = silicon_experiment.Upper_Limit(dm, shm, 0.95);

   	std::cout << "Argon experiment: \tsigma_e < " << In_Units(limit_Ar, cm * cm) << " cm^2 (95%CL)" << std::endl;
   	std::cout << "Silicon experiment: \tsigma_e < " << In_Units(limit_Si, cm * cm) << " cm^2 (95%CL)" << std::endl;

   	return 0;
   }

.. raw:: html

	</details>

.. raw:: html

	<details>
	<summary><a>The terminal output</a></summary>
 
.. code-block::
  
   ----------------------------------------
   DM particle summary:
           Mass:                   100 MeV
           Spin:                   0.5
           Low mass:               [ ]

           Interaction:            Spin-Independent (SI)

           Coupling ratio fixed:   [x]
           Isospin conservation:   [x]
           Coupling ratio:         fn/fp = 1

           Sigma_P[cm^2]:          1e-40
           Sigma_N[cm^2]:          1e-40
           Sigma_E[cm^2]:          1e-40

           Interaction type:       Contact
   ----------------------------------------
   Dark matter distribution - Summary
           Standard halo model (SHM)

           Local DM density[GeV/cm^3]:     0.4
           Speed domain [km/sec]:          [0,777]
           Average DM velocity [km/sec]:   (-11.1 , -232 , -7.3)
           Average DM speed [km/sec]:      330

           Speed dispersion v_0[km/sec]:   220
           Gal. escape velocity [km/sec]:  544
           Observer's velocity [km/sec]:   (11.1 , 232 , 7.3)
           Observer's speed [km/sec]:      233


   ----------------------------------------
   Experiment summary:     Argon toy experiment
           Target particles:       Electrons
           Exposure [kg year]:     100
           Flat efficiency [%]:    100
           Observed events:        0
           Expected background:    0
           Statistical analysis:   Poisson


           Electron recoil experiment (ionization).
           Target(s):
                           Ar      (100%)
           Electron bins:          [ ]
           PE (S2) bins:           [ ]
                   Ne threshold:   4
                   Ne max:         15
   ----------------------------------------


   ----------------------------------------
   Experiment summary:     Silicon toy experiment
           Target particles:       Electrons
           Exposure [kg year]:     0.01
           Flat efficiency [%]:    100
           Observed events:        0
           Expected background:    0
           Statistical analysis:   Poisson


           Electron recoil experiment (semiconductor).
           Target:                 Si semiconductor
           eh pair threshold:      2
   ----------------------------------------

   Argon experiment:       sigma_e < 1.67038e-41 cm^2 (95%CL)
   Silicon experiment:     sigma_e < 1.1756e-39 cm^2 (95%CL)

.. raw:: html

	</details>