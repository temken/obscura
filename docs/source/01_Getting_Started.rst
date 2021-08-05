==================
1. Getting started
==================

------------
Installation
------------

Before building *obscura*, there are a few libraries that need to be installed.

^^^^^^^^^^^^
Dependencies
^^^^^^^^^^^^

""""""""""""""""""""""""""""""""""""
1. `boost <https://www.boost.org/>`_
""""""""""""""""""""""""""""""""""""

To install *boost* on a Mac, we can use `homebrew <https://brew.sh/>`_ ::

	brew install boost

On Linux machines, run::

   sudo apt-get update && sudo apt-get install -yq libboost-all-dev


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
2. `libconfig <https://hyperrealm.github.io/libconfig/>`_
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To install *boost* on a Mac, we can use `homebrew <https://brew.sh/>`_ ::

	brew install libconfig

On Linux machines, you can build `libconfig` via::

	wget https://hyperrealm.github.io/libconfig/dist/libconfig-1.7.2.tar.gz
	tar -xvzf libconfig-1.7.2.tar.gz
	pushd libconfig-1.7.2
	./configure
	make
	sudo make install
	popd

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
3. `libphysica <https://github.com/temken/libphysica>`_
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

*libphysica* does not need to be installed. It will be downloaded and compiled during the CMake build.

^^^^^^^^^^^^^^^^
Download & Build
^^^^^^^^^^^^^^^^

The `obscura` source code can be downloaded by cloning this `git repository <https://github.com/temken/obscura>`_: ::

   git clone https://github.com/temken/obscura.git
   cd obscura

The code is compiled and the executable and library is built by `CMake <https://cmake.org/>`_. To build run the following commands from the repository's root folder.::

	cmake -E make_directory build
	cd build
	cmake -DCMAKE_BUILD_TYPE=Release -DCODE_COVERAGE=OFF ..
	cmake --build . --config Release
	cmake --install .

If everything worked well, the executable and library file are created as::

	bin/obscura
	lib/libobscura.a


-------------------------
Using *obscura* as a tool
-------------------------

*Obscura* can be used as a tool and builds an executable which can be run from */bin/* via::

./obscura config.cfg

As can be seen in the `/src/main.cpp <https://github.com/temken/obscura/blob/master/src/main.cpp>`_ file, this script computes direct detection limits and saves them in the */results/* folder.
The specifications of the exclusion limits (DM physics and halo model, statistics, experiment, mass range,...) are defined in a configuration file, in this case *config.cfg*.
For the handling of configuration files, *obscura* relies on `libconfig <https://hyperrealm.github.io/libconfig/>`_. 

^^^^^^^^^^^^^^^^^^^^^^
The configuration file
^^^^^^^^^^^^^^^^^^^^^^

The configuration file contains all input parameters necessary to define the various *obscura* models.

.. warning::

	The import of these parameters via libconfig is very case-sensitive. A float parameter has to be set to e.g. *1.0*, and **not** just *1*.

.. raw:: html

	<details>
	<summary><a>The full configuration file</a></summary>
 
.. code-block:: c++

   //obscura - Configuration File

   //ID
   	ID		=	"test";

   //Dark matter particle
   	DM_mass		  	=	0.1;		// in GeV
   	DM_spin		  	=	0.5;
   	DM_fraction		=	1.0;		// the DM particle's fractional abundance (set to 1.0 for 100%)
   	DM_light		=	false;		// Options: true or false. low mass mode

   	DM_interaction		=	"SI";		// Options: "SI" or "SD"

   	DM_isospin_conserved		=	true; 		// only relevant for SI and SD
   	DM_relative_couplings		=	(1.0, 0.0); //relation between proton (left) and neutron (right) couplings.
   												//only relevant if 'DM_isospin_conserved' is false.
   	DM_cross_section_nucleon	=	1.0e-36;	//in cm^2
   	DM_cross_section_electron	=	1.0e-36;	//in cm^2 (only relevant for SI and SD)
   	DM_form_factor		=	"Contact";	// Options: "Contact", "Electric-Dipole", "Long-Range", "General"
   												//(only relevant for SI)
   	DM_mediator_mass	=	0.0;		// in MeV (only relevant if 'DM_form_factor' is "General")

   //Dark matter distribution
   	DM_distribution 	=	"SHM";		//Options: "SHM", "SHM++", "File"
   	DM_local_density	=	0.4;		//in GeV / cm^3
   	
   	//Options for "SHM" and "SHM++"
   		SHM_v0		=	220.0;				//in km/sec
   		SHM_vObserver	=	(0.0, 232.0, 0.0);	//in km/sec
   		SHM_vEscape	=	544.0;				//in km/sec
   	//Options for "SHM++"
   		SHMpp_eta	=	0.2;
   		SHMpp_beta	=	0.9;
   	//Options for "File" (The file has to be a 2-column table of format v[km/sec] :: f(v) [sec/km])
   		file_path  = "DM_Speed_PDF.txt";

   //Dark matter detection experiment
   	DD_experiment	=	"Electron recoil";	//Options for nuclear recoils: "Nuclear recoil", "DAMIC_N_2011", "XENON1T_N_2017", "CRESST-II","CRESST-III", "CRESST-surface"
							//Options for electron recoils: "Semiconductor","protoSENSEI@MINOS","protoSENSEI@surface", "SENSEI@MINOS", "CDMS-HVeV_2018", "CDMS-HVeV_2020", "Electron recoil", "XENON10_S2", "XENON100_S2", "XENON1T_S2", "DarkSide-50_S2"

   	//Options for user-defined experiments ("Nuclear recoil", "Electron recoil", and "Semiconductor")
	  //General
	  DD_exposure 		=	1.0;	//in kg years
	  DD_efficiency 		=	1.0;	//flat efficiency
	  DD_observed_events 	=	0;		//observed signal events
	  DD_expected_background 	=	0.0;	//expected background events

	  //Specific options for "Nuclear recoil"
	  DD_targets_nuclear	=	(
	  				(4.0, 8),
	  				(1.0, 20),
	  				(1.0, 74)
	  			);				// Nuclear targets defined by atom ratio/abundances and Z
	  DD_threshold_nuclear    =	4.0;    //in keV
	  DD_Emax_nuclear         =	40.0;	//in keV
	  DD_energy_resolution    =	0.0;    //in keV

	  //Specific options for "Electron recoil" and "Semiconductor:
	  DD_target_electron	=	"Xe";	//Options for "Electron recoil": 	"Xe", "Ar"
	  								//Options for "Semiconductor":	"Si", "Ge"
	  DD_threshold_electron	=	4;		//In number of electrons or electron hole pairs.

   //Computation of exclusion limits
   	constraints_certainty	=	0.95;	//Certainty level
   	constraints_mass_min	=	0.02;	//in GeV										
   	constraints_mass_max	=	1.0;	//in GeV
   	constraints_masses	=	10;										
 
.. raw:: html

	</details>

----------------------------
Using *obscura* as a library
----------------------------

If we want to use *obscura* functions in an external code, we can do so and import it as a library.
We recommend to do this inside your CMake build, where *obscura* can be downloaded, built, included, and linked automatically during the build of your code.


As an instructional example `this repository <https://github.com/temken/template_cpp_cmake_obscura>`_ contains a C++ project template built with CMake that imports and uses the *obscura* library.