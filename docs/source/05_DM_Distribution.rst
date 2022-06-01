==================================
5. The ``DM_Distribution`` classes
==================================

In order to make predictions for direct detection experiments, the statistical properties of the incoming DM flux need to be specified. In particular, we need to know how many DM particles pass through the detector and with what energy.
In other words, we need to know the DM particle flux, or alternatively the local DM density and the energy distribution.



--------------------------
The interface / base class
--------------------------


The class ``DM_Distribution``, that is declared in `/include/obscura/DM_Distribution.hpp <https://github.com/temken/obscura/blob/main/include/obscura/DM_Distribution.hpp>`_, is an abstract base class or interface that defines all the functions we require to characterize a distribution and flux of DM particles.

Most importantly, the class provides interfaces to probability density functions (PDFs) for the DM particles' velocity or speed, their local energy density, differential particle flux, etc.

-----------------------------
The standard halo model (SHM)
-----------------------------

The conventional assumptions on the halo DM particles' properties is the Standard Halo Model (SHM).
The SHM describes the galactic DM by a truncated Maxwell-Boltzmann distribution. It is characterized by the following 4 parameters (for details see e.g. chapter 3.2 of [Emken2019]_

.. math::
	\rho_\chi, v_0, v_\mathrm{esc}, \mathbf{v}_\mathrm{obs}

In `/include/obscura/DM_Halo_Models.hpp <https://github.com/temken/obscura/blob/main/include/obscura/DM_Halo_Models.hpp>`_ we define the ``Standard_Halo_Model`` class which is an implemenation of this model.
It is a derived class of ``DM_Distribution``.

We can construct the SHM model by the default constructor, which assumes default values for the 4 parameters.

.. code-block:: c++

	#include "obscura/DM_Halo_Models.hpp"

	// ...

	obscura::Standard_Halo_Model shm;

Or we define the parameters explicitly.

.. code-block:: c++

	#include "libphysica/Natural_Units.hpp"

	#include "obscura/DM_Halo_Models.hpp"

	using namespace libphysica::natural_units;

	// ...
	double rho = 0.4 * GeV / cm / cm / cm;
	double v_0 = 230.0 * km / sec;
	double v_esc = 600 * km / sec;
	double v_obs = 232.0 * km / sec;
	obscura::Standard_Halo_Model shm(rho, v_0, v_obs, v_esc);

---------
The SHM++
---------

As a second example for a DM halo model, *obscura* also implements the SHM++ as proposed in [Evans2019]_.

Since it extends the SHM, the corresponding class ``SHM_Plus_Plus`` is a derived class of ``Standard_Halo_Model`` which is in turn derived from ``DM_Distribution``.
The class is also declared in `/include/obscura/DM_Halo_Models.hpp <https://github.com/temken/obscura/blob/main/include/obscura/DM_Halo_Models.hpp>`_.

This halo model can be constructed and used essentially identically to the SHM.

-------------------------
Imported DM distributions
-------------------------

It is also possible to import a DM distribution from a file.
This is the purpose of the ``Imported_DM_Distribution`` class, another derived class of ``DM_Distribution`` which can be found in `/include/obscura/DM_Distribution.hpp <https://github.com/temken/obscura/blob/main/include/obscura/DM_Distribution.hpp>`_.

As input file, we need a two-column table of the DM speed PDF using the format (v[km/sec] :: f(v) [sec/km]).
Additionally we need to specify the local DM density.

Here is an example of using this class assuming a tabulated speed pdf given in the file *DM_Speed_PDF.txt*.


.. code-block:: c++

	#include "libphysica/Natural_Units.hpp"

	#include "obscura/DM_Distribution.hpp"

	using namespace libphysica::natural_units;

	// ...
	double rho = 0.4 * GeV / cm / cm / cm;
	obscura::Imported_DM_Distribution dm_distribution(rho, "DM_Speed_PDF.txt");