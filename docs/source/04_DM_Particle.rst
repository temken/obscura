.. _Section_DM_Particle:

==============================
4. The ``DM_Particle`` classes
==============================

The ``DM_Particle`` class and its derived classes are responsible for the particle physics aspects of direct detection.
In particular, an instance of ``DM_Particle`` entails the particle properties of a DM candidate particle, such as its mass, spin, and its differential and total interaction cross sections with nuclei or electrons.
The base class's functions provide an interface that is sufficient for the calculation of e.g. event rates at direct DM search experiments.
Furthermore, it contains a number of functions regarding scattering angles, their distributions and sampling, which might not be relevant for direct detection, but can be used in the context of e.g. MC simulations.

--------------------------
The interface / base class
--------------------------

The abstract base class is defined in `/include/obscura/DM_Particle.hpp <https://github.com/temken/obscura/blob/main/include/obscura/DM_Particle.hpp>`_ and all the member functions and parameters can be seen there.

The most important (virtual) functions for direct detection specific calculations are the differential cross sections.

.. code-block:: c++

   //Differential cross sections for nuclear targets
   virtual double dSigma_dq2_Nucleus(double q, const Isotope& target, double vDM, double param = -1.0) const { return 0.0; };
   double dSigma_dER_Nucleus(double ER, const Isotope& target, double vDM, double param = -1.0) const;
   double d2Sigma_dER_dEe_Migdal(double ER, double Ee, double vDM, const Isotope& isotope, Atomic_Electron& shell) const;

   // Differential cross section for electron targets
   virtual double dSigma_dq2_Electron(double q, double vDM, double param = -1.0) const { return 0.0; };
   virtual double d2Sigma_dq2_dEe_Ionization(double q, double Ee, double vDM, Atomic_Electron& shell) const { return 0.0; };
   virtual double d2Sigma_dq2_dEe_Crystal(double q, double Ee, double vDM, Crystal& crystal) const { return 0.0; };

We point out that here we have to pass instances of the target classes discussed in the previous section (i.e. nuclear isotopes, atomic electrons, and electrons in crystals).
Also included is a simple implementation of Migdal scatterings with atomic targets based on [Essig2020]_.

The most standard DM candidate considered in the direct detection literature is a WIMP with SI or SD interactions. *obscura* contains derived classes for each of these scenarios, which are declared in `/include/obscura/DM_Particle_Standard.hpp <https://github.com/temken/obscura/blob/main/include/obscura/DM_Particle_Standard.hpp>`_.

----------------------------------
Spin-Independent (SI) interactions
----------------------------------

The differential cross section for SI nuclear interactions is given by

.. math::
	\frac{\mathrm{d}\sigma^{\rm SI}_{N}}{\mathrm{d} E_R} =\frac{m_N}{2\pi v_\chi^2}\left[f_p Z+f_n(A-Z) \right]^2 \left|F^{\rm SI}_N\left(E_R\right)\right|^2\, .

For details, we refere to e.g. chapter 3.4 of [Emken2019]_.

The class ``DM_Particle_SI`` is derived from ``DM_Particle`` and evaluates the cross sections of SI interactions with nuclei (and electrons).

The following example demonstrates how to

* construct an instance of ``DM_Particle_SI`` that describes a DM particle of 10 GeV mass.
* set the SI proton cross section to :math:`\sigma_p=10^{-40}\mathrm{cm}^2`, the electron cross section of :math:`\sigma_e=10^{-36}\mathrm{cm}^2`.
* to evaluate the differential and total scattering cross section with argon nuclei.

.. code-block:: c++

  #include "libphysica/Natural_Units.hpp"

  #include "obscura/DM_Particle_Standard.hpp"

  using namespace libphysica::natural_units;

  // ...

  // Declare the DM particle
  obscura::DM_Particle_SI dm(10.0 * GeV);
  dm.Set_Sigma_Proton(1.0e-40 * cm * cm);
  dm.Set_Sigma_Electron(1.0e-36 * cm * cm);

  // Define the target
  obscura::Isotope argon = obscura::Get_Isotope(18, 40);

  // Evaluate cross sections
  double E_R = 1.0 * keV;
  double v_DM = 300.0 * km/sec;
  double diff_cross_section = dm.dSigma_dER_Nucleus(E_R, argon, v_DM);
  double tot_cross_section = dm.Sigma_Total_Nucleus(argon, v_DM);

  // Convert to other units
  std::cout <<In_Units(diff_cross_section, cm * cm / keV)<<std::endl;
  std::cout <<In_Units(tot_cross_section, cm * cm)<<std::endl;

--------------------------------
Spin-Dependent (SD) interactions
--------------------------------

The differential cross section for SD nuclear interactions is given by

.. math::
	\frac{\mathrm{d} \sigma_N^{\rm SD}}{\mathrm{d} E_R} = \frac{2m_N}{\pi v_\chi^2}\frac{J+1}{J}\left(f_p \langle S_p\rangle +f_n \langle S_N\rangle\right)^2 \left.F_N^{\rm SD}(E_R)^2\right|

Similarly to ``DM_Particle_SI``, we also define a ``DM_Particle_SD`` class, which evaluates this cross section for nuclear targets with spin :math:`S\neq0`.