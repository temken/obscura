=====================
3. The target classes
=====================

The basic hope of direct detection experiments is that the DM particles from the galactic halo occasionally collide with ordinary particles, see e.g. [Nobile2021]_.
The original target of direct DM searches were nuclear recoils.
Later on also electron targets gained more and more attention in the context of sub-GeV dark matter [Essig2012]_.

In *obscura* each target type is represented by a class, that can e.g. be passed to the cross section functions of the ``DM_Particle`` class, see :ref:`Section_DM_Particle`.

---------------
Nuclear targets
---------------

For nuclear recoil experiments, we define two target classes, ``Isotope`` and ``Nucleus`` that are declared in `/include/obscura/Target_Nucleus.hpp <https://github.com/temken/obscura/blob/master/include/obscura/Target_Nucleus.hpp>`_.

^^^^^^^^^^^^^^^^^^^^^
The ``Isotope`` class
^^^^^^^^^^^^^^^^^^^^^

A nuclear isotope is characterized by the number Z of protons and A of nucleons (protons *and* neutrons), its mass, spin, and average spin contribution for protons and neutrons as required e.g. in the context of spin-dependent nuclear interactions.

^^^^^^^^^^^^^^^^^^^^^
The ``Nucleus`` class
^^^^^^^^^^^^^^^^^^^^^

In addition to ``Isotope``, we also define a ``Nucleus`` class which mainly consists of a number of isotopes with given relative abundances.


"""""""""""""""""""""""""""""""
Construction of nuclear targets
"""""""""""""""""""""""""""""""

There are different ways to construct instances of ``Isotope`` and ``Nucleus``.

**Example:** Assume we are interested in oxygen as a target, either the isotope O-16 or the element of various isotopes.

.. code-block:: c++

    #include "obscura/Target_Nucleus.hpp"

    // ...

    // We can define O-16 via the constructor.
    Isotope oxygen_16(8,16);

This instance of an oxygen isotope however has no knowledge of e.g. its spin or relative abundance in nature.

For this purpose, *obscura* contains a nuclear data set, see `/data/Nuclear_Data.txt <https://github.com/temken/obscura/blob/master/data/Nuclear_Data.txt>`_ ([Bednyakov2005]_ [Klos2013]_), which can be accessed through the following function defined in *Target_Nucleus.hpp*.

.. code-block:: c++

   extern Isotope Get_Isotope(unsigned int Z, unsigned int A);
   extern Nucleus Get_Nucleus(unsigned int Z);
   extern Nucleus Get_Nucleus(std::string name);

Using these functions, we can construct isotopes and nuclei simply as

.. code-block:: c++

    #include "obscura/Target_Nucleus.hpp"

    // ...

    Isotope oxygen_16 = Get_Isotope(8,16);
    Nucleus oxygen = Get_Nucleus(8);
    Nucleus oxygen_alternative = Get_Nucleus("O");


The last two lines construct an instance of the ``Nucleus`` class containing all isotopes of oxygen including their relative abundance, spin, and average spin contribution of protons and neutrons.

-------------------------
Electron targets in atoms
-------------------------

For sub-GeV DM searches, an important target are electrons bound in atoms [Essig2012]_.
To take into account the fact that electrons are bound states, we need to evaluate the *ionization form factor* or *atomic response function* for each electronic orbital [Catena2019]_.

The target classes for atomic electrons are declared in `/include/obscura/Target_Atom.hpp <https://github.com/temken/obscura/blob/master/include/obscura/Target_Atom.hpp>`_.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The ``Atomic_Electron`` class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first target class in this context is ``Atomic_Electron``.

By constructing an instance of this class, the tabulated ionization form factor is imported from `/data/Form_Factors_Ionization/ <https://github.com/temken/obscura/tree/master/data/Form_Factors_Ionization>`_.


^^^^^^^^^^^^^^^^^^
The ``Atom`` class
^^^^^^^^^^^^^^^^^^

Having target classes for nuclei and bound electrons, we can combine them into a single atomic target, consisting of a nucleus and a number of bound electrons.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Included ionization form factors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
At this point, *obscura* comes with the ionization form factors of

* Xenon (5p, 5s, 4d, 4p, 4s)
* Argon (3p, 3s, 2p, 2s, 1s)

The tables can be found under `/data/Form_Factors_Ionization/ <https://github.com/temken/obscura/tree/master/data/Form_Factors_Ionization>`_.
They have been tabulated using the `DarkARC <https://github.com/temken/DarkARC>`_ code as described in detail in [Catena2019]_.

The easiest way to access the ionization form factors is by constructing an instance of ``Atom``, as seen in this example.

.. code-block:: c++

    #include "libphysica/Natural_Units.hpp"

    #include "obscura/Target_Atom.hpp"

    using namespace libphysica::natural_units;
    // ... 

    Atom xenon("Xe");
    Atom argon("Ar");

    // For example, to access the ionization form factor of xenon's 5s (quantum numbers n=5, l=0) orbital for a given momentum transfer q and energy E_e:
    int n = 5; int l = 0;
    double q = 0.5 * keV;
    double E_e = 10.0 * eV;

    std::cout << xenon.Electron(n, l).Ionization_Form_Factor(q, E_e) << std::endl;


----------------------------
Electron targets in crystals
----------------------------

One of the most important targets for sub-GeV DM detectors are crystals, such as e.g. semiconductors [Essig2016]_.
The electronic properties of the target material is encapsulated in the crystal form factor which is tabulated and can be found in `/data/Semiconductors/ <https://github.com/temken/obscura/tree/master/data/Form_Factors_Ionization>`_.
The included crystals are

* Silicon semiconductors
* Germanium semiconductors

The tables have been generated using `QEdark <http://ddldm.physics.sunysb.edu/ddlDM/>`_, a module of `Quantum ESPRESSO <https://www.quantum-espresso.org/>`_.

Also for crystals, *obscura* contains a target class ``Crystal`` declared in `/include/obscura/Target_Crystal.hpp <https://github.com/temken/obscura/blob/master/include/obscura/Target_Crystal.hpp>`_.

The crystal form factor, similarly to the ionization form factors, are imported by the class constructor. Here is an example of how to access the crystal form factor.

.. code-block:: c++

  #include "libphysica/Natural_Units.hpp"

  #include "obscura/Target_Crystal.hpp"

  using namespace libphysica::natural_units;

  // ...

  Crystal silicon("Si");
  Crystal germanium("Ge");

  double q = 0.5 * keV;
  double E_e = 10.0 * eV;

  std::cout << silicon.Crystal_Form_Factor(q, E_e) << std::endl;
  std::cout << germanium.Crystal_Form_Factor(q, E_e) << std::endl;