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

For nuclear recoil experiments, we define two target classes, ``Isotope`` and ``Nucleus``.

^^^^^^^^^^^
``Isotope``
^^^^^^^^^^^

A nuclear isotope is characterized by the number Z of protons and A of nucleons (protons *and* neutrons), its mass, spin, and average spin contribution for protons and neutrons as required e.g. in the context of spin-dependent nuclear interactions.

.. raw:: html

	<details>
	<summary><a>The class</a></summary>
 
.. code-block:: c++

   struct Isotope
   {
   	unsigned int Z, A;
   	double abundance;
   	double spin;
   	double sp, sn;

   	std::string name;
   	double mass;

   	Isotope();
   	Isotope(unsigned int z, unsigned int a, double abund = 1.0, double Spin = 0.0, double Sp = 0.0, double Sn = 0.0);

   	double Thomas_Fermi_Radius() const;

   	//Nuclear form factor for SI interactions
   	double Helm_Form_Factor(double q) const;

   	//Nuclear form factor for SD interactions
   	// to do

   	void Print_Summary(unsigned int MPI_rank) const;
   };

.. raw:: html

	</details>

^^^^^^^^^^^
``Nucleus``
^^^^^^^^^^^

In addition to ``Isotope``, we also define a ``Nucleus`` class which mainly consists of a number of isotopes with given relative abundances.

.. raw:: html

	<details>
	<summary><a>The class</a></summary>
 
.. code-block:: c++

   struct Nucleus
   {
   	int Z;
   	std::vector<Isotope> isotopes;
   	std::string name;

   	Nucleus();
   	Nucleus(const std::vector<Isotope>& iso);
   	Nucleus(const Isotope& iso);

   	unsigned int Number_of_Isotopes() const;

   	Isotope Get_Isotope(unsigned int A) const;

   	Isotope& operator[](int i)
   	{
   		return isotopes[i];
   	}
   	const Isotope& operator[](int i) const
   	{
   		return isotopes[i];
   	}

   	double Average_Nuclear_Mass() const;

   	void Print_Summary(unsigned int MPI_rank = 0) const;
   };

.. raw:: html

	</details>


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


----------------------------
Electron targets in crystals
----------------------------