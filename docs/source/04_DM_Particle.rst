.. _Section_DM_Particle:

==============================
4. The ``DM_Particle`` classes
==============================

--------------------------
The interface / base class
--------------------------


.. raw:: html

	<details>
	<summary><a>The DM_Particle base class</a></summary>
 
.. code-block:: c++

   class DM_Particle
   {
     protected:
   	bool low_mass, using_cross_section;

   	// Base class implementations
   	void Print_Summary_Base(int MPI_rank = 0) const;

   	// Some function have an additional function argument 'param' which is not used by the classes included in obscura.
   	// It might be relevant for more complex derived classes to have an additional argument.
   	double Sigma_Total_Nucleus_Base(const Isotope& target, double vDM, double param = -1.0) const;
   	double Sigma_Total_Electron_Base(double vDM, double param = -1.0) const;

   	double PDF_Scattering_Angle_Nucleus_Base(double cos_alpha, const Isotope& target, double vDM, double param = -1.0);
   	double PDF_Scattering_Angle_Electron_Base(double cos_alpha, double vDM, double param = -1.0);
   	double CDF_Scattering_Angle_Nucleus_Base(double cos_alpha, const Isotope& target, double vDM, double param = -1.0);
   	double CDF_Scattering_Angle_Electron_Base(double cos_alpha, double vDM, double param = -1.0);
   	double Sample_Scattering_Angle_Nucleus_Base(std::mt19937& PRNG, const Isotope& target, double vDM, double param = -1.0);
   	double Sample_Scattering_Angle_Electron_Base(std::mt19937& PRNG, double vDM, double param = -1.0);

     public:
   	double mass, spin, fractional_density;
   	bool DD_use_eta_function;

   	//Constructors:
   	DM_Particle();
   	explicit DM_Particle(double m, double s = 1.0 / 2.0);

   	virtual void Set_Mass(double mDM);
   	void Set_Spin(double s);
   	void Set_Low_Mass_Mode(bool ldm);
   	void Set_Fractional_Density(double f);

   	//Primary interaction parameter, such as a coupling constant or cross section
   	virtual double Get_Interaction_Parameter(std::string target) const
   	{
   		return 0.0;
   	};
   	virtual void Set_Interaction_Parameter(double par, std::string target) {};
   	bool Interaction_Parameter_Is_Cross_Section() const;

   	// Set reference cross sections
   	virtual void Set_Sigma_Proton(double sigma) {};
   	virtual void Set_Sigma_Neutron(double sigma) {};
   	virtual void Set_Sigma_Electron(double sigma) {};

   	//Differential cross sections for nuclear targets
   	virtual double dSigma_dq2_Nucleus(double q, const Isotope& target, double vDM, double param = -1.0) const { return 0.0; };
   	double dSigma_dER_Nucleus(double ER, const Isotope& target, double vDM, double param = -1.0) const;
   	double d2Sigma_dER_dEe_Migdal(double ER, double Ee, double vDM, const Isotope& isotope, Atomic_Electron& shell) const;

   	// Differential cross section for electron targets
   	virtual double dSigma_dq2_Electron(double q, double vDM, double param = -1.0) const { return 0.0; };
   	virtual double d2Sigma_dq2_dEe_Ionization(double q, double Ee, double vDM, Atomic_Electron& shell) const { return 0.0; };
   	virtual double d2Sigma_dq2_dEe_Crystal(double q, double Ee, double vDM, Crystal& crystal) const { return 0.0; };

   	// Reference cross sections
   	virtual double Sigma_Proton() const { return 0.0; };
   	virtual double Sigma_Neutron() const { return 0.0; };
   	virtual double Sigma_Electron() const { return 0.0; };

   	virtual double Sigma_Total_Nucleus(const Isotope& target, double vDM, double param = -1.0) const;
   	virtual double Sigma_Total_Electron(double vDM, double param = -1.0) const;

   	virtual void Print_Summary(int MPI_rank = 0) const;

   	// Scattering angle functions
   	virtual double PDF_Scattering_Angle_Nucleus(double cos_alpha, const Isotope& target, double vDM, double param = -1.0);
   	virtual double PDF_Scattering_Angle_Electron(double cos_alpha, double vDM, double param = -1.0);
   	virtual double CDF_Scattering_Angle_Nucleus(double cos_alpha, const Isotope& target, double vDM, double param = -1.0);
   	virtual double CDF_Scattering_Angle_Electron(double cos_alpha, double vDM, double param = -1.0);
   	virtual double Sample_Scattering_Angle_Nucleus(std::mt19937& PRNG, const Isotope& target, double vDM, double param = -1.0);
   	virtual double Sample_Scattering_Angle_Electron(std::mt19937& PRNG, double vDM, double param = -1.0);
   };

.. raw:: html

	</details>

----------------------------------
Spin-Independent (SI) interactions
----------------------------------


--------------------------------
Spin-Dependent (SD) interactions
--------------------------------