#include "DM_Particle_Standard.hpp"

#include <cmath>

//Headers from libphysica library
#include "Numerics.hpp"

namespace obscura
{
using namespace libphysica::natural_units;

//1. Abstract parent class for SI and SD interactions
DM_Particle_Standard::DM_Particle_Standard()
: DM_Particle(), prefactor(1.0), fixed_coupling_relation(true), fp_relative(0.5), fn_relative(0.5), sigma_electron(0.0)
{
	DD_use_eta_function = true;
}

DM_Particle_Standard::DM_Particle_Standard(double mDM, double pre)
: DM_Particle(mDM), prefactor(pre), fixed_coupling_relation(true), fp_relative(0.5), fn_relative(0.5), sigma_electron(0.0)
{
	DD_use_eta_function = true;
}

void DM_Particle_Standard::Set_Sigma_Proton(double sigma)
{
	fp = sqrt(M_PI * sigma / prefactor) / libphysica::Reduced_Mass(mass, mProton);
	if(fixed_coupling_relation)
	{
		if(sigma == 0.0 && fp_relative != 0.0)
		{
			std::cerr << "Error in obscura::DM_Particle_Standard::Set_Sigma_Proton(): fp is set to zero while fixed fn/fp." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(sigma != 0.0 && fp_relative == 0.0)
		{
			std::cerr << "Error in obscura::DM_Particle_Standard::Set_Sigma_Proton(): fp is set while fp/fn is fixed to zero." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
			fn = fn_relative / fp_relative * fp;
	}
}

//Primary interaction parameter, in this case the proton or neutron cross section
double DM_Particle_Standard::Get_Interaction_Parameter(std::string target) const
{
	if(target == "Nuclei")
	{
		if(fp_relative > 0.0)
			return Sigma_Proton();
		else
			return Sigma_Neutron();
	}
	else if(target == "Electrons")
		return Sigma_Electron();
	else
	{
		std::cerr << "Error in obscura::DM_Particle_Standard::Get_Interaction_Parameter(std::string): Target " << target << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void DM_Particle_Standard::Set_Interaction_Parameter(double par, std::string target)
{
	if(target == "Nuclei")
	{
		if(fp_relative > 0.0)
			Set_Sigma_Proton(par);
		else
			Set_Sigma_Neutron(par);
	}
	else if(target == "Electrons")
		Set_Sigma_Electron(par);
	else
	{
		std::cerr << "Error in obscura::DM_Particle_Standard::Get_Interaction_Parameter(std::string): Target " << target << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void DM_Particle_Standard::Set_Sigma_Neutron(double sigma)
{
	fn = sqrt(M_PI * sigma / prefactor) / libphysica::Reduced_Mass(mass, mProton);
	if(fixed_coupling_relation)
	{
		if(sigma == 0.0 && fn_relative != 0.0)
		{
			std::cerr << "Error in obscura::DM_Particle_Standard::Set_Sigma_Neutron(): fn is set to zero while fixed fn/fp." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(sigma != 0.0 && fn_relative == 0.0)
		{
			std::cerr << "Error in obscura::DM_Particle_Standard::Set_Sigma_Neutron(): fn is set while fn/fp is fixed to zero." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
			fp = fp_relative / fn_relative * fn;
	}
}

void DM_Particle_Standard::Set_Sigma_Electron(double sigma)
{
	sigma_electron = sigma;
}

void DM_Particle_Standard::Fix_Coupling_Ratio(double fp_rel, double fn_rel)
{
	fixed_coupling_relation = true;

	double tot	= fp_rel + fn_rel;
	fp_relative = fp_rel / tot;
	fn_relative = fn_rel / tot;
	if(fp != 0.0)
		fn = fn_relative / fp_relative * fp;
	else if(fn != 0.0)
		fp = fp_relative / fn_relative * fp;
	else
	{
		std::cerr << "Error in obscura::DM_Particle_Standard::Fix_Coupling_Ratio(double, double): Both couplings zero." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void DM_Particle_Standard::Fix_fn_over_fp(double ratio)
{
	fixed_coupling_relation = true;

	fp_relative = 1.0 / (1.0 + ratio);
	fn_relative = ratio / (1.0 + ratio);
	if(fp != 0.0)
		fn = fn_relative / fp_relative * fp;
	else if(fn != 0.0)
		fp = fp_relative / fn_relative * fp;
	else
	{
		std::cerr << "Error in obscura::DM_Particle_Standard::Fix_fn_over_fp(double): Both couplings zero." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void DM_Particle_Standard::Fix_fp_over_fn(double ratio)
{
	fixed_coupling_relation = true;

	fp_relative = ratio;
	fn_relative = 1.0;
	if(fp != 0.0)
		fn = fn_relative / fp_relative * fp;
	else if(fn != 0.0)
		fp = fp_relative / fn_relative * fp;
	else
	{
		std::cerr << "Error in obscura::DM_Particle_Standard::Fix_fp_over_fn(double): Both couplings zero." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void DM_Particle_Standard::Unfix_Coupling_Ratios()
{
	fixed_coupling_relation = false;
}

//Reference cross sections
double DM_Particle_Standard::Sigma_Proton() const
{
	return prefactor * fp * fp * libphysica::Reduced_Mass(mass, mProton) * libphysica::Reduced_Mass(mass, mProton) / M_PI;
}

double DM_Particle_Standard::Sigma_Neutron() const
{
	return prefactor * fn * fn * libphysica::Reduced_Mass(mass, mProton) * libphysica::Reduced_Mass(mass, mProton) / M_PI;
}

double DM_Particle_Standard::Sigma_Electron() const
{
	return sigma_electron;
}

void DM_Particle_Standard::Print_Summary_Standard(int MPI_rank) const
{
	if(MPI_rank == 0)
	{
		std::cout << "\tCoupling ratio fixed:\t" << (fixed_coupling_relation ? "[x]" : "[ ]") << std::endl
				  << "\tIsospin conservation:\t" << ((fn == fp) ? "[x]" : "[ ]") << std::endl;
		if(fixed_coupling_relation)
		{
			std::cout << "\tCoupling ratio:\t\t" << ((fp != 0.0) ? "fn/fp = " : "fp/fn = ") << ((fp != 0.0) ? libphysica::Round(fn / fp) : libphysica::Round(fp / fn)) << std::endl
					  << std::endl;
		}
		std::cout << "\tSigma_P[cm^2]:\t\t" << libphysica::Round(In_Units(Sigma_Proton(), cm * cm)) << std::endl
				  << "\tSigma_N[cm^2]:\t\t" << libphysica::Round(In_Units(Sigma_Neutron(), cm * cm)) << std::endl
				  << "\tSigma_E[cm^2]:\t\t" << libphysica::Round(In_Units(Sigma_Electron(), cm * cm)) << std::endl
				  << std::endl;
	}
}

//2. Spin-independent (SI) interactions
//Constructors:
DM_Particle_SI::DM_Particle_SI()
: DM_Particle_Standard(), FF_DM("Contact"), mMediator(0.0)
{
	qRef = aEM * mElectron;
	Set_Sigma_Proton(1e-40 * cm * cm);
}

DM_Particle_SI::DM_Particle_SI(double mDM)
: DM_Particle_Standard(mDM, 1.0), FF_DM("Contact"), mMediator(0.0)
{
	qRef = aEM * mElectron;
	Set_Sigma_Proton(1e-40 * cm * cm);
}

DM_Particle_SI::DM_Particle_SI(double mDM, double sigmaP)
: DM_Particle_Standard(mDM, 1.0), FF_DM("Contact"), mMediator(0.0)
{
	qRef = aEM * mElectron;
	Set_Sigma_Proton(sigmaP);
}

void DM_Particle_SI::Set_FormFactor_DM(std::string ff, double mMed)
{
	if(ff == "Contact" || ff == "Electric-Dipole" || ff == "Long-Range" || ff == "General")
		FF_DM = ff;
	else
	{
		std::cerr << "Error in obscura::DM_Particle_SI::Set_FormFactor_DM(): Form factor " << ff << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if(FF_DM == "General" && mMed > 0.0)
		mMediator = mMed;
}

void DM_Particle_SI::Set_Mediator_Mass(double m)
{
	mMediator = m;
}

//DM form factir
double DM_Particle_SI::FormFactor2_DM(double q) const
{
	double FF;
	if(FF_DM == "Contact")
		FF = 1.0;
	else if(FF_DM == "General")
		FF = (qRef * qRef + mMediator * mMediator) / (q * q + mMediator * mMediator);
	else if(FF_DM == "Long-Range")
		FF = qRef * qRef / q / q;
	else if(FF_DM == "Electric-Dipole")
		FF = qRef / q;
	else
	{
		std::cerr << "Error in obscura::DM_Particle_SI::FormFactor2_DM(): Form factor " << FF_DM << "not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	return FF * FF;
}

//Differential Cross Sections
double DM_Particle_SI::dSigma_dq2_Nucleus(double q, const Isotope& target, double vDM) const
{
	double nuclear_form_factor = (low_mass) ? 1.0 : target.Helm_Form_Factor(q);
	return 1.0 / 4.0 / M_PI / vDM / vDM * pow((fp * target.Z + fn * (target.A - target.Z)), 2.0) * FormFactor2_DM(q) * nuclear_form_factor * nuclear_form_factor;
}

double DM_Particle_SI::dSigma_dq2_Electron(double q, double vDM) const
{
	return sigma_electron / pow(2.0 * libphysica::Reduced_Mass(mass, mElectron) * vDM, 2.0) * FormFactor2_DM(q);
}

//Total cross sections
double DM_Particle_SI::Sigma_Nucleus(const Isotope& isotope, double vDM) const
{
	double sigmatot = 0.0;
	if(FF_DM != "Contact" && FF_DM != "General")
	{
		std::cerr << "Error in obscura::DM_Particle_SI::Sigma_Nucleus(): Divergence in the IR." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(!low_mass)
		sigmatot = Sigma_Nucleus_Base(isotope, vDM);
	else
	{
		sigmatot = pow(libphysica::Reduced_Mass(mass, isotope.mass), 2.0) / M_PI * pow(fp * isotope.Z + fn * (isotope.A - isotope.Z), 2.0);
		if(FF_DM == "General")
		{
			double q2max = 4.0 * pow(libphysica::Reduced_Mass(mass, isotope.mass) * vDM, 2.0);
			sigmatot *= pow(qRef * qRef + mMediator * mMediator, 2.0) / mMediator / mMediator / (mMediator * mMediator + q2max);
		}
	}
	return sigmatot;
}

void DM_Particle_SI::Print_Summary(int MPI_rank) const
{
	if(MPI_rank == 0)
	{
		Print_Summary_Base();
		std::cout << std::endl
				  << "\tInteraction:\t\tSpin-Independent (SI)" << std::endl
				  << std::endl;
		Print_Summary_Standard();
		std::cout << "\tInteraction type:\t" << FF_DM << std::endl;
		if(FF_DM == "General")
		{
			double massunit			= (mMediator < keV) ? eV : ((mMediator < MeV) ? keV : ((mMediator < GeV) ? MeV : GeV));
			std::string massunitstr = (mMediator < keV) ? "eV" : ((mMediator < MeV) ? "keV" : ((mMediator < GeV) ? "MeV" : "GeV"));
			std::cout << "\t\tMediator mass:\t" << In_Units(mMediator, massunit) << " " << massunitstr << std::endl;
		}
		std::cout << "----------------------------------------" << std::endl;
	}
}

//3. Spin-dependent (SD) interactions
//Constructors:
DM_Particle_SD::DM_Particle_SD()
: DM_Particle_Standard(1.0 * GeV, 3.0)
{
	Set_Sigma_Proton(1e-40 * cm * cm);
}

DM_Particle_SD::DM_Particle_SD(double mDM)
: DM_Particle_Standard(mDM, 3.0)
{
	Set_Sigma_Proton(1e-40 * cm * cm);
}

DM_Particle_SD::DM_Particle_SD(double mDM, double sigmaP)
: DM_Particle_Standard(mDM, 3.0)
{
	Set_Sigma_Proton(sigmaP);
}

//Differential Cross Sections
double DM_Particle_SD::dSigma_dq2_Nucleus(double q, const Isotope& target, double vDM) const
{
	return (target.spin == 0) ? 0.0 : 1.0 / M_PI / vDM / vDM * (target.spin + 1) / target.spin * pow((fp * target.sp + fn * target.sn), 2);
}

double DM_Particle_SD::dSigma_dq2_Electron(double q, double vDM) const
{
	return sigma_electron / pow(2.0 * libphysica::Reduced_Mass(mass, mElectron) * vDM, 2.0);
}

//Total cross sections with nuclear isotopes, elements, and electrons
double DM_Particle_SD::Sigma_Nucleus(const Isotope& isotope, double vDM) const
{
	return (isotope.spin != 0) ? 4.0 * pow(libphysica::Reduced_Mass(mass, isotope.mass), 2.0) / M_PI * (isotope.spin + 1.0) / isotope.spin * pow(fp * isotope.sp + fn * isotope.sn, 2.0) : 0.0;
}

void DM_Particle_SD::Print_Summary(int MPI_rank) const
{
	if(MPI_rank == 0)
	{
		Print_Summary_Base(MPI_rank);
		std::cout << "Interaction:\t\tSpin-Dependent (SD)" << std::endl
				  << std::endl;
		Print_Summary_Standard();
		std::cout << "----------------------------------------" << std::endl
				  << std::endl;
	}
}

// 4. Dark photon model with kinetic mixing
// Constructors
DM_Particle_DP::DM_Particle_DP()
: DM_Particle(), alpha_dark(aEM), q_reference(aEM * mElectron), FF_DM("Contact"), m_dark_photon(GeV)
{
	DD_use_eta_function = true;
	Set_Sigma_Proton(1.0e-40 * cm * cm);
}

DM_Particle_DP::DM_Particle_DP(double mDM)
: DM_Particle(mDM), alpha_dark(aEM), q_reference(aEM * mElectron), FF_DM("Contact"), m_dark_photon(GeV)
{
	DD_use_eta_function = true;
	Set_Sigma_Proton(1.0e-40 * cm * cm);
}

DM_Particle_DP::DM_Particle_DP(double mDM, double sigma_p)
: DM_Particle(mDM), alpha_dark(aEM), q_reference(aEM * mElectron), FF_DM("Contact"), m_dark_photon(GeV)
{
	DD_use_eta_function = true;
	Set_Sigma_Proton(sigma_p);
}

double DM_Particle_DP::FormFactor2_DM(double q) const
{
	double FF;
	if(FF_DM == "Contact")
		FF = 1.0;
	else if(FF_DM == "General")
		FF = (q_reference * q_reference + m_dark_photon * m_dark_photon) / (q * q + m_dark_photon * m_dark_photon);
	else if(FF_DM == "Long-Range")
		FF = q_reference * q_reference / q / q;
	else if(FF_DM == "Electric-Dipole")
		FF = q_reference / q;
	else
	{
		std::cerr << "Error in obscura::DM_Particle_DP::FormFactor2_DM(): Form factor " << FF_DM << "not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	return FF * FF;
}

void DM_Particle_DP::Set_Mass(double mDM)
{
	double sigma_p = Sigma_Proton();
	mass		   = mDM;
	Set_Sigma_Proton(sigma_p);
}

void DM_Particle_DP::Set_Epsilon(double e)
{
	epsilon = e;
}

double DM_Particle_DP::Get_Epsilon()
{
	return epsilon;
}

// Dark matter form factor
void DM_Particle_DP::Set_FormFactor_DM(std::string ff, double mMed)
{
	if(ff == "Contact" || ff == "Electric-Dipole" || ff == "Long-Range" || ff == "General")
		FF_DM = ff;

	else
	{
		std::cerr << "Error in obscura::DM_Particle_DP::Set_FormFactor_DM(): Form factor " << ff << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if(FF_DM == "General" && mMed > 0.0)
		m_dark_photon = mMed;
	else if(FF_DM == "Long-Range")
		m_dark_photon = 0.0;
}

void DM_Particle_DP::Set_Dark_Photon_Mass(double m)
{
	if(FF_DM != "Long-Range")
		m_dark_photon = m;
}

// Primary interaction parameter, such as a coupling constant or cross section
double DM_Particle_DP::Get_Interaction_Parameter(std::string target) const
{
	if(target == "Nuclei")
		return Sigma_Proton();
	else if(target == "Electrons")
		return Sigma_Electron();
	else
	{
		std::cerr << "Error in obscura::DM_Particle_DP::Get_Interaction_Parameter(std::string): Target " << target << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void DM_Particle_DP::Set_Interaction_Parameter(double par, std::string target)
{
	if(target == "Nuclei")
		Set_Sigma_Proton(par);
	else if(target == "Electrons")
		Set_Sigma_Electron(par);
	else
	{
		std::cerr << "Error in obscura::DM_Particle_DP::Get_Interaction_Parameter(std::string): Target " << target << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void DM_Particle_DP::Set_Sigma_Proton(double sigma)
{
	epsilon = (q_reference * q_reference + m_dark_photon * m_dark_photon) / 4.0 / libphysica::Reduced_Mass(mass, mProton) * sqrt(sigma / M_PI / aEM / alpha_dark);
}

void DM_Particle_DP::Set_Sigma_Electron(double sigma)
{
	epsilon = (q_reference * q_reference + m_dark_photon * m_dark_photon) / 4.0 / libphysica::Reduced_Mass(mass, mElectron) * sqrt(sigma / M_PI / aEM / alpha_dark);
}

// Differential cross sections with nuclear isotopes, elements, and electrons
double DM_Particle_DP::dSigma_dq2_Nucleus(double q, const Isotope& target, double vDM) const
{
	double nuclear_form_factor = (low_mass) ? 1.0 : target.Helm_Form_Factor(q);
	double mu				   = libphysica::Reduced_Mass(mass, mProton);
	return Sigma_Proton() / 4.0 / mu / mu / vDM / vDM * FormFactor2_DM(q) * nuclear_form_factor * target.Z * target.Z;
}

double DM_Particle_DP::dSigma_dq2_Electron(double q, double vDM) const
{
	double mu = libphysica::Reduced_Mass(mass, mElectron);
	return Sigma_Electron() / 4.0 / mu / mu / vDM / vDM * FormFactor2_DM(q);
}

// Total cross sections with nuclear isotopes, elements, and electrons
double DM_Particle_DP::Sigma_Proton() const
{
	double mu = libphysica::Reduced_Mass(mass, mProton);
	return 16.0 * M_PI * aEM * alpha_dark * epsilon * epsilon * mu * mu / pow((q_reference * q_reference + m_dark_photon * m_dark_photon), 2.0);
}

double DM_Particle_DP::Sigma_Electron() const
{
	double mu = libphysica::Reduced_Mass(mass, mElectron);
	return 16.0 * M_PI * aEM * alpha_dark * epsilon * epsilon * mu * mu / pow((q_reference * q_reference + m_dark_photon * m_dark_photon), 2.0);
}

double DM_Particle_DP::Sigma_Nucleus(const Isotope& target, double vDM) const
{
	double sigmatot = 0.0;
	if(FF_DM != "Contact" && FF_DM != "General")
	{
		std::cerr << "Error in obscura::DM_Particle_DP::Sigma_Nucleus(): Divergence in the IR." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(!low_mass)
		sigmatot = Sigma_Nucleus_Base(target, vDM);
	else
	{
		double mu_p = libphysica::Reduced_Mass(mass, mProton);
		double mu_N = libphysica::Reduced_Mass(mass, target.mass);
		sigmatot	= Sigma_Proton() * mu_N * mu_N / mu_p / mu_p * target.Z * target.Z;
		if(FF_DM == "General")
		{
			double q2max = 4.0 * pow(mu_N * vDM, 2.0);
			sigmatot *= pow(q_reference * q_reference + m_dark_photon * m_dark_photon, 2.0) / m_dark_photon / m_dark_photon / (m_dark_photon * m_dark_photon + q2max);
		}
	}
	return sigmatot;
}

void DM_Particle_DP::Print_Summary(int MPI_rank) const
{
	if(MPI_rank == 0)
	{
		Print_Summary_Base();
		double massunit			= (m_dark_photon < keV) ? eV : ((m_dark_photon < MeV) ? keV : ((m_dark_photon < GeV) ? MeV : GeV));
		std::string massunitstr = (m_dark_photon < keV) ? "[eV]" : ((m_dark_photon < MeV) ? "[keV]" : ((m_dark_photon < GeV) ? "[MeV]" : "[GeV]"));
		std::cout << std::endl
				  << "\tInteraction:\t\tDark photon (DP)" << std::endl
				  << std::endl;
		std::cout << "\tInteraction type:\t" << FF_DM << std::endl
				  << "\tKinetic mixing:\t\t" << libphysica::Round(epsilon) << std::endl
				  << "\tGauge coupling a_D:\t" << libphysica::Round(alpha_dark) << std::endl
				  << "\tDark photon mass" << massunitstr << ":\t" << libphysica::Round(In_Units(m_dark_photon, massunit)) << std::endl
				  << "\tSigma_P[cm^2]:\t\t" << libphysica::Round(In_Units(Sigma_Proton(), cm * cm)) << std::endl
				  << "\tSigma_N[cm^2]:\t\t" << libphysica::Round(In_Units(Sigma_Neutron(), cm * cm)) << std::endl
				  << "\tSigma_E[cm^2]:\t\t" << libphysica::Round(In_Units(Sigma_Electron(), cm * cm)) << std::endl
				  << std::endl;

		std::cout
			<< "----------------------------------------" << std::endl;
	}
}

}	// namespace obscura