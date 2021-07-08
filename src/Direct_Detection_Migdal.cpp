#include "obscura/Direct_Detection_Migdal.hpp"

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"
#include "libphysica/Utilities.hpp"

namespace obscura
{
using namespace libphysica::natural_units;

double vMinimal_Migdal(double ER, double Ee, double mDM, double mN)
{
	//eq 5 of https://arxiv.org/pdf/1711.09906.pdf
	double mu = libphysica::Reduced_Mass(mDM, mN);
	return sqrt(mN * ER / 2.0 / mu / mu) + Ee / sqrt(2.0 * mN * ER);
}

double dRdEe_Ionization_Migdal(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, const Isotope& isotope, Atomic_Electron& shell)
{
	double NT = 1.0 / isotope.mass;

	std::function<double(double)> ER_integrand = [Ee, &DM, &DM_distr, &isotope, &shell](double ER) {
		double vMin = (ER > 0) ? vMinimal_Migdal(ER, Ee, DM.mass, isotope.mass) : 1.0e-20 * cm / sec;
		if(vMin >= DM_distr.Maximum_DM_Speed())
			return 0.0;
		if(DM.DD_use_eta_function && DM_distr.DD_use_eta_function)
		{
			double vDM = 1.0e-3;   // cancels
			return DM_distr.DM_density / DM.mass * DM_distr.Eta_Function(vMin) * DM.d2Sigma_dER_dEe_Migdal(ER, Ee, vDM, isotope, shell) * vDM * vDM;
		}
		else
		{
			std::function<double(double)> v_integrand = [ER, Ee, &DM_distr, &DM, &isotope, &shell](double v) {
				return DM_distr.Differential_DM_Flux(v, DM.mass) * DM.d2Sigma_dER_dEe_Migdal(ER, Ee, v, isotope, shell);
			};
			double v_integral = libphysica::Integrate(v_integrand, vMin, DM_distr.Maximum_DM_Speed());
			return v_integral;
		}
	};

	// Integral over nuclear recoil energies
	double vMax	  = DM_distr.Maximum_DM_Speed();
	double qMin	  = shell.binding_energy / vMax;
	double qMax	  = 2.0 * libphysica::Reduced_Mass(DM.mass, isotope.mass) * vMax;
	double ER_min = qMin * qMin / 2.0 / isotope.mass;
	double ER_max = qMax * qMax / 2.0 / isotope.mass;

	return NT * libphysica::Integrate(ER_integrand, ER_min, ER_max);
}

extern double dRdEe_Ionization_Migdal(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, const Nucleus& nucleus, Atomic_Electron& shell)
{
	double result = 0.0;
	for(auto& isotope : nucleus.isotopes)
		result += isotope.abundance * dRdEe_Ionization_Migdal(Ee, DM, DM_distr, isotope, shell);
	return result;
}

double dRdEe_Ionization_Migdal(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, Atom& atom)
{
	double result = 0.0;
	for(auto& electron : atom.electrons)
		result += dRdEe_Ionization_Migdal(Ee, DM, DM_distr, atom.nucleus, electron);
	return result;
}

DM_Detector_Ionization_Migdal::DM_Detector_Ionization_Migdal()
: DM_Detector_Ionization("Migdal scattering experiment", kg * day, "Nuclei", "Xe") {}
DM_Detector_Ionization_Migdal::DM_Detector_Ionization_Migdal(std::string label, double expo, std::string atom)
: DM_Detector_Ionization(label, expo, "Nuclei", atom) {}
DM_Detector_Ionization_Migdal::DM_Detector_Ionization_Migdal(std::string label, double expo, std::vector<std::string> atoms, std::vector<double> mass_fractions)
: DM_Detector_Ionization(label, expo, "Nuclei", atoms, mass_fractions) {}

double DM_Detector_Ionization_Migdal::dRdE_Ionization(double E, const DM_Particle& DM, DM_Distribution& DM_distr, const Nucleus& nucleus, Atomic_Electron& shell)
{
	return flat_efficiency * dRdEe_Ionization_Migdal(E, DM, DM_distr, nucleus, shell);
}

}	// namespace obscura
