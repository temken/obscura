#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>

#include "libphysica/Integration.hpp"
#include "libphysica/Linear_Algebra.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/Configuration.hpp"
#include "obscura/DM_Particle_Standard.hpp"
#include "obscura/Direct_Detection_ER.hpp"
#include "version.hpp"

using namespace libphysica::natural_units;
using namespace obscura;

int main(int argc, char* argv[])
{
	// Initial terminal output
	auto time_start	  = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	auto* ctime_start = ctime(&time_start_t);
	if(ctime_start[std::strlen(ctime_start) - 1] == '\n')
		ctime_start[std::strlen(ctime_start) - 1] = '\0';
	std::cout << "[Started on " << ctime_start << "]" << std::endl;
	std::cout << PROJECT_NAME << "-v" << PROJECT_VERSION << "\tgit:" << GIT_BRANCH << "/" << GIT_COMMIT_HASH << std::endl
			  << OBSCURA_LOGO
			  << std::endl;

	// Import configuration file
	obscura::Configuration cfg(argv[1]);
	cfg.Print_Summary();
	////////////////////////////////////////////////////////////////////////

	// Import response matrix
	std::vector<std::vector<double>> data = libphysica::Import_Table(TOP_LEVEL_DIR "data/XENON1Te/s2_response_er.dat");
	libphysica::Matrix matrix(data);
	matrix = matrix.Transpose();
	std::cout << matrix.Rows() << " x " << matrix.Columns() << std::endl;

	std::vector<std::vector<double>> energy_bins = libphysica::Import_Table(TOP_LEVEL_DIR "data/XENON1Te/s2_response_er_energybins.dat", {keV, keV, keV});
	std::cout << energy_bins.size() << std::endl;

	// Compute energy spectrum vector
	Atom Xe("Xe");
	std::vector<double> dRdE_list_1, dRdE_list_2;
	for(unsigned int i = 0; i < energy_bins.size(); i++)
	{
		double E  = energy_bins[i][0];
		double dE = energy_bins[i][2] - energy_bins[i][1];
		std::cout << "dE " << dE << std::endl;
		dRdE_list_1.push_back(dE * dRdEe_Ionization_ER(E, *(cfg.DM), *(cfg.DM_distr), Xe));

		double E1						   = energy_bins[i][1];
		double E2						   = energy_bins[i][2];
		std::function<double(double)> drde = [&cfg, &Xe](double e) {
			return dRdEe_Ionization_ER(e, *(cfg.DM), *(cfg.DM_distr), Xe);
		};
		dRdE_list_2.push_back(libphysica::Integrate(drde, E1, E2));
	}
	for(unsigned int i = 0; i < dRdE_list_1.size(); i++)
		std::cout << i << "\t" << dRdE_list_1[i] << "\t" << dRdE_list_2[i] << "\tratio = " << dRdE_list_1[i] / dRdE_list_2[i] << std::endl;

	libphysica::Vector dRdE(dRdE_list_2);
	std::cout << dRdE.Size() << std::endl;

	libphysica::Vector RS2 = matrix * dRdE;
	std::cout << RS2.Size() << std::endl;
	for(unsigned int i = 0; i < RS2.Size(); i++)
		std::cout << i << "\t" << RS2[i] * kg * year << std::endl;
	////////////////////////////////////////////////////////////////////////
	// Atom Xe("Xe");
	// double q = 23.0 * keV;
	// double E = 19.0 * eV;
	// std::cout << Xe.Electron(5, 1).Atomic_Response_Function(1, q, E) << std::endl;
	// std::cout << Xe.Electron(5, 1).Ionization_Form_Factor(q, E) << std::endl;
	// std::vector<double> DM_masses = libphysica::Log_Space(cfg.constraints_mass_min, cfg.constraints_mass_max, cfg.constraints_masses);

	// std::vector<std::vector<double>> exclusion_limits = cfg.DM_detector->Upper_Limit_Curve(*(cfg.DM), *(cfg.DM_distr), DM_masses, cfg.constraints_certainty);
	// for(unsigned int i = 0; i < exclusion_limits.size(); i++)
	// 	std::cout << i + 1 << "/" << exclusion_limits.size()
	// 			  << "\tmDM = " << libphysica::Round(In_Units(exclusion_limits[i][0], (exclusion_limits[i][0] < GeV) ? MeV : GeV)) << ((exclusion_limits[i][0] < GeV) ? " MeV" : " GeV")
	// 			  << "\tUpper Bound:\t" << libphysica::Round(In_Units(exclusion_limits[i][1], cm * cm)) << std::endl;

	// int CL = std::round(100.0 * cfg.constraints_certainty);
	// libphysica::Export_Table(TOP_LEVEL_DIR "results/" + cfg.ID + "/DD_Constraints_" + std::to_string(CL) + ".txt", exclusion_limits, {GeV, cm * cm});

	// DM_Particle_SI DM(100.0 * MeV);
	// DM.Set_Sigma_Electron(pb);
	// auto energies = libphysica::Log_Space(0.1 * eV, 400 * eV, 50);
	// std::vector<std::vector<double>> spectrum;
	// for(auto E : energies)
	// 	spectrum.push_back({E, E * dRdEe_Ionization_ER(E, DM, *cfg.DM_distr, Xe)});
	// libphysica::Export_Table("dRdlnEe_mX_100MeV_sigma_1pb.dat", spectrum, {eV, 1.0 / kg / year}, "# E [eV]\tdR/dln(Ee) [1/kg/year]");
	// extern double dRdEe_Ionization_ER(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, Atom& atom);

	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << std::round(1000. * durationTotal) / 1000. << "s";
	if(durationTotal > 60.0)
		std::cout << " (" << floor(durationTotal / 3600.0) << ":" << floor(fmod(durationTotal / 60.0, 60.0)) << ":" << floor(fmod(durationTotal, 60.0)) << ")]." << std::endl;
	else
		std::cout << "]" << std::endl;

	return 0;
}