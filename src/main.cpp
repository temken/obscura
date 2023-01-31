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

	std::vector<double> DM_masses = libphysica::Log_Space(cfg.constraints_mass_min, cfg.constraints_mass_max, cfg.constraints_masses);

	std::vector<std::vector<double>> exclusion_limits = cfg.DM_detector->Upper_Limit_Curve(*(cfg.DM), *(cfg.DM_distr), DM_masses, cfg.constraints_certainty);
	for(unsigned int i = 0; i < exclusion_limits.size(); i++)
		std::cout << i + 1 << "/" << exclusion_limits.size()
				  << "\tmDM = " << libphysica::Round(In_Units(exclusion_limits[i][0], (exclusion_limits[i][0] < GeV) ? MeV : GeV)) << ((exclusion_limits[i][0] < GeV) ? " MeV" : " GeV")
				  << "\tUpper Bound:\t" << libphysica::Round(In_Units(exclusion_limits[i][1], cm * cm)) << std::endl;

	int CL = std::round(100.0 * cfg.constraints_certainty);
	libphysica::Export_Table(TOP_LEVEL_DIR "results/" + cfg.ID + "/DD_Constraints_" + std::to_string(CL) + "_Response_Matrix_4bins+background.txt", exclusion_limits, {GeV, cm * cm});

	// // Import response matrix
	// std::vector<std::vector<double>> data = libphysica::Import_Table(TOP_LEVEL_DIR "data/XENON1T_S2/s2_response_er.dat");
	// libphysica::Matrix matrix(data);
	// matrix = matrix.Transpose();
	// std::cout << matrix.Rows() << " x " << matrix.Columns() << std::endl;

	// std::vector<std::vector<double>> energy_bins = libphysica::Import_Table(TOP_LEVEL_DIR "data/XENON1T_S2/s2_response_er_energybins.dat", {keV, keV, keV});
	// std::cout << energy_bins.size() << std::endl;
	// std::vector<std::vector<double>> s2_bin_info = libphysica::Import_Table(TOP_LEVEL_DIR "data/XENON1T_S2/s2_binning_info.dat", {}, 1);

	// // Compute energy spectrum vector
	// Atom Xe("Xe");
	// std::vector<double> dRdE_list(energy_bins.size(), 0.0);

	// std::function<double(double)> drde = [&cfg, &Xe](double e) {
	// 	return dRdEe_Ionization_ER(e, *(cfg.DM), *(cfg.DM_distr), Xe);
	// };

	// std::vector<double> E_list = libphysica::Log_Space(0.01 * eV, energy_bins[energy_bins.size() - 1][2], 500);
	// libphysica::Export_Function("dRdE.dat", drde, E_list, {keV, 1.0 / keV / kg / year}, "# E [keV]\tdR/dE [keV^-1 kg^-1 yr^-1]");

	// for(unsigned int i = 0; i < energy_bins.size(); i++)
	// {
	// 	double E_min = energy_bins[i][1];
	// 	double E_max = energy_bins[i][2];
	// 	double dR	 = libphysica::Integrate(drde, E_min, E_max);
	// 	if(dR > 0.0)
	// 		dRdE_list[i] = dR;
	// 	else
	// 		break;
	// }
	// for(unsigned int i = 0; i < dRdE_list.size(); i++)
	// 	std::cout << i << "\t" << dRdE_list[i] << std::endl;

	// libphysica::Vector dRdE(dRdE_list);

	// libphysica::Vector RS2 = matrix * dRdE;
	// std::cout << RS2.Size() << std::endl;

	// std::vector<unsigned int> s2_bins = {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2050, 2100, 2150, 2200, 2250, 2300, 2350, 2400, 2450, 2500, 2550, 2600, 2650, 2700, 2750, 2800, 2850, 2900, 2950, 3000};
	// unsigned int n_bins				  = s2_bins.size() - 1;
	// std::vector<double> R_S2_binned(n_bins, 0.0);
	// for(unsigned int i = 0; i < n_bins; i++)
	// {
	// 	double s2_min = s2_bins[i];
	// 	double s2_max = s2_bins[i + 1];
	// 	for(unsigned int j = 0; j < s2_bin_info.size(); j++)
	// 		if(s2_bin_info[j][1] >= s2_min && s2_bin_info[j][2] < s2_max)
	// 		{
	// 			R_S2_binned[i] += RS2[j];
	// 			continue;
	// 		}
	// }
	// std::cout << "\nS2 bins of width 50 PE" << std::endl;
	// for(unsigned int i = 0; i < n_bins; i++)
	// 	std::cout << s2_bins[i] << "\t" << s2_bins[i + 1] << "\t" << R_S2_binned[i] << std::endl;

	// std::ofstream f;
	// f.open("R_S2_unbinned.dat");
	// for(unsigned int i = 0; i < RS2.Size(); i++)
	// 	f << s2_bin_info[i][1] << "\t" << RS2[i] * kg * year << std::endl;
	// f.close();
	// f.open("R_S2_binned.dat");
	// for(unsigned int i = 0; i < n_bins; i++)
	// 	f << s2_bins[i] << "\t" << R_S2_binned[i] * kg * year << std::endl;
	// f.close();

	// DM_Detector_Ionization_ER detector("XENON1T_S2", 1.0, "Xe");
	// detector.Initialize_S2_Spectrum("Response matrix");
	// detector.Use_PE_Bins(s2_bins);
	// std::vector<double> spectrum = detector.DM_Signals_PE_Bins(*(cfg.DM), *(cfg.DM_distr));

	// for(int i = 0; i < spectrum.size(); i++)
	// 	std::cout << spectrum[i] << "\t" << R_S2_binned[i] << std::endl;

	// ////////////////////////////////////////////////////////////////////////
	// Atom Xe("Xe");
	// // double q = 23.0 * keV;
	// // double E = 19.0 * eV;
	// // std::cout << Xe.Electron(5, 1).Atomic_Response_Function(1, q, E) << std::endl;
	// // std::cout << Xe.Electron(5, 1).Ionization_Form_Factor(q, E) << std::endl;
	// // std::vector<double> DM_masses = libphysica::Log_Space(cfg.constraints_mass_min, cfg.constraints_mass_max, cfg.constraints_masses);

	// // std::vector<std::vector<double>> exclusion_limits = cfg.DM_detector->Upper_Limit_Curve(*(cfg.DM), *(cfg.DM_distr), DM_masses, cfg.constraints_certainty);
	// // for(unsigned int i = 0; i < exclusion_limits.size(); i++)
	// // 	std::cout << i + 1 << "/" << exclusion_limits.size()
	// // 			  << "\tmDM = " << libphysica::Round(In_Units(exclusion_limits[i][0], (exclusion_limits[i][0] < GeV) ? MeV : GeV)) << ((exclusion_limits[i][0] < GeV) ? " MeV" : " GeV")
	// // 			  << "\tUpper Bound:\t" << libphysica::Round(In_Units(exclusion_limits[i][1], cm * cm)) << std::endl;

	// // int CL = std::round(100.0 * cfg.constraints_certainty);
	// // libphysica::Export_Table(TOP_LEVEL_DIR "results/" + cfg.ID + "/DD_Constraints_" + std::to_string(CL) + ".txt", exclusion_limits, {GeV, cm * cm});

	// DM_Particle_SI DM(35.0 * MeV);
	// DM.Set_Sigma_Electron(pb);
	// auto energies = libphysica::Log_Space(0.1 * eV, 200 * eV, 100);
	// std::vector<std::vector<double>> spectrum;
	// for(auto E : energies)
	// 	spectrum.push_back({E, dRdEe_Ionization_ER(E, DM, *cfg.DM_distr, Xe)});
	// libphysica::Export_Table("dRdEe_mX_35MeV_sigma_1pb.dat", spectrum, {eV, 1.0 / keV / kg / year}, "# E [eV]\tdR/dEe [1/keV/kg/year]");
	// // extern double dRdEe_Ionization_ER(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, Atom& atom);

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