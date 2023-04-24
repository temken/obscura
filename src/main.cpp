#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/Configuration.hpp"
#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"
#include "obscura/Direct_Detection_Crystal.hpp"
#include "obscura/Target_Crystal.hpp"
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
	////////////////////////////////////////////////////////////////////////

	// // Import configuration file
	// obscura::Configuration cfg(argv[1]);
	// cfg.Print_Summary();

	// std::vector<double> DM_masses = libphysica::Log_Space(cfg.constraints_mass_min, cfg.constraints_mass_max, cfg.constraints_masses);

	// std::vector<std::vector<double>> exclusion_limits = cfg.DM_detector->Upper_Limit_Curve(*(cfg.DM), *(cfg.DM_distr), DM_masses, cfg.constraints_certainty);
	// for(unsigned int i = 0; i < exclusion_limits.size(); i++)
	// 	std::cout << i + 1 << "/" << exclusion_limits.size()
	// 			  << "\tmDM = " << libphysica::Round(In_Units(exclusion_limits[i][0], (exclusion_limits[i][0] < GeV) ? MeV : GeV)) << ((exclusion_limits[i][0] < GeV) ? " MeV" : " GeV")
	// 			  << "\tUpper Bound:\t" << libphysica::Round(In_Units(exclusion_limits[i][1], cm * cm)) << std::endl;

	// int CL = std::round(100.0 * cfg.constraints_certainty);
	// libphysica::Export_Table(TOP_LEVEL_DIR "results/" + cfg.ID + "/DD_Constraints_" + std::to_string(CL) + ".txt", exclusion_limits, {GeV, cm * cm});

	DM_Particle_SI DM(100 * MeV);
	DM.Set_Sigma_Electron(pb);
	Standard_Halo_Model SHM;
	Crystal germanium("Ge");
	Crystal silicon("Si");

	int Q	  = 6;
	double Ee = 18 * eV;
	std::cout << silicon.Ionization_Yield(Ee, Q) << std::endl;

	double sum = 0.0;
	for(int Q = 1; Q <= 20; Q++)
		sum += silicon.Ionization_Yield(Ee, Q);
	std::cout << "Sum " << sum << std::endl;

	std::ofstream f;

	f.open("RQ_Ge_new_2.txt");
	std::cout << "Germanium" << std::endl
			  << "Q\tR(Q) [per gram yr]" << std::endl;
	for(int Q = 1; Q < 11; Q++)
	{
		double RQ = R_Q_Crystal(Q, DM, SHM, germanium);
		std::cout << Q << "\t" << In_Units(RQ, 1.0 / gram / year) << std::endl;
		f << Q << "\t" << In_Units(RQ, 1.0 / gram / year) << std::endl;
	}
	f.close();
	f.open("RQ_Si_new_2.txt");
	std::cout << "Silicon" << std::endl
			  << "Q\tR(Q) [per gram yr]" << std::endl;
	for(int Q = 1; Q < 11; Q++)
	{
		double RQ = R_Q_Crystal(Q, DM, SHM, silicon);
		std::cout << Q << "\t" << In_Units(RQ, 1.0 / gram / year) << std::endl;
		f << Q << "\t" << In_Units(RQ, 1.0 / gram / year) << std::endl;
	}
	f.close();

	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	auto time_end			 = std::chrono::system_clock::now();
	double durationTotal	 = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::string final_output = "[Finished " + libphysica::Time_Display(durationTotal) + "]";
	std::cout << libphysica::Formatted_String(final_output, "Green", true) << std::endl;
	return 0;
}