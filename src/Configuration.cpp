#include "Configuration.hpp"

#include <fstream>
#include <libconfig.h++>
// #include <sys/types.h> // required for stat.h
#include <sys/stat.h>	//required to create a folder

#include "DM_Particle_Standard.hpp"
#include "Direct_Detection_Nucleus.hpp"

using namespace libconfig;
// Read in the configuration file
	Configuration::Configuration()
	{
		Configuration("config.cfg");
	}
	Configuration::Configuration(std::string cfg_filename, int MPI_rank)
	: cfg_file(cfg_filename), results_path("./")
	{
		Config cfg;

		//1. Read the file. If there is an error, report it and std::exit.
		try
		{		
			cfg.readFile(cfg_filename.c_str());
		}
		catch(const FileIOException &fioex)
		{
			std::cerr << "Error in Configuration::Configuration(std::string): I/O error while reading configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		catch(const ParseException &pex)
		{
			std::cerr << "Error in Configuration::Configuration(std::string): Configurate file parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << std::endl;
			std::exit(EXIT_FAILURE);
		}

		//2. Identifier for the program's run.
		try
		{
			ID = cfg.lookup("ID").c_str();
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "Error in Configuration::Configuration(std::string): No 'ID' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}

		//3. DM particle
		double DM_mass, DM_spin, DM_fraction;
		bool DM_light;
		//3.1 General properties
		try
		{
			DM_mass = cfg.lookup("DM_mass");
			DM_mass *= GeV;
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'DM_mass' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		try
		{
			DM_spin = cfg.lookup("DM_spin");
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'DM_spin' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		try
		{
			DM_fraction = cfg.lookup("DM_fraction");
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'DM_fraction' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		try
		{
			DM_light = cfg.lookup("DM_light");
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'DM_light' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}

		//3.2 DM interactions

		std::string DM_interaction;
		try
		{
			DM_interaction = cfg.lookup("DM_interaction").c_str();
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'DM_interaction' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}

		//3.2.1 SI and SD
		if(DM_interaction == "SI" || DM_interaction == "SD")
		{
			//SI
			if(DM_interaction == "SI") 
			{
				DM = new DM_Particle_SI();

				//DM form factor
				std::string DM_form_factor;
				double DM_mediator_mass = -1.0;
				try
				{
					DM_form_factor = cfg.lookup("DM_form_factor").c_str();
				}
				catch(const SettingNotFoundException &nfex)
				{
					std::cerr << "No 'DM_form_factor' setting in configuration file." << std::endl;
					std::exit(EXIT_FAILURE);
				}
				if(DM_form_factor == "General")
				{
					try
					{
						DM_mediator_mass = cfg.lookup("DM_mediator_mass");
						DM_mediator_mass *= MeV;
					}
					catch(const SettingNotFoundException &nfex)
					{
						std::cerr << "No 'DM_mediator_mass' setting in configuration file." << std::endl;
						std::exit(EXIT_FAILURE);
					}
				}
				dynamic_cast<DM_Particle_SI*>(DM)->Set_FormFactor_DM(DM_form_factor, DM_mediator_mass);
			}

			//SD
			else if (DM_interaction == "SD") 
			{
				DM = new DM_Particle_SD();
			}

			//SI and SD
			bool DM_isospin_conserved;
			try
			{
				DM_isospin_conserved = cfg.lookup("DM_isospin_conserved");
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DM_isospin_conserved' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			double fp_rel, fn_rel;
			if(DM_isospin_conserved)
			{
				fp_rel = 1.0;
				fn_rel = 1.0;
			}
			else
			{
				try
				{
					fp_rel = cfg.lookup("DM_relative_couplings")[0];
					fn_rel = cfg.lookup("DM_relative_couplings")[1];
				}
				catch(const SettingNotFoundException &nfex)
				{
					std::cerr << "No 'DM_relative_couplings' setting in configuration file." << std::endl;
					std::exit(EXIT_FAILURE);
				}
			}
			dynamic_cast<DM_Particle_Standard*>(DM)->Fix_Coupling_Ratio(fp_rel,fn_rel);

			double DM_cross_section_nucleon;
			try
			{
				DM_cross_section_nucleon = cfg.lookup("DM_cross_section_nucleon");
				DM_cross_section_nucleon *= cm*cm;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DM_cross_section_nucleon' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			DM->Set_Interaction_Parameter(DM_cross_section_nucleon);

			double DM_cross_section_electron;
			try
			{
				DM_cross_section_electron = cfg.lookup("DM_cross_section_electron");
				DM_cross_section_electron *= cm*cm;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DM_cross_section_electron' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			DM->Set_Sigma_Electron(DM_cross_section_electron);
		}
		else
		{
			std::cerr << "Error: 'DM_interaction' setting "<<DM_interaction <<" in configuration file not recognized." << std::endl;
			std::exit(EXIT_FAILURE);
		}

		DM->Set_Mass(DM_mass);
		DM->Set_Spin(DM_spin);
		DM->Set_Fractional_Density(DM_fraction);
		DM->Set_Low_Mass_Mode(DM_light);
		
		//4. DM Distribution
		std::string DM_distribution;
		double DM_local_density;
		try
		{
			DM_distribution = cfg.lookup("DM_distribution").c_str();
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'DM_distribution' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		try
		{
			DM_local_density = cfg.lookup("DM_local_density");
			DM_local_density *= GeV /cm/cm/cm;
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'DM_local_density' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}

		//4.1 Standard halo model
		if(DM_distribution == "SHM")
		{
			double SHM_v0, SHM_vEarth, SHM_vEscape;
			try
			{
				SHM_v0 = cfg.lookup("SHM_v0");
				SHM_v0 *= km/sec;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'SHM_v0' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				SHM_vEarth = cfg.lookup("SHM_vEarth");
				SHM_vEarth *= km/sec;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'SHM_vEarth' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				SHM_vEscape = cfg.lookup("SHM_vEscape");
				SHM_vEscape *= km/sec;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'SHM_vEscape' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			DM_distr = new Standard_Halo_Model(DM_local_density,SHM_v0, SHM_vEarth, SHM_vEscape);
		}
		else
		{
			std::cerr << "Error: 'DM_distribution' setting "<<DM_distribution <<" in configuration file not recognized." << std::endl;
			std::exit(EXIT_FAILURE);
		}

		//5. DM-detection experiment
		std::string DD_experiment;
		try
		{
			DD_experiment = cfg.lookup("DD_experiment").c_str();
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'DD_experiment' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		if(DD_experiment == "Nuclear recoil")
		{
			std::vector<Element> DD_targets_nuclear;
			std::vector<double>  DD_targets_nuclear_abundances;
			try
			{
				int element_count=cfg.lookup("DD_targets_nuclear").getLength();
				for(int j=0;j<element_count;j++)
				{
					double abund = cfg.lookup("DD_targets_nuclear")[j][0];
					int Z = cfg.lookup("DD_targets_nuclear")[j][1];
					DD_targets_nuclear_abundances.push_back(abund);
					DD_targets_nuclear.push_back( Get_Element(Z) );
				}
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_targets_nuclear' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}

			double DD_threshold_nuclear, DD_Emax_nuclear, DD_exposure_nuclear,  DD_efficiency_nuclear;
			unsigned int DD_background_nuclear;
			try
			{
				DD_threshold_nuclear = cfg.lookup("DD_threshold_nuclear");
				DD_threshold_nuclear *= keV;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_threshold_nuclear' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				DD_Emax_nuclear = cfg.lookup("DD_Emax_nuclear");
				DD_Emax_nuclear *= keV;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_Emax_nuclear' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				DD_exposure_nuclear = cfg.lookup("DD_exposure_nuclear");
				DD_exposure_nuclear *= kg*yr;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_exposure_nuclear' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				DD_efficiency_nuclear = cfg.lookup("DD_efficiency_nuclear");
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_efficiency_nuclear' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				DD_background_nuclear = cfg.lookup("DD_background_nuclear");
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_background_nuclear' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			// DM_detector = new DM_Detector_Nucleus(DD_targets_nuclear, DD_exposure_nuclear, DD_threshold_nuclear, DD_Emax_nuclear,DD_targets_nuclear_abundances);
			DM_detector = new DM_Detector_Nucleus("Nuclear Recoil", DD_exposure_nuclear, DD_targets_nuclear, DD_threshold_nuclear, DD_Emax_nuclear,DD_targets_nuclear_abundances);
			DM_detector->Set_Flat_Efficiency(DD_efficiency_nuclear);
			DM_detector->Set_Background(DD_background_nuclear);
		}
		else if(DD_experiment == "DAMIC")
		{
			double DAMIC_exposure = 0.107*kg*day;
			std::vector<Element> DAMIC_targets = {Get_Element(14)};
			double DAMIC_threshold = 0.55*keV;
			double DAMIC_Emax = 7.0*keV;
			unsigned int DAMIC_background = 106;
			
			DM_detector = new DM_Detector_Nucleus(DD_experiment, DAMIC_exposure, DAMIC_targets, DAMIC_threshold, DAMIC_Emax);
			DM_detector->Set_Background(DAMIC_background);
		}
		else if(DD_experiment == "XENON1T")
		{
			double XENON1T_exposure = 34.2*day*1042*kg;
			std::vector<Element> XENON1T_targets = {Get_Element(54)};
			double XENON1T_threshold = 5.0*keV;
			double XENON1T_Emax = 40.0*keV;
			double XENON1T_efficiency = 0.82;
			unsigned int XENON1T_background = 0;
			
			DM_detector = new DM_Detector_Nucleus(DD_experiment, XENON1T_exposure, XENON1T_targets, XENON1T_threshold, XENON1T_Emax);
			DM_detector->Set_Background(XENON1T_background);
			DM_detector->Set_Flat_Efficiency(XENON1T_efficiency);
		}
		else if(DD_experiment == "CRESST-II")
		{
			double CRESST_II_exposure = 52.15*kg*day;
			std::vector<Element> CRESST_II_targets = {Get_Element(8),Get_Element(20),Get_Element(74)}; //CaOW
			std::vector<double> CRESST_II_target_ratios = {4,1,1};
			double CRESST_II_threshold = 307*eV;
			double CRESST_II_Emax = 40.0*keV;
			double CRESST_II_resolution = CRESST_II_threshold / 5.0;
			std::vector<std::string> efficiency_files = {"../data/CRESST-II/Lise_eff_AR_O.dat","../data/CRESST-II/Lise_eff_AR_Ca.dat","../data/CRESST-II/Lise_eff_AR_W.dat"};

			DM_detector = new DM_Detector_Nucleus(DD_experiment, CRESST_II_exposure, CRESST_II_targets, CRESST_II_threshold, CRESST_II_Emax, CRESST_II_target_ratios);
			
			DM_detector->Use_Maximum_Gap("../data/CRESST-II/Lise_AR.dat");
			dynamic_cast<DM_Detector_Nucleus*>(DM_detector)->Set_Resolution(CRESST_II_resolution);
			dynamic_cast<DM_Detector_Nucleus*>(DM_detector)->Import_Efficiency(efficiency_files);
		}
		else if(DD_experiment == "CRESST-surface")
		{
			double CRESST_surface_exposure = 0.046*gram*day;
			std::vector<Element> CRESST_surface_targets = {Get_Element(8),Get_Element(13)};
			std::vector<double> CRESST_surface_target_ratios = {3,2};
			double CRESST_surface_threshold = 19.7*eV;
			double CRESST_surface_Emax = 600*eV;
			double CRESST_surface_resolution = 3.74*eV;

			DM_detector = new DM_Detector_Nucleus(DD_experiment, CRESST_surface_exposure, CRESST_surface_targets, CRESST_surface_threshold, CRESST_surface_Emax, CRESST_surface_target_ratios);
			
			DM_detector->Use_Maximum_Gap("../data/CRESST-surface/data.txt");
			dynamic_cast<DM_Detector_Nucleus*>(DM_detector)->Set_Resolution(CRESST_surface_resolution);
		}
		else if(DD_experiment == "CRESST-III")
		{
			double CRESST_III_exposure = 5.594*kg*day;
			std::vector<Element> CRESST_III_targets = {Get_Element(8),Get_Element(20),Get_Element(74)}; //CaOW
			std::vector<double> CRESST_III_target_ratios = {4,1,1};
			double CRESST_III_threshold = 30.1*eV;
			double CRESST_III_Emax = 16*keV;
			double CRESST_III_resolution = 4.6*eV;
			double CRESST_III_efficiency = 0.5;
			std::vector<std::string> efficiency_files = {"../data/CRESST-III/C3P1_DetA_eff_AR_O.dat","../data/CRESST-III/C3P1_DetA_eff_AR_Ca.dat","../data/CRESST-III/C3P1_DetA_eff_AR_W.dat"};

			DM_detector = new DM_Detector_Nucleus(DD_experiment, CRESST_III_exposure, CRESST_III_targets, CRESST_III_threshold, CRESST_III_Emax, CRESST_III_target_ratios);
			
			DM_detector->Use_Maximum_Gap("../data/CRESST-III/C3P1_DetA_AR.dat");
			DM_detector->Set_Flat_Efficiency(CRESST_III_efficiency);
			dynamic_cast<DM_Detector_Nucleus*>(DM_detector)->Set_Resolution(CRESST_III_resolution);
			dynamic_cast<DM_Detector_Nucleus*>(DM_detector)->Import_Efficiency(efficiency_files);
		}
		else if(DD_experiment == "XENON10e")
		{
		// 		detector = new Detector_Ionization("Xe",15*kg*day);
		// 		detector->Set_Flat_Efficiency(0.92);
		// 		std::vector<unsigned long int> data={126, 60, 12, 3, 2, 0, 2};
		// 		dynamic_cast<Detector_Ionization*> (detector)->Set_Binned_Events(data);
		// 		double muPE=27.0;
		// 		double sigPE=6.7;
		// 		std::vector<int> binsizes ={14,41,68,95,122,149,176,203};
		// 		dynamic_cast<Detector_Ionization*> (detector)->Set_PE_Distribution(muPE,sigPE,binsizes);
		// 		dynamic_cast<Detector_Ionization*> (detector)->Import_Trigger_Efficiency_PE("../data/XENON10e/PE_Trigger_Efficiency.txt");
		}
		else if(DD_experiment == "XENON100e")
		{
		// 		detector = new Detector_Ionization("Xe",30*kg*year);
		// 		std::vector<unsigned long int> data100={794, 1218, 924, 776, 669, 630, 528, 488, 433, 387};
		// 		dynamic_cast<Detector_Ionization*> (detector)->Set_Binned_Events(data100);
		// 		double muPE100=19.7;
		// 		double sigPE100=6.2;
		// 		std::vector<int> binsizes100 ={80, 90, 110,130,150,170,190,210,230,250,270};
		// 		dynamic_cast<Detector_Ionization*> (detector)->Set_PE_Distribution(muPE100,sigPE100,binsizes100);
		// 		dynamic_cast<Detector_Ionization*> (detector)->Import_Trigger_Efficiency_PE("../data/XENON100e/PE_Trigger_Efficiency.txt");
		// 		dynamic_cast<Detector_Ionization*> (detector)->Import_Acceptance_Efficiency_PE("../data/XENON100e/PE_Acceptance_Efficiency.txt");
		}
		else if(DD_experiment == "XENON1Te")
		{
		// 		// detector = new Detector_Ionization("Xe",22000*kg*day);
		// 		detector = new Detector_Ionization("Xe",80755.2*kg*day);
		// 	  	detector->Set_Name("XENON1T");
		// 		std::vector<unsigned long int> data_X1T={8, 7, 2, 1};
		// 		dynamic_cast<Detector_Ionization*> (detector)->Set_Binned_Events(data_X1T);
		// 		double muPE_X1T=33.0;
		// 		double sigPE_X1T=7.0;
		// 		std::vector<int> binsizes_X1T ={150,200,250,300,350};
		// 		dynamic_cast<Detector_Ionization*> (detector)->Set_PE_Distribution(muPE_X1T,sigPE_X1T,binsizes_X1T);
		// 		dynamic_cast<Detector_Ionization*> (detector)->Import_Trigger_Efficiency_PE("../data/XENON1Te/XENON1T_TotalEfficiency.txt");
		}
		else if(DD_experiment == "DarkSide-50")
		{
		// 		detector = new Detector_Ionization("Ar",6786.0*kg*day);
		// 		dynamic_cast<Detector_Ionization*> (detector)->Set_Electron_Threshold(3);
		// 		std::vector<unsigned long int>data={118643, 219893, 6131, 673, 252, 227, 198, 199, 189, 247, 230, 261,249, 329, 336, 349, 351, 352, 384, 411, 405, 461, 460, 436, 500, 546, 538, 536, 556, 583, 573, 630, 603, 635, 639, 682, 736, 755, 804, 811, 809, 882, 934, 935, 871, 965, 946, 1072, 997, 1060};
		// 		dynamic_cast<Detector_Ionization*> (detector)->Set_Binned_Events(data);
		}
		else if(DD_experiment == "Semiconductor")
		{
			// std::string DD_target_semiconductor;
			// unsigned int DD_threshold_semiconductor, DD_background_semiconductor;
			// double DD_exposure_semiconductor, DD_efficiency_semiconductor;
			// try
			// {
			// 	DD_target_semiconductor = cfg.lookup("DD_target_semiconductor").c_str();
			// }
			// catch(const SettingNotFoundException &nfex)
			// {
			// 	std::cerr << "No 'DD_target_semiconductor' setting in configuration file." << std::endl;
			// 	std::exit(EXIT_FAILURE);
			// }
			// try
			// {
			// 	DD_threshold_semiconductor = cfg.lookup("DD_threshold_semiconductor");
			// }
			// catch(const SettingNotFoundException &nfex)
			// {
			// 	std::cerr << "No 'DD_threshold_semiconductor' setting in configuration file." << std::endl;
			// 	std::exit(EXIT_FAILURE);
			// }
			// try
			// {
			// 	DD_exposure_semiconductor = cfg.lookup("DD_exposure_semiconductor");
			// 	DD_exposure_semiconductor *= gram*yr;
			// }
			// catch(const SettingNotFoundException &nfex)
			// {
			// 	std::cerr << "No 'DD_exposure_semiconductor' setting in configuration file." << std::endl;
			// 	std::exit(EXIT_FAILURE);
			// }
			// try
			// {
			// 	DD_efficiency_semiconductor = cfg.lookup("DD_efficiency_semiconductor");
			// }
			// catch(const SettingNotFoundException &nfex)
			// {
			// 	std::cerr << "No 'DD_efficiency_semiconductor' setting in configuration file." << std::endl;
			// 	std::exit(EXIT_FAILURE);
			// }
			// try
			// {
			// 	DD_background_semiconductor = cfg.lookup("DD_background_semiconductor");
			// }
			// catch(const SettingNotFoundException &nfex)
			// {
			// 	std::cerr << "No 'DD_background_semiconductor' setting in configuration file." << std::endl;
			// 	std::exit(EXIT_FAILURE);
			// }
			
			// DM_detector = new Detector_Semiconductor(DD_experiment, DD_exposure_semiconductor, DD_target_semiconductor, DD_threshold_semiconductor);
			// DM_detector->Set_Flat_Efficiency(DD_efficiency_semiconductor);
			// DM_detector->Set_Background(DD_background_semiconductor);
		}
		else if(DD_experiment == "SENSEI-surface")
		{
		// 		detector = new Detector_Semiconductor("Si",0.07*gram*456*minute,1);
		// 		std::vector<double> eff={0.668,0.41,0.32,0.27,0.24};
		// 		std::vector<unsigned long int> data ={140302,4676,131,1,0};
		// 		dynamic_cast<Detector_Semiconductor*>(detector)->Set_Binned_Events(data,eff);
		}
		else if(DD_experiment == "SuperCDMS")
		{
		// 		detector = new Detector_Semiconductor("Si",0.487*gram*day,1);
		// 		std::vector<double>  eff={0.88,0.91,0.91,0.91,0.91,0.91};
		// 		std::vector<unsigned long int> data ={53000, 400, 74, 18, 7, 14};
		// 		dynamic_cast<Detector_Semiconductor*>(detector)->Set_Binned_Events(data,eff);
		// 		detector->Set_Flat_Efficiency(0.9545);
		}
		else
		{
			std::cerr << "Error: Experiment " <<DD_experiment<<" not recognized." << std::endl;
			std::exit(EXIT_FAILURE);
		}

		// 6. Computation of exclusion limits
		try
		{
			constraints_certainty = cfg.lookup("constraints_certainty");
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'constraints_certainty' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		try
		{
			constraints_mass_min = cfg.lookup("constraints_mass_min");
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'constraints_mass_min' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		try
		{
			constraints_mass_max = cfg.lookup("constraints_mass_max");
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'constraints_mass_max' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		try
		{
			constraints_masses = cfg.lookup("constraints_masses");
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'constraints_masses' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}

		Create_Result_Folder(MPI_rank);
		Copy_Config_File(MPI_rank);
		Print_Summary(MPI_rank);
	}

	void Configuration::Create_Result_Folder(int MPI_rank)
	{
		if(MPI_rank == 0) 
		{
			//1. Create the /results/ folder if necessary
			std::string results_folder = "../results";
			mode_t nMode = 0733; // UNIX style permissions
			int nError_1 = 0;
			#if defined(_WIN32)
				nError_1 = _mkdir(results_folder.c_str()); // can be used on Windows
			#else
				nError_1 = mkdir(results_folder.c_str(),nMode); // can be used on non-Windows
			#endif

			//2. Create a /result/<ID>/ folder for result files.
			results_path = results_folder + "/" + ID + "/";
			int nError = 0;
			#if defined(_WIN32)
				nError = _mkdir(results_path.c_str()); // can be used on Windows
			#else
			  nError = mkdir(results_path.c_str(),nMode); // can be used on non-Windows
			#endif
			if (nError != 0) 
			{
				std::cerr <<"\nWarning in Configuration::Create_Result_Folder(int): The folder exists already, data will be overwritten."<< std::endl;
			}
		}
	}
	void Configuration::Copy_Config_File(int MPI_rank)
	{
		if(MPI_rank == 0)
		{
			std::ifstream inFile;
			std::ofstream outFile;
			inFile.open(cfg_file);
	    	outFile.open("../results/"+ID+"/"+ID+".cfg");
			outFile << inFile.rdbuf();
			inFile.close();
			outFile.close();
		}
	}

	void Configuration::Print_Summary(int MPI_rank)
	{
		if(MPI_rank == 0)
		{
			std::string line="----------------------------------------";
			std::cout <<std::endl<<line<<std::endl <<"Config file:\t\t"<<cfg_file<< std::endl<< std::endl;
			std::cout <<"\tID:\t\t" <<ID << std::endl;
			DM->Print_Summary(MPI_rank);
			DM_distr->Print_Summary(MPI_rank);
			DM_detector->Print_Summary(MPI_rank);
			std::cout <<line <<std::endl;
		}
	}