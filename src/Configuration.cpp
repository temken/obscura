#include "Configuration.hpp"

#include <fstream>
// #include <sys/types.h> // required for stat.h
#include <sys/stat.h>	//required to create a folder

#include "DM_Particle_Standard.hpp"
#include "Direct_Detection_Nucleus.hpp"
#include "Direct_Detection_Ionization.hpp"
#include "Direct_Detection_Semiconductor.hpp"

using namespace libconfig;
// Read in the configuration file
	Configuration::Configuration()
	{
		Configuration("config.cfg");
	}
	Configuration::Configuration(std::string cfg_filename, int MPI_rank)
	: cfg_file(cfg_filename), results_path("./")
	{
		Read_Config_File();

		//2. Identifier for the program's run.
		try
		{
			ID = config.lookup("ID").c_str();
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "Error in Configuration::Configuration(std::string): No 'ID' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}

		//3. DM particle
		Construct_DM_Particle();
		
		//4. DM Distribution
		Construct_DM_Distribution();

		//5. DM-detection experiment
		Construct_DM_Detector();

		// 6. Computation of exclusion limits
		try
		{
			constraints_certainty = config.lookup("constraints_certainty");
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'constraints_certainty' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		
		try
		{
			constraints_mass_min = config.lookup("constraints_mass_min");
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'constraints_mass_min' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		
		try
		{
			constraints_mass_max = config.lookup("constraints_mass_max");
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'constraints_mass_max' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		
		try
		{
			constraints_masses = config.lookup("constraints_masses");
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'constraints_masses' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}

		Create_Result_Folder(MPI_rank);
		Copy_Config_File(MPI_rank);
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
				std::cerr <<"\nWarning in Configuration::Create_Result_Folder(int): The folder exists already, data will be overwritten."<< std::endl <<std::endl;
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
			std::cout<<"Direct detection constraints"<<std::endl
						<<"\tCertainty level [%]:\t" <<100.0*constraints_certainty <<std::endl
						<<"\tMass range [GeV]:\t[" <<constraints_mass_min<<","<<constraints_mass_max<<"]" <<std::endl
						<<"\tMass steps:\t\t" <<constraints_masses <<std::endl
						<<line<<std::endl<<std::endl;
		}	
	}

	void Configuration::Read_Config_File()
	{
		try
		{		
			config.readFile(cfg_file.c_str());
		}
		catch(const FileIOException &fioex)
		{
			std::cerr << "Error in Configuration::Read_Config_File(): I/O error while reading configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		catch(const ParseException &pex)
		{
			std::cerr << "Error in Configuration::Read_Config_File(): Configurate file parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	void Configuration::Construct_DM_Particle()
	{
		double DM_mass, DM_spin, DM_fraction;
		bool DM_light;
		//3.1 General properties
		try
		{
			DM_mass = config.lookup("DM_mass");
			DM_mass *= GeV;
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'DM_mass' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		
		try
		{
			DM_spin = config.lookup("DM_spin");
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'DM_spin' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		
		try
		{
			DM_fraction = config.lookup("DM_fraction");
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'DM_fraction' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		
		try
		{
			DM_light = config.lookup("DM_light");
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
			DM_interaction = config.lookup("DM_interaction").c_str();
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
					DM_form_factor = config.lookup("DM_form_factor").c_str();
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
						DM_mediator_mass = config.lookup("DM_mediator_mass");
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
				DM_isospin_conserved = config.lookup("DM_isospin_conserved");
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
					fp_rel = config.lookup("DM_relative_couplings")[0];
					fn_rel = config.lookup("DM_relative_couplings")[1];
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
				DM_cross_section_nucleon = config.lookup("DM_cross_section_nucleon");
				DM_cross_section_nucleon *= cm*cm;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DM_cross_section_nucleon' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			DM->Set_Interaction_Parameter(DM_cross_section_nucleon, "Nuclei");

			double DM_cross_section_electron;
			try
			{
				DM_cross_section_electron = config.lookup("DM_cross_section_electron");
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
			std::cerr << "Error in Configuration::Construct_DM_Particle(): 'DM_interaction' setting "<<DM_interaction <<" in configuration file not recognized." << std::endl;
			std::exit(EXIT_FAILURE);
		}

		DM->Set_Mass(DM_mass);
		DM->Set_Spin(DM_spin);
		DM->Set_Fractional_Density(DM_fraction);
		DM->Set_Low_Mass_Mode(DM_light);
	}

	void Configuration::Construct_DM_Distribution()
	{
		std::string DM_distribution;
		double DM_local_density;
		try
		{
			DM_distribution = config.lookup("DM_distribution").c_str();
		}
		catch(const SettingNotFoundException &nfex)
		{
			std::cerr << "No 'DM_distribution' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		
		try
		{
			DM_local_density = config.lookup("DM_local_density");
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
				SHM_v0 = config.lookup("SHM_v0");
				SHM_v0 *= km/sec;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'SHM_v0' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				SHM_vEarth = config.lookup("SHM_vEarth");
				SHM_vEarth *= km/sec;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'SHM_vEarth' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				SHM_vEscape = config.lookup("SHM_vEscape");
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
			std::cerr << "Error in Configuration::Construct_DM_Distribution(): 'DM_distribution' setting "<<DM_distribution <<" in configuration file not recognized." << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	void Configuration::Construct_DM_Detector()
	{
		std::string DD_experiment;
		try
		{
			DD_experiment = config.lookup("DD_experiment").c_str();
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
				int element_count=config.lookup("DD_targets_nuclear").getLength();
				for(int j=0;j<element_count;j++)
				{
					double abund = config.lookup("DD_targets_nuclear")[j][0];
					int Z = config.lookup("DD_targets_nuclear")[j][1];
					DD_targets_nuclear_abundances.push_back(abund);
					DD_targets_nuclear.push_back( Get_Element(Z) );
				}
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_targets_nuclear' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}

			double DD_threshold_nuclear, DD_Emax_nuclear, DD_exposure_nuclear,  DD_efficiency_nuclear, DD_expected_background_nuclear;
			unsigned int DD_observed_events_nuclear;
			try
			{
				DD_threshold_nuclear = config.lookup("DD_threshold_nuclear");
				DD_threshold_nuclear *= keV;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_threshold_nuclear' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				DD_Emax_nuclear = config.lookup("DD_Emax_nuclear");
				DD_Emax_nuclear *= keV;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_Emax_nuclear' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				DD_exposure_nuclear = config.lookup("DD_exposure_nuclear");
				DD_exposure_nuclear *= kg*yr;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_exposure_nuclear' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				DD_efficiency_nuclear = config.lookup("DD_efficiency_nuclear");
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_efficiency_nuclear' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				DD_observed_events_nuclear = config.lookup("DD_observed_events_nuclear");
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_observed_events_nuclear' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				DD_expected_background_nuclear = config.lookup("DD_expected_background_nuclear");
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_expected_background_nuclear' setting in configuration file." << std::endl;
			}
			DM_detector = new DM_Detector_Nucleus(DD_experiment, DD_exposure_nuclear, DD_targets_nuclear, DD_threshold_nuclear, DD_Emax_nuclear,DD_targets_nuclear_abundances);
			DM_detector->Set_Flat_Efficiency(DD_efficiency_nuclear);
			DM_detector->Set_Observed_Events(DD_observed_events_nuclear);
			DM_detector->Set_Expected_Background(DD_expected_background_nuclear);
		}
		else if(DD_experiment == "DAMIC")
		{
			double DAMIC_exposure = 0.107*kg*day;
			std::vector<Element> DAMIC_targets = {Get_Element(14)};
			double DAMIC_threshold = 0.55*keV;
			double DAMIC_Emax = 7.0*keV;
			unsigned int DAMIC_observed_events = 106;
			
			DM_detector = new DM_Detector_Nucleus(DD_experiment, DAMIC_exposure, DAMIC_targets, DAMIC_threshold, DAMIC_Emax);
			DM_detector->Set_Observed_Events(DAMIC_observed_events);
		}
		else if(DD_experiment == "XENON1T")
		{
			double XENON1T_exposure = 34.2*day*1042*kg;
			std::vector<Element> XENON1T_targets = {Get_Element(54)};
			double XENON1T_threshold = 5.0*keV;
			double XENON1T_Emax = 40.0*keV;
			double XENON1T_efficiency = 0.82;
			unsigned int XENON1T_observed_events = 0;
			
			DM_detector = new DM_Detector_Nucleus(DD_experiment, XENON1T_exposure, XENON1T_targets, XENON1T_threshold, XENON1T_Emax);
			DM_detector->Set_Observed_Events(XENON1T_observed_events);
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
		else if(DD_experiment == "Ionization")
		{
			// std::string DD_target_ionization;
			// unsigned int DD_threshold_ionization, DD_observed_events_ionization;
			// double DD_exposure_ionization, DD_efficiency_ionization, DD_expected_background_ionization;
			// try
			// {
			// 	DD_target_ionization = config.lookup("DD_target_ionization").c_str();
			// }
			// catch(const SettingNotFoundException &nfex)
			// {
			// 	std::cerr << "No 'DD_target_ionization' setting in configuration file." << std::endl;
			// 	std::exit(EXIT_FAILURE);
			// }
			// try
			// {
			// 	DD_threshold_ionization = config.lookup("DD_threshold_ionization");
			// }
			// catch(const SettingNotFoundException &nfex)
			// {
			// 	std::cerr << "No 'DD_threshold_ionization' setting in configuration file." << std::endl;
			// 	std::exit(EXIT_FAILURE);
			// }
			// try
			// {
			// 	DD_exposure_ionization = config.lookup("DD_exposure_ionization");
			// 	DD_exposure_ionization *= kg*yr;
			// }
			// catch(const SettingNotFoundException &nfex)
			// {
			// 	std::cerr << "No 'DD_exposure_ionization' setting in configuration file." << std::endl;
			// 	std::exit(EXIT_FAILURE);
			// }
			// try
			// {
			// 	DD_efficiency_ionization = config.lookup("DD_efficiency_ionization");
			// }
			// catch(const SettingNotFoundException &nfex)
			// {
			// 	std::cerr << "No 'DD_efficiency_ionization' setting in configuration file." << std::endl;
			// 	std::exit(EXIT_FAILURE);
			// }
			// try
			// {
			// 	DD_observed_events_ionization = config.lookup("DD_observed_events_ionization");
			// }
			// catch(const SettingNotFoundException &nfex)
			// {
			// 	std::cerr << "No 'DD_observed_events_ionization' setting in configuration file." << std::endl;
			// 	std::exit(EXIT_FAILURE);
			// }
			// try
			// {
			// 	DD_expected_background_ionization = config.lookup("DD_expected_background_ionization");
			// }
			// catch(const SettingNotFoundException &nfex)
			// {
			// 	std::cerr << "No 'DD_expected_background_ionization' setting in configuration file." << std::endl;
			// 	std::exit(EXIT_FAILURE);
			// }
			
			// DM_detector = new DM_Detector_Ionization(DD_experiment, DD_exposure_ionization, DD_target_ionization);
			// DM_detector->Use_Electron_Bins(DD_threshold_ionization);
			// DM_detector->Set_Flat_Efficiency(DD_efficiency_ionization);
			// DM_detector->Set_Observed_Events(DD_observed_events_ionization);
			// DM_detector->Set_Expected_Background(DD_expected_background_ionization);
		}
		else if(DD_experiment == "XENON10e")
		{
			std::string target_name = "Xenon";
			Atom target = Import_Ionization_Form_Factors(target_name);
			double exposure = 15*kg*day;
			double flat_efficiency = 0.92;
			std::vector<unsigned long int> observed_event_bins = {126, 60, 12, 3, 2, 0, 2};
			double muPE = 27.0;
			double sigPE = 6.7;
			std::vector<int> S2_bin_ranges = {14,41,68,95,122,149,176,203};
			std::string trigger_efficiency = "../data/XENON10e/PE_Trigger_Efficiency.txt";

			DM_detector = new DM_Detector_Ionization(DD_experiment, exposure, target);
			DM_detector->Set_Flat_Efficiency(flat_efficiency);
			DM_detector->Set_Observed_Events(observed_event_bins);
			dynamic_cast<DM_Detector_Ionization*>(DM_detector)->Use_PE_Bins(muPE, sigPE, S2_bin_ranges);
			dynamic_cast<DM_Detector_Ionization*>(DM_detector)->Import_Trigger_Efficiency_PE(trigger_efficiency);
		}
		else if(DD_experiment == "XENON100e")
		{
			std::string target_name = "Xenon";
			Atom target = Import_Ionization_Form_Factors(target_name);
			double exposure = 30*kg*yr;
			std::vector<unsigned long int> observed_event_bins = {794, 1218, 924, 776, 669, 630, 528, 488, 433, 387};
			double muPE = 19.7;
			double sigPE = 6.2;
			std::vector<int> S2_bin_ranges = {80, 90, 110,130,150,170,190,210,230,250,270};
			std::string trigger_efficiency = "../data/XENON100e/PE_Trigger_Efficiency.txt";
			std::string acceptance_efficiency = "../data/XENON100e/PE_Acceptance_Efficiency.txt";

			DM_detector = new DM_Detector_Ionization(DD_experiment, exposure, target);
			DM_detector->Set_Observed_Events(observed_event_bins);
			dynamic_cast<DM_Detector_Ionization*>(DM_detector)->Use_PE_Bins(muPE, sigPE, S2_bin_ranges);
			dynamic_cast<DM_Detector_Ionization*>(DM_detector)->Import_Trigger_Efficiency_PE(trigger_efficiency);
			dynamic_cast<DM_Detector_Ionization*>(DM_detector)->Import_Acceptance_Efficiency_PE(acceptance_efficiency);
		}
		else if(DD_experiment == "XENON1Te")
		{
			std::string target_name = "Xenon";
			Atom target = Import_Ionization_Form_Factors(target_name);
			double exposure = 80755.2*kg*day;
			std::vector<unsigned long int> observed_event_bins = {8, 7, 2, 1};
			double muPE = 33.0;
			double sigPE = 7.0;
			std::vector<int> S2_bin_ranges = {150,200,250,300,350};
			std::string trigger_efficiency = "../data/XENON1Te/XENON1T_TotalEfficiency.txt";

			DM_detector = new DM_Detector_Ionization(DD_experiment, exposure, target);
			DM_detector->Set_Observed_Events(observed_event_bins);
			dynamic_cast<DM_Detector_Ionization*>(DM_detector)->Use_PE_Bins(muPE, sigPE, S2_bin_ranges);
			dynamic_cast<DM_Detector_Ionization*>(DM_detector)->Import_Trigger_Efficiency_PE(trigger_efficiency);
		}
		else if(DD_experiment == "DarkSide-50")
		{
			std::string target_name = "Argon";
			Atom target = Import_Ionization_Form_Factors(target_name);
			double exposure = 6786.0*kg*day;
			unsigned int ne_threshold = 3;
			std::vector<unsigned long int> observed_event_bins = {118643, 219893, 6131, 673, 252, 227, 198, 199, 189, 247, 230, 261,249, 329, 336, 349, 351, 352, 384, 411, 405, 461, 460, 436, 500, 546, 538, 536, 556, 583, 573, 630, 603, 635, 639, 682, 736, 755, 804, 811, 809, 882, 934, 935, 871, 965, 946, 1072, 997, 1060};

			DM_detector = new DM_Detector_Ionization(DD_experiment, exposure, target);
			DM_detector->Set_Observed_Events(observed_event_bins);
			dynamic_cast<DM_Detector_Ionization*>(DM_detector)->Use_Electron_Bins(ne_threshold,15);
		}
		else if(DD_experiment == "Semiconductor")
		{
			std::string DD_target_semiconductor;
			unsigned int DD_threshold_semiconductor, DD_observed_events_semiconductor;
			double DD_exposure_semiconductor, DD_efficiency_semiconductor, DD_expected_background_semiconductor;
			try
			{
				DD_target_semiconductor = config.lookup("DD_target_semiconductor").c_str();
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_target_semiconductor' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				DD_threshold_semiconductor = config.lookup("DD_threshold_semiconductor");
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_threshold_semiconductor' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				DD_exposure_semiconductor = config.lookup("DD_exposure_semiconductor");
				DD_exposure_semiconductor *= gram*yr;
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_exposure_semiconductor' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				DD_efficiency_semiconductor = config.lookup("DD_efficiency_semiconductor");
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_efficiency_semiconductor' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				DD_observed_events_semiconductor = config.lookup("DD_observed_events_semiconductor");
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_observed_events_semiconductor' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			try
			{
				DD_expected_background_semiconductor = config.lookup("DD_expected_background_semiconductor");
			}
			catch(const SettingNotFoundException &nfex)
			{
				std::cerr << "No 'DD_expected_background_semiconductor' setting in configuration file." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			
			DM_detector = new DM_Detector_Semiconductor(DD_experiment, DD_exposure_semiconductor, DD_target_semiconductor, DD_threshold_semiconductor);
			DM_detector->Set_Flat_Efficiency(DD_efficiency_semiconductor);
			DM_detector->Set_Observed_Events(DD_observed_events_semiconductor);
			DM_detector->Set_Expected_Background(DD_expected_background_semiconductor);
		}
		else if(DD_experiment == "SENSEI-surface")
		{
				double SENSEI_surface_exposure = 0.07*gram*456*minute;
				unsigned int SENSEI_surface_Q_threshold = 1;
				unsigned int SENSEI_surface_N_bins = 5;
				std::vector<double> SENSEI_surface_efficiencies = {0.668,0.41,0.32,0.27,0.24};
				std::vector<unsigned long int> SENSEI_surface_observed_events = {140302,4676,131,1,0};

				DM_detector = new DM_Detector_Semiconductor(DD_experiment,SENSEI_surface_exposure, "Si" ,SENSEI_surface_Q_threshold);
				dynamic_cast<DM_Detector_Semiconductor*>(DM_detector)->Use_Q_Bins(SENSEI_surface_Q_threshold, SENSEI_surface_N_bins);
				DM_detector->Set_Observed_Events(SENSEI_surface_observed_events);
				DM_detector->Set_Bin_Efficiencies(SENSEI_surface_efficiencies);
		}
		else if(DD_experiment == "SENSEI")
		{
				double SENSEI_exposure = 0.246*gram*day;
				unsigned int SENSEI_Q_threshold = 1;
				unsigned int SENSEI_N_bins = 3;
				std::vector<double> SENSEI_efficiencies = {1.0,0.62,0.48};
				std::vector<unsigned long int> SENSEI_observed_events = {8516,87,0};

				DM_detector = new DM_Detector_Semiconductor(DD_experiment,SENSEI_exposure, "Si" ,SENSEI_Q_threshold);
				dynamic_cast<DM_Detector_Semiconductor*>(DM_detector)->Use_Q_Bins(SENSEI_Q_threshold, SENSEI_N_bins);
				DM_detector->Set_Observed_Events(SENSEI_observed_events);
				DM_detector->Set_Bin_Efficiencies(SENSEI_efficiencies);
		}
		else if(DD_experiment == "SuperCDMS")
		{
			double SuperCDMS_exposure = 0.487*gram*day;
			double SuperCDMS_flat_efficiency = 0.9545;
			unsigned int SuperCDMS_Q_threshold = 1;
			unsigned int SuperCDMS_N_bins = 6;
			std::vector<double> SuperCDMS_efficiencies = {0.88,0.91,0.91,0.91,0.91,0.91};
			std::vector<unsigned long int> SuperCDMS_observed_events = {53000, 400, 74, 18, 7, 14};

			DM_detector = new DM_Detector_Semiconductor(DD_experiment,SuperCDMS_exposure, "Si" ,SuperCDMS_Q_threshold);
			DM_detector->Set_Flat_Efficiency(SuperCDMS_flat_efficiency);
			dynamic_cast<DM_Detector_Semiconductor*>(DM_detector)->Use_Q_Bins(SuperCDMS_Q_threshold, SuperCDMS_N_bins);
			DM_detector->Set_Observed_Events(SuperCDMS_observed_events);
			DM_detector->Set_Bin_Efficiencies(SuperCDMS_efficiencies);
		}
		else
		{
			std::cerr << "Error in Configuration::Construct_DM_Detector(): Experiment " <<DD_experiment<<" not recognized." << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

