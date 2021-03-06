#ifndef __Configuration__hpp_
#define __Configuration__hpp_

#include <libconfig.h++>
#include <string>

#include "obscura/DM_Distribution.hpp"
#include "obscura/DM_Particle.hpp"
#include "obscura/Direct_Detection.hpp"

namespace obscura
{

class Configuration
{
  protected:
	libconfig::Config config;
	std::string cfg_file;

	void Read_Config_File();

	void Initialize_Result_Folder(int MPI_rank = 0);
	void Create_Result_Folder(int MPI_rank = 0);
	void Copy_Config_File(int MPI_rank = 0);

	virtual void Construct_DM_Particle();
	void Construct_DM_Particle_Standard(std::string DM_interaction);

	virtual void Construct_DM_Distribution();
	void Construct_DM_Halo_Model(std::string model_label);
	void Construct_Imported_Distribution();

	virtual void Construct_DM_Detector();
	void Construct_DM_Detector_Nuclear();
	void Construct_DM_Detector_Ionization(std::string ER_or_Migdal);
	void Construct_DM_Detector_Crystal();

	void Initialize_Parameters();

	void Print_Summary_Base(int MPI_rank);

  public:
	std::string ID;
	std::string results_path;

	//Direct detection constraints
	double constraints_mass_min, constraints_mass_max;
	unsigned int constraints_masses;

	double constraints_certainty;

	DM_Particle* DM			  = {nullptr};
	DM_Distribution* DM_distr = {nullptr};
	DM_Detector* DM_detector  = {nullptr};

	//Constructors
	Configuration();
	explicit Configuration(std::string cfg_filename, int MPI_rank = 0);

	virtual void Print_Summary(int MPI_rank = 0) { Print_Summary_Base(MPI_rank); };
};

}	// namespace obscura

#endif
