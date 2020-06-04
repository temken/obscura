#ifndef __Configuration__hpp_
#define __Configuration__hpp_

#include <libconfig.h++>
#include <string>

#include "DM_Distribution.hpp"
#include "DM_Particle.hpp"
#include "Direct_Detection.hpp"

namespace obscura
{

class Configuration
{
  protected:
	libconfig::Config config;
	std::string cfg_file;

	virtual void Read_Config_File();

	void Initialize_Result_Folder(int MPI_rank = 0);
	void Create_Result_Folder(int MPI_rank = 0);
	void Copy_Config_File(int MPI_rank = 0);

	virtual void Construct_DM_Particle();
	virtual void Construct_DM_Particle_Standard(std::string DM_interaction);

	virtual void Construct_DM_Distribution();
	virtual void Construct_DM_Distribution_SHM();

	virtual void Construct_DM_Detector();
	virtual void Construct_DM_Detector_Nuclear();
	virtual void Construct_DM_Detector_Ionization();
	virtual void Construct_DM_Detector_Semiconductor();

	virtual void Initialize_Parameters();

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

	virtual void Print_Summary(int MPI_rank = 0);
};

}	// namespace obscura

#endif