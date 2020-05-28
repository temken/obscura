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
  private:
	libconfig::Config config;
	void Read_Config_File();
	void Construct_DM_Particle();
	void Construct_DM_Distribution();
	void Construct_DM_Detector();

	void Create_Result_Folder(int MPI_rank = 0);
	void Copy_Config_File(int MPI_rank = 0);
	std::string cfg_file;

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
	Configuration(std::string cfg_filename, int MPI_rank = 0);

	void Print_Summary(int MPI_rank = 0);
};

}	// namespace obscura

#endif