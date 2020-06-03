#ifndef __Configuration__hpp_
#define __Configuration__hpp_

#include <string>
#include <libconfig.h++>

#include "Direct_Detection.hpp"
#include "DM_Distribution.hpp"
#include "DM_Particle.hpp"

namespace obscura
{

class Configuration
{
	private:
		libconfig::Config config;
		virtual void Read_Config_File();
		virtual void Construct_DM_Particle();
		virtual void Construct_DM_Distribution();
		virtual void Construct_DM_Detector();

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

		DM_Particle *DM =	{ nullptr };
		DM_Distribution *DM_distr =	{ nullptr };
		DM_Detector *DM_detector =	{ nullptr };

		//Constructors
		Configuration();
		Configuration(std::string cfg_filename, int MPI_rank = 0);
 
		virtual void Print_Summary(int MPI_rank = 0);
};

}	// namespace obscura

#endif