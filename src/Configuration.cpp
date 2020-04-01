#include "Configuration.hpp"

#include <libconfig.h++>

// Read in the configuration file
	Configuration::Configuration(std::string cfg_filename, int MPI_rank)
	{
		libconfig::Config cfg;
	}