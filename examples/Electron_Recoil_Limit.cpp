#include <iostream>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"
#include "obscura/Direct_Detection_Crystal.hpp"
#include "obscura/Direct_Detection_ER.hpp"
#include "obscura/Target_Atom.hpp"
#include "obscura/Target_Crystal.hpp"

using namespace libphysica::natural_units;

int main()
{

	// 1. DM particle
	obscura::DM_Particle_SI dm(100.0 * MeV);
	dm.Print_Summary();

	// 2. DM distribution
	obscura::Standard_Halo_Model shm;
	shm.Print_Summary();

	// 3. Argon target experiment
	obscura::DM_Detector_Ionization_ER argon_experiment("Argon toy experiment", 100.0 * kg * year, "Ar");
	argon_experiment.Use_Electron_Threshold(4);
	argon_experiment.Print_Summary();

	// 4. Si target experiment
	obscura::DM_Detector_Crystal silicon_experiment("Silicon toy experiment", 10.0 * gram * year, "Si");
	silicon_experiment.Use_Q_Threshold(2);
	silicon_experiment.Print_Summary();

	// 5. Compute the 95% CL exclusion limits for m = 100.0 MeV
	double limit_Ar = argon_experiment.Upper_Limit(dm, shm, 0.95);
	double limit_Si = silicon_experiment.Upper_Limit(dm, shm, 0.95);

	std::cout << "Argon experiment: \tsigma_e < " << In_Units(limit_Ar, cm * cm) << " cm^2 (95%CL)" << std::endl;
	std::cout << "Silicon experiment: \tsigma_e < " << In_Units(limit_Si, cm * cm) << " cm^2 (95%CL)" << std::endl;

	return 0;
}