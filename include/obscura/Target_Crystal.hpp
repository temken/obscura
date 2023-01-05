#ifndef __Target_Crystal_hpp_
#define __Target_Crystal_hpp_

#include "libphysica/Numerics.hpp"

namespace obscura
{

class Crystal
{
  private:
	libphysica::Interpolation_2D form_factor_interpolation;

  public:
	int N_E, N_q;
	std::string name;
	double dE, dq, E_max, q_max;
	double M_cell;
	double energy_gap, epsilon;
	unsigned int Q_max;

	explicit Crystal(std::string target);

	double Crystal_Form_Factor(double q, double E);
};
}	// namespace obscura

#endif