#include "obscura/Astronomy.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"

namespace obscura
{

using namespace libphysica::natural_units;

//1. Transform between astronomical coordinate systems.
//Transformation matrices
libphysica::Matrix Transformation_Matrix_P(double T)
{
	double zeta_A  = 2306.083227 * arcsec * T + 0.298850 * arcsec * T * T;
	double z_A	   = 2306.077181 * arcsec * T + 1.092735 * arcsec * T * T;
	double theta_A = 2004.191903 * arcsec * T - 0.429493 * arcsec * T * T;
	libphysica::Matrix P(
		{{cos(zeta_A) * cos(theta_A) * cos(z_A) - sin(zeta_A) * sin(z_A), -sin(zeta_A) * cos(theta_A) * cos(z_A) - cos(zeta_A) * sin(z_A), -sin(theta_A) * cos(z_A)},
		 {cos(zeta_A) * cos(theta_A) * sin(z_A) + sin(zeta_A) * cos(z_A), -sin(zeta_A) * cos(theta_A) * sin(z_A) + cos(zeta_A) * cos(z_A), -sin(theta_A) * sin(z_A)},
		 {cos(zeta_A) * sin(theta_A), -sin(zeta_A) * sin(theta_A), cos(theta_A)}});
	return P;
}

libphysica::Matrix Transformation_Matrix_M()
{
	double l_CP		= 122.932 * deg;
	double alpha_GP = 192.85948 * deg;
	double delta_GP = 27.12825 * deg;
	libphysica::Matrix M(
		{{-sin(l_CP) * sin(alpha_GP) - cos(l_CP) * cos(alpha_GP) * sin(delta_GP), sin(l_CP) * cos(alpha_GP) - cos(l_CP) * sin(alpha_GP) * sin(delta_GP), cos(l_CP) * cos(delta_GP)},
		 {cos(l_CP) * sin(alpha_GP) - sin(l_CP) * cos(alpha_GP) * sin(delta_GP), -cos(l_CP) * cos(alpha_GP) - sin(l_CP) * sin(alpha_GP) * sin(delta_GP), sin(l_CP) * cos(delta_GP)},
		 {cos(alpha_GP) * cos(delta_GP), sin(alpha_GP) * cos(delta_GP), sin(delta_GP)}});
	return M;
}

libphysica::Matrix Transformation_Matrix_R(double T)
{
	double eps			 = 23.4393 * deg - 0.0130 * deg * T;
	libphysica::Matrix R = Rotation_Matrix(eps, 3, libphysica::Vector({1, 0, 0}));
	return R;
}

//Coordinate transformations
libphysica::Vector Transform_Equat_to_Gal(const libphysica::Vector& v_Equat, double T)
{
	return Transformation_Matrix_M() * Transformation_Matrix_P(T).Inverse() * v_Equat;
}

libphysica::Vector Transform_GeoEcl_to_Gal(const libphysica::Vector& v_GeoEcl, double T)
{
	return Transformation_Matrix_M() * Transformation_Matrix_P(T).Inverse() * Transformation_Matrix_R(T) * v_GeoEcl;
}

libphysica::Vector Transform_HelEcl_to_Gal(const libphysica::Vector& v_HelEcl, double T)
{
	return (-1.0) * Transformation_Matrix_M() * Transformation_Matrix_P(T).Inverse() * Transformation_Matrix_R(T) * v_HelEcl;
}

libphysica::Vector Transform_Gal_to_Equat(const libphysica::Vector& v_Gal, double T)
{
	return Transformation_Matrix_P(T) * Transformation_Matrix_M().Inverse() * v_Gal;
}

//2. Times
//Fractional days since epoch J2000.0.
double Fractional_Days_since_J2000(int day, int month, int year, double hour, double minute, double second)
{
	double n = 0.0;
	if(month == 1 || month == 2)
	{
		n += floor(365.25 * (year - 1)) + floor(30.61 * (month + 13));
	}
	else
	{
		n += floor(365.25 * year) + floor(30.61 * (month + 1));
	}
	n += day - 730563.5 + hour / 24. + minute / 1440. + second / 86400.;
	return n;
}

double Local_Apparent_Sidereal_Time(double nJ2000, double longitude)
{
	double T = nJ2000 / 36525.;
	//GMST:
	double t		= 86400. * (0.779057273264 + 0.00273781191135448 * nJ2000 + fmod(nJ2000, 1)) + 0.00096707 + 307.47710227 * T + 0.092772113 * pow(T, 2);
	double epsilonA = 23.439279444 * deg - 0.01301021361 * deg * T;
	double Omega	= 125.04455501 * deg - 0.05295376 * deg * nJ2000;
	double L		= 280.47 * deg - 0.98565 * deg * nJ2000;
	double DeltaPsi = -1.1484 * sin(Omega) - 0.0864 * cos(2 * L);
	//GAST:
	t += DeltaPsi * cos(epsilonA) + 0.000176 * sin(Omega) + 0.000004 * sin(2 * Omega);
	//LAST:
	t += longitude / 2.0 / M_PI * 86400.;
	return fmod(t, 86400) * sec;
}

//3. Sun's and earth's velocity in the galactic frame
libphysica::Vector Sun_Velocity()
{
	//1. Sun's rotation around galactic center:
	libphysica::Vector vGal({0, 220 * km / sec, 0});
	//2. Sun's peculiar motion:
	libphysica::Vector vPec({11.1 * km / sec, 12.2 * km / sec, 7.3 * km / sec});

	return vGal + vPec;
}

libphysica::Vector Earth_Velocity(double nJ2000)
{
	//1. Sun's velocity
	libphysica::Vector vSun = Sun_Velocity();

	//2. Earth's rotation around galactic center:
	double e	 = 0.01671;	  //Ellipticity of earth's orbit
	double L	 = fmod(280.46 * deg + nJ2000 * 0.9856474 * deg, 2 * M_PI);
	double omega = fmod(282.932 * deg + nJ2000 * 0.0000471 * deg, 2 * M_PI);
	double T	 = nJ2000 / 36525.0;
	libphysica::Vector ex({1, 0, 0});
	libphysica::Vector ey({0, 1, 0});

	//Basis vectors in heliocentric ecliptic coordinate system
	libphysica::Vector exEcl = Transform_HelEcl_to_Gal(ex, T);	 //(0.054876-0.024232 * T,-0.494109-0.002689 * T,0.867666 + 1.546e-6 * T);
	libphysica::Vector eyEcl = Transform_HelEcl_to_Gal(ey, T);	 //(0.993824 + 0.001316 * T,0.110992-0.011851 * T,0.000352 + 0.021267 * T);
	double ve				 = 29.79 * km / sec;
	libphysica::Vector uE	 = -ve * (sin(L) + e * sin(2 * L - omega)) * exEcl + ve * (cos(L) + e * cos(2 * L - omega)) * eyEcl;

	return vSun + uE;
}

}	// namespace obscura
