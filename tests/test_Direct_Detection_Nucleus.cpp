#include "gtest/gtest.h"

#include "obscura/Direct_Detection_Nucleus.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Particle_Standard.hpp"
#include "obscura/Target_Nucleus.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

TEST(TestDirectDetectionNucleus, TestdRdERNucleusIsotope)
{
	// ARRANGE
	double mDM	   = 10.0 * GeV;
	double sigma_p = 1.0 * pb;
	DM_Particle_SI DM(mDM);
	DM.Set_Sigma_Proton(sigma_p);
	DM.Set_Low_Mass_Mode(true);
	Standard_Halo_Model SHM;
	Isotope hydrogen(1, 1);
	double ER	  = 2.0 * keV;
	double result = SHM.DM_density / mDM * sigma_p / 2.0 / libphysica::Reduced_Mass(mDM, mProton) / libphysica::Reduced_Mass(mDM, mProton) * SHM.Eta_Function(vMinimal_Nucleus(ER, mDM, mProton));

	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(dRdER_Nucleus(ER, DM, SHM, hydrogen), result);
}

TEST(TestDirectDetectionNucleus, TestdRdERNucleusElement)
{
	// ARRANGE
	double mDM	   = 10.0 * GeV;
	double sigma_p = 1.0 * pb;
	DM_Particle_SI DM(mDM);
	DM.Set_Sigma_Proton(sigma_p);
	DM.Set_Low_Mass_Mode(true);
	Standard_Halo_Model SHM;
	std::vector<Isotope> hydrogen = {Isotope(1, 1, 0.5), Isotope(1, 2, 0.5)};
	double ER					  = 2.0 * keV;
	double result				  = SHM.DM_density / mDM * sigma_p / 2.0 / libphysica::Reduced_Mass(mDM, mProton) / libphysica::Reduced_Mass(mDM, mProton) * 0.5 * (SHM.Eta_Function(vMinimal_Nucleus(ER, mDM, mProton)) + 4.0 * SHM.Eta_Function(vMinimal_Nucleus(ER, mDM, 2 * mNucleon)));

	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(dRdER_Nucleus(ER, DM, SHM, Element(hydrogen)), result);
}

TEST(TestDirectDetectionNucleus, TestDefaultConstructor)
{
	// ARRANGE
	DM_Detector_Nucleus detector;
	// ACT & ASSERT
	ASSERT_EQ(detector.name, "Nuclear recoil experiment");
}

// TEST(TestDirectDetectionNucleus, TestMaximumEnergyDeposit)
// {
// 	// ARRANGE
// 	double exposure = 1.0*kg*day;
// 	std::vector<Element> elements = {Element(Isotope(8,16))};
// 	double mT = 16*mNucleon;
// 	DM_Detector_Nucleus detector("Test", exposure, elements);

// 	double mDM = 10.0 * GeV;
// 	DM_Particle_SI DM(mDM);
// 	Standard_Halo_Model SHM;
// 	// ACT & ASSERT
// 	ASSERT_DOUBLE_EQ(detector.Maximum_Energy_Deposit(DM,SHM), 2.0 * pow(libphysica::Reduced_Mass(mDM,mT)*SHM.Maximum_DM_Speed(),2.0)/mT);
// }

// TEST(TestDirectDetectionNucleus, TestMinimumDMMass)
// {
// 	// ARRANGE

// 	// ACT & ASSERT
// }

TEST(TestDirectDetectionNucleus, TestSimpleLimit)
{
	// ARRANGE
	double mDM = 10.0 * GeV;
	DM_Particle_SI DM(mDM);

	Standard_Halo_Model SHM;

	Element target(Isotope(54, 131));
	DM_Detector_Nucleus detector("Test", kg * day, {target});
	detector.Use_Energy_Threshold(3 * keV, 30 * keV);

	double tol = 1e-6;
	// ACT
	double upper_limit_90 = detector.Upper_Limit(DM, SHM, 0.9);
	double upper_limit_95 = detector.Upper_Limit(DM, SHM, 0.95);
	// ASSERT
	DM.Set_Sigma_Proton(upper_limit_90);
	EXPECT_NEAR(detector.DM_Signals_Total(DM, SHM), log(10), tol);
	DM.Set_Sigma_Proton(upper_limit_95);
	EXPECT_NEAR(detector.DM_Signals_Total(DM, SHM), log(20), tol);
}

// TEST(TestDirectDetectionNucleus, ImportEfficiency1)
// {
// 	// ARRANGE

// 	// ACT & ASSERT
// }

// TEST(TestDirectDetectionNucleus, ImportEfficiency2)
// {
// 	// ARRANGE

// 	// ACT & ASSERT
// }

// TEST(TestDirectDetectionNucleus, Minimum_DM_Speed)
// {
// 	// ARRANGE

// 	// ACT & ASSERT
// }

// TEST(TestDirectDetectionNucleus, dRdE)
// {
// 	// ARRANGE

// 	// ACT & ASSERT
// }

// TEST(TestDirectDetectionNucleus, PrintSummary)
// {
// 	// ARRANGE

// 	// ACT & ASSERT

// }
