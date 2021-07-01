#include "obscura/DM_Particle.hpp"
#include "gtest/gtest.h"

#include <random>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Particle_Standard.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

TEST(TestDMParticle, TestDefaultConstructor)
{
	// ARRANGE
	DM_Particle dm;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(dm.mass, 10.0);
	EXPECT_DOUBLE_EQ(dm.spin, 0.5);
	EXPECT_DOUBLE_EQ(dm.fractional_density, 1.0);
	EXPECT_FALSE(dm.DD_use_eta_function);
	EXPECT_FALSE(dm.Interaction_Parameter_Is_Cross_Section());
}

TEST(TestDMParticle, TestConstructor)
{
	// ARRANGE
	double mDM	= MeV;
	double spin = 0.0;
	DM_Particle dm(mDM, spin);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(dm.mass, mDM);
	EXPECT_DOUBLE_EQ(dm.spin, spin);
	EXPECT_DOUBLE_EQ(dm.fractional_density, 1.0);
	EXPECT_FALSE(dm.DD_use_eta_function);
	EXPECT_FALSE(dm.Interaction_Parameter_Is_Cross_Section());
}

TEST(TestDMParticle, TestSetFunctions)
{
	// ARRANGE
	double mDM	= MeV;
	double spin = 0.0;
	double f	= 0.4;
	DM_Particle dm;
	// ACT
	dm.Set_Mass(mDM);
	dm.Set_Spin(spin);
	dm.Set_Fractional_Density(f);
	// ASSERT
	EXPECT_DOUBLE_EQ(dm.mass, mDM);
	EXPECT_DOUBLE_EQ(dm.spin, spin);
	EXPECT_DOUBLE_EQ(dm.fractional_density, f);
}

TEST(TestDMParticle, TestSetLowMassMode)
{
	// ARRANGE
	DM_Particle_SI dm;
	dm.Set_Interaction_Parameter(pb, "Nuclei");
	double q	   = MeV;
	Isotope target = Get_Isotope(54, 131);
	double vDM	   = 1e-3;
	double dsdq2   = dm.dSigma_dq2_Nucleus(q, target, vDM);
	// ACT
	dm.Set_Low_Mass_Mode(true);
	// ASSERT
	EXPECT_TRUE(dm.Interaction_Parameter_Is_Cross_Section());
	EXPECT_LT(dsdq2, dm.dSigma_dq2_Nucleus(q, target, vDM));
}

TEST(TestDMParticle, TestInteractionParameter)
{
	// ARRANGE
	DM_Particle_SI dm;
	// ACT
	dm.Set_Interaction_Parameter(pb, "Nuclei");
	dm.Set_Interaction_Parameter(0.5 * pb, "Electrons");
	// ASSERT
	EXPECT_DOUBLE_EQ(dm.Get_Interaction_Parameter("Nuclei"), pb);
	EXPECT_DOUBLE_EQ(dm.Get_Interaction_Parameter("Electrons"), 0.5 * pb);
}

TEST(TestDMParticle, TestReferenceCrossSections)
{
	// ARRANGE
	DM_Particle_SI dm;
	dm.Unfix_Coupling_Ratios();
	// ACT
	dm.Set_Sigma_Proton(0.2 * pb);
	dm.Set_Sigma_Neutron(0.3 * pb);
	dm.Set_Sigma_Electron(0.4 * pb);
	// ASSERT
	EXPECT_DOUBLE_EQ(dm.Sigma_Proton(), 0.2 * pb);
	EXPECT_DOUBLE_EQ(dm.Sigma_Neutron(), 0.3 * pb);
	EXPECT_DOUBLE_EQ(dm.Sigma_Electron(), 0.4 * pb);
}

TEST(TestDMParticle, TestDifferentialCrossSectionsNucleus)
{
	// ARRANGE
	DM_Particle_SI dm;
	dm.Set_Sigma_Proton(pb);
	double q	   = MeV;
	Isotope target = Get_Isotope(54, 131);
	double ER	   = q * q / 2.0 / target.mass;
	double vDM	   = 1e-3;
	// ACT & ASSERT
	EXPECT_LT(dm.dSigma_dq2_Nucleus(2.0 * q, target, vDM), dm.dSigma_dq2_Nucleus(q, target, vDM));
	dm.Set_Low_Mass_Mode(true);
	EXPECT_DOUBLE_EQ(dm.dSigma_dq2_Nucleus(2.0 * q, target, vDM), dm.dSigma_dq2_Nucleus(q, target, vDM));
	EXPECT_DOUBLE_EQ(dm.dSigma_dER_Nucleus(ER, target, vDM), 2.0 * target.mass * dm.dSigma_dq2_Nucleus(q, target, vDM));
}

TEST(TestDMParticle, TestDifferentialCrossSectionsElectron)
{
	// ARRANGE
	DM_Particle_SI dm;
	dm.Set_Sigma_Electron(pb);
	double q   = MeV;
	double vDM = 1e-3;
	// ACT
	double dsdq2 = dm.dSigma_dq2_Electron(q, vDM);
	//ACT & ASSERT
	EXPECT_DOUBLE_EQ(dm.dSigma_dq2_Electron(2.0 * q, vDM), dsdq2);
	dm.Set_Sigma_Electron(0.5 * pb);
	EXPECT_DOUBLE_EQ(dm.dSigma_dq2_Electron(2.0 * q, vDM), 0.5 * dsdq2);
}

TEST(TestDMParticle, TestTotalCrossSections)
{
	// ARRANGE
	DM_Particle_SI dm;
	dm.Set_Low_Mass_Mode(true);
	dm.Set_Sigma_Proton(pb);
	dm.Set_Sigma_Electron(0.9 * pb);
	Isotope hydrogen = Get_Isotope(1, 1);
	double vDM		 = 1e-3;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(dm.Sigma_Total_Nucleus(hydrogen, vDM), dm.Sigma_Proton());
	EXPECT_DOUBLE_EQ(dm.Sigma_Total_Electron(vDM), dm.Sigma_Electron());
}

TEST(TestDMParticle, TestPrintSummary)
{
	// ARRANGE
	DM_Particle dm;
	// ACT & ASSERT
	dm.Print_Summary();
}

TEST(TestDMParticle, TestScatteringAnglePDF)
{
	// ARRANGE
	DM_Particle_SI dm;
	dm.Set_Sigma_Proton(pb);
	dm.Set_Sigma_Electron(0.9 * pb);
	Isotope target = Get_Isotope(54, 131);
	double vDM	   = 1e-3;

	// ACT & ASSERT
	EXPECT_GT(dm.PDF_Scattering_Angle_Nucleus(+0.5, target, vDM), dm.PDF_Scattering_Angle_Nucleus(-0.5, target, vDM));
	dm.Set_Low_Mass_Mode(true);
	EXPECT_DOUBLE_EQ(dm.PDF_Scattering_Angle_Nucleus(+0.5, target, vDM), dm.PDF_Scattering_Angle_Nucleus(-0.5, target, vDM));
	EXPECT_DOUBLE_EQ(dm.PDF_Scattering_Angle_Electron(+0.5, vDM), dm.PDF_Scattering_Angle_Electron(-0.5, vDM));
}

TEST(TestDMParticle, TestScatteringAngleCDF)
{
	// ARRANGE
	DM_Particle_SI dm;
	dm.Set_Sigma_Proton(pb);
	dm.Set_Sigma_Electron(0.9 * pb);
	Isotope target = Get_Isotope(54, 131);
	double vDM	   = 1e-3;

	// ACT & ASSERT
	EXPECT_GT(dm.CDF_Scattering_Angle_Nucleus(+0.5, target, vDM), dm.CDF_Scattering_Angle_Nucleus(-0.5, target, vDM));
	EXPECT_GT(dm.CDF_Scattering_Angle_Electron(+0.5, vDM), dm.CDF_Scattering_Angle_Electron(-0.5, vDM));

	EXPECT_DOUBLE_EQ(dm.CDF_Scattering_Angle_Nucleus(-1.0, target, vDM), 0.0);
	EXPECT_DOUBLE_EQ(dm.CDF_Scattering_Angle_Nucleus(+1.0, target, vDM), 1.0);
	EXPECT_DOUBLE_EQ(dm.CDF_Scattering_Angle_Electron(-1.0, vDM), 0.0);
	EXPECT_DOUBLE_EQ(dm.CDF_Scattering_Angle_Electron(1.0, vDM), +1.0);
}

TEST(TestDMParticle, TestScatteringAngleSampling)
{
	// ARRANGE
	std::random_device rd;
	std::mt19937 PRNG(rd());
	DM_Particle_SI dm;
	dm.Set_Sigma_Proton(pb);
	dm.Set_Sigma_Electron(0.9 * pb);
	Isotope target = Get_Isotope(54, 131);
	double vDM	   = 1e-3;

	// ACT & ASSERT
	double costheta_n = 0.0;
	double costheta_e = 0.0;
	for(int i = 0; i < 100; i++)
	{
		double ct_n = dm.Sample_Scattering_Angle_Nucleus(PRNG, target, vDM);
		EXPECT_LT(ct_n, 1.0);
		EXPECT_GT(ct_n, -1.0);
		EXPECT_NE(ct_n, costheta_n);
		costheta_n = ct_n;

		double ct_e = dm.Sample_Scattering_Angle_Electron(PRNG, vDM);
		EXPECT_LT(ct_e, 1.0);
		EXPECT_GT(ct_e, -1.0);
		EXPECT_NE(ct_e, costheta_e);
		costheta_e = ct_e;
	}
}