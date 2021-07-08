#include "gtest/gtest.h"

#include "libphysica/Natural_Units.hpp"
#include "obscura/DM_Particle_Standard.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

TEST(TestDMParticleStandard, TestDefaultConstructor)
{
	// ARRANGE
	DM_Particle_Standard dm;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(dm.mass, 10.0);
}

TEST(TestDMParticleStandard, TestFixedRatios1)
{
	// ARRANGE
	DM_Particle_SI dm;
	dm.Set_Sigma_Proton(pb);
	double ratio = 0.25;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(dm.Sigma_Proton(), dm.Sigma_Neutron());
	dm.Fix_Coupling_Ratio(1.0, 8.0);
	EXPECT_DOUBLE_EQ(dm.Sigma_Proton(), 1.0 / 64. * dm.Sigma_Neutron());
	dm.Fix_fn_over_fp(ratio);
	EXPECT_DOUBLE_EQ(dm.Sigma_Neutron(), ratio * ratio * dm.Sigma_Proton());
	dm.Fix_fp_over_fn(ratio);
	EXPECT_DOUBLE_EQ(dm.Sigma_Proton(), ratio * ratio * dm.Sigma_Neutron());
	dm.Fix_Coupling_Ratio(1.0, 0.0);
	EXPECT_DOUBLE_EQ(dm.Sigma_Neutron(), 0.0);
	EXPECT_NE(dm.Sigma_Proton(), 0.0);
	dm.Fix_Coupling_Ratio(0.0, 1.0);
	EXPECT_DOUBLE_EQ(dm.Sigma_Proton(), 0.0);
	EXPECT_NE(dm.Sigma_Neutron(), 0.0);
}

TEST(TestDMParticleStandard, TestFixedRatios2)
{
	// ARRANGE
	DM_Particle_SI dm;
	dm.Set_Sigma_Neutron(pb);
	double ratio = 0.25;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(dm.Sigma_Proton(), dm.Sigma_Neutron());
	dm.Fix_Coupling_Ratio(1.0, 8.0);
	EXPECT_DOUBLE_EQ(dm.Sigma_Proton(), 1.0 / 64. * dm.Sigma_Neutron());
	dm.Fix_fn_over_fp(ratio);
	EXPECT_DOUBLE_EQ(dm.Sigma_Neutron(), ratio * ratio * dm.Sigma_Proton());
	dm.Fix_fp_over_fn(ratio);
	EXPECT_DOUBLE_EQ(dm.Sigma_Proton(), ratio * ratio * dm.Sigma_Neutron());
	dm.Fix_Coupling_Ratio(1.0, 0.0);
	EXPECT_DOUBLE_EQ(dm.Sigma_Neutron(), 0.0);
	EXPECT_NE(dm.Sigma_Proton(), 0.0);
}

TEST(TestDMParticleStandard, TestFnOverFp)
{
	// ARRANGE
	DM_Particle_SI dm;
	dm.Set_Sigma_Neutron(pb);
	double ratio = 0.25;
	// ACT & ASSERT
	dm.Fix_fn_over_fp(ratio);
	EXPECT_DOUBLE_EQ(dm.Sigma_Neutron(), ratio * ratio * dm.Sigma_Proton());
	dm.Fix_fn_over_fp(0.0);
	EXPECT_DOUBLE_EQ(dm.Sigma_Neutron(), 0.0);
	EXPECT_NE(dm.Sigma_Proton(), 0.0);
}

TEST(TestDMParticleStandard, TestFpOverFn)
{
	// ARRANGE
	DM_Particle_SI dm;
	dm.Set_Sigma_Neutron(pb);
	double ratio = 0.25;
	// ACT & ASSERT
	dm.Fix_fp_over_fn(ratio);
	EXPECT_DOUBLE_EQ(dm.Sigma_Proton(), ratio * ratio * dm.Sigma_Neutron());
	dm.Fix_fp_over_fn(0.0);
	EXPECT_DOUBLE_EQ(dm.Sigma_Proton(), 0.0);
	EXPECT_NE(dm.Sigma_Neutron(), 0.0);
}

TEST(TestDMParticleSI, TestDefaultConstructor)
{
	// ARRANGE
	DM_Particle_SI dm;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(dm.mass, 10.0);
	ASSERT_DOUBLE_EQ(dm.Sigma_Proton(), 1e-40 * cm * cm);
	ASSERT_DOUBLE_EQ(dm.Sigma_Neutron(), 1e-40 * cm * cm);
}

TEST(TestDMParticleSI, TestConstructor1)
{
	// ARRANGE
	DM_Particle_SI dm(0.5);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(dm.mass, 0.5);
	ASSERT_DOUBLE_EQ(dm.Sigma_Proton(), 1e-40 * cm * cm);
	ASSERT_DOUBLE_EQ(dm.Sigma_Neutron(), 1e-40 * cm * cm);
}

TEST(TestDMParticleSI, TestConstructor2)
{
	// ARRANGE
	DM_Particle_SI dm(0.5, pb);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(dm.mass, 0.5);
	ASSERT_DOUBLE_EQ(dm.Sigma_Proton(), pb);
	ASSERT_DOUBLE_EQ(dm.Sigma_Neutron(), pb);
}

TEST(TestDMParticleSI, TestSetFormFactorDM)
{
	// ARRANGE
	DM_Particle_SI dm(0.5, pb);
	double q		 = MeV;
	double vDM		 = 1e-3;
	double q0		 = aEM * mElectron;
	double mMediator = GeV;
	Isotope target	 = Get_Isotope(54, 131);
	double dsdq2_e	 = dm.dSigma_dq2_Electron(q, vDM);
	double dsdq2_n	 = dm.dSigma_dq2_Nucleus(q, target, vDM);
	// ACT & ASSERT
	dm.Set_FormFactor_DM("Electric-Dipole");
	EXPECT_DOUBLE_EQ(dm.dSigma_dq2_Electron(q, vDM), pow(q0 / q, 2) * dsdq2_e);
	dm.Set_FormFactor_DM("Long-Range");
	EXPECT_DOUBLE_EQ(dm.dSigma_dq2_Electron(q, vDM), pow(q0 / q, 4) * dsdq2_e);
	dm.Set_FormFactor_DM("General", mMediator);
	EXPECT_DOUBLE_EQ(dm.dSigma_dq2_Electron(q, vDM), pow((q0 * q0 + mMediator * mMediator) / (q * q + mMediator * mMediator), 2) * dsdq2_e);

	dm.Set_FormFactor_DM("Electric-Dipole");
	EXPECT_DOUBLE_EQ(dm.dSigma_dq2_Nucleus(q, target, vDM), pow(q0 / q, 2) * dsdq2_n);
	dm.Set_FormFactor_DM("Long-Range");
	EXPECT_DOUBLE_EQ(dm.dSigma_dq2_Nucleus(q, target, vDM), pow(q0 / q, 4) * dsdq2_n);
	dm.Set_FormFactor_DM("General", mMediator);
	EXPECT_DOUBLE_EQ(dm.dSigma_dq2_Nucleus(q, target, vDM), pow((q0 * q0 + mMediator * mMediator) / (q * q + mMediator * mMediator), 2) * dsdq2_n);
}

TEST(TestDMParticleSI, TestSetMass)
{
	// ARRANGE
	DM_Particle_SI dm;
	double sigma_n = pb;
	double mDM	   = 100.0;
	dm.Unfix_Coupling_Ratios();
	dm.Set_Sigma_Proton(0.0);
	dm.Set_Sigma_Neutron(sigma_n);
	// ACT
	dm.Set_Mass(mDM);
	// ASSERT
	EXPECT_DOUBLE_EQ(dm.mass, mDM);
	EXPECT_DOUBLE_EQ(dm.Sigma_Proton(), 0.0);
	EXPECT_DOUBLE_EQ(dm.Sigma_Neutron(), sigma_n);
}

TEST(TestDMParticleSI, TestPrintSummary)
{
	// ARRANGE
	DM_Particle_SI dm;
	// ACT & ASSERT
	dm.Print_Summary();
}

// 3. Spin-dependent (SD) interactions
TEST(TestDMParticleSD, TestDefaultConstructor)
{
	// ARRANGE
	DM_Particle_SD dm;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(dm.mass, 10.0);
	ASSERT_DOUBLE_EQ(dm.Sigma_Proton(), 1e-40 * cm * cm);
	ASSERT_DOUBLE_EQ(dm.Sigma_Neutron(), 1e-40 * cm * cm);
}

TEST(TestDMParticleSD, TestConstructor1)
{
	// ARRANGE
	DM_Particle_SD dm(0.5);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(dm.mass, 0.5);
	ASSERT_DOUBLE_EQ(dm.Sigma_Proton(), 1e-40 * cm * cm);
	ASSERT_DOUBLE_EQ(dm.Sigma_Neutron(), 1e-40 * cm * cm);
}

TEST(TestDMParticleSD, TestConstructor2)
{
	// ARRANGE
	DM_Particle_SD dm(0.5, pb);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(dm.mass, 0.5);
	ASSERT_DOUBLE_EQ(dm.Sigma_Proton(), pb);
	ASSERT_DOUBLE_EQ(dm.Sigma_Neutron(), pb);
}

TEST(TestDMParticleSD, TestSDInteractions)
{
	// ARRANGE
	DM_Particle_SD dm(0.5, pb);
	double vDM				 = 1e-3;
	double q				 = MeV;
	Isotope spinless_isotope = Get_Isotope(54, 130);
	Isotope spin_isotope	 = Get_Isotope(54, 129);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(dm.Sigma_Total_Nucleus(Get_Isotope(1, 1), vDM), pb);
	EXPECT_DOUBLE_EQ(dm.Sigma_Total_Nucleus(spinless_isotope, vDM), 0.0);
	EXPECT_GT(dm.Sigma_Total_Nucleus(spin_isotope, vDM), 0.0);
	for(int Z = 1; Z < 93; Z++)
	{
		Nucleus nucleus = Get_Nucleus(Z);
		for(auto& isotope : nucleus.isotopes)
			if(isotope.spin == 0)
			{
				EXPECT_DOUBLE_EQ(dm.dSigma_dq2_Nucleus(q, isotope, vDM), 0.0);
				EXPECT_DOUBLE_EQ(dm.Sigma_Total_Nucleus(isotope, vDM), 0.0);
			}
	}
}

TEST(TestDMParticleSD, TestPrintSummary)
{
	// ARRANGE
	DM_Particle_SD dm;
	// ACT & ASSERT
	dm.Print_Summary();
}

TEST(TestDMParticleSD, TestScatteringAnglePDF)
{
	// ARRANGE
	DM_Particle_SD dm;
	dm.Set_Sigma_Proton(pb);
	dm.Set_Sigma_Electron(0.9 * pb);
	Isotope target = Get_Isotope(54, 131);
	double vDM	   = 1e-3;

	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(dm.PDF_Scattering_Angle_Nucleus(+0.5, target, vDM), dm.PDF_Scattering_Angle_Nucleus(-0.5, target, vDM));
	EXPECT_DOUBLE_EQ(dm.PDF_Scattering_Angle_Electron(+0.5, vDM), dm.PDF_Scattering_Angle_Electron(-0.5, vDM));
}

TEST(TestDMParticleSD, TestScatteringAngleCDF)
{
	// ARRANGE
	DM_Particle_SD dm;
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

TEST(TestDMParticleSD, TestScatteringAngleSampling)
{
	// ARRANGE
	std::random_device rd;
	std::mt19937 PRNG(rd());
	DM_Particle_SD dm;
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