#include <gmock/gmock.h>
#include <neb/neb_fortran_interface.h>

#include "LepsPotential.h"

using namespace pele;
using namespace testing;


template<typename PotentialType>
double potential_callback(int ncoords, double *coords, double *grad, void *userdata)
{
	if (ncoords != ncoords)
	{
		throw std::logic_error("inconsistent number of coordinates in potential call");
	}

	if (coords == nullptr) 
	{
		throw std::logic_error("potential call with null pointer for coordinates");
	}

	if (grad == nullptr)
	{
		throw std::logic_error("potential call with null pointer for gradient");
	}

	if (userdata == nullptr)
	{
		throw std::logic_error("userdata in potential call is zero");
	}

	PotentialType *pot = static_cast<PotentialType *>(userdata);

	return pot->get_energy_gradient(
		Array<double>(coords, ncoords),
		Array<double>(grad, ncoords)
		);
}

class MockPotential : public BasePotential
{
public:
	MOCK_METHOD1(get_energy, double(Array<double> x));
	MOCK_METHOD2(get_energy_gradient, double(Array<double> x, Array<double> grad));
};


TEST(FortranInterfaceTest, SetupAndCleanup_GivesNoError)
{
	LepsPotential pot;
	neb_setup(&potential_callback<LepsPotential>, static_cast<void *>(&pot));
	neb_cleanup();
}

TEST(FortranInterfaceTest, InitializePath_GivesNoErrors)
{
	LepsPotential pot;
	neb_setup(&potential_callback<LepsPotential>, static_cast<void *>(&pot));
	neb_initialize_path(3, 2);
	neb_cleanup();
}

TEST(FortranInterfaceTest, SingleStep_CallsPotential)
{
	MockPotential pot;

	EXPECT_CALL(pot, get_energy_gradient(_, _))
		.Times(AtLeast(3));

	neb_setup(&potential_callback<MockPotential>, static_cast<void *>(&pot));
	neb_initialize_path(3, 1);
	double x = 0;
	neb_set_image_coords(0, 1, &x);
	neb_set_image_coords(1, 1, &x);
	neb_set_image_coords(2, 1, &x);
	neb_start();
	neb_step();
	neb_cleanup();
}

TEST(FortranInterfaceTest, RunNebTillEnd_WithoutErro)
{
	LepsPotential pot;
	neb_setup(&potential_callback<LepsPotential>, static_cast<void *>(&pot));
	neb_initialize_path(3, 2);

	double x_in[2] = { 1.23, 4.56 };
	neb_set_image_coords(0, 2, x_in);
	double x_out[2];
	neb_get_image_coords(0, 2, x_out);
	neb_cleanup();

	ASSERT_DOUBLE_EQ(x_in[0], x_out[0]);
	ASSERT_DOUBLE_EQ(x_in[1], x_out[1]);
}