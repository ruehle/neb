#include <gtest/gtest.h>
#include <neb/neb_fortran_interface.h>

#include "LepsPotential.h"

using namespace pele;

double leps_potential(int ncoords, double *coords, double *grad, void *userdata)
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

	LepsPotential *pot = static_cast<LepsPotential *>(userdata);

	return pot->get_energy_gradient(
		Array<double>(coords, ncoords),
		Array<double>(grad, ncoords)
		);
}

TEST(FortranInterfaceTest, SetupAndCleanup_GivesNoError)
{
	LepsPotential pot;
	neb_setup(&leps_potential, static_cast<void *>(&pot));
	neb_cleanup();
}

TEST(FortranInterfaceTest, InitializePath_GivesNoErrors)
{
	LepsPotential pot;
	neb_setup(&leps_potential, static_cast<void *>(&pot));
	neb_initialize_path(3, 2);
	neb_cleanup();
}

TEST(FortranInterfaceTest, SetAndRetrieveImageCoordinates_Agree)
{
	LepsPotential pot;
	neb_setup(&leps_potential, static_cast<void *>(&pot));
	neb_initialize_path(3, 2);

	double x_in[2] = { 1.23, 4.56 };
	neb_set_image_coords(0, 2, x_in);
	double x_out[2];
	neb_get_image_coords(0, 2, x_out);
	neb_cleanup();

	ASSERT_DOUBLE_EQ(x_in[0], x_out[0]);
	ASSERT_DOUBLE_EQ(x_in[1], x_out[1]);
}