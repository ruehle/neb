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