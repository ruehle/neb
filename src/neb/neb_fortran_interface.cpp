#include "neb.h"
#include "neb_fortran_interface.h"

using namespace pele;

namespace {
	NEB *g_neb = nullptr;

	class PotentialWrapper : public BasePotential
	{
	public:
		PotentialWrapper(potential_callback_t *potentialCallback, void *userdata)
			: _potentialCallback(potentialCallback), _userdata(userdata) {}

		double get_energy(Array<double> x) override
		{
			throw std::logic_error("The method or operation is not implemented.");
		}

		double get_energy_gradient(Array<double> x, Array<double> grad) override
		{
			return (*_potentialCallback)(x.size(), x.data(), grad.data(), _userdata);
		}

		potential_callback_t *_potentialCallback;
		void *_userdata;
	};
}

void neb_setup(potential_callback_t *potential, void *userdata)
{
	if (g_neb)
		neb_cleanup();
	g_neb = new NEB(new PotentialWrapper(potential, userdata));
}

void neb_cleanup()
{
	if (!g_neb)
		throw std::logic_error("trying to cleanup neb which was not initialized");
	delete g_neb;
	g_neb = nullptr;
}

void neb_initialize_path(int nimages, int num_coords_per_image)
{
	vector < Array<double> > path;
	Array<double> image(num_coords_per_image, 0);
	for (int i = 0; i < nimages; ++i) {
		path.push_back(image.copy());
	}
	g_neb->set_path(path);
}

void neb_set_image_coords(int image, int ncoords, const double *coords)
{
	auto image_coords = g_neb->images()[image];
	for (int i = 0; i < ncoords; ++i)
	{
		image_coords[i] = coords[i];
	}
}

void neb_get_image_coords(int image, int ncoords, double *coords)
{
	auto image_coords = g_neb->images()[image];
	for (int i = 0; i < ncoords; ++i)
	{
		coords[i] = image_coords[i];
	}
}


void neb_start()
{
	g_neb->start();
}

bool neb_step()
{
	return g_neb->step();
}



