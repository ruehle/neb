#include <fenv.h>
#include <signal.h>
#include "neb.h"
#include "neb_fortran_interface.h"

using namespace cpp_neb;

// nullptr is used, but it is only defined in c++11
#define nullptr NULL;


namespace {
	NEB *g_neb = nullptr;

	class PotentialWrapper : public BasePotential
	{
	public:
		PotentialWrapper(potential_callback_t *potentialCallback, void *userdata)
			: _potentialCallback(potentialCallback), _userdata(userdata) {}

		double get_energy(Array<double> x)
		{
			throw std::logic_error("The method or operation is not implemented.");
		}

		double get_energy_gradient(Array<double> x, Array<double> grad)
		{
			return (*_potentialCallback)(x.size(), x.data(), grad.data(), _userdata);
		}

		potential_callback_t *_potentialCallback;
		void *_userdata;
	};

}

void neb_setup(potential_callback_t *potential, distance_callback_t *distance, void *userdata)
{
	if (g_neb)
		neb_cleanup();
	g_neb = new NEB(new PotentialWrapper(potential, userdata), new DistanceWrapper(distance));
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
//	std::cout << "nimages in initialise_path: " << nimages << std::endl;
//	std::cout << "number of coords per image: " << num_coords_per_image << std::endl;
	vector < Array<double> > path;
	Array<double> image(num_coords_per_image, 0);
	for (int i = 0; i < nimages; ++i) {
		path.push_back(image.copy());
	}
	g_neb->set_path(path);
}

void neb_set_image_coords(int image, int ncoords, const double *coords)
{
	Array<double> image_coords = g_neb->images()[image];
	for (int i = 0; i < ncoords; ++i)
	{
		image_coords[i] = coords[i];
	}
}

void neb_parameters(int double_nudging, double rmstol, double k_initial, double adjust_k_tol,
		double adjust_k_factor, double maxstep, int maxiter, int iprint, int verbosity)
{
	if (double_nudging==0) {
		g_neb->set_double_nudging(false);
	} else {
		g_neb->set_double_nudging(true);
	}

	g_neb->set_k(k_initial);
	g_neb->set_k_tol(adjust_k_tol);
	g_neb->set_k_factor(adjust_k_factor);
	g_neb->set_verbosity(verbosity);

	g_neb->get_optimizer()->set_tol(rmstol);
	g_neb->get_optimizer()->set_max_iter(maxiter);
	g_neb->get_optimizer()->set_maxstep(maxstep);
	g_neb->get_optimizer()->set_iprint(iprint);
}

void neb_get_image_coords(int image, int ncoords, double *coords)
{
    Array<double> image_coords = g_neb->images()[image];
	for (int i = 0; i < ncoords; ++i)
	{
		coords[i] = image_coords[i];
	}
}

void neb_get_image_energies(int nimages, double *energies)  // sn402: added
{
    Array<double> image_energies = g_neb->energies();
	for (int i=0; i<nimages; ++i)
		energies[i] = image_energies[i];
}

void neb_start()  // sn402: this should be superseded by a start function that specifies
				  // the desired optimizer (currently only lbfgs is coded).
{
	g_neb->start();
}

void neb_start_with_lbfgs(double rmstol, int setM, double max_f_rise, double H0)  // sn402: added
{
	g_neb->start_with_lbfgs(rmstol, setM, max_f_rise, H0);
}


bool neb_step()
{
	return g_neb->step();
}

void neb_adjust_k()
{
	g_neb->adjust_k();
}

double neb_rms()
{
        return g_neb->get_rms();
}
