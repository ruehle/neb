#pragma once

extern "C" {

	typedef double potential_callback_t(int ncoords, double *coords, double *grad, void 
		*userdata);

	typedef double distance_callback_t(int ncoords, double *left, double *right,
			double *gradient_left, double *gradient_right);

	void neb_setup(potential_callback_t *potential, distance_callback_t *distance, void *userdata);
	void neb_distance(distance_callback_t *distance);
	void neb_cleanup();

	void neb_initialize_path(int nimages, int num_coords_per_image);
	void neb_set_image_coords(int image, int ncoords, const double *coords);
	void neb_get_image_coords(int image, int ncoords, double *coords);

	// sn402
	void neb_get_image_energies(int nimages, double *energies);
	void neb_parameters(int double_nudging, double rmstol, double k_initial, double adjust_k_tol,
			double adjust_k_factor, double maxstep, int maxiter, int iprint, int verbosity);

	void neb_start();  // sn402: this should be superseded by a start function that specifies
	// the desired optimizer (currently only lbfgs is coded).
	void neb_start_with_lbfgs(double rmstol, int setM, double H0); // sn402 added
	bool neb_step();
	void neb_adjust_k();
	double neb_rms();
}
