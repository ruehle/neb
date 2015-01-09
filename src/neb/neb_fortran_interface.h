#pragma once

extern "C" {

	typedef double potential_callback_t(int ncoords, double *coords, double *grad, void 
		*userdata);

	void neb_setup(potential_callback_t *potential, void *userdata);
	void neb_cleanup();

	void neb_initialize_path(int nimages, int num_coords_per_image);
	void neb_set_image_coords(int image, int ncoords, const double *coords);
	void neb_get_image_coords(int image, int ncoords, double *coords);

	void neb_start();
	bool neb_step();
	double neb_rms();
}
