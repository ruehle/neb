#pragma once

extern "C" {

	typedef double potential_callback_t(int ncoords, double *coords, double *grad, void 
		*userdata);

	void neb_setup(potential_callback_t *potential, void *userdata);
	void neb_cleanup();

	void neb_start();
	bool neb_step();
}
