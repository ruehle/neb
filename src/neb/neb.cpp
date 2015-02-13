#include <exception>
#include "neb.h"
#include <algorithm>
#include "lbfgs.h"
#include <cmath>

namespace pele {

void NEB::set_path(std::vector< Array<double> > path)
{
	// check if there are images in the path
	if(path.size() == 0)
		throw std::runtime_error("cannot initialize neb with empty path");

	int nimages = _nimages = path.size();
	int N = _N = path[0].size();

	// check if all images have the same number of coordinates
	for(int i=0; i<nimages; ++i)
		if(path[i].size() != N) 
			throw std::runtime_error("number of coordinates in path images differs");

	// reset the old image shortcuts and allocate memory to store new path
	_images.clear();
	_coords.resize(N*nimages);

	// generate array views for quick access
	for(int i=0, j=0; i<nimages; ++i, j+=N) {
		// generate a view
		Array<double> image(_coords.view(j, j+N));
		// copy data
		image.assign(path[i]);
		_images.push_back(image);
	}

	// update all variables where the work is performed in
	// e.g. distances, tangent vectors, ...
	adjust_worker_variables();
	_nimages = nimages;
	_N = N;
//	std::cout << "Path set up\n";
//	for (int i=0; i<nimages;i++) {
//		for(int j=0; j<N;j+=3) {
//			std::cout << _images[i][j] << "\t" << _images[i][j+1] << "\t" _images[i][j+2] << std::endl;
//		}
//	}
}

void NEB::resize_array_collection(std::vector< Array<double> > &items, size_t size, size_t nelements)
{
	if(items.size() > size) {
		items.resize(size);
	} else while(items.size() < size-2) {
		items.push_back(Array<double>(nelements));
	}
}

void resize_array_vector(vector< Array<double> > &x, int nimages, int n)
{
	if(x.size() > nimages)
		x.resize(nimages);
	else while(x.size() < nimages) {
		x.push_back(Array<double>(n));
	}
}

void NEB::adjust_worker_variables()
{
	_energies.resize(_nimages);
	_distances.resize(_nimages-1);

	resize_array_vector(_tangents, _nimages-2, _N);
	resize_array_vector(_true_gradients, _nimages, _N);
	resize_array_vector(_tau_left, _nimages, _N);
	resize_array_vector(_tau_right, _nimages, _N);
}

void NEB::adjust_k()  // sn402
{
	// Compute the mean distance between consecutive images
	std::cout << "Distances:" << std::endl;
	double average_d = 0;
	for (int i=0;i<_nimages-1;i++)
	{
		std::cout<<_distances[i]<<std::endl;
		average_d += _distances[i];
	}
	average_d /= (_nimages-1);
	std::cout << "Average distance: " << average_d << std::endl;
	// Compute the average deviation of the distances (normalised to the average distance)
	std::cout << "Deviations:\n";
	std::vector<double> deviations;
	double ave_dev = 0;
	for (int i=0;i<_nimages-1;i++)
	{
		deviations.push_back(fabs((_distances[i]-average_d)/average_d));
		ave_dev += deviations[i];
		std::cout<<deviations[i]<<std::endl;
	}
	ave_dev /= (_nimages-1);

	// If this average deviation is larger than a specified tolerance, the band is too loose
	// so we increase the force constant. If the deviation is smaller, we relax the force constant.
	if (ave_dev > _adjust_k_tol)
	{
		_k *= _adjust_k_factor;
		if (_verbosity > 0)
			std::cout << "Increasing NEB force constant to " << _k <<
					" (average deviation is " << ave_dev << ")" << std::endl;
	}
	else
	{
		_k /= _adjust_k_factor;
		if (_verbosity > 0)
			std::cout << "Decreasing NEB force constant to " << _k <<
					" (average deviation is " << ave_dev << ")" << std::endl;
	}

}

std::vector< Array<double> > NEB::generate_image_views(Array<double> coords)
{
	std::vector< Array<double> > images;
	// generate array views for quick access
	for(size_t i=0; i<_nimages; ++i)
		images.push_back(coords.view(i*_N, (i+1)*_N));
	return images;
}


double NEB::get_energy(Array<double> coords)
{	
	// first wrap coordinates for convenient access
	std::vector< Array<double> > images
		= generate_image_views(coords);
	
	// now loop over all images and sum up the energy
	double energy = 0.0;
	for(size_t i=0; i<_images.size(); ++i)
	{ // sn402
		double old_energy=energy; // sn402
		energy += _potential->get_energy(images[i]);
		std::cout << "Energy for image " << i << " is " << energy-old_energy << std::endl;
	} // sn402
	return energy;
}

double NEB::get_energy_gradient(Array<double> coords, Array<double> grad)
{
	if (grad.size() != coords.size()) {
		throw std::runtime_error("coords and grad should have the same size");
	}
	if (coords.size() != _N * _nimages) {
		throw std::runtime_error("coords has the wrong size");
	}
	std::cout << "cpp energy gradient function\n";  // sn402

	// first wrap coordinates for convenient access
	std::vector< Array<double> > images
		= generate_image_views(coords);
	// same for the gradients
	std::vector< Array<double> > image_gradients
		= generate_image_views(grad);

	// calculate the true energy and gradient
	grad.assign(0);
	double energy = 0.0;
	for(size_t i=0; i<_images.size(); ++i) {
//		std::cout << "Computing energy for image " << i <<"\n";  // sn402
		//std::cout << "image coordinates";  // sn402
		//std::cout << images[i];  // sn402
		// calculate the true gradient
		_energies[i] = _potential->get_energy_gradient(images[i], _true_gradients[i]);
		energy += _energies[i];
		//jdf43
//		std::cout << "energy " << i << " " << _energies[i] << std::endl;  // sn402
	}

	// update distances and tangents
	update_distances(images, true);

	for(size_t i=1; i<_images.size()-1; ++i) {

		Array<double> & tangent = _tangents[i-1];

		// perpendicular part of true gradient
		double project = dot(_true_gradients[i], _tangents[i-1]);
		for(size_t j=0; j<_N; ++j)
			image_gradients[i][j] = _true_gradients[i][j] - project*_tangents[i-1][j];

		// double nudging
		// we first do the double nudging since in this case image_gradient still contains
		// the perpendicular part of the true gradient which saves an extra storage variable
		Array<double> spring(_N);
		if(_double_nudging) {
			    for(size_t j=0; j<_N; ++j)
					spring[j] = _k*(_tau_left[i-1][j] + _tau_right[i-1][j]);
				std::cout << "spring " << norm(spring) << std::endl;
				// first project out parallel part since this this is treated separately
				// in the normal nudging
				double project1 = dot(spring, _tangents[i-1]);
				for(size_t j=0; j<_N; ++j)
					spring[j] -= project1 * _tangents[i-1][j];

				// project out the part which goes along the direction of the true gradient
				double project2 = dot(spring,image_gradients[i])/dot(image_gradients[i],image_gradients[i]);
				for(size_t j=0; j<_N; ++j)
					spring[j] -= project2*image_gradients[i][j];

                // add the spring force
				for(size_t j=0; j<_N; ++j)
					image_gradients[i][j] += spring[j];
		}

		// spring force
		double d = _k * (_distances[i-1] - _distances[i]);

		double temp1 = 0; // rms real grad
		double temp2 = 0; // rms spring force
		for(size_t j=0; j<_N; ++j) {
			temp1 += _true_gradients[i][j]* _true_gradients[i][j];
			temp2 += (d*_tangents[i-1][j])*(d*_tangents[i-1][j]);
		}
		temp1 = sqrt(temp1);
		temp2 = sqrt(temp2);

		for(size_t j=0; j<_N; ++j)
			image_gradients[i][j] += d*_tangents[i-1][j];

		std::cout << _k
				<< " " << _distances[i-1]
				<< " " << _distances[i]
				<< " " << _distances[i-1] - _distances[i]
				<< " " << d
				<< " " << temp1
				<< " " << temp2
				<< std::endl;

	}

	// the gradient of first and last image is zero
	image_gradients[0].assign(0.0);
	image_gradients[_nimages-1].assign(0.0);

//	double temp3 = 0;
//	for(int j=1; j<_images.size(); j++)
//		temp3 += 0.5*_k*_distances[j-1]*_distances[j-1];
//	std::cout << "Spring energy 1: " << temp3 <<"\n";
//
//	double temp4 = 0;
//	for(int i=0; i<_images.size(); i++) {
//		for(int j=0; j<_N; j++)
//			temp4 += 0.5*_k*(_tau_left[i][j]+_tau_right[i][j])*(_tau_left[i][j]+_tau_right[i][j]);
//	}
//
//	std::cout << "Spring energy 2: " << temp4 <<"\n";

	std::cout << "Gradients\n";
	for(size_t i=0; i<_images.size()-1; i++) {
//		std::cout << "Image " << i << std::endl;
		for (size_t j=0; j<_N; j+=3)
			std::cout << grad[i*_N+j] << "\t" << grad[i*_N+j+1] << "\t" <<grad[i*_N+j+2] << std::endl;
	}
	//jdf43 std::cout << "rms: " << get_rms() << std::endl;
	return energy;//+temp3;
}

void NEB::update_distances(std::vector< Array<double> > images, bool update_tangent)
{
	// adjust array to store distances + tangents
	_distances.resize(images.size()-1);
	// TODO: this probably won't work
	//_tangents.resize(images.size()-2);

	// these arrays store the left and right distance gradient
	Array<double> tau_left(_N);
	Array<double> tau_right(_N);
	
	// these variables store the previous distance calculation to perform the tangent interpolation
	Array<double> tau_save(_N);

	// do first distance calculation, we cannot do tangent interpolation yet since
	// we only have distances between one pair
	_distances[0] = _distance->get_distance(images[0], images[1], tau_left, tau_right);
	// store distance gradient from previous calculation
	tau_save.assign(tau_right);

	// calculate the distances
	for(size_t i=1; i<_nimages-1; ++i) {
		Array<double> tau(_N);
		// calculate distance to next image
//		std::cout << "Calculating distances between images " << i << " and " << i+1 << std::endl;  // sn402
//		std::cout << "tau_save" << tau_save << std::endl;
//		std::cout << "tau_left" << tau_left << std::endl;

		_distances[i] = _distance->get_distance(images[i], images[i+1], tau_left, tau_right);

		_tau_left[i-1].assign(tau_save);
		_tau_right[i-1].assign(tau_left);
//std::cout << "tau_save" << tau_save << std::endl;
//std::cout << "tau_left" << tau_left << std::endl;
		// interpolate tangent based on previous step
		interpolate_tangent(_tangents[i-1], _energies[i], _energies[i-1],
							_energies[i+1], tau_save, tau_left);
//		interpolate_tangent(_tangents[i-1], _energies[i], _energies[i-1],
//							_energies[i+1], tau_save, tau_right);

		// save latest distance calculation
//		std::cout << "tau_overlap" << dot(tau_save,tau_right)/norm(tau_left)/norm(tau_right)<< std::endl;
		tau_save.assign(tau_right);
	}
}


/*
    New uphill tangent formulation

	The method was  described in
	"Improved tangent estimate in the nudged elastic band method for finding
	minimum energy paths and saddle points"
	Graeme Henkelman and Hannes Jonsson
	J. Chem. Phys 113 (22), 9978 (2000)
*/
void NEB::interpolate_tangent(Array<double> tau, double energy, double energy_left, double energy_right, 
							  Array<double> tau_left, Array<double> tau_right)
{      
    // if central point is a maximum or minimum, then interpolate gradients
	if((energy >= energy_left && energy >= energy_right) || (energy <= energy_left && energy <= energy_right)) {
		double vmax = std::max(std::abs(energy - energy_left), std::abs(energy - energy_right));
		double vmin = std::min(std::abs(energy - energy_left), std::abs(energy - energy_right));

		if(energy_left > energy_right)
			for( size_t k=0; k<tau.size(); ++k)				
                tau[k] = vmax * tau_left[k] - vmin*tau_right[k];
		else
			for( size_t k=0; k<tau.size(); ++k)				
                tau[k] = vmin * tau_left[k] - vmax*tau_right[k];
	}
	// otherwise take the tangent which points to the image higher in energy
	else if (energy_left > energy_right) {
		tau.assign(tau_left);
	}
	else {
		tau.assign(tau_right);
		tau *= -1.0;
	}

	// normalize tangent vector
	double inv_norm = 1./norm(tau);
	tau *= inv_norm;
//	for(size_t i=0; i<tau.size(); ++i)
//		tau[i] *= inv_norm;
}

void NEB::start()
{
	LBFGS *lbfgs = new LBFGS(this, this->_coords);
//	lbfgs->set_max_f_rise(0.1);

	_optimizer = lbfgs;

}

// Currently, lbfgs is the only optimiser implemented and is hard-coded into the NEB.
// So we always call this function.
// When a new optimiser is implemented, we will need to call a different parameters function
// depending on which optimiser is being used.
void NEB::set_lbfgs_parameters(int setM, double max_f_rise, double H0)
{
	LBFGS *lbfgs = (LBFGS*) _optimizer;
	lbfgs->set_max_f_rise(max_f_rise);
	lbfgs->set_H0(H0);
	lbfgs->set_M(setM);
}

void NEB::start_with_lbfgs(double rmstol, int setM, double max_f_rise, double H0)
{
	LBFGS *lbfgs = new LBFGS(this, this->_coords, rmstol, setM);
	lbfgs->set_max_f_rise(max_f_rise);
	lbfgs->set_H0(H0);
	_optimizer = lbfgs;
}

void NEB::set_parameters(int double_nudging, double rmstol, double k_initial, double adjust_k_tol,
		double adjust_k_factor, double maxstep, int maxiter, int iprint, int verbosity)
{
	_optimizer->set_tol(rmstol);
	_optimizer->set_maxstep(maxstep);
	_optimizer->set_max_iter(maxiter);
	_optimizer->set_iprint(iprint);
	_optimizer->set_verbosity(verbosity);
	_verbosity = verbosity;
	_k = k_initial;
	_adjust_k_tol = adjust_k_tol;
	_adjust_k_factor = adjust_k_factor;
	if(double_nudging>0) _double_nudging = true;

}

bool NEB::step()
{
//	std::cout << "Got to step call\n";  // sn402
	_optimizer->one_iteration();
	_coords.assign(_optimizer->get_x());
	//jdf43	std::cout << "energy: " << _optimizer->get_f() << std::endl;
	return _optimizer->success();
}

double NEB::get_rms()
{
        return _optimizer->get_rms();
}

}
