#include <exception>
#include "neb.h"
#include <algorithm>
#include "lbfgs.h"
#include <cmath>

namespace cpp_neb {

void NEB::set_path(std::vector< Array<double> > path)
{
	// check if there are images in the path
	if(path.size() == 0) {
		throw std::runtime_error("cannot initialize neb with empty path");
	}

	int nimages = _nimages = path.size();
	size_t N = _N = path[0].size();

	// check if all images have the same number of coordinates
	for(int i=0; i < path.size(); ++i) {
		if(path[i].size() != N) {
			throw std::runtime_error("number of coordinates in path images differs");
		}
	}

	// reset the old image shortcuts and allocate memory to store new path
	_images.clear();
	_coords = Array<double>(N*nimages);

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
}

void resize_array_vector(vector< Array<double> > &x, size_t nimages, size_t n)
{
	if(x.size() > nimages) {
		x.resize(nimages);
	} else {
	    while(x.size() < nimages) {
            x.push_back(Array<double>(n));
        }
	}
}

void NEB::adjust_worker_variables()
{
	_energies = Array<double>(_nimages);
	_distances = Array<double>(_nimages-1);

	resize_array_vector(_tangents, _nimages-2, _N);
	resize_array_vector(_true_gradients, _nimages, _N);
	resize_array_vector(_tau_left, _nimages, _N);
	resize_array_vector(_tau_right, _nimages, _N);
}

/**
 * Compute the average deviation of the distances between consecutive images in the path.
 * If this deviation is larger than a specified tolerance, k is increased to tighten the band and
 * encourage an even distribution of image. Otherwise, k is decreased to loosen the band and allow exploration
 * of longer paths.
 */

void NEB::adjust_k()  // sn402
{
	// Compute the mean distance between consecutive images
	double average_d = _distances.sum() / _distances.size();

	// Compute the average deviation of the distances (normalised to the average distance)
	double ave_dev = 0;
	for (size_t i=0; i < _distances.size(); i++) {
		double deviation = fabs((_distances[i]-average_d)/average_d);
		ave_dev += deviation;
//		std::cout << deviation << std::endl;
	}
	ave_dev /= _distances.size();

	// If this average deviation is larger than a specified tolerance, the band is too loose
	// so we increase the force constant. If the deviation is smaller, we relax the force constant.
	if (ave_dev > _adjust_k_tol) {
		_k *= _adjust_k_factor;
		if (_verbosity > 0) {
			std::cout << "Increasing NEB force constant to " << _k <<
					" (average deviation is " << ave_dev << ")" << std::endl;
		}
	} else {
		_k /= _adjust_k_factor;
		if (_verbosity > 0) {
			std::cout << "Decreasing NEB force constant to " << _k <<
					" (average deviation is " << ave_dev << ")" << std::endl;
		}
	}

}

std::vector< Array<double> > NEB::generate_image_views(Array<double> coords)
{
	std::vector< Array<double> > images;
	// generate array views for quick access
	for(size_t i=0; i<_nimages; ++i) {
		images.push_back(coords.view(i*_N, (i+1)*_N));
	}
	return images;
}


double NEB::get_energy(Array<double> coords)
{	
	// first wrap coordinates for convenient access
	std::vector< Array<double> > images
		= generate_image_views(coords);
	
	// now loop over all images and sum up the energy
	double energy = 0.0;
	for(size_t i=0; i<_images.size(); ++i) {
		energy += _potential->get_energy(images[i]);
	}
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
//	std::cout << "cpp energy gradient function\n";  // sn402

	// first wrap coordinates for convenient access
	std::vector< Array<double> > images
		= generate_image_views(coords);
	// same for the gradients
	std::vector< Array<double> > image_gradients
		= generate_image_views(grad);

	// calculate the true energy and gradient
	grad.assign(0.);
	double energy = 0.0;
	for(size_t i=0; i<_images.size(); ++i) {
	    // js850> NOTE: image 0 and N-1 never change, so we only need
	    // to compute the energy and gradient once.
		// calculate the true gradient
		_energies[i] = _potential->get_energy_gradient(images[i], _true_gradients[i]);
		energy += _energies[i];
	}

	// update distances and tangents
	update_distances(images, true);

    Array<double> spring(_N);
//    double E_neb = 0;   // sn402: this will hold the spring contribution to the energy

	for(size_t i=1; i<_images.size()-1; ++i) {

	    // define some aliases for convenient access
		Array<double> & tangent = _tangents[i-1];
        Array<double> & image_gradient = image_gradients[i]; // this is what we're computing
        Array<double> & true_gradient = _true_gradients[i];


		// perpendicular part of true gradient
		double project = dot(true_gradient, tangent);
		for(size_t j=0; j<_N; ++j) {
		    image_gradient[j] = true_gradient[j] - project * tangent[j];
		}

		// double nudging
		// we first do the double nudging since in this case image_gradient still contains
		// the perpendicular part of the true gradient which saves an extra storage variable
		if(_double_nudging) {
            spring.assign(0.);
            for(size_t j=0; j<_N; ++j) {
                spring[j] = _k * (_tau_left[i-1][j] + _tau_right[i-1][j]);
            }
//            E_neb += dot(spring, spring)/(2*_k);  // sn402: spring energy  - this doesn't seem to work

            // first project out parallel part since this this is treated separately
            // in the normal nudging
            double project1 = dot(spring, tangent);
            for(size_t j=0; j<_N; ++j) {
                spring[j] -= project1 * tangent[j];
            }

            // project out the part which goes along the direction of the true gradient
            double project2 = dot(spring, image_gradient)
                    / dot(image_gradient, image_gradient);
            for(size_t j=0; j<_N; ++j) {
                spring[j] -= project2 * image_gradient[j];
            }

            // add the spring force
            image_gradient += spring;
		}

		// spring force
		double d = _k * (_distances[i-1] - _distances[i]);

		for(size_t j=0; j<_N; ++j) {
		    image_gradient[j] += d * tangent[j];
		}

//		double temp1 = sqrt(dot(true_gradient, true_gradient));
//		double temp2 = sqrt(dot(tangent, tangent)) * d;
//		std::cout << _k
//				<< " " << _distances[i-1]
//				<< " " << _distances[i]
//				<< " " << _distances[i-1] - _distances[i]
//				<< " " << d
//				<< " " << temp1
//				<< " " << temp2
//				<< std::endl;

	}

	// the gradient of first and last image is zero
	image_gradients[0].assign(0.0);
	image_gradients[_nimages-1].assign(0.0);

	//jdf43 std::cout << "rms: " << get_rms() << std::endl;

	// sn402: The lbfgs algorithm uses the "energy" returned from this function to reject steps that would
	// increase the total energy by more than a specified amount. However here the energy would have to take
	// into account the spring force as well. In practise it becomes simpler to return a nominal "energy"
	// of 0, so that all steps are accepted (in effect, we switch off the linesearch in lbfgs).
	//	return energy + (_k/2)*dot(_distances,_distances);
	return 0;
}

void NEB::update_distances(std::vector< Array<double> > images, bool update_tangent)
{
	// adjust array to store distances + tangents
	_distances = Array<double>(images.size()-1);
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
		_distances[i] = _distance->get_distance(images[i], images[i+1], tau_left, tau_right);

		_tau_left[i-1].assign(tau_save);
		_tau_right[i-1].assign(tau_left);

		// interpolate tangent based on previous step
		interpolate_tangent(_tangents[i-1], _energies[i], _energies[i-1],
							_energies[i+1], tau_save, tau_left);

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
}

void NEB::start()
{
	LBFGS *lbfgs = new LBFGS(this, this->_coords);

	_optimizer = lbfgs;

}

// Currently, lbfgs is the only optimiser implemented, so this is the only option to initialise the NEB.
void NEB::start_with_lbfgs(double rmstol, int setM, double max_f_rise, double H0)
{
//	std::cout << "lbfgs parameters passed in:" << setM << "  " << max_f_rise << "  " << H0;
	LBFGS *lbfgs = new LBFGS(this, this->_coords, rmstol, setM);
	lbfgs->set_max_f_rise(max_f_rise);
	lbfgs->set_H0(H0);
	_optimizer = lbfgs;
}


bool NEB::step()
{
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
