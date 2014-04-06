#ifndef PELE_NEB_H
#define PELE_NEB_H

#include <vector>
#include "array.h"
#include "base_potential.h"
#include "optimizer.h"

namespace pele {

class NEB : public BasePotential
{
public:
	NEB(BasePotential *potential)
		: _potential(potential), _k(100.0) { _distance = new CartesianDistance(); }

	virtual ~NEB() {};

	void set_path(std::vector< Array<double> > path);

	double get_energy(Array<double> coords);
	double get_energy_gradient(Array<double> coords, Array<double> grad);


	std::vector< Array<double> > &images() { return _images; }

	void start(void);
	bool step();

	class INEBDistance {
	public:
		// this defines the interface for distance calculations in the NEB
		// the routine is called with the coordinates of the left and right image
		// and should return the distance between the left and right image.
		// In addition, it must write its derivatives with respect to the coordinates
		// of the left and right image in gradient_left and gradient_right, respectively
		virtual double get_distance(Array<double> left, Array<double> right,
							Array<double> gradient_left, Array<double> gradient_right) = 0;
	};

	class CartesianDistance : public INEBDistance{
	public:
		double get_distance(Array<double> left, Array<double> right,
							Array<double> gradient_left, Array<double> gradient_right) {
			double d=0;
			for(size_t i=0; i<left.size(); ++i) {
				d += (right[i] - left[i])*(right[i] - left[i]);
				gradient_left[i] = -(right[i] - left[i]);
				gradient_right[i] = right[i] - left[i];
			}
			return sqrt(d);
		}
	};

protected:
	// TODO: this function involves the copy constructor of std::vector, is this acceptable?
	// but it is not too bad since only array views need to be copied, not the actual data
	std::vector< Array<double> > generate_image_views(Array<double> coords);

	void resize_array_collection(std::vector< Array<double> > &items, size_t size, size_t nelements);
	void adjust_worker_variables();
	void update_distances(std::vector< Array<double> > images, bool update_tangents=true);
	void interpolate_tangent(Array<double> tau, double energy, double energy_left, double energy_right, 
							  Array<double> tau_left, Array<double> tau_right);

	BasePotential *_potential;
	INEBDistance *_distance;
	double _k;
	bool _double_nudging;

	Array<double> _coords;
	std::vector< Array<double> > _images;	
	size_t _nimages;
	size_t _N;

	Array<double> _energies;
	Array<double> _distances;
	std::vector< Array<double> > _tangents;
	std::vector< Array<double> > _tau_left, _tau_right;
	std::vector< Array<double> > _true_gradients;

	Optimizer *_optimizer;
};

}
#endif











