#ifndef PELE_NEB_H
#define PELE_NEB_H

#include <vector>
#include "array.h"
#include "base_potential.h"
#include "base_distance.h"
#include "optimizer.h"

namespace pele {

class NEB : public BasePotential
{
public:
	NEB(BasePotential *potential, DistanceWrapper *distance)
		: _potential(potential), _k(100.0)  { _distance = distance; }

	virtual ~NEB() {};

	void set_path(std::vector< Array<double> > path);

	double get_energy(Array<double> coords);
	double get_energy_gradient(Array<double> coords, Array<double> grad);


	std::vector< Array<double> > &images() { return _images; }

	void start(void);
	bool step();
	double get_rms();

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
	DistanceWrapper *_distance;
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











