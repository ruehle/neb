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
	NEB(BasePotential *potential, BaseDistance *distance)
		: _potential(potential),
	      _distance(distance),
		  _k(100.0),
		  _adjust_k_tol(.1),
		  _adjust_k_factor(1.05),
		  _double_nudging(false),
		  _nimages(2), // unknown until set_path is called
		  _N(1), // unknown until set_path is called
	      _verbosity(0),
	      _optimizer(NULL)
    {}

	virtual ~NEB() {};

	void set_path(std::vector< Array<double> > path);
	void set_parameters(int double_nudging, double rmstol, double k_initial, double adjust_k_tol,
			double adjust_k_factor, double maxstep, int maxiter, int iprint, int verbosity);
	void set_lbfgs_parameters(int setM, double max_f_rise, double H0);

	double get_energy(Array<double> coords);
	double get_energy_gradient(Array<double> coords, Array<double> grad);


	std::vector< Array<double> > &images() { return _images; }
	Array<double> energies() { return _energies; }

	void start(void);  // sn402: deprecated
	void start_with_lbfgs(double rmstol, int setM, double max_f_rise, double H0);
	bool step();
	void adjust_k();

	void set_k(double k) { _k = k; }
	void set_double_nudging(bool on_off) { _double_nudging = on_off; }
	void set_verbosity(int verbosity) { _verbosity = verbosity; }

    Array<double> get_true_energies() { return _energies; }
    Array<double> get_distancese() { return _distances; }
	double get_rms();

protected:
	// TODO: this function involves the copy constructor of std::vector, is this acceptable?
	// but it is not too bad since only array views need to be copied, not the actual data.
    // js850> it should be fine.  See http://stackoverflow.com/questions/19454068/c-return-value-optimization
	std::vector< Array<double> > generate_image_views(Array<double> coords);

	void adjust_worker_variables();
	void update_distances(std::vector< Array<double> > images, bool update_tangents=true);
	void interpolate_tangent(Array<double> tau, double energy, double energy_left, double energy_right, 
							  Array<double> tau_left, Array<double> tau_right);

	BasePotential *_potential;
	BaseDistance *_distance;
	double _k;

	double _adjust_k_tol;
	double _adjust_k_factor;

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

	int _verbosity;

	GradientOptimizer *_optimizer;
};

}
#endif











