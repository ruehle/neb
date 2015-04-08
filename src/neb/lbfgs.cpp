#include "lbfgs.h"

namespace cpp_neb {

	LBFGS::LBFGS(
		cpp_neb::BasePotential * potential,
		const cpp_neb::Array<double> x0,
		double tol,
		int M)
		:
		GradientOptimizer(potential, x0, tol),
		M_(M),
		max_f_rise_(1e-4),
		rho_(M_),
		H0_(0.1),
		k_(0)
	{
		// set the precision of the printing
		cout << std::setprecision(12);
		// allocate space for s_ and y_
		for (size_t i = 0; i < M_; ++i){
			s_.push_back(Array<double>(x_.size()));
			y_.push_back(Array<double>(x_.size()));
		}
	}

	/**
	* Do one iteration iteration of the optimization algorithm
	*/
	void LBFGS::one_iteration()
	{
		if (!func_initialized_)
			initialize_func_gradient();

		// make a copy of the position and gradient
		Array<double> xold(x_.copy());
		Array<double> gold(g_.copy());

		// get the stepsize and direction from the LBFGS algorithm
		Array<double> step(x_.size());
		compute_lbfgs_step(step);

		// reduce the stepsize if necessary
		double stepsize = backtracking_linesearch(step);
		// update the LBFGS memeory
		update_memory(xold, gold, x_, g_);
		// print some status information
		if ((iprint_ > 0) && (iter_number_ % iprint_ == 0)){
			cout << "lbgs: " << iter_number_
//				<< " E " << f_       // sn402: This isn't a real energy - see note in neb.cpp
				<< " rms " << rms_
				<< " nfev " << nfev_
				<< " stepsize " << stepsize << std::endl;
		}
		iter_number_ += 1;
	}

	void LBFGS::update_memory(
		Array<double> xold,
		Array<double> gold,
		Array<double> xnew,
		Array<double> gnew)
	{
		// update the lbfgs memory
		// This updates s_, y_, rho_, and H0_, and k_
		int klocal = k_ % M_;
		for (size_t j2 = 0; j2 < x_.size(); ++j2){
			y_[klocal][j2] = gnew[j2] - gold[j2];
			s_[klocal][j2] = xnew[j2] - xold[j2];
		}
		double ys = dot(y_[klocal], s_[klocal]);
		if (ys == 0.) {
			if (verbosity_ > 0) {
				cout << "warning: resetting YS to 1.\n";
			}
			ys = 1.;
		}
		rho_[klocal] = 1. / ys;

		double yy = dot(y_[klocal], y_[klocal]);
		if (yy == 0.) {
			if (verbosity_ > 0) {
				cout << "warning: resetting YY to 1.\n";
			}
			yy = 1.;
		}
		H0_ = ys / yy;
		// increment k
		k_ += 1;
	}

	void LBFGS::compute_lbfgs_step(Array<double> step)
	{
		if (k_ == 0){
			// take a conservative first step
			double gnorm = norm(g_);
			if (gnorm > 1.) gnorm = 1. / gnorm;
			for (size_t j2 = 0; j2 < x_.size(); ++j2){
				step[j2] = -gnorm * H0_ * g_[j2];
			}
			return;
		}

		// copy the gradient into step
		step.assign(g_);
		int jmin = std::max(0, k_ - M_);
		int jmax = k_;
		int i;
		double beta;
		Array<double> alpha(M_);

		// loop backwards through the memory
		for (int j = jmax - 1; j >= jmin; --j){
			i = j % M_;
			alpha[i] = rho_[i] * dot(s_[i], step);
			for (size_t j2 = 0; j2 < step.size(); ++j2){
				step[j2] -= alpha[i] * y_[i][j2];
			}
		}
		// scale the step size by H0
		step *= H0_;

		// loop forwards through the memory
		for (int j = jmin; j < jmax; ++j){
			i = j % M_;
			beta = rho_[i] * dot(y_[i], step);
			for (size_t j2 = 0; j2 < step.size(); ++j2){
				step[j2] += s_[i][j2] * (alpha[i] - beta);
			}
		}

		// invert the step to point downhill
		step *= -1;

	}

	double LBFGS::backtracking_linesearch(Array<double> step)
	{
		Array<double> xnew(x_.size());
		Array<double> gnew(x_.size());
		double fnew;

		// if the step is pointing uphill, invert it
		if (dot(step, g_) > 0.){
			if (verbosity_ > 1) {
				cout << "warning: step direction was uphill.  inverting\n";
			}
			for (size_t j2 = 0; j2 < step.size(); ++j2){
				step[j2] *= -1;
			}
		}

		double factor = 1.;
		double stepsize = norm(step);

		// make sure the step is no larger than maxstep_
		if (factor * stepsize > maxstep_){
			factor = maxstep_ / stepsize;
		}

		int nred;
		int nred_max = 10;
		for (nred = 0; nred < nred_max; ++nred){
			for (size_t j2 = 0; j2 < xnew.size(); ++j2){
				xnew[j2] = x_[j2] + factor * step[j2];
			}
			compute_func_gradient(xnew, fnew, gnew);

			double df = fnew - f_;
			if (df < max_f_rise_){
				break;
			}
			else {
				factor /= 10.;
				if (verbosity_ > 2) {
					cout
						<< "energy increased by " << df
						<< " to " << fnew
						<< " from " << f_
						<< " reducing step size to " << factor * stepsize
						<< " H0 " << H0_ << "\n";
				}
			}
		}

		if (nred >= nred_max){
			// possibly raise an error here
			if (verbosity_ > 0) {
				cout << "warning: the line search backtracked too many times\n";
			}
		}

		x_.assign(xnew);
		g_.assign(gnew);
		f_ = fnew;
		rms_ = norm(gnew) / sqrt(gnew.size());

		return stepsize * factor;
	}
}
