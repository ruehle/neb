#ifndef PYGMIN_DISTANCE_H
#define PYGMIN_DISTANCE_H

#include <assert.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include "array.h"
#include "neb_fortran_interface.h"

namespace cpp_neb {
	class BaseDistance {
	public:
		virtual ~BaseDistance() {}
		// this defines the interface for distance calculations in the NEB
		// the routine is called with the coordinates of the left and right image
		// and should return the distance between the left and right image.
		// In addition, it must write its derivatives with respect to the coordinates
		// of the left and right image in gradient_left and gradient_right, respectively
		virtual double get_distance(Array<double> left, Array<double> right,
							Array<double> gradient_left, Array<double> gradient_right) = 0;
	};

	class DistanceWrapper : public BaseDistance
	{
	public:
		DistanceWrapper(distance_callback_t *distanceCallback)
			: _distanceCallback(distanceCallback) {}

		double get_distance(Array<double> left, Array<double> right,
                                                        Array<double> gradient_left, Array<double> gradient_right)
		{
			return (*_distanceCallback)(left.size(), left.data(), right.data(), gradient_left.data(), gradient_right.data());
		}

		distance_callback_t *_distanceCallback;
	};
}

#endif
