#include <fstream>
#include <gtest/gtest.h>
#include "LepsPotential.h"

using namespace std;
using namespace pele;

struct LepsDataPoint {
	Array<double> coords;
	double energy;
	Array<double> grad;
};

vector<LepsDataPoint> LepsReference(string filename)
{
	vector<LepsDataPoint> points;

	ifstream ref;
	ref.open("../tests/AcceptanceTests//data/leps_reference.txt");
	if (!ref.is_open()) throw std::runtime_error("could not open reference");
	
	while (true) {
		LepsDataPoint p;
		p.grad.resize(2);
		p.coords.resize(2);
		ref >> p.coords[0];
		ref >> p.coords[1];
		ref >> p.energy;
		ref >> p.grad[0];
		ref >> p.grad[1];
		if (!ref) break;
		points.push_back(p);
	}
	return points;
}



TEST(LepsPoential, EnergiesOnGrid_AgreeWithReference)
{
	LepsPotential pot;
	for (auto ref : LepsReference("data/leps_reference.txt")) {
		double energy;
		energy = pot.get_energy(ref.coords);
		ASSERT_NEAR(ref.energy, energy, 1e-10);	
	}
}

TEST(LepsPoential, DerivativesOnGrid_AgreeWithReference)
{
	LepsPotential pot;
	for (auto ref : LepsReference("data/leps_reference.txt")) {
		double energy;
		Array<double> grad(2);
		energy = pot.get_energy_gradient(ref.coords, grad);

		ASSERT_NEAR(ref.grad[0], grad[0], 1e-10);
		ASSERT_NEAR(ref.grad[1], grad[1], 1e-10);		
	}
}

TEST(LepsPotentialTests, MiniminzedPath_AgreesWithReference)
{
	FAIL();
}

TEST(LepsPotentialTests, MinimizedPath_GradientPerpendicularToPathVanishes)
{
	FAIL();
}