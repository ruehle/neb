#include <fstream>
#include <gtest/gtest.h>
#include "LepsPotential.h"

using namespace std;

struct LepsDataPoint {
	double coords[2];
	double energy;
	double grad[2];
};

vector<LepsDataPoint> LepsReference(string filename)
{
	vector<LepsDataPoint> points;

	ifstream ref;
	ref.open("../tests/AcceptanceTests//data/leps_reference.txt");
	if (!ref.is_open()) throw std::runtime_error("could not open reference");
	
	while (ref) {
		LepsDataPoint p;
		ref >> p.coords[0];
		ref >> p.coords[1];
		ref >> p.energy;
		ref >> p.grad[0];
		ref >> p.grad[1];
		points.push_back(p);
	}
	return points;
}



TEST(LepsPoential, EnergiesOnGrid_AgreeWithReference)
{
	LepsPotential pot;
	for (auto ref : LepsReference("data/leps_reference.txt")) {
		double energy, grad[2];
		energy = pot.getEnergy(ref.coords);
		
		ASSERT_NEAR(ref.energy, energy, 1e-10);		
	}
}

TEST(LepsPoential, DerivativesOnGrid_AgreeWithReference)
{
	LepsPotential pot;
	for (auto ref : LepsReference("data/leps_reference.txt")) {
		double energy, grad[2];
		energy = pot.getEnergyGradient(ref.coords, grad);

		ASSERT_NEAR(ref.grad[0], grad[0], 1e-10);
		ASSERT_NEAR(ref.grad[1], grad[1], 1e-10);		
	}
}

TEST(LepsPotentialTests, MiniminzedPath_AgreesWithReference)
{
	FAIL();
}