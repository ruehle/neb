class LepsPotential
{
public:
	const double a = 0.05;
	const double b = 0.3;
	const double c = 0.05;
	const double alpha = 1.942;
	const double r0 = 0.742;
	const double dAB = 4.746;
	const double dBC = 4.746;
	const double dAC = 3.445;

	double getEnergy(double *r)
	{
		double rAB = r[0];
		double rBC = r[1];
		double rAC = rAB + rBC;

		double JABred = J(dAB, rAB) / (1 + a);
		double JBCred = J(dBC, rBC) / (1 + b);
		double JACred = J(dAC, rAC) / (1 + c);

		return Q(dAB, rAB) / (1 + a) +
			Q(dBC, rBC) / (1 + b) +
			Q(dAC, rAC) / (1 + c) -
			sqrt(JABred*JABred +
			JBCred*JBCred +
			JACred*JACred -
			JABred*JBCred -
			JBCred*JACred -
			JABred*JACred);
	}

	double getEnergyGradient(double *r, double *grad)
	{
		double rAB = r[0];
		double rBC = r[1];
		double rAC = rAB + rBC;

		double JABred = J(dAB, rAB) / (1 + a);
		double JBCred = J(dBC, rBC) / (1 + b);
		double JACred = J(dAC, rAC) / (1 + c);

		double dJABred = dJ(dAB, rAB) / (1 + a);
		double dJBCred = dJ(dBC, rBC) / (1 + b);
		double dJACred = dJ(dAC, rAC) / (1 + c);
		
		double Fx = dQ(dAB, rAB) / (1 + a) +
			dQ(dAC, rAC) / (1 + c) -
			(2 * JABred*dJABred +
			2 * JACred*dJACred -
			dJABred*JBCred -
			JBCred*dJACred -
			dJABred*JACred -
			JABred*dJACred) / 
			(2 * sqrt(JABred*JABred + 
			JBCred*JBCred + 
			JACred*JACred - 
			JABred*JBCred - 
			JBCred*JACred - 
			JABred*JACred));
		
		double Fy = dQ(dBC, rBC) / (1 + b) +
			dQ(dAC, rAC) / (1 + c) -
			(2 * JBCred*dJBCred +
			2 * JACred*dJACred -
			JABred*dJBCred -
			dJBCred*JACred -
			JBCred*dJACred -
			JABred*dJACred) /
			(2 * sqrt(JABred*JABred +
			JBCred*JBCred +
			JACred*JACred -
			JABred*JBCred -
			JBCred*JACred -
			JABred*JACred));
		
		grad[0] = Fx;
		grad[1] = Fy;
		return getEnergy(r);
	}

protected:
	double Q(double d, double r) const
	{
		return d*(3 * exp(-2 * alpha*(r - r0)) / 2 - exp(-alpha*(r - r0))) / 2.0;
	}

	double J(double d, double r) const
	{
		return d*(exp(-2 * alpha*(r - r0)) - 6 * exp(-alpha*(r - r0))) / 4.0;
	}

	double dQ(double d, double r) const 
	{
		return alpha*d*(-3 * exp(-2 * alpha*(r - r0)) + exp(-alpha*(r - r0))) / 2;
	}
	
	double dJ(double d, double r) const 
	{
		return alpha*d*(-2 * exp(-2 * alpha*(r - r0)) + 6 * exp(-alpha*(r - r0))) / 4;
	}

};
