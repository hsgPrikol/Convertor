#include "WGS84toPZ90_Gaus.h"

WGS84toPZ90_Gaus::WGS84toPZ90_Gaus(double _B, double _L)
{
	this->B_degree = _B;
	this->L_degree = _L;
}

double WGS84toPZ90_Gaus::getX_Pz90()
{
	return this->X;
}

double WGS84toPZ90_Gaus::getY_Pz90()
{
	return this->Y;
}

void WGS84toPZ90_Gaus::calcPZ90_X_Y()
{
	B_radian = B_degree * Pi / 180;
	L_radian = L_degree * Pi / 180;

	this->nu = e_ * cos(B_radian);

	this->E = (6 + L_degree) / 6;

	this->n = E;
	std::cout << "n: " << this->n << "\n";

	this->L0 = (6 * n - 3) * Pi / 180;
	this->l = L_radian - L0;

	float tmp = 0.00002480158;

	double N = a * pow((1 - e2 * pow((sin(B_radian)), 2)), -0.5);

	double a2 = 0.5 * N * sin(B_radian) * cos(B_radian);

	double a4 = 0.04166 * N * sin(B_radian) * pow(cos(B_radian), 3) * (5 - pow(tan(B_radian), 2) + 9 * pow(nu, 2) + 4 * pow(nu, 4));

	double a6 = 1.0/720.0 * N * sin(B_radian) * pow(cos(B_radian), 5) * (61.0 - 58.0 * pow(tan(B_radian), 2) + pow(tan(B_radian), 4) + 270.0 * pow(nu, 2) - 330.0 * pow(nu, 2) * pow(tan(B_radian), 2));

	double a8 = 1.0/40320.0 * N * sin(B_radian) * pow(cos(B_radian), 7) * (1385.0 - 3111.0 * pow(tan(B_radian), 2) + 543.0 * pow(tan(B_radian), 4) - pow(tan(B_radian), 6));

	double b1 = N * cos(B_radian);

	double b3 = 1.0/6.0 * N * pow(cos(B_radian), 3) * (1.0 - pow(tan(B_radian), 2) + pow(nu, 2));

	double b5 = 1.0 / 120.0 * N * pow(cos(B_radian), 5) * (5.0 - 18.0 * pow(tan(B_radian), 2) + pow(tan(B_radian), 4) + 14.0 * pow(nu, 2) - 58.0 * pow(nu, 2) * pow(tan(B_radian), 2));

	double b7 = 1.0 / 5040.0 * N * pow(cos(B_radian), 7) * (61 - 479 * pow(tan(B_radian), 2) + 179 * pow(tan(B_radian), 4) - pow(tan(B_radian), 6));

	this->S = C0_WGS84 * B_radian - C2_WGS84 * sin(2 * B_radian) + C4_WGS84 * sin(4 * B_radian) - C6_WGS84 * sin(6 * B_radian);

	this->X = S + a2 * pow(l, 2) + a4 * pow(l, 4) + a6 * pow(l, 6) + a8 * pow(l, 8);

	double yTmp = b1 * l + b3 * pow(l, 3) + b5 * pow(l, 5) + b7 * pow(l, 7);

	this->Y = yTmp + (5 + 10 * n) * pow(10, 5);
}
