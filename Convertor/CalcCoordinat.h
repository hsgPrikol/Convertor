#pragma once
#include "Oko.h"
#include <math.h>

class CalcCoordinat
{
private:
	double base;
	double distanceSide;
	double anglePurpose;
	double coordinatePurposeX;
	double coordinatePurposeY;
	Oko *okoRight = nullptr;
	Oko *okoLeft = nullptr;
	const double radian = 57.29578;

public:
	CalcCoordinat();
	CalcCoordinat(Oko* okoRight, Oko* okoLeft);

	double getCalcBase();
	double getDistanceSide();
	double getAnglePurpose();
	double getNewCoordinateX();
	double getNewCoordinateY();

};

