#include "CalcCoordinat.h"
#include <iostream>

CalcCoordinat::CalcCoordinat()
{
	this->anglePurpose = 0;
	this->base = 0;
	this->coordinatePurposeX = 0;
	this->coordinatePurposeY = 0;
	this->distanceSide = 0;
	this->okoLeft = nullptr;
	this->okoRight = nullptr;
}

CalcCoordinat::CalcCoordinat(Oko* okoRight, Oko* okoLeft)
{
	this->anglePurpose = 0;
	this->base = 0;
	this->coordinatePurposeX = 0;
	this->coordinatePurposeY = 0;
	this->distanceSide = 0;
	this->okoLeft = okoLeft;
	this->okoRight = okoRight;
}

double CalcCoordinat::getCalcBase()
{
	this->base = sqrt((pow(okoRight->getCoordinateXGsKr() - okoLeft->getCoordinateXGsKr(), 2) + (pow(okoRight->getCoordinateYGsKr() - okoLeft->getCoordinateYGsKr(), 2))));
	return this->base;
}

double CalcCoordinat::getDistanceSide()
{
	this->distanceSide = (getCalcBase() / sin(this->getAnglePurpose() / this->radian)) * sin(this->okoLeft->getAngleOko() / this->radian);
	return this->distanceSide;
}

double CalcCoordinat::getAnglePurpose()
{
	this->anglePurpose = 180 - this->okoLeft->getAngleOko() - this->okoRight->getAngleOko();
	return this->anglePurpose;
}

double CalcCoordinat::getNewCoordinateX()
{
	double deltaX = cos(this->okoRight->getAngleOko() / this->radian) * this->getDistanceSide();
	this->coordinatePurposeX = this->okoRight->getCoordinateXGsKr() - deltaX;
	return this->coordinatePurposeX;
}

double CalcCoordinat::getNewCoordinateY()
{

	double deltaY = sin(this->okoRight->getAngleOko() / this->radian) * this->getDistanceSide();
	this->coordinatePurposeY = this->okoRight->getCoordinateYGsKr() + deltaY;
	return this->coordinatePurposeY;
}



