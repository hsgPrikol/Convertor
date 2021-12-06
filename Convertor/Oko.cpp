#include "Oko.h"

Oko::Oko()
{
	this->angleNordPolus = 0;
	this->angleOko = 0;
	this->coordinateXGsKr = 0;
	this->coordinateYGsKr = 0;
	this->height = 0;
}

Oko::Oko(double coordinateXGsKr, double coordinateYGsKr, double height, double angleNordPolus, double angleOko)
{
	this->angleNordPolus = angleNordPolus;
	this->angleOko = angleOko;
	this->coordinateXGsKr = coordinateXGsKr;
	this->coordinateYGsKr = coordinateYGsKr;
	this->height = height;
}

void Oko::setCoordinateXGsKr(double lattitude)
{
	this->coordinateXGsKr = coordinateXGsKr;
}

void Oko::setCoordinateYGsKr(double longtitude)
{
	this->coordinateYGsKr = coordinateYGsKr;
}

void Oko::setAngleNordPolus(double angleNordPolus)
{
	this->angleNordPolus = angleNordPolus;
}

void Oko::setAngleOko(double angleOko)
{
	this->angleOko = angleOko;
}

void Oko::setHeight(double height)
{
	this->height = height;
}

double Oko::getCoordinateXGsKr()
{
	return this->coordinateXGsKr;
}

double Oko::getCoordinateYGsKr()
{
	return this->coordinateYGsKr;
}

double Oko::getAngleNordPolus()
{
	return this->angleNordPolus;
}

double Oko::getAngleOko()
{
	return this->angleOko;
}

double Oko::getHeight()
{
	return this->height;
}
