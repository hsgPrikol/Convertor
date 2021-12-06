#pragma once
class Oko
{
private:
	double coordinateXGsKr;
	double coordinateYGsKr;
	double height;
	double angleNordPolus;
	double angleOko;

public:
	Oko();
	Oko(double coordinateXGsKr = 0, double coordinateYGsKr = 0, double height = 0, double angleNordPolus = 0, double angleOko = 0);

	void setCoordinateXGsKr(double coordinateXGsKr);
	void setCoordinateYGsKr(double coordinateYGsKr);
	void setAngleNordPolus(double angleNordPolus);
	void setAngleOko(double angleOko);
	void setHeight(double height);

	double getCoordinateXGsKr();
	double getCoordinateYGsKr();
	double getAngleNordPolus();
	double getAngleOko();
	double getHeight();
};

