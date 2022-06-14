#pragma once
#include <math.h>
#include <iostream>

class WGS84toPZ90_Gaus
{
private:
	const double Pi = 3.14159265358979; // Число Пи

	const double compression = 298.257223563; // Сжатие эллипса
	const double f = 1 / compression;

	const double a = 6378137; // Большая полуось
	const double b = a * (1 - f); // Малая полуось

	const double e2 = 1 - ((b * b) / (a * a)); // Квадрат эксонтрицитета
	const double e_ = ((a * a) / (b * b)) - 1; // Второй эксонтрицитет
	
	double nu;

	double X; // X по гауссу-крюгеру пз90.11
	double Y; // Y по гауссу-крюгеру пз90.11

	double B_degree; //широта в градусах
	double L_degree; //долгота в градусах

	double B_radian; //широта в радианах
	double L_radian; //долгота в радианах

	double E;

	int n; // Номер зоны

	double L0; // Середина текущей зоны

	double l; // Разность текушей долготы от средней 
	const double k = (a - b) / (a + b);

	const double C0_WGS84 = 6367449.1458;
	const double C2_WGS84 = 16038.5086;
	const double C4_WGS84 = 16.8326;
	const double C6_WGS84 = 0.0220;

	double S; // Длина дуги мереданы

public:
	WGS84toPZ90_Gaus(double _B, double _L);

	void calcPZ90_X_Y();

	double getX_Pz90();
	double getY_Pz90();



};

