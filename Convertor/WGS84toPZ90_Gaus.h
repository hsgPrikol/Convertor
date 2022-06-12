#pragma once
#include <math.h>
#include <iostream>

class WGS84toPZ90_Gaus
{
private:
	const double Pi = 3.14159265358979; // ����� ��

	const double compression = 298.257223563; // ������ �������
	const double f = 1 / compression;

	const int a = 6378137; // ������� �������
	const int b = a * (1 - f); // ����� �������

	const double e2 = 1 - ((b * b) / (a * a)); // ������� ��������������
	const double e_ = ((b * b) / (a * a)) - 1; // ������ �������������
	
	double nu;

	double X; // X �� ������-������� ��90.11
	double Y; // Y �� ������-������� ��90.11

	double B_degree; //������ � ��������
	double L_degree; //������� � ��������

	double B_radian; //������ � ��������
	double L_radian; //������� � ��������

	double E;

	int n; // ����� ����

	double L0; // �������� ������� ����

	double l; // �������� ������� ������� �� ������� 
	const double k = (a - b) / (a + b);

	const double C0_WGS84 = 6367449.1458;
	const double C2_WGS84 = 16038.5086;
	const double C4_WGS84 = 16.8326;
	const double C6_WGS84 = 0.0220;

	double S; // ����� ���� ��������

public:
	WGS84toPZ90_Gaus(double _B, double _L);

	void calcPZ90_X_Y();

	double getX_Pz90();
	double getY_Pz90();



};

