// Convertor.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

//#define _USE_MATH_DEFINES

#include <iostream>

#include <math.h>
#include "Oko.h"
#include "CalcCoordinat.h"
#include "WGS84toPZ90_Gaus.h"


class CK42
{
public:
    CK42() : longt(0), latt(0), h(0) {};

    CK42(double longt, double latt, double h)
    {
        this->longt = longt;
        this->latt = latt;
        this->h = h;
    }
    double longt, latt, h;
};

class GsKr
{
public:
    double x, y, h;
      
    GsKr() : x(0), y(0), h(0) {};
};

class WGS84
{
public:
    WGS84() : longt(0), latt(0), h(0) {};
    WGS84(double longt, double latt, double h)
    {
        this->longt = longt;
        this->latt = latt;
        this->h = h;
    }

    double X_WGS84 = 0;
    double Y_WGS84 = 0;
    double Z_WGS84 = 0;

    double longt, latt, h;
};



class PZ90
{
public:
    PZ90() : longt(0), latt(0), h(0) {};
    PZ90(double longt, double latt, double h)
    {
        this->longt = longt;
        this->latt = latt;
        this->h = h;
    }
    double longt, latt, h;


    double X_PZ90= 0;
    double Y_PZ90 = 0;
    double Z_PZ90 = 0;
};


class Convertor
{
public:
    Convertor(void);
    ~Convertor(void);
    CK42 WGS84ToCK42(WGS84 wgs84);
    GsKr CK42ToGsKr(CK42 ck42);

    GsKr CK42ToGsKr(PZ90 ck42);

    const double Pi = 3.14159265358979;

    //parameters

    const double radian = Pi / 180;

    double deltaX = -0.013; //m
    double deltaY = 0.106; //m 
    double deltaZ = 0.022; //m

    double omegaX = -0.00230; // radian
    double omegaY = 0.00354;// radian
    double omegaZ = 0.00421;// radian

    double scaleM = -0.000000008; // масштабный коэффициент

    double alpha = 1 / 298.257223563;

    double bigHalfOS_A = 6378137;

    double koefficient = 0.999999992;

    double firstMatrix[3][3] = {//+
         {1 * koefficient, -2.041066 * pow(10, -8) * koefficient, -1.716240 * pow(10, -8) * koefficient},
         {2.041066 * pow(10, -8) * koefficient, 1 * koefficient, -1.115071 * pow(10, -8) * koefficient},
        {1.716240 * pow(10, -8) * koefficient, 1.115071 * pow(10, -8) * koefficient, 1 * koefficient}
    };

    double secondMatrix[3][1] = {
        {-0.013},
        {0.106},
        {0.022}
    };

    double secondM[3] = { -0.013 , 0.106, 0.022 };
    
    WGS84 WGS84ToWGS84_XYZ(WGS84 wgs84);//+

    PZ90 WGS84ToPZ90(WGS84 wgs84);

};


const double Pi = 3.14159265358979; // Число Пи
const double ro = 206264.8062; // Число угловых секунд в радиане

// Эллипсоид Красовского
const double aP = 6378245; // Большая полуось
const double alP = 1 / 298.3; // Сжатие
const double e2P = 2 * alP - pow(alP, 2); // Квадрат эксцентриситета

// Эллипсоид WGS84 (GRS80, эти два эллипсоида сходны по большинству параметров)
const double aW = 6378137; // Большая полуось
const double alW = 1 / 298.257223563; // Сжатие
const double e2W = 2 * alW - pow(alW, 2); // Квадрат эксцентриситета

// Вспомогательные значения для преобразования эллипсоидов
const double a = (aP + aW) / 2;
const double e2 = (e2P + e2W) / 2;
const double da = aW - aP;
const double de2 = e2W - e2P;

// Линейные элементы трансформирования, в метрах
const double dx = 23.92;
const double dy = -141.27;
const double dz = -80.9;

// Угловые элементы трансформирования, в секундах
const double wx = 0;
const double wy = 0;
const double wz = 0;

// Дифференциальное различие масштабов
const double ms = 0;

//*********************************************************************************************************************************/
// Функции преобразования между WGS84 и CK42
double dB(double Bd, double Ld, double H)
{
    double B, L, M, N;
    B = Bd * Pi / 180;
    L = Ld * Pi / 180;
    M = a * (1 - e2) / pow((1 - e2 * pow(sin(B), 2)), 1.5);
    N = a * pow((1 - e2 * pow(sin(B), 2)), -0.5);
    return ro / (M + H) * (N / a * e2 * sin(B) * cos(B) * da + ((N * N) / (a * a) + 1) * N * sin(B) * cos(B) * de2 / 2 - (dx * cos(L) + dy * sin(L)) * sin(B) + dz * cos(B)) - wx * sin(L) * (1 + e2 * cos(2 * B)) + wy * cos(L) * (1 + e2 * cos(2 * B)) - ro * ms * e2 * sin(B) * cos(B);
}

double dL(double Bd, double Ld, double H)
{
    double  B, L, N;
    B = Bd * Pi / 180;
    L = Ld * Pi / 180;
    N = a * pow((1 - e2 * pow(sin(B), 2)), -0.5);
    return ro / ((N + H) * cos(B)) * (-dx * sin(L) + dy * cos(L)) + tan(B) * (1 - e2) * (wx * cos(L) + wy * sin(L)) - wz;
}

double WGS84Alt(double Bd, double Ld, double H)
{
    double B, L, N, dH;
    B = Bd * Pi / 180;
    L = Ld * Pi / 180;
    N = a * pow((1 - e2 * pow(sin(B), 2)), -0.5);
    dH = -a / N * da + N * pow(sin(B), 2) * de2 / 2 + (dx * cos(L) + dy * sin(L)) * cos(B) + dz * sin(B) - N * e2 * sin(B) * cos(B) * (wx / ro * sin(L) - wy / ro * cos(L)) + ((a * a) / N + H) * ms;
    return H + dH;
}

double WGS84_SK42_Lat(double Bd, double Ld, double H)
{
    return Bd - dB(Bd, Ld, H) / 3600;
}

double SK42_WGS84_Lat(double Bd, double Ld, double H)
{
    return Bd + dB(Bd, Ld, H) / 3600;
}

double WGS84_SK42_Long(double Bd, double Ld, double H)
{
    return Ld - dL(Bd, Ld, H) / 3600;
}

double SK42_WGS84_Long(double Bd, double Ld, double H)
{
    return Ld + dL(Bd, Ld, H) / 3600;
}
//*********************************************************************************************************************************/

//*********************************************************************************************************************************/
// Функции преобразования в координатную проекцию Гаусса-Крюгера
double SK42BTOX(double B, double L, double H)
{
    double No = (6 + L) / 6;
    double Lo = (L - (3 + 6 * (No - 1))) / 57.29577951;
    double Bo = B * Pi / 180;
    double Xa = pow(Lo, 2) * (109500 - 574700 * pow(sin(Bo), 2) + 863700 * pow(sin(Bo), 4) - 398600 * pow(sin(Bo), 6));
    double Xb = pow(Lo, 2) * (278194 - 830174 * pow(sin(Bo), 2) + 572434 * pow(sin(Bo), 4) - 16010 * pow(sin(Bo), 6) + Xa);
    double Xc = pow(Lo, 2) * (672483.4 - 811219.9 * pow(sin(Bo), 2) + 5420 * pow(sin(Bo), 4) - 10.6 * pow(sin(Bo), 6) + Xb);
    double Xd = pow(Lo, 2) * (1594561.25 + 5336.535 * pow(sin(Bo), 2) + 26.79 * pow(sin(Bo), 4) + 0.149 * pow(sin(Bo), 6) + Xc);
    return 6367558.4968 * Bo - sin(Bo * 2) * (16002.89 + 66.9607 * pow(sin(Bo), 2) + 0.3515 * pow(sin(Bo), 4) - Xd);
}
double SK42LTOY(double B, double L, double H)
{
    double No = (6 + L) / 6;
    double Lo = (L - (3 + 6 * (No - 1))) / 57.29577951;
    double Bo = B * Pi / 180;
    double Ya = pow(Lo, 2) * (79690 - 866190 * pow(sin(Bo), 2) + 1730360 * pow(sin(Bo), 4) - 945460 * pow(sin(Bo), 6));
    double Yb = pow(Lo, 2) * (270806 - 1523417 * pow(sin(Bo), 2) + 1327645 * pow(sin(Bo), 4) - 21701 * pow(sin(Bo), 6) + Ya);
    double Yc = pow(Lo, 2) * (1070204.16 - 2136826.66 * pow(sin(Bo), 2) + 17.98 * pow(sin(Bo), 4) - 11.99 * pow(sin(Bo), 6) + Yb);  
    return (5 + 10 * No) * 100000 + Lo * cos(Bo) * (6378245 + 21346.1415 * pow(sin(Bo), 2) + 107.159 * pow(sin(Bo), 4) + 0.5977 * pow(sin(Bo), 6) + Yc);
}

double SK42XTOB(double X, double Y, double Z)
{ 
    double No = pow(Y * 10, -6);
    double Bi = X / 6367558.4968;
    double Bo = Bi + sin(Bi * 2) * (0.00252588685 - 0.0000149186 * pow(sin(Bi), 2) + 0.00000011904 * pow(sin(Bi), 4));
    double Zo = (Y - (10 * No + 5) * 100000) / (6378245 * cos(Bo));
    double Ba = Zo * Zo * (0.01672 - 0.0063 * pow(sin(Bo), 2) + 0.01188 * pow(sin(Bo), 4) - 0.00328 * pow(sin(Bo), 6));
    double Bb = Zo * Zo * (0.042858 - 0.025318 * pow(sin(Bo), 2) + 0.014346 * pow(sin(Bo), 4) - 0.001264 * pow(sin(Bo), 6) - Ba);
    double Bc = Zo * Zo * (0.10500614 - 0.04559916 * pow(sin(Bo), 2) + 0.00228901 * pow(sin(Bo), 4) - 0.00002987 * pow(sin(Bo), 6) - Bb);
    double dB = Zo * Zo * sin(Bo * 2) * (0.251684631 - 0.003369263 * pow(sin(Bo), 2) + 0.000011276 * pow(sin(Bo), 4) - Bc);
    return (Bo - dB) * 180 / Pi;
}

double SK42YTOL(double X, double Y, double Z)
{
    double No = pow(Y * 10, -6);
    double Bi = X / 6367558.4968;
    double Bo = Bi + sin(Bi * 2) * (0.00252588685 - 0.0000149186 * pow(sin(Bi), 2) + 0.00000011904 * pow(sin(Bi), 4));
    double Zo = (Y - (10 * No + 5) * 100000) / (6378245 * cos(Bo));
    double La = Zo * Zo * (0.0038 + 0.0524 * pow(sin(Bo), 2) + 0.0482 * pow(sin(Bo), 4) + 0.0032 * pow(sin(Bo), 6));
    double Lb = Zo * Zo * (0.01225 + 0.09477 * pow(sin(Bo), 2) + 0.03282 * pow(sin(Bo), 4) - 0.00034 * pow(sin(Bo), 6) - La);
    double Lc = Zo * Zo * (0.0420025 + 0.1487407 * pow(sin(Bo), 2) + 0.005942 * pow(sin(Bo), 4) - 0.000015 * pow(sin(Bo), 6) - Lb);
    double Ld = Zo * Zo * (0.16778975 + 0.16273586 * pow(sin(Bo), 2) - 0.0005249 * pow(sin(Bo), 4) - 0.00000846 * pow(sin(Bo), 6) - Lc);
    double dL = Zo * (1 - 0.0033467108 * pow(sin(Bo), 2) - 0.0000056002 * pow(sin(Bo), 4) - 0.0000000187 * pow(sin(Bo), 6) - Ld);
    return (6 * (No - 0.5) / 57.29577951 + dL) * 180 / Pi;
}
//*********************************************************************************************************************************/
Convertor::Convertor(void)
{
}

Convertor::~Convertor(void)
{
}

CK42 Convertor::WGS84ToCK42(WGS84 wgs84)
{
    CK42 ck42;
    ck42.latt = WGS84_SK42_Lat(wgs84.latt, wgs84.longt, wgs84.h);
    ck42.longt = WGS84_SK42_Long(wgs84.latt, wgs84.longt, wgs84.h);
    ck42.h = wgs84.h;
    return ck42;
}

GsKr Convertor::CK42ToGsKr(CK42 ck42)
{
    GsKr gk;
    gk.x = SK42BTOX(ck42.latt, ck42.longt, ck42.h);
    //gk.x = gk.x - 500000;
    gk.y = SK42LTOY(ck42.latt, ck42.longt, ck42.h);
    //gk.y = gk.y - 500000;
    gk.h = ck42.h;
    return gk;
}

GsKr Convertor::CK42ToGsKr(PZ90 ck42)
{
    GsKr gk;
    gk.x = SK42BTOX(ck42.latt, ck42.longt, ck42.h);
    //gk.x = gk.x - 500000;
    gk.y = SK42LTOY(ck42.latt, ck42.longt, ck42.h);
    //gk.y = gk.y - 500000;
    gk.h = ck42.h;
    return gk;


    //GsKr gk;
    //gk.x = SK42BTOX(ck42.latt, ck42.longt, 0);
    ////gk.x = gk.x - 500000;
    //gk.y = SK42LTOY(ck42.latt, ck42.longt, 0);
    ////gk.y = gk.y - 500000;
    //gk.h = ck42.h;
    //return gk;
}

WGS84 Convertor::WGS84ToWGS84_XYZ(WGS84 wgs84)
{
   

    double e = (2 * alpha) - (alpha * alpha); //+

    double N = bigHalfOS_A / (sqrt((1 - (e) * pow(sin(wgs84.longt * radian), 2))));

    double sinw = sqrt((1 - (e)*pow(sin(wgs84.longt * radian), 2)));

    double tempCosLat = cos(wgs84.longt * radian);
    double tempCoslong = cos(wgs84.latt * radian) ;
    double tempSinLat = sin(wgs84.longt * radian) ;
    double tempSinLongt = sin(wgs84.latt* radian) ;


    wgs84.X_WGS84 = (N + wgs84.h) * tempCosLat * tempCoslong;
    wgs84.Y_WGS84 = (N + wgs84.h) * tempCosLat * tempSinLongt;
    wgs84.Z_WGS84 = ((1 - e) * N + wgs84.h) * tempSinLat;

    return wgs84;
}

PZ90 Convertor::WGS84ToPZ90(WGS84 wgs84)
{
    PZ90 pz90;

    double resultMatrix[3][1] = { {0}, {0}, {0} };

    double mainMatrix[3] = { wgs84.X_WGS84, wgs84.Y_WGS84, wgs84.Z_WGS84 };

    double mainMatrixXYZ[3][1] = {
        {wgs84.X_WGS84},
        {wgs84.Y_WGS84},
        {wgs84.Z_WGS84}
    };

    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 1; col++) {
            // Multiply the row of A by the column of B to get the row, column of product.
            for (int inner = 0; inner < 3; inner++) {
                resultMatrix[row][col] += firstMatrix[row][inner] * mainMatrixXYZ[inner][col];
            }
        }
    }

    //for (int row = 0; row < 3; row++) {
    //    
    //    for (int col = 0; col < 3; col++) {
    //            resultMatrix[row][col] += firstMatrix[row][col] * mainMatrix[col];

    //        //std::cout << product[row][col] << "  ";
    //    }
    //    //std::cout << "\n";
    //}

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 1; j++)
            resultMatrix[i][j] = resultMatrix[i][j] + secondM[j];

    pz90.X_PZ90 = resultMatrix[0][0];
    pz90.Y_PZ90 = resultMatrix[1][0];
    pz90.Z_PZ90 = resultMatrix[2][0];


    return pz90;
}


    //pz90.X_PZ90 = wgs84.X_WGS84 + this->deltaX - this->omegaY * wgs84.Z_WGS84 + this->omegaZ * wgs84.Y_WGS84 + this->scaleM * wgs84.X_WGS84;



    //pz90.Y_PZ90 = wgs84.Y_WGS84 + this->deltaY + this->omegaX * wgs84.Z_WGS84 - this->omegaZ * wgs84.X_WGS84 + this->scaleM * wgs84.Y_WGS84;



    //pz90.Z_PZ90 = wgs84.Z_WGS84 + this->deltaZ - this->omegaX * wgs84.Y_WGS84 + this->omegaY * wgs84.X_WGS84 + this->scaleM * wgs84.Z_WGS84;


int main()
{    
    std::cout << std::fixed;
    std::cout.precision(10);
    WGS84toPZ90_Gaus wgs(50.45519, 30.52973);

    wgs.calcPZ90_X_Y();
    std::cout << "X: " << wgs.getX_Pz90() << "\n";
    std::cout << "Y: " << wgs.getY_Pz90() << "\n";



    return 0;

    //Convertor convertor;

    ////WGS84 wgs84right(37.350719, 44.899007, 0);
    ////WGS84 wgs84left(37.349291, 44.899298, 0);
    ////44.899829, 37.358395
    //WGS84 wgs84right(44.899829, 37.358395, 0);
    //WGS84 wgs84left(44.899833, 37.359399, 0);

    //wgs84right = convertor.WGS84ToWGS84_XYZ(wgs84right);
    //wgs84left = convertor.WGS84ToWGS84_XYZ(wgs84left);

    //PZ90 pz90right = convertor.WGS84ToPZ90(wgs84right);
    //PZ90 pz90left = convertor.WGS84ToPZ90(wgs84left);

    //Oko* okoRight = new Oko(pz90right.X_PZ90, pz90right.Y_PZ90, pz90right.Z_PZ90, 0, 30);
    //Oko* okoLeft = new Oko(pz90left.X_PZ90, pz90left.Y_PZ90, pz90left.Z_PZ90, 0, 30);

    //CalcCoordinat* calcCoordinat = new CalcCoordinat(okoRight, okoLeft);

    //GsKr gs1;

    //gs1 = convertor.CK42ToGsKr(pz90right);


    //std::cout << "Base\t" << calcCoordinat->getCalcBase() << "\n";
    //std::cout << "Side\t" << calcCoordinat->getDistanceSide() << "\n";
    //std::cout << "Angle\t" << calcCoordinat->getAnglePurpose() << "\n";
    //std::cout << "New X\t" << calcCoordinat->getNewCoordinateX() << "\n";
    //std::cout << "New Y\t" << calcCoordinat->getNewCoordinateY() << "\n";

    

    //WGS84 wgs1(37.618, 55.752, 0);
    //WGS84 wgs2(44.89930111, 37.35263611, 0);

    //CK42 ck1;
    //ck1 = convertor.WGS84ToCK42(wgs1);
    //CK42 ck2;
    //ck2 = convertor.WGS84ToCK42(wgs2);
    //CK42 test(44.938004,37.309905, 0);
    //GsKr gskr1;
    //gskr1 = convertor.CK42ToGsKr(ck1);
    //GsKr gskr2;
    //gskr2 = convertor.CK42ToGsKr(ck2);
    // 
    //std::cout << gskr1.x << " " << gskr2.x << "\n";
    //std::cout << gskr1.y << " " << gskr2.y << "\n";
    

 }
