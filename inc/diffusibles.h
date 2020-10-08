#include "CellGrids.h"

#ifndef MTC_1D_V2_DIFFUSIBLES_H
#define MTC_1D_V2_DIFFUSIBLES_H

class Diffusibles {
public:
    Diffusibles(double DX, double Dconst, double KIL4Tumor, double KIL4M2, double KIFNG, double KM1F);
    void diffusion(CellGrids cg, double tstep);
    double M1f[100][100];
    double IL4[100][100];
    double IFNG[100][100];
private:
    double dx;
    double D;
    double kIL4Tumor;
    double kIL4M2;
    double kIFNG;
    double kM1f;
};


#endif //MTC_1D_V2_DIFFUSIBLES_H
