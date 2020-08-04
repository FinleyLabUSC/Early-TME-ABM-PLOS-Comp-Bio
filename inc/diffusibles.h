#include "CellGrids.h"

#ifndef MTC_1D_V2_DIFFUSIBLES_H
#define MTC_1D_V2_DIFFUSIBLES_H

class Diffusibles {
public:
    Diffusibles(double DX, double Dconst, double KIL4Tumor, double KIL4M2, double KIFNG, double KM1F);
    void diffusion(CellGrids cg, double tstep);

    double activationFactor[100][100];
    double IL4[100][100];
    double IFNG[100][100];

private:
    double dx; // spatial step
    double D; // diffusion constant
    double kIL4Tumor; // tumor IL-4 secretion rate
    double kIL4M2; // M2 IL-4 secretion rate
    double kIFNG; // T cell IFN-G secretion rate
    double kAct; // tumor macrophage activation factor secretion rate
};


#endif //MTC_1D_V2_DIFFUSIBLES_H
