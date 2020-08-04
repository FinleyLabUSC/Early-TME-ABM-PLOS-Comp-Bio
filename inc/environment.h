#include <vector>
#include <string>
#include "macrophage.h"
#include "CellGrids.h"
#include "CD8.h"
#include "cancer.h"
#include <random>
#include <ctime>
#include <array>
#include "diffusibles.h"

#ifndef MTC_ENVIRONMENT_H
#define MTC_ENVIRONMENT_H

class Environment{
public:
    Environment(double stepSize, std::string folder, int set, double ifng, double il4, int perturb, double perturbTime, double perturbLvl);
    void simulate(double days, double treatmentModulation, double timeOn);

private:
    void run_Cancer(CellGrids &cg);
    void run_Macrophage(CellGrids &cg, Diffusibles &diff, double pi3k, double depletion);
    void recruitMacrophage(CellGrids &cg, double rec);
    void recruitTcells(int simStep, CellGrids &cg);
    void run_CD8(CellGrids &cg);
    void plot(CellGrids &cg, Diffusibles &diff, int s);
    void save();
    void clean(CellGrids &cg);
    void loadNNParameters();
    void loadParameters(int perturb, double time, double level);
    void removeDead(CellGrids &cg);
    void initializeCells(CellGrids &cg);
    void updateTimeCourses(int s, CellGrids &cg, Diffusibles &diff);
    void checkError(int s, CellGrids &cg);
    void printStep();
    void clearGrids(CellGrids &cg, Diffusibles &diff);
    void printFinalInfo(CellGrids &cg);
    double treatmentLevel(double treatmentModulation, double timeOn);

    std::string saveDir;

    // treatment params
    double PI3KTime;
    double PI3Klvl;
    double depTime;
    double deplvl;
    double recTime;
    double reclvl;
    double carBin;

    // cycling  treatment
    double treatmentOn;
    double treatmentOff;

    std::vector<int> treatmentOnOff;

    // cell lists
    std::vector<Cancer> cc_list;
    std::vector<Macrophage> mp_list;
    std::vector<CD8> c8_list;

    // diffusion parameters
    double D;
    double dx;
    double k_Act; // tumor macrophage activation factor secretion rate
    double k_IFNG; // t cell ifng secretion rate
    double k_IL4; // tumor il4 secretion rate

    // time courses
    std::vector<double> times;
    std::vector<int> cc_num;
    std::vector<int> m0_num;
    std::vector<int> m1_num;
    std::vector<int> m2_num;
    std::vector<int> c8_num;
    std::vector<int> actT_num;
    std::vector<double> M1f_max;
    std::vector<double> M1f_avg;
    std::vector<double> IL4_max;
    std::vector<double> IL4_avg;
    std::vector<double> TLS_max;
    std::vector<double> TLS_avg;
    std::vector<double> IFNG_max;
    std::vector<double> IFNG_avg;
    std::vector<double> volume;

    // T cell recruitment params (taken from Gong 2017)
    std::vector<int> cancer_deaths;
    double tdelay;
    double twindow;
    double ka;
    double ki;
    double r1;

    double tstep; // time in hours of a time step
    double macProb; // probability of recruiting a macrophage to a lattice site

    double endTime;

    std::vector<int> attackedCancer; // cancer cell targeted by T cells at each simulation step

    std::vector<std::vector<std::vector<double>>> macrophageWeights;
    std::vector<std::vector<std::vector<double>>> macrophageBiases;
};

#endif //MTC_ENVIRONMENT_H
