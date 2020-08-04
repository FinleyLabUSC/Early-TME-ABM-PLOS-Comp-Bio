#include <array>
#include <string>
#include <random>
#include <numeric>
#include <tuple>
#include "CellGrids.h"

#ifndef MTC_CANCER_H
#define MTC_CANCER_H

class Cancer{
public:
    Cancer(std::array<int, 2> loc, int index, int initial);
    void proliferate(CellGrids &cg, std::vector<Cancer> &cc_list);
    void simulate(std::vector<int> attackedCancer, double tsim, CellGrids &cg, std::vector<Cancer> &cc_list);

    std::string state; // cell state
    int idx; // index in cell list
    std::array<int,2> location; // grid coordinates

private:
    double cellCycle; // length of cell cycle in hourse
    double growth; // time spent in current cell cycle
    int maxDivisions; // maximum number of divisions
    int nDivisions; // current number of divisions
    double life_span; // lifespan hours
    double age; // current age in hours
    double engagement; // total hours to be engaged  with a T cell before death
    double dying; // hours of T cell engagement remaining before death

    double kIL4c;
    double kIFNGc;
};

#endif //MTC_CANCER_H
