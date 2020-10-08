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
    Cancer(std::array<int, 2> loc, double prolTime, int index, int initial);
    void proliferate(CellGrids &cg, std::vector<Cancer> &cc_list);
    void simulate(std::vector<int> attackedCancer, double tsim, CellGrids &cg, std::vector<Cancer> &cc_list);
    std::string state;
    int idx;
    std::array<int,2> location;

private:
    double div;
    double growth;
    int maxDiv;
    int nDiv;
    double life_span;
    double age;
    double engagement;
    double dying;
};

#endif //MTC_CANCER_H
