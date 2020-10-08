#include <array>
#include <string>
#include <vector>
#include "CellGrids.h"

#ifndef MTC_CD8_H
#define MTC_CD8_H

class CD8{
public:
    CD8(std::array<int, 2> loc);
    void migrate(CellGrids &cg);
    void proliferate(CellGrids &cg, std::vector<CD8> &c8_list);
    int kill(CellGrids &cg);
    void simulate(double tstep, CellGrids &cg, std::vector<CD8> &c8_list);
    std::string state;
    int cc;
    std::array<int, 2> location;

private:
    double div;           // hr, division time
    double growth;        // hr, time since last divide
    int maxDiv;          // maximum number of divisions
    int nDiv;             // current number of divisions
    double life_span;  // hr, max age
    double age;           // hr, current age

    double killing;
    int kills;  // max number of tumor cells that can be killed
    double engagementTime;
};

#endif //MTC_CD8_H
