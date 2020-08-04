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

    std::string state; // cell state
    int targetCancerCell; // cancer cell to kill
    std::array<int, 2> location; // grid coordinates

private:
    double div;           // hr, division time
    double growth;        // hr, time since last divide
    int maxDiv;          // maximum number of divisions
    int nDiv;             // current number of divisions
    double life_span;  // hr, max age
    double age;           // hr, current age

    double killing; // time left to kill a cancer cell
    int kills;  // max number of tumor cells that can be killed
    double engagementTime; // total time needed to kill a cancer cell
};

#endif //MTC_CD8_H
