#ifndef MTC_1D_V2_CELLGRIDS_H
#define MTC_1D_V2_CELLGRIDS_H


class CellGrids {
public:
    CellGrids();
    int ccg[100][100]; // cancer cells
    int ccid[100][100]; // cancer cell list indices
    int mpg[100][100]; // all macrophages
    int m0g[100][100]; // naive macropahges
    int m1g[100][100]; //  M1 macrophages
    int m2g[100][100]; // M2 macrophages
    int c8g[100][100]; // T  cells
    int actT[100][100]; //  active T cells
    int allCells[100][100]; // all cells
    int vas[100][100]; // vasculature (only used as entry points for immune cells. excludes the center of the environment)
};


#endif //MTC_1D_V2_CELLGRIDS_H
