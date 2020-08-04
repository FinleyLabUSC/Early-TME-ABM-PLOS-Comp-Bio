#include "cancer.h"
#include <random>
#include <numeric>
#include <iostream>
#include <algorithm>

extern std::random_device rd;

Cancer::Cancer(std::array<int, 2> loc, int index, int initial){
    location = loc;
    idx = index;

    state = "alive";
    cellCycle = 30;
    growth = 0;
    maxDivisions = 8;
    nDivisions = 0;
    life_span = 24 * 5;
    age = 0;

    // set initial age and cellCycle to random values for the initial cells
    // to avoid death and division happening at the same time
    std::mt19937 g(rd());
    std::uniform_real_distribution<> growthStart(0.0,cellCycle);
    std::uniform_real_distribution<>  ageStart(0.0,life_span);
    if(initial == 1){
        growth = growthStart(g);
        age = ageStart(g);
    }

    engagement = 6; // hrs of engagement to CTL before death
    dying = 0; // time engaged to CTL
}

void Cancer::proliferate(CellGrids &cg, std::vector<Cancer> &cc_list){
    // spawn a new cell if space is available
    
    int i = location[0];
    int j = location[1];

    // moore neighborhood
    int z = 0;
    std::vector<int> ix;
    std::vector<int> jx;
    for(int il=-1; il<2; il++){
        for(int jl=-1; jl<2; jl++){
            ix.push_back(il);
            jx.push_back(jl);
            z++;
        }
    }
    double probs[z];
    double sum = 0;

    // probability of proliferating into each site
    // based on if the site is free or not
    for(int q=0; q<z; q++){
        probs[q] = (1 - cg.allCells[i+ix[q]][j+jx[q]]);
        sum += probs[q];
    }

    // if sum == 0, there are no free sites
    if(sum == 0){return;}

    // normalize probabilities so that they add to 1
    double norm_probs[z];
    for(int q=0;q<z;q++){norm_probs[q] = probs[q]/(sum);}

    // choose a site at random based on normalized probabilities
    for(int q=1; q<z; q++){norm_probs[q] = norm_probs[q] + norm_probs[q-1];}
    std::uniform_real_distribution<> dis(0.0,1.0);
    double p = dis(rd);
    int choice = 0;
    for(double norm_prob : norm_probs){
        if(p > norm_prob){choice++;}
    }

    // get coordinates of new site
    int ni = i + ix[choice];
    int nj = j + jx[choice];

    // reset cell cycle and increase the number of cell division
    growth = 0;
    nDivisions++;

    // add new cell to the grids
    cg.allCells[ni][nj] = 1 ;
    cg.ccg[ni][nj] = 1;
    // add new cell to the cell list
    cc_list.push_back(Cancer({ni,nj}, cc_list.size(), 0));
    cg.ccid[ni][nj] = cc_list.size()-1;
}

void Cancer::simulate(std::vector<int> attackedCancer, double tstep, CellGrids &cg, std::vector<Cancer> &cc_list){
    if(state == "dead"){return;}

    // if being killed by T cell, reduce the timer until death
    if(dying > 0){
        dying -= tstep;
        if(dying <= 0){
            state = "dead";
            return;
        }
    }

    // die of age
    age = age + tstep;
    if(age >= life_span){
        state = "dead";
        return;
    }

    // see if attacked by T cell
    if(std::find(attackedCancer.begin(), attackedCancer.end(), idx) != attackedCancer.end()){
        dying = engagement; // t cell engagement
        return;
    }

    // see if the cell can proliferate
    growth = growth + tstep;
    if(growth >= cellCycle && nDivisions < maxDivisions){
        proliferate(cg, cc_list);
    }
}

