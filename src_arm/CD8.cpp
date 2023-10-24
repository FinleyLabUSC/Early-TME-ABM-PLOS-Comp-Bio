#include "CD8.h"
#include "cancer.h"
#include <string>
#include <random>
#include <math.h>
#include <iostream>

extern std::random_device rd;

CD8::CD8(std::array<int, 2> loc){
    location = loc;
    state = "inactive";

    div = 8;           // hr, division time
    growth = 0;        // hr, time since last divide
    maxDiv = 5;          // maximum number of divisions
    nDiv = 0;             // current number of divisions
    life_span = 41;  // hr, max age
    age = 0;           // hr, current age

    cc = -1;          // idx of cancer cell to kill

    killing = 0;
    kills = 5; // number of tumor cells it  can kill
    engagementTime = 6; // hrs
}

void CD8::migrate(CellGrids &cg){
    // migrate towards the center of the environment
    // this is them essentially chemotaxing towards the tumor

    int i = location[0];
    int j = location[1];

    int ix[] = {-1,1,0,0};
    int jx[] = {0,0,-1,1};
    double probs[4];
    double sum = 0;
    double distance = 0;

    for(int q=0; q<4; q++){
        distance = std::pow(2,std::sqrt(std::pow((i+ix[q])-50,2) + std::pow((j+jx[q])-50,2)));
        probs[q] = (1 - cg.allCells[i+ix[q]][j+jx[q]])/distance;
        sum += probs[q];
    }

    if(sum == 0){
        return;
    }

    double norm_probs[4];
    for(int q=0; q<4; q++){
        norm_probs[q] = probs[q]/sum;
    }
    for(int q=1; q<4; q++){
        norm_probs[q] = norm_probs[q] + norm_probs[q-1];
    }

    std::uniform_real_distribution<> dis(0.0,1.0);
    double p = dis(rd);

    int choice = 0;

    for(double norm_prob : norm_probs){
        if(p > norm_prob){
            choice++;
        }
    }

    int ni = i + ix[choice];
    int nj = j + jx[choice];

    // update stored location
    location[0] = ni;
    location[1] = nj;

    cg.allCells[i][j] = 0;
    cg.c8g[i][j] = 0;
    cg.actT[i][j] = 0;

    cg.allCells[ni][nj] = 1;
    cg.c8g[ni][nj] = 1;
    if(state == "active"){
        cg.actT[ni][nj] = 1;
    }
}

void CD8::proliferate(CellGrids &cg, std::vector<CD8> &c8_list){
    // spawn a new cell if there is space

    int i = location[0];
    int j = location[1];

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

    for(int q=0; q<z; q++){
        probs[q] = (1 - cg.allCells[i+ix[q]][j+jx[q]]);
        sum += probs[q];
    }

    if(sum == 0){
        return;
    }

    double norm_probs[z];
    for(int q=0; q<z; q++){
        norm_probs[q] = probs[q]/sum;
    }
    for(int q=1; q<z; q++){
        norm_probs[q] += norm_probs[q-1];
    }

    std::uniform_real_distribution<> dis(0.0,1.0);
    double p = dis(rd);
    int choice = 0;

    for(double prob : norm_probs){
        if(p > prob){
            choice++;
        }
    }

    int ni = i + ix[choice];
    int nj = j + jx[choice];

    growth = 0;
    nDiv++;

    cg.allCells[ni][nj] = 1;
    cg.c8g[ni][nj] = 1;
    c8_list.push_back(CD8({ni, nj}));
}

int CD8::kill(CellGrids &cg){
    // kill a cancer cell if one is adjacent

    int i = location[0];
    int j = location[1];

    int ix[] = {-1,1,0,0};
    int jx[] = {0,0,-1,1};
    double probs[4];
    double sum = 0;

    for(int q=0; q<4; q++){
        probs[q] = cg.ccg[i+ix[q]][j+jx[q]];
        sum += probs[q];
    }

    if(sum == 0){
        return -1;
    }

    double norm_probs[4];
    for(int q=0; q<4; q++){
        norm_probs[q] = probs[q]/sum;
    }
    for(int q=1; q<4; q++){
        norm_probs[q] += norm_probs[q-1];
    }

    std::uniform_real_distribution<> dis(0.0,1.0);
    double p = dis(rd);

    int choice = 0;

    for(double norm_prob : norm_probs){
        if(p > norm_prob){
            choice++;
        }
    }

    return cg.ccid[i+ix[choice]][j+jx[choice]];
}

void CD8::simulate(double tstep, CellGrids &cg, std::vector<CD8> &c8_list){
    cc = -1; // reset cancer cell to kill

    if(state == "dead"){return;}

    if(killing > 0){
        killing -= tstep;
        age += tstep;
        return;
    }

    int i = location[0];
    int j = location[1];

    age = age + tstep;
    if(age >= life_span){
        state = "dead";
        return;
    }

    if(kills == 0){
        state = "exhausted";
        cg.actT[i][j] = 0;
    }

    // proliferate
    growth = growth + tstep;
    if(growth>=div && nDiv<maxDiv  && state=="active"){
        proliferate(cg, c8_list);
        return;
    }

    // if active, try to kill
    if(state == "active" and kills > 0){
        cc = kill(cg);
        if(cc != -1){
            kills -= 1;
            killing = engagementTime;
            return;
        }
    }

    // try to migrate
    migrate(cg);

    i = location[0];
    j = location[1];

    // Try to activate
    int tumorPresence = 0;
    if((cg.ccg[i+1][j] + cg.ccg[i-1][j] + cg.ccg[i][j+1] + cg.ccg[i][j-1]) > 0){tumorPresence = 1;}
    int m1Presence = 0;
    if((cg.m1g[i+1][j] + cg.m1g[i-1][j] + cg.m1g[i][j+1] + cg.m1g[i][j-1]) > 0){m1Presence = 1;}

    int activatorPresence = std::min(1,tumorPresence+m1Presence);

    int numM1 = cg.m1g[i+1][j+1] + cg.m1g[i][j+1] + cg.m1g[i-1][j+1]
                    +cg.m1g[i+1][j] + cg.m1g[i-1][j]
                    +cg.m1g[i+1][j-1] + cg.m1g[i][j-1] + cg.m1g[i-1][j-1] + 1;

    int numM2 = cg.m2g[i+1][j+1] + cg.m2g[i][j+1] + cg.m2g[i-1][j+1]
                     +cg.m2g[i+1][j] + cg.m2g[i-1][j]
                     +cg.m2g[i+1][j-1] + cg.m2g[i][j-1] + cg.m2g[i-1][j-1];

    double actProb = 1/(1 + exp(-3*((activatorPresence - numM2/numM1) - 0.5)));

    std::uniform_real_distribution<double> dis(0.0,1.0);
    if(activatorPresence == 1 && dis(rd) <= actProb && state != "active" && state != "exhausted"){
        state = "active";
        cg.actT[i][j] = 1;
    }
}