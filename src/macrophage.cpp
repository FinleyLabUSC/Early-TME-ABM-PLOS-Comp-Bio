#include "macrophage.h"
#include <string>
#include <random>
#include <cmath>
#include <algorithm>

extern std::random_device rd;

Macrophage::Macrophage(std::array<int, 2> loc, std::vector<std::vector<std::vector<double>>> nnWeights, std::vector<std::vector<std::vector<double>>> nnBiases, int initial){
    location = loc;
    state = "M0";

    reDiff = 0;

    std::uniform_real_distribution<double> life(0.0,1.0);
    life_span = 24*30;
    // randomize age of initial macrophages so  they don't all die at once
    if(initial == 1) {
        age = life_span*life(rd);           // hr, current age
    } else {age=0;}

    // from Wells 2015
    activationThreshold = 8e-6;

    // neural network parameters
    weights = nnWeights;
    biases = nnBiases;

    // NN inputs
    k15 = 1.16; 
}

void Macrophage::migrate(CellGrids &cg, Diffusibles &diff){
    // chemotax towards IL4
    // similar process as in CD8.cpp

    int i = location[0];
    int j = location[1];

    int ix[] = {-1,1,0,0};
    int jx[] = {0,0,-1,1};
    double probs[4];
    double sum = 0;

    for(int q=0; q<4; q++){
        probs[q] = (1 - cg.allCells[i+ix[q]][j+jx[q]])
                   *diff.IL4[i+ix[q]][j+jx[q]];
        sum = sum + probs[q];
    }

    if(sum == 0){
        return;
    }

    int maxIdx = 0;
    double maxProb = 0;

    for(int q=0; q<4; q++){
        if(probs[q] > maxProb){
            maxProb = probs[q];
            maxIdx = q;
        }
    }

    // weight the highest site a little bit
    probs[maxIdx] = 3*probs[maxIdx];

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
    cg.mpg[i][j] = 0;
    cg.m0g[i][j] = 0;
    cg.m1g[i][j] = 0;
    cg.m2g[i][j] = 0;

    cg.allCells[ni][nj] = 1;
    cg.mpg[ni][nj] = 1;
    if(state == "M0"){
        cg.m0g[ni][nj] = 1;
    }
    if(state == "M1"){
        cg.m1g[ni][nj] = 1;
    }
    if(state == "M2"){
        cg.m2g[ni][nj] = 1;
    }
}

int Macrophage::activate(std::vector<std::vector<double> > outside) {
    // use the neural network to determine phenotype upon activation

    // scaling parameters from Monte Carlo simulations and NN training
    std::vector<double> maxes = {5e4, 2e4, 1.16};
    std::vector<double> mins = {0.0, 0.0, 0.0};

    // scale the inputs
    for(int i=0; i<outside.size(); ++i){
        for(int j=0; j<outside[0].size(); ++j){
            outside[i][j] = (outside[i][j] - mins[j])/(maxes[j] - mins[j]);
        }
    }

    // predict differentiation
    int diffState = neuralNetwork(outside);

    return diffState;
}

int Macrophage::neuralNetwork(std::vector<std::vector<double> > outside) {
    // ANN via matrix multiplication
    // one hidden layer with sigmoid activation
    // output layer is a single neuron with sigmoid activation

    std::vector<std::vector<double>> currentLayer = outside;

    for(int l=0; l<weights.size(); ++l){
        std::vector<std::vector<double>> neuronLayer = matmul(currentLayer, weights[l]);
        for(int i=0; i<neuronLayer.size(); ++i){
            for(int j=0; j<neuronLayer[0].size(); ++j){
                neuronLayer[i][j] += biases[l][j][0];
            }
        }
        currentLayer = sigmoid(neuronLayer);
    }

    double output = currentLayer[0][0];
    if(output >= 0.5){return 1;}
    else{return 0;}
}

std::vector<std::vector<double>> Macrophage::matmul(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b){
    // matrix multiplication
    // inputs (l x M) multiplied by weights (m x n)

    std::vector<std::vector<double>> output;

    int l = a.size();
    int m = a[0].size();
    int n = b[0].size();

    for(int i=0; i<l; ++i){
        std::vector<double> blah;
        for(int j=0; j<n; ++j){
            blah.push_back(0);
        }
        output.push_back(blah);
    }

    for(int i=0; i<l; ++i){
        for(int j=0; j<n; ++j){
            for(int k=0; k<m; ++k){
                output[i][j] += a[i][k]*b[k][j];
            }
        }
    }

    return output;
}

std::vector<std::vector<double>> Macrophage::sigmoid(std::vector<std::vector<double>> x){
    // sigmoid function used as activation function for hidden layer
    // sig(x) = 1/(1 + exp(-x))

    std::vector<std::vector<double>> y;

    for(int i=0; i<x.size(); ++i){
        std::vector<double> someVector;
        for(int j=0; j<x[0].size(); j++){
            someVector.push_back(0);
        }
        y.push_back(someVector);
    }

    for(int i=0; i<x.size(); ++i){
        for(int j=0; j<x[0].size(); ++j){
            y[i][j] = 1/(1 + exp(-x[i][j]));
        }
    }

    return y;
}

void Macrophage::simulate(double tstep, CellGrids &cg, Diffusibles &diff, double depletion){
    if(state == "dead"){
        return;
    }

    // if macrophage depletion immunotherapy is present, see if macrophage dies
    std::uniform_real_distribution<double> dis(0.0,1.0);
    if(dis(rd) < depletion){state = "dead"; return;}

    int i = location[0];
    int j = location[1];

    double actSig = diff.activationFactor[i][j];

    // die of age
    age = age + tstep;
    if(age >= life_span){
        state = "dead";
        return;
    }

    reDiff -= tstep;
    if(reDiff < 0){
        reDiff = 0;
    }

    // activation from Wells 2015 (threshold of activation signal)
    // differentiation based on neural network
    // neural network trained on ODE model from Zhao 2019
    if((actSig >= activationThreshold && state == "M0") || (actSig>=activationThreshold && state!="M0" && reDiff==0)){
        std::vector<std::vector<double>> nnInputs = {{diff.IFNG[i][j],
                                                      diff.IL4[i][j],
                                                      k15}};

        int differentiate = activate(nnInputs);
        reDiff = 24;
        if(differentiate == 1){
            state = "M1";
            cg.m0g[i][j] = 0;
            cg.m1g[i][j] = 1;
            return;
        } else {
            state = "M2";
            cg.m0g[i][j] = 0;
            cg.m2g[i][j] = 1;
            return;
        }
    }

    migrate(cg, diff);
}

/* REFERENCES
 * ----------
 * Zhao et al. "A mechanistic integrative computational model of macrophage polarization: implications in human pathophysiology." PLOS Comp Bio 2019
 */
