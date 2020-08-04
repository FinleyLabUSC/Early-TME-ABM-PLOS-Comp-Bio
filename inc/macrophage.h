#include <array>
#include <string>
#include <vector>
#include "CellGrids.h"
#include "diffusibles.h"

#ifndef MTC_MACROPHAGE_H
#define MTC_MACROPHAGE_H

class Macrophage{
public:
    Macrophage(std::array<int, 2> loc,
               std::vector<std::vector<std::vector<double>>> nnWeights,
               std::vector<std::vector<std::vector<double>>> nnBiases,
               int initial);
    void migrate(CellGrids &cg, Diffusibles &diff);
    void simulate(double tstep, CellGrids &cg, Diffusibles &diff, double depletion);
    int activate(std::vector<std::vector<double>> outside);
    int neuralNetwork(std::vector<std::vector<double>> outside);
    std::vector<std::vector<double>> sigmoid(std::vector<std::vector<double>> x);
    std::vector<std::vector<double>> matmul(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b);

    std::string state; // cell state
    std::array<int, 2> location; // grid coordinates
    double k15; // phosphorylation of AKT, represents PI3K | NN training: 0 - 1.16

private:
    double life_span;  // hr, max age
    double age;           // hr, current age

    double activationThreshold;

    std::vector<std::vector<std::vector<double>>> weights; // NN weights
    std::vector<std::vector<std::vector<double>>> biases; // NN biases

    double reDiff; // time until redifferentiation

};

#endif //MTC_MACROPHAGE_H
