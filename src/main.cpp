#include "environment.h"
#include <iostream>
#include <omp.h>
#include <fstream>
#include <random>
#include <thread>
#include <chrono>
#include "macrophage.h"
#include <sstream>

/* Input arguments
 * ---------------
 * 1 - folder to save outputs to
 * 2 - which set of simulations this is (subfolder to arg 1)
 * 3 - number of simulations to run (usually do 100 then average them)
 * 4 - T cell ifng secretion rate
 * 5 - tumor il4 secretion rate
 * 6 - which immunotherapy to use
 * 7 - when to start immunotherapy (usually at 100 days)
 * 8 - therapy strength (0 if no treatment)
 * 9 - length of treatment cycle
 * 10 - days treat is on for within a cycle
 */

int main(int argc, char **argv){
    std::string folder = argv[1];
    std::string set = argv[2];
    int N = std::stoi(argv[3]);

    double ifng = std::stod(argv[4]);
    double il4 = std::stod(argv[5]);

    double simTime = 200;

    int perturb = std::stoi(argv[6]);
    double perturbTime = std::stod(argv[7]);
    double perturbLvl = std::stod(argv[8]);

    double treatmentModulation = std::stod(argv[9]);
    double timeOn = std::stod(argv[10]);

    std::string str = "mkdir -p ./"+folder+"/set_" + set;
    const char *command = str.c_str();
    std::system(command);

    double start = omp_get_wtime();
    //std::cout << "Running...\n";
    for (int i = 0; i < N; i++) {
        std::cout << i << std::endl;
        Environment model(0.5, folder, std::stoi(set), ifng, il4, perturb, perturbTime, perturbLvl);
        model.simulate(simTime, treatmentModulation, timeOn);
    }
    double stop = omp_get_wtime();
    std::cout << "Duration: " << (stop-start)/(60) << std::endl;

    return 0;
}
