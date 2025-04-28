#include <iostream>
#include <string> //stoi
#include "simulator.hpp"


// 7 arguments must be given to main program in order to be able to run sine cosine algorithm with the 6 criterion functions
// arguments: configuration_optimizer_file num_runs function_id dimensions range_min range_max
// config_opt_file = settings for sine cosine algorithm
// num_runs = number of iterations of sine cosine algorithm and its criterion function 
// function_id = id of criterion function
// dimensions = dimensions of each vector (size of each  array element)
// range_min = range lower limit
// range_max = range upper limit
// -------------------

int main(int argc, char** argv){

    int function_id;
    int num_runs;
    int dimensions;
    double range_min;
    double range_max;
    std::string config_opt_file;

    if(argc ==7 ) {
        config_opt_file = argv[1];
        num_runs = std::stoi(argv[2]);
        function_id = std::stoi(argv[3]);
        dimensions = std::stoi(argv[4]);
        range_min = std::stod(argv[5]);
        range_max = std::stod(argv[6]);
    }else {
        printf("arguments: config_opt_file num_runs function_id dimensions range_min range_max\n");
        return 1;
    }

    std::cout << config_opt_file << " " << num_runs << " " << function_id << " " << dimensions << " " << range_min << " " << range_max << " " << std::endl;
    Simulator simulator(num_runs);
    simulator.Optimizer( range_min, range_max , function_id, dimensions , config_opt_file);

    return 0;
}