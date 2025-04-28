#ifndef _SIMULATOR_H_
#define _SIMULATOR_H_

#define MACRO2

#include <random> //mt19937
#include <iostream>

#ifdef MACRO1 
    #include "sinecosine.hpp"
#endif
#ifdef MACRO2 
    #include "modifiedsca.hpp"
#endif
#ifdef MACRO3
    #include "chaoticsca.hpp"
#endif
#ifdef MACRO4
    #include "modifiedexpsca.hpp"
#endif
#ifdef MACRO5
    #include "greedylevysca.hpp"
#endif
#ifdef MACRO6
    #include "cobsca.hpp"
#endif
#ifdef MACRO7
    #include "oppositionsca.hpp"
#endif
#ifdef MACRO8
    #include "oblaquilasca.hpp"
#endif
#ifdef MACRO9
    #include "springproblemsca.hpp"
#endif
#ifdef MACRO10
    #include "springproblemscachaos.hpp"
#endif
#ifdef MACRO11
    #include "springproblemscalevy.hpp"
#endif
#ifdef MACRO12
    #include "springproblemscaaquila.hpp"
#endif

#include <string>
#include <algorithm>    //std::sort
#include <fstream>  //file
#include <vector>    // process csv file
#include <sstream>   // read csv file

// class which runs sine cosine algorithm, and save the statistics in a file

class Simulator{
    private:

    //struct of statistics

    struct statistic{
        double standard_deviation;
        double range;
        double time_milsec;
        double median_value;
        double mean_value;
        double min_range;
        double max_range;
    };

    statistic statistical_analysis; // instance of statistic struct [compute_stats() result]
    int dimension; //  size of the array pointed by each element of vectors  [ parameter for optimizer()]
    int num_runs; // number of runs
    double min_range_rand; // lowest value that the random generator should provide for each sample  [ parameter for optimizer()]
    double max_range_rand; // highgest value that the random generator should provide for each sample  [ parameter for optimizer()]
    double* output; // result of the criterion function by pointer
    int function_type; // type of the criterion function  [ parameter for optimizer()]

    /// @brief compute the result and save it in #statistical_analysis
    /// @param time_ms // running time in milliseconds
    void compute_stats(double time_ms);

    /// @brief save result in csv file
    void save_stats(const std::string& input);
   
    SineCosineClass *sca_algo;

    SCAParamClass sca_params;

    public:
    //constructor of the class Simulator
    // allocate space for vectors (array of array) 
    /// @param num_runs // number of runs 
    Simulator(int num_runs); //constructor 

    //destructor of the class Simulator
    /// free memory
    ~Simulator(); //destructor

    void Optimizer(double min_range, double max_range,  int function_type, int dim , std::string conf_file);
    void printoutput();
    void SCAParserFile(std::string conf_file);  //parameter = address of the configuration file

};

#endif