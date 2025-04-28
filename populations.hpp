#ifndef _POPULATION_H_
#define _POPULATION_H_

#include <time.h> 
#include <iostream>
#include "utility.hpp"
#include <random>

/// @brief Template Populations Class
///
/// The populations class in used to store the populations 

class Populations
{
    private:
    
    std::mt19937 rng; // mersenne twister random generator

    protected:
    double** dataset; // 2d array
    int dimension; //  size of the 2nd array
    int num_items; // total number of items in dataset or size of the 1st array

    public:
    Populations(int num_items, int dimension);
    ~Populations();
    void assign_random(double range_low, double range_high);
    void assign_random2(double range_low, double range_high,int offset);
    void print_population();
    double get_data(int i, int j);
    double* get_data(int i);
    void set_data(double *item, int i);
    void set_data(double data_j, int i, int j);
};

#endif