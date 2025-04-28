#ifndef _POPULATIONCRITERION_H_
#define _POPULATIONCRITERION_H_

#include "populations.hpp" // class
#include "functions.hpp"
#include <algorithm> //sort


/// @brief  PopulationCriterion Class
// The PopulationCriterion class inherit all features of PopulationCriterion base class
// this class is useful for different populationcriterion based optimization algorithms

class PopulationCriterion: public Populations 
{
    private:

    double* costval;  // array that contains the cost of every item
    double min_costval; // minimum cost of the array
    int min_costval_i; // index of the mininum cost

    public:
    /// @brief function to print array of double type
    /// @param separatorch separator character after  printing 

    void print_array(double* array, int num, char separatorch){
        for (int i = 0; i < num; i++)   
            std::cout<< array[i] << separatorch;
        std::cout <<std::endl;
    }

    ~PopulationCriterion();
    PopulationCriterion(int num_items, int dimension);
    void evaluate_costval(int function_id);
    void print_costval();
    void calculate_min_costval();
    double get_min_costval();
    double get_costval(int i);
    double* get_min_costval_data();
     
};

#endif