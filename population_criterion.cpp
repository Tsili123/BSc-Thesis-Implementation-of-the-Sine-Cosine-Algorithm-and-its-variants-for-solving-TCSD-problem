#include "population_criterion.hpp"


PopulationCriterion::PopulationCriterion(int num_items, int dimension): Populations(num_items,dimension){
    // allocate memory for array for costval
    costval = new double[this->num_items];
}

PopulationCriterion::~PopulationCriterion(){
    // free array of costval
    if(costval)
        delete[] costval;
}

void PopulationCriterion::evaluate_costval(int function_type){
    Functions Criterionfunctions;
    Functions::functionctr_pointer funcp = Criterionfunctions.get_function_type(function_type);
    for (int i = 0; i<this->num_items; i++)
        costval[i] = (Criterionfunctions.*funcp)(this->dataset[i],this->dimension); 
}

void PopulationCriterion::print_costval(){
    std::cout<<" ___ costval ____"<<std::endl;
    print_array(costval, this->num_items,'\n');
}

void PopulationCriterion::calculate_min_costval(){
    min_costval = this->costval[0];
    min_costval_i=0;

    for (int i=1; i<this->num_items; i++){
        if(this->costval[i]<min_costval){
            min_costval = this->costval[i];
            min_costval_i = i;
        }
    }
}

double PopulationCriterion::get_min_costval(){
    return min_costval;
}

double PopulationCriterion::get_costval(int i){
    return costval[i];
}

double* PopulationCriterion::get_min_costval_data()
{
    return this->dataset[min_costval_i];
}