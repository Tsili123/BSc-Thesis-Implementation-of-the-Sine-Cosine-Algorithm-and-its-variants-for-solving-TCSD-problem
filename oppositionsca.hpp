
#ifndef _SINECOSINECLASS_H_
#define _SINECOSINECLASS_H_

#include "population_criterion.hpp"
#include <math.h>  
#include <fstream>  //file
#include <vector>    // process csv file
#include <sstream>   // read csv file

struct SCAParamClass{
    
    struct Limits
    {
        double upper; 
        double lower; 
    };

    Limits limits; // instance of struct #Limits
    int population_size; 
    int dimension; 
    double a; //  r1 calculation
    int iter_max; // maximum number of iterations
};

struct points {
    double fitness;
    int id;
    static bool compare  (const points& a , const points& b) {
            return a.fitness < b.fitness;
    }
};

class SineCosineClass
{
    private:
    
    std::mt19937 rng; // mersenne twister random generator
    SCAParamClass param; // SCA parameters

    double global_optimal_cost; //  global optimal cost
    double* global_optimal_data; // global optimal data

    void stay_in_range(double & val){
        if(val < param.limits.lower)
            val = param.limits.lower;
        else if (val > param.limits.upper)
            val = param.limits.upper;
    }
    
    void global_optimal_set_data(double* data){
        for (int i = 0; i<param.dimension; i++)
            global_optimal_data[i] = data[i];
    }

    public:
    PopulationCriterion *population; // population
    PopulationCriterion *op_population; // opposite population
    PopulationCriterion *best_population; // opposite population
    PopulationCriterion *variable_population;
    PopulationCriterion *opvariable_population;

    points *arrayp;

    void save_result(double optimal_cost, double* optimal_vector, std::string output_file){

        // file pointer 
        std::ifstream instream;
        std::ofstream outstream; 
    
        // opens a created csv file or creates a new csv file in the same folder
        std::string file_name = output_file;
        outstream.open(file_name, std::ios::out | std::ios::app); 
    
        outstream << optimal_cost << "," ;

        for(int i =0; i<param.dimension; i++)
            outstream << optimal_vector[i] <<",";

        outstream<<"\n";

        outstream.close();
    }

    SineCosineClass(SCAParamClass param)
    {
        this->param = param;

        //print parameters
        std::cout << "a: " << param.a << std::endl;
        std::cout << "limits.lower: " << param.limits.lower << std::endl;
        std::cout << "limits.uppper: " << param.limits.upper << std::endl;
        std::cout << "population size: " << param.population_size<< std::endl;
        std::cout << "maximum number of iterations: " << param.iter_max << std::endl;
        std::cout << "dimension: " << param.dimension << std::endl;
        
        //allocate memory for population citerion class
        population = new PopulationCriterion(param.population_size, param.dimension);

        //allocate memory for opposite population citerion class
        op_population = new PopulationCriterion(param.population_size, param.dimension);

        //allocate memory for opposite population citerion class
        best_population = new PopulationCriterion(param.population_size, param.dimension);

        //allocate memory for opposite population citerion class
        variable_population = new PopulationCriterion(1, param.dimension);

        //allocate memory for opposite population citerion class
        opvariable_population = new PopulationCriterion(1, param.dimension);

        //data: allocate memory for optimal data vector
        global_optimal_data = new double[param.dimension];

        arrayp = new points[2*param.population_size];

        //initialize seed of random generator
        rng.seed(time(NULL));
    }

    ~SineCosineClass(){
        if (global_optimal_data)
            delete[] global_optimal_data;

        if (population)
            delete population;

        if (op_population)
            delete op_population;
        
        if (best_population)
            delete best_population;

        if (variable_population)
            delete variable_population;

        if (opvariable_population)
            delete opvariable_population;
        
        if (arrayp)
            delete arrayp;
    }

    double execute(int function_id){

        // initialize distribution in order to fetch random numbers
        std::uniform_real_distribution<double> mstw(0.0,2*M_PI);
        std::uniform_real_distribution<double> mstw1(0,2.0);
        std::uniform_real_distribution<double> mstw2(0.0,1.0);

        population->assign_random(param.limits.lower, param.limits.upper);
        debug1(population->print_population());
        population->evaluate_costval(function_id);

        // create opposite points population
        for (int i = 0; i<param.population_size; i++){
            // update data of every dimension
            for(int j = 0; j<param.dimension; j++){
                double op_var = population->get_data(i,j);
                op_var = param.limits.upper + param.limits.lower - op_var;
                op_population->set_data(op_var,i,j);
            }
        }
        //compute fitness of new population
        op_population->evaluate_costval(function_id);

        //pick the 100 greatest values of fitness of both sets
        for (int i = 0; i< param.population_size; i++){
            arrayp[i].id = i;
            arrayp[i].fitness = population->get_costval(i);
        }

        for (int i = param.population_size; i< 2*param.population_size; i++){
            arrayp[i].id = i;
            arrayp[i].fitness = op_population->get_costval(i-param.population_size);
        }

        //sort in ascending order
        //std::sort(std::begin(arrayp), std::end(arrayp), [](const points &a, const points &b){return a.fitness < b.fitness});
        std::sort(arrayp,arrayp+2*param.population_size,points::compare);

        // for (int i = 0; i< 2*param.population_size; i++){
        //     std::cout << arrayp[i].fitness << " ";
        // }
        // std::cout << std::endl;

        //get the first 100 points
        for (int i = 0; i< param.population_size; i++){
            for(int j = 0; j<param.dimension; j++){
                if(arrayp[i].id >= param.population_size ){
                    arrayp[i].id = arrayp[i].id - param.population_size;
                    best_population->set_data(op_population->get_data(arrayp[i].id,j),i,j);
                }
                else{
                    best_population->set_data(population->get_data(arrayp[i].id,j),i,j);
                }
            }
        }

        best_population->evaluate_costval(function_id);
        // save first global optimal value
        best_population->calculate_min_costval();
        global_optimal_cost = best_population->get_min_costval();
        global_optimal_set_data(best_population->get_min_costval_data());

        debug3(global_optimal_cost);

//EDW0
        double new_xvar2;

        for (int tvar = 0; tvar<param.iter_max; tvar++){

            double r1 = param.a - tvar*((double)param.a/(double)param.iter_max);
            for (int i = 0; i<param.population_size; i++){

                // update data of every dimension
                for(int j =0; j<param.dimension; j++){

                    double xvar = best_population->get_data(i,j);
                    double new_xvar;
                    double r2 = mstw(rng);
                    double r3 = mstw1(rng);
                    double r4 = mstw2(rng);

                    if (r4 < 0.5) //sine
                        new_xvar = xvar + r1*sin(r2)*fabs(r3*global_optimal_data[j]-xvar);
                    else //cosine
                        new_xvar = xvar + r1*cos(r2)*fabs(r3*global_optimal_data[j]-xvar);

                    stay_in_range(new_xvar);
                    best_population->set_data(new_xvar,i,j);
                }

            }

            // update the new optimal cost
            best_population->evaluate_costval(function_id);
            best_population->calculate_min_costval();
            if(best_population->get_min_costval() < global_optimal_cost){
                global_optimal_cost = best_population->get_min_costval();
                global_optimal_set_data(best_population->get_min_costval_data());
            }

            for(int j =0; j<param.dimension; j++){
                new_xvar2 =  param.limits.upper + param.limits.lower - global_optimal_data[j];
                opvariable_population->set_data(new_xvar2,0,j);
            }
            
            opvariable_population->evaluate_costval(function_id);
            if(opvariable_population->get_costval(0) < global_optimal_cost ){
                global_optimal_cost = opvariable_population->get_min_costval();
                global_optimal_set_data(opvariable_population->get_min_costval_data());
            }

            debug2(save_result(global_optimal_cost,global_optimal_data, "sinecosine_f"+std::to_string(function_id)+"_iter.csv"));
        } 

// //EDW */

//         double new_xvar;
//         for (int tvar = 0; tvar<param.iter_max; tvar++){

//             double r1 = param.a - tvar*((double)param.a/(double)param.iter_max);
//             for (int i = 0; i<param.population_size; i++){

//                 // update data of every dimension
//                 for(int j=0; j<param.dimension; j++){

//                     double xvar = best_population->get_data(i,j);

//                     double r2 = mstw(rng);
//                     double r3 = mstw1(rng);
//                     double r4 = mstw2(rng);

//                     if (r4 < 0.5) //sine
//                         new_xvar = xvar + r1*sin(r2)*fabs(r3*global_optimal_data[j]-xvar);
//                     else //cosine
//                         new_xvar = xvar + r1*cos(r2)*fabs(r3*global_optimal_data[j]-xvar);

//                     stay_in_range(new_xvar);
//                     population->set_data(new_xvar, i, j);
//                     //best_population->set_data(new_xvar, i, j);

//                     double op_var ;//= op_population->get_data(i,j);

//                     op_var = param.limits.upper + param.limits.lower - xvar;

//                     // if (r4 < 0.5) //sine
//                     //     new_xvar = op_var + r1*sin(r2)*fabs(r3*global_optimal_data[j]-op_var);
//                     // else //cosine
//                     //     new_xvar = op_var + r1*cos(r2)*fabs(r3*global_optimal_data[j]-op_var);

//                     stay_in_range(op_var);
//                     op_population->set_data(op_var,i,j);

//                     // op_population->set_data(new_xvar, i, j);
//                     //best_population->set_data(new_xvar, i, j);
//                 }
//             }
        
//         //compute fitness of new population
//         op_population->evaluate_costval(function_id);
//         //compute fitness of new population
//         population->evaluate_costval(function_id);

//             //pick the 100 greatest values of fitness of both sets
//         for (int i = 0; i< param.population_size; i++){
//             arrayp[i].id = i;
//             arrayp[i].fitness = population->get_costval(i);
//         }

//         for (int i = param.population_size; i< 2*param.population_size; i++){
//             arrayp[i].id = i;
//             arrayp[i].fitness = op_population->get_costval(i-param.population_size);
//         }

//         //sort in ascending order
//         //std::sort(std::begin(arrayp), std::end(arrayp), [](const points &a, const points &b){return a.fitness < b.fitness});
//         std::sort(arrayp,arrayp+2*param.population_size,points::compare);

//         // for (int i = 0; i< 2*param.population_size; i++){
//         //     std::cout << arrayp[i].fitness << " ";
//         // }
//         // std::cout << std::endl;

//         //get the first 100 points
//         for (int i = 0; i< param.population_size; i++) {
//             for(int j = 0; j<param.dimension; j++){
//                 if(arrayp[i].id >= param.population_size ){
//                     arrayp[i].id = arrayp[i].id - param.population_size;
//                     best_population->set_data(op_population->get_data(arrayp[i].id,j),i,j);
//                 }
//                 else{
//                     best_population->set_data(population->get_data(arrayp[i].id,j),i,j);
//                 }
//             }
//         }

//             // op_population->evaluate_costval(function_id);
//             //     population->evaluate_costval(function_id);

//             //     variable_population->evaluate_costval(function_id);
//             //     opvariable_population->evaluate_costval(function_id);

//             //     if(opvariable_population->get_costval(0) < variable_population->get_costval(0) ){
//             //         best_population->set_data(opvariable_population->get_data(0),i);
//             //     }else{
//             //         best_population->set_data(variable_population->get_data(0),i);
//             //     }

//             // update the new optimal cost
//             best_population->evaluate_costval(function_id);
//             best_population->calculate_min_costval();
//             if(best_population->get_min_costval() < global_optimal_cost){
//                 global_optimal_cost = best_population->get_min_costval();
//                 global_optimal_set_data(best_population->get_min_costval_data());
//             }

//             debug2(save_result(global_optimal_cost,global_optimal_data, "sinecosine_f"+std::to_string(function_id)+"_iter.csv"));
//         }   
        
        debug3(global_optimal_cost);
        return global_optimal_cost;
    }
};

#endif