
#ifndef _SINECOSINECLASS_H_
#define _SINECOSINECLASS_H_

#include "population_criterion.hpp"
#include <math.h>  
#include <fstream>  //file
#include <vector>    // process csv file
#include <sstream>   // read csv file
#include <algorithm>   

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
    PopulationCriterion *population2; // population2
    double* x_mut;
    int* r_mut;

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

        //allocate memory for population citerion class
        population2 = new PopulationCriterion(param.population_size, param.dimension);

        //data: allocate memory for optimal data vector
        global_optimal_data = new double[param.dimension];

        //data: allocate memory for optimal data vector
        r_mut = new int[param.population_size];

         //data: allocate memory for optimal data vector
        x_mut = new double[param.population_size];

        //initialize seed of random generator
        rng.seed(time(NULL));
    }

    ~SineCosineClass(){
        if (global_optimal_data)
            delete[] global_optimal_data;

        if (population)
            delete population;

        if (population2)
            delete population2;

        if (r_mut)
            delete r_mut;

        if (x_mut)
            delete x_mut;
    }

    double execute(int function_id){

        // initialize distribution in order to fetch random numbers
        std::uniform_real_distribution<double> mstw(0.0,2*M_PI);
        std::uniform_real_distribution<double> mstw1(0,2.0);
        std::uniform_real_distribution<double> mstw2(0.0,1.0);
        std::uniform_int_distribution<int> msca(0,param.dimension-1);
        std::uniform_real_distribution<double> unipop2(param.limits.lower,param.limits.upper);

        population->assign_random(param.limits.lower, param.limits.upper);
        debug1(population->print_population());
        population->evaluate_costval(function_id);

        // save first global optimal value
        population->calculate_min_costval();
        global_optimal_cost = population->get_min_costval();
        global_optimal_set_data(population->get_min_costval_data());

        debug3(global_optimal_cost);

        for (int tvar = 0; tvar<param.iter_max; tvar++){

            double r1 = param.a - tvar*((double)param.a/(double)param.iter_max);
            for (int i = 0; i<param.population_size; i++){

                // update data of every dimension
                for(int j =0; j<param.dimension; j++){

                    double xvar = population->get_data(i,j);
                    double new_xvar;
                    double r2 = mstw(rng);
                    double r3 = mstw1(rng);
                    double r4 = mstw2(rng);

                    if (r4 < 0.5) //sine
                        new_xvar = xvar + r1*sin(r2)*fabs(r3*global_optimal_data[j]-xvar);
                    else //cosine
                        new_xvar = xvar + r1*cos(r2)*fabs(r3*global_optimal_data[j]-xvar);

                    stay_in_range(new_xvar);
                    population->set_data(new_xvar, i, j);
                }
            }

            population->evaluate_costval(function_id);

            for(int i=0; i<param.population_size; i++){
                r_mut[i] = msca(rng);
            }

            for(int i=0; i<param.population_size; i++){
                x_mut[i] = unipop2(rng);
                stay_in_range(x_mut[i]);
            }

            //copy old population to new
            for(int i=0; i<param.population_size; i++){
                for(int j=0; j<param.dimension; j++){
                    population2->set_data( population->get_data(i,j) ,i ,j );
                }
            }

            //update new population
            for(int i=0; i<param.population_size; i++){
                    population2->set_data(x_mut[i],i, r_mut[i]);
            }

            // update the new optimal cost
            population2->evaluate_costval(function_id);
            for(int i=0; i<param.population_size; i++){
                // if fitness x_new1 < fitness x_new
                // update element of old population with the one from new population
                if(population2->get_costval(i) < population->get_costval(i) ){
                    population->set_data(population2->get_data(i,r_mut[i]), i, r_mut[i]);
                }
            }

            population->evaluate_costval(function_id);
            //calculate new  min cost
            population->calculate_min_costval();
            // if necessary , update the global min
            if(population->get_min_costval() < global_optimal_cost){
                global_optimal_cost = population->get_min_costval();
                global_optimal_set_data(population->get_min_costval_data());
            }
            
            debug2(save_result(global_optimal_cost,global_optimal_data, "sinecosine_f"+std::to_string(function_id)+"_iter.csv"));

        }   
        
        debug3(global_optimal_cost);
        return global_optimal_cost;
    }
};

#endif