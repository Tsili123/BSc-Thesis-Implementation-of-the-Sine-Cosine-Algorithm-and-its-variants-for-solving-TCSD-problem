
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

        // search space for spring problem 
        double upper_x = 2.0; 
        double lower_x = 0.05; 
        double upper_y = 1.3; 
        double lower_y = 0.25; 
        double upper_z = 15.0; 
        double lower_z = 2.0; 
    };

    Limits limits; // instance of struct #Limits
    int population_size; 
    int dimension = 3; 
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

    void stay_in_range(double & val , double parameterdown , double parameterup){
        if(val < parameterdown)
            val = parameterdown;
        else if (val > parameterup)
            val = parameterup;
    }
    
    void global_optimal_set_data(double* data){
        for (int i = 0; i<param.dimension; i++)
            global_optimal_data[i] = data[i];
    }

    public:
    PopulationCriterion *population; // population
    PopulationCriterion *V; // V

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
        std::cout << "limits.lowerx: " << param.limits.lower_x << std::endl;
        std::cout << "limits.uppperx: " << param.limits.upper_x << std::endl;
        std::cout << "limits.lowery: " << param.limits.lower_y << std::endl;
        std::cout << "limits.upppery: " << param.limits.upper_y << std::endl;
        std::cout << "limits.lowerz: " << param.limits.lower_z << std::endl;
        std::cout << "limits.uppperz: " << param.limits.upper_z << std::endl;
        std::cout << "population size: " << param.population_size<< std::endl;
        std::cout << "maximum number of iterations: " << param.iter_max << std::endl;
        std::cout << "dimension: " << param.dimension << std::endl;
        
        //allocate memory for population citerion class
        population = new PopulationCriterion(param.population_size, param.dimension);

        //data: allocate memory for optimal data vector
        global_optimal_data = new double[param.dimension];

        //allocate memory for population citerion class
        V = new PopulationCriterion(1, param.dimension);

        //initialize seed of random generator
        rng.seed(time(NULL));
    }

    ~SineCosineClass(){
        if (global_optimal_data)
            delete[] global_optimal_data;

        if (population)
            delete population;
    }

    double execute(int function_id){

        // initialize distribution in order to fetch random numbers
        std::uniform_real_distribution<double> mstw(0.0,2*M_PI);
        std::uniform_real_distribution<double> mstw1(0,2.0);
        std::uniform_real_distribution<double> mstw2(0.0,1.0);

        for (int i = 0; i< param.dimension; i++){
            if(i==0)
                population->assign_random2(param.limits.lower_x, param.limits.upper_x,i);
            if(i==1)
                population->assign_random2(param.limits.lower_y, param.limits.upper_y,i);
            if(i==2)
                population->assign_random2(param.limits.lower_z, param.limits.upper_z,i);
        }

        debug1(population->print_population());
        population->evaluate_costval(function_id);

        // save first global optimal value
        population->calculate_min_costval();
        global_optimal_cost = population->get_min_costval();
        global_optimal_set_data(population->get_min_costval_data());

        debug3(global_optimal_cost);

        double lamda ,y_k;
        for (int tvar = 0; tvar<param.iter_max; tvar++){

            double r1 = param.a - tvar*((double)param.a/(double)param.iter_max);
            //double r1 = 4*(1-(double)tvar/(double)param.iter_max)*(1-pow(2, ((double)tvar/(double)param.iter_max) -1));

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

                    if(j==0)
                        stay_in_range(new_xvar,param.limits.lower_x, param.limits.upper_x);
                    if(j==1)
                        stay_in_range(new_xvar,param.limits.lower_y, param.limits.upper_y);
                    if(j==2)
                        stay_in_range(new_xvar,param.limits.lower_z, param.limits.upper_z);

                    population->set_data(new_xvar, i, j);
                }
            }

            // update the new optimal cost
            population->evaluate_costval(function_id);
            population->calculate_min_costval();
            if(population->get_min_costval() < global_optimal_cost){
                global_optimal_cost = population->get_min_costval();
                global_optimal_set_data(population->get_min_costval_data());
            }

            lamda = (double)( param.iter_max - tvar + 1 )/  (double)param.iter_max;
            
            if(tvar==0)
                y_k = 0.87;//if I pick a value from the set{0.25,0.5,0.75,1} global minimum 
                //increases significantly as it is stated on the paper
            else
                y_k = 4 * y_k *(1-y_k);

            
            //for (int k = 0; k<param.population_size; k++){
                // update data of every dimension
                for(int i=0; i<param.dimension; i++){

                    double var = (1-lamda)*global_optimal_data[i] + lamda * (param.limits.lower + y_k* (param.limits.upper - param.limits.lower) );
                    if(i==0)
                        stay_in_range(var,param.limits.lower_x, param.limits.upper_x);
                    if(i==1)
                        stay_in_range(var,param.limits.lower_y, param.limits.upper_y);
                    if(i==2)
                        stay_in_range(var,param.limits.lower_z, param.limits.upper_z);

                    V->set_data(var, 0, i);
                }
               // y_k = 4 * y_k *(1-y_k);
            //}

            // update the new optimal cost
            V->evaluate_costval(function_id);
            V->calculate_min_costval();
            // if necessary , update the global min
            if(V->get_min_costval() < global_optimal_cost){
                global_optimal_cost = V->get_min_costval();
                global_optimal_set_data(V->get_min_costval_data());
            }
            //debug2(save_result(global_optimal_cost,global_optimal_data, "sinecosine_f"+std::to_string(function_id)+"_iter.csv"));
        }   
        
        debug2(save_result(global_optimal_cost,global_optimal_data, "sinecosine_f"+std::to_string(function_id)+"_iter.csv"));
        debug3(global_optimal_cost);
        return global_optimal_cost;
    }
};

#endif