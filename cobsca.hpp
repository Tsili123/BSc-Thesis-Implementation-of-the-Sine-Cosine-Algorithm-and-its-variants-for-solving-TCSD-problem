
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


class SineCosineClass
{
    private:
    
    std::mt19937 rng; // mersenne twister random generator
    SCAParamClass param; // SCA parameters

    double global_optimal_cost; //  global optimal cost
    double* global_optimal_data; // global optimal data
    double* centroid; // global optimal data
    double* Minj; // Minj
    double* Maxj; // Maxj

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
    PopulationCriterion *var_population; // opposite population
    PopulationCriterion *variable_population;
    PopulationCriterion *opvariable_population;

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
        var_population = new PopulationCriterion(param.population_size, param.dimension);

        //allocate memory for opposite population citerion class
        variable_population = new PopulationCriterion(1, param.dimension);

        //allocate memory for opposite population citerion class
        opvariable_population = new PopulationCriterion(1, param.dimension);

        //global_optimal_data : allocate memory for optimal data vector
        global_optimal_data = new double[param.dimension];

        //centroid: allocate memory for centroid vector
        centroid = new double[param.dimension];

        //Minj : allocate memory for Minj vector
        Minj = new double[param.dimension];

        //Maxj : allocate memory for Maxj vector
        Maxj = new double[param.dimension];

        //initialize seed of random generator
        rng.seed(time(NULL));
    }

    ~SineCosineClass(){
        // if (global_optimal_data)
        //     delete[] global_optimal_data;
        
        // if (centroid)
        //     delete[] centroid;

        // if (Maxj)
        //     delete[] Maxj;

        // if (Minj)
        //     delete[] Minj;

        // if (population)
        //     delete population;

        // if (op_population)
        //     delete op_population;
        
        // if (var_population)
        //     delete var_population;

        // if (variable_population)
        //     delete variable_population;

        // if (opvariable_population)
        //     delete opvariable_population;
    }

    double execute(int function_id){

        // initialize distribution in order to fetch random numbers
        std::uniform_real_distribution<double> mstw(0.0,2*M_PI);
        std::uniform_real_distribution<double> mstw1(0,2.0);
        std::uniform_real_distribution<double> mstw2(0.0,1.0);
        std::uniform_real_distribution<double> rand(0.0,1.0);
        std::uniform_real_distribution<double> rand2(0.0,1.0);

        population->assign_random(param.limits.lower, param.limits.upper);
        debug1(population->print_population());
        //population->evaluate_costval(function_id);

        double sum;
        double min = param.limits.upper + 1;
        double max = param.limits.lower - 1;
        double xij;
        double op_var;

        for(int j = 0; j < param.dimension; j++){
            sum=0.0;
            // create opposite points population
            for (int i = 0; i < param.population_size; i++){

                xij = population->get_data(i,j);
                // update data of every dimension
                sum += xij;

                if(xij < min){
                    min = xij;
                    Minj[j] = min;
                }

                if(xij > max){
                    max = xij;
                    Maxj[j] = max;
                }
            }
            centroid[j]=sum;
        }

        for(int j = 0; j < param.dimension; j++){
            
             centroid[j] = (double)centroid[j]/(double)param.population_size;
        }    

        // create opposite points population
        for (int i = 0; i<param.population_size; i++){
            // update data of every dimension
            for(int j = 0; j<param.dimension; j++){

                op_var = population->get_data(i,j);
                op_var = 2.0*centroid[j] - op_var;

                if(op_var < param.limits.lower){
                    op_var = Minj[j] + ( centroid[j] - Minj[j] )*rand(rng);
                    op_population->set_data(op_var,i,j);
                }else if(op_var >= param.limits.upper){
                    op_var = centroid[j] + ( Maxj[j] - centroid[j] )*rand(rng);
                    op_population->set_data(op_var,i,j);
                }else{
                    op_population->set_data(op_var,i,j);
                }
            }
        }


        op_population->evaluate_costval(function_id);
        // save first global optimal value
        op_population->calculate_min_costval();
        global_optimal_cost = op_population->get_min_costval();
        global_optimal_set_data(op_population->get_min_costval_data());

        debug3(global_optimal_cost);

        double Pgj,Pgj_init=0.8,Pgj_min=0.01;
        

        for (int tvar = 0; tvar<param.iter_max; tvar++){

            //Pgj = Pgj_init * exp(-0.001*((double)tvar/(double)param.iter_max));

            Pgj = Pgj_init - (Pgj_init - Pgj_min) * ((double)tvar/(double)param.iter_max);
            //Pgj = (1 - exp(-0.001 * tvar)) * Pgj_init;
            // double valp = Pgj_init - ((double)tvar/(double)param.iter_max)*Pgj_init;

            // if(valp > Pgj_min)
            //     Pgj = valp;
            // else   
            //     Pgj = Pgj_min;

            if(rand(rng) < Pgj ){

                for (int i = 0; i<param.population_size; i++){

                    // update data of every dimension
                    for(int j =0; j<param.dimension; j++){

                        op_var = population->get_data(i,j);
                        op_var = 2.0*centroid[j] - op_var;

                        if(op_var < param.limits.lower){
                            op_var = Minj[j] + ( centroid[j] - Minj[j] )*rand2(rng);
                            var_population->set_data(op_var,i,j);
                        }else if(op_var >= param.limits.upper){
                            op_var = centroid[j] + ( Maxj[j] - centroid[j] )*rand2(rng);
                            var_population->set_data(op_var,i,j);
                        }else{
                            var_population->set_data(op_var,i,j);
                        }
                    }
                }

                var_population->evaluate_costval(function_id);

                for (int i = 0; i<param.population_size; i++){
                      if(var_population->get_costval(i) < op_population->get_costval(i) ){
                            for(int j=0; j<param.dimension; j++){
                                double var = var_population->get_data(i,j);
                                op_population->set_data(var,i,j);
                            }
                      }
                }


            }else{

                double r1 = param.a - tvar*((double)param.a/(double)param.iter_max);
                //double r1 = 4*(1-(double)tvar/(double)param.iter_max)*(1-pow(2, ((double)tvar/(double)param.iter_max) -1));

                for (int i = 0; i<param.population_size; i++){

                    // update data of every dimension
                    for(int j =0; j<param.dimension; j++){

                        double xvar = op_population->get_data(i,j);
                        double new_xvar;
                        double r2 = mstw(rng);
                        double r3 = mstw1(rng);
                        double r4 = mstw2(rng);

                        if (r4 < 0.5) //sine
                            new_xvar = xvar + r1*sin(r2)*fabs(r3*global_optimal_data[j]-xvar);
                        else //cosine
                            new_xvar = xvar + r1*cos(r2)*fabs(r3*global_optimal_data[j]-xvar);

                        stay_in_range(new_xvar);
                        op_population->set_data(new_xvar,i,j);
                    }
                }

            }

            // update the new optimal cost
            op_population->evaluate_costval(function_id);
            op_population->calculate_min_costval();
            if(op_population->get_min_costval() < global_optimal_cost){
                global_optimal_cost = op_population->get_min_costval();
                global_optimal_set_data(op_population->get_min_costval_data());
            }

             debug2(save_result(global_optimal_cost,global_optimal_data, "sinecosine_f"+std::to_string(function_id)+"_iter.csv"));
         } 

         debug3(global_optimal_cost);
         return global_optimal_cost;
    }
};

#endif