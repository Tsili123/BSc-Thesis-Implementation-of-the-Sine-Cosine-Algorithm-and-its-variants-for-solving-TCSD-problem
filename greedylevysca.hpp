
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
    double* r_j; 
    double* rmax_j; 
    double* theta_j; 

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
    PopulationCriterion *gs_best; // population
    double * rmax;
    double * r;
    double * theta;

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
        gs_best = new PopulationCriterion(1, param.dimension);

        //data: allocate memory for optimal data vector
        global_optimal_data = new double[param.dimension];

        //data: allocate memory for r_j
        r = new double[param.dimension];

        //data: allocate memory for theta_j
        theta = new double[param.dimension];

        //data: allocate memory for rmax_j
        rmax = new double[param.dimension];

        //initialize seed of random generator
        rng.seed(time(NULL));
    }

    ~SineCosineClass(){
        if (global_optimal_data)
            delete[] global_optimal_data;

        if (r)
            delete r;

        if (theta)
            delete theta;

        if (rmax)
            delete rmax;

    }

    double execute(int function_id){

        // initialize distribution in order to fetch random numbers
        std::uniform_real_distribution<double> mstw(0.0,2*M_PI);
        std::uniform_real_distribution<double> mstw1(0,2.0);
        std::uniform_real_distribution<double> mstw2(0.0,1.0);
        std::uniform_real_distribution<double> mstw3(-1.0,1.0);

        double mean = 0.0;
        double stddev  = 1.0;
        std::normal_distribution<double> u(mean, stddev);
        double beta = 3.0/2.0;
        stddev = pow( ( tgamma(1.0+beta) * sin( M_PI*beta/2.0) )/( tgamma((1.0+beta)/2.0) * beta * pow(2,(beta-1.0)/2.0) ) , 1.0/beta );
        std::normal_distribution<double> v(mean, stddev);

        population->assign_random(param.limits.lower, param.limits.upper);
        debug1(population->print_population());
        population->evaluate_costval(function_id);

        // save first global optimal value
        population->calculate_min_costval();
        global_optimal_cost = population->get_min_costval();
        global_optimal_set_data(population->get_min_costval_data());

        debug3(global_optimal_cost);

        for (int tvar = 0; tvar<param.iter_max; tvar++){

            double r1 = param.a * exp((double)param.a/(double)param.iter_max);
            double omega_max = 0.9;
            double omega_min = 0.4;
            double omega = omega_max - (omega_max - omega_min) * ((double)param.a/(double)param.iter_max);

            for (int i = 0; i<param.population_size; i++){

                // update data of every dimension
                for(int j =0; j<param.dimension; j++){

                    double xvar = population->get_data(i,j);
                    double new_xvar;
                    double r2 = mstw(rng);
                    double r3 = mstw1(rng);
                    double r4 = mstw2(rng);
                    double urd = mstw3(rng);
                    double lambda = 0.01;

                    if (r4 < 0.5) //sine
                        new_xvar = omega * xvar + r1*sin(r2)*fabs(r3* global_optimal_data[j] * (1+lambda*urd) - xvar);
                    else //cosine
                        new_xvar = omega * xvar + r1*cos(r2)*fabs(r3* global_optimal_data[j] * (1+lambda*urd) - xvar);

                    stay_in_range(new_xvar);
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

            double sum , max , min, levy;
            for(int j=0;j<param.dimension;j++){

                sum = 0.0;
                //compute series
                for(int i=0;i<param.population_size;i++){
                    sum += population->get_data(i, j);
                }
                sum = sum/param.population_size;

                r[j] = abs( global_optimal_data[j] - sum );
            
                max = param.limits.lower;
                min = param.limits.upper;
                for(int i=0;i<param.population_size;i++){
                    double var = population->get_data(i, j);
                    if( var > max)
                        max =  var;
                    
                    if( var < min)
                        min =  var;
                }

                rmax[j] = max - min ; 

                theta[j] = exp( ( -30.0 * tvar / (double)param.iter_max) * ( 1.0- ((double)r[j]/(double)rmax[j]) ) );

                levy = u(rng) / pow (abs(v(rng)) , 1.0/beta );

                double xvarl = (global_optimal_data[j] + theta[j] * levy * global_optimal_data[j]);
                stay_in_range(xvarl);
                gs_best->set_data( xvarl , 0, j ); 
            }

            // update the new optimal cost
            gs_best->evaluate_costval(function_id);
            gs_best->calculate_min_costval();
            // if necessary , update the global min
            if(gs_best->get_min_costval() < global_optimal_cost){
                global_optimal_cost = gs_best ->get_min_costval();
                global_optimal_set_data(gs_best->get_min_costval_data());
            }

            debug2(save_result(global_optimal_cost,global_optimal_data, "sinecosine_f"+std::to_string(function_id)+"_iter.csv"));
        }   
        
        debug3(global_optimal_cost);
        return global_optimal_cost;
    }
};

#endif