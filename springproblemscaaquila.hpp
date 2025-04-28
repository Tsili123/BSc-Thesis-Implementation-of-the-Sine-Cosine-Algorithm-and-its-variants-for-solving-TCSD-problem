
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
    PopulationCriterion *op_population; // opposite population
    PopulationCriterion *best_population; // opposite population

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
        
        if (arrayp)
            delete arrayp;
    }

    double execute(int function_id){

        // initialize distribution in order to fetch random numbers
        std::uniform_real_distribution<double> mstw(0.0,2*M_PI);
        std::uniform_real_distribution<double> mstw2(0.0,1.0);
        std::uniform_int_distribution<int> r1(0,20);
        std::uniform_int_distribution<int> D1(1,param.dimension);
        std::uniform_int_distribution<int> x_r(0,param.population_size-1);

        //Initializing Aquila Optimizer
        double r = mstw2(rng);
        double xvar;
        double dij;

        for (int i = 0; i< param.dimension; i++){
            if(i==0)
                population->assign_random2(param.limits.lower_x, param.limits.upper_x,i);
            if(i==1)
                population->assign_random2(param.limits.lower_y, param.limits.upper_y,i);
            if(i==2)
                population->assign_random2(param.limits.lower_z, param.limits.upper_z,i);
        }

        debug1(population->print_population());
        // population->evaluate_costval(function_id);

        // // save first global optimal value
        // population->calculate_min_costval();
        // global_optimal_cost = population->get_min_costval();
        // global_optimal_set_data(population->get_min_costval_data());

        // create opposite points population
        for (int i = 0; i<param.population_size; i++){
            // update data of every dimension
            for(int j = 0; j<param.dimension; j++){
        
                xvar = population->get_data(i,j);
                double newxvar;
                
                if(j==0){
                    newxvar = param.limits.upper_x + param.limits.lower_x - xvar;
                    dij = (param.limits.upper_x + param.limits.lower_x)/2.0;
                }
                if(j==1){
                    newxvar = param.limits.upper_y + param.limits.lower_y - xvar;
                    dij = (param.limits.upper_y + param.limits.lower_y)/2.0;
                }
                if(j==2){
                    newxvar = param.limits.upper_z + param.limits.lower_z - xvar;
                    dij = (param.limits.upper_z + param.limits.lower_z)/2.0;
                }
                
                if(xvar < dij){
                    newxvar = dij + (newxvar - dij) * r; 
                }else{
                    newxvar = newxvar +( dij - newxvar)* r;
                }
                if(j==0)
                    stay_in_range(newxvar,param.limits.lower_x, param.limits.upper_x);
                if(j==1)
                    stay_in_range(newxvar,param.limits.lower_y, param.limits.upper_y);
                if(j==2)
                    stay_in_range(newxvar,param.limits.lower_z, param.limits.upper_z);
                op_population->set_data(newxvar, i, j);
            }
        }

        //compute levy for x3
        double mean = 0.0;
        double stddev  = 1.0;
        std::normal_distribution<double> u(mean, stddev);
        double beta = 3.0/2.0;
        stddev = pow( ( tgamma(1.0+beta) * sin( M_PI*beta/2.0) )/( tgamma((1+beta)/2.0) * beta * pow(2,(beta-1.0)/2.0) ) , 1.0/beta );
        std::normal_distribution<double> v(mean, stddev);
        double levy;
        
        //compute constrants for x2
        double p,y,x,theta;

        debug1(population->print_population());
        //compute fitness of new population
        op_population->evaluate_costval(function_id);

        // save first global optimal value
        op_population->calculate_min_costval();
        global_optimal_cost = op_population->get_min_costval();
        global_optimal_set_data(op_population->get_min_costval_data());

        if (std::isnan(global_optimal_cost))
            return global_optimal_cost;
        debug3(global_optimal_cost);

        double sum,newxvar2,xvar2;

        for (int tvar = 0; tvar<param.iter_max; tvar++){
    
            double d1var;
            d1var = D1(rng);
            levy = 0.01 * u(rng) / pow (abs(v(rng)) , 1.0/beta ) ;
            theta = 0.005 * d1var + 3.0 * M_PI / 2.0;
            p = (double)r1(rng) +0.00565 * d1var;
            x = p * cos(theta);
            y = p * sin(theta);

            for (int i = 0; i<param.population_size; i++){
                // update data of every dimension
                for(int j=0; j<param.dimension; j++){

                        if(tvar< 2.0*param.iter_max/3.0){

                            //Expanded exploration
                            if( r < 0.5 ){
                                sum = 0.0;
                                //compute series
                                for(int i=0;i<param.population_size;i++){
                                    sum += op_population->get_data(i, j);
                                }
                                sum = sum/param.population_size;
                           
                                xvar2 = global_optimal_data[j] * (1 - (double)tvar/(double)param.iter_max) + ( sum - global_optimal_data[j] * r );

                                if(j==0){
                                    stay_in_range(xvar2,param.limits.lower_x, param.limits.upper_x);
                                    dij = (param.limits.upper_x + param.limits.lower_x)/2.0;
                                    newxvar2 = param.limits.upper_x + param.limits.lower_x - xvar2;
                                }
                                if(j==1){
                                    stay_in_range(xvar2,param.limits.lower_y, param.limits.upper_y);
                                    dij = (param.limits.upper_y + param.limits.lower_y)/2.0;
                                    newxvar2 = param.limits.upper_y + param.limits.lower_y - xvar2;
                                }
                                if(j==2){
                                     stay_in_range(xvar2,param.limits.lower_z, param.limits.upper_z);
                                     dij = (param.limits.upper_z + param.limits.lower_z)/2.0;
                                     newxvar2 = param.limits.upper_z + param.limits.lower_z - xvar2;
                                }
         

                                if(xvar2 < dij){
                                    newxvar2 = dij + (newxvar2 - dij) * r; 
                                }else{
                                    newxvar2 = newxvar2 +( dij - newxvar2)* r;
                                }

                                if(j==0)
                                    stay_in_range(newxvar2,param.limits.lower_x, param.limits.upper_x);
                                if(j==1)
                                    stay_in_range(newxvar2,param.limits.lower_y, param.limits.upper_y);
                                if(j==2)
                                    stay_in_range(newxvar2,param.limits.lower_z, param.limits.upper_z);

                                op_population->set_data(newxvar2, i, j);
                            }
                            //Narrowed exploration
                            else{
                                    int popvar = x_r(rng);
                                    double* x_rvar = op_population->get_data(popvar);
                                    xvar2 = global_optimal_data[j] * x_rvar[j] * levy * (y - x) * r;
                                    
                                    if(j==0){
                                        stay_in_range(xvar2,param.limits.lower_x, param.limits.upper_x);
                                        newxvar2 = param.limits.upper_x + param.limits.lower_x - xvar2;
                                        dij = (param.limits.upper_x + param.limits.lower_x)/2.0;
                                    }
                                    if(j==1){
                                        stay_in_range(xvar2,param.limits.lower_y, param.limits.upper_y);
                                        newxvar2 = param.limits.upper_y + param.limits.lower_y - xvar2;
                                        dij = (param.limits.upper_y + param.limits.lower_y)/2.0;
                                    }
                                    if(j==2){
                                        stay_in_range(xvar2,param.limits.lower_z, param.limits.upper_z);
                                        newxvar2 = param.limits.upper_z + param.limits.lower_z - xvar2;
                                        dij = (param.limits.upper_z + param.limits.lower_z)/2.0;
                                    }

                                    if(xvar2 < dij){
                                        newxvar2 = dij + (newxvar2 - dij) * r; 
                                    }else{
                                        newxvar2 = newxvar2 +( dij - newxvar2)* r;
                                    }

                                    if(j==0)
                                        stay_in_range(newxvar2,param.limits.lower_x, param.limits.upper_x);
                                    if(j==1)
                                        stay_in_range(newxvar2,param.limits.lower_y, param.limits.upper_y);
                                    if(j==2)
                                        stay_in_range(newxvar2,param.limits.lower_z, param.limits.upper_z);

                                    op_population->set_data(newxvar2, i, j);

                            }
                        }else{
                            if( r < 0.5){
                                
                                    sum = 0.0;
                                    //compute series
                                    for(int i=0;i<param.population_size;i++){
                                        sum += op_population->get_data(i, j);
                                    }
                                    sum = sum/param.population_size;

                                    
                                    if(j==0){
                                        xvar2 = ( global_optimal_data[j] - sum ) * 0.1 - r + ((param.limits.upper_x - param.limits.lower_x) * r + param.limits.lower_x) * 0.1;
                                        stay_in_range(xvar2,param.limits.lower_x, param.limits.upper_x);
                                         newxvar2 = param.limits.upper_x + param.limits.lower_x - xvar2;
                                         dij = (param.limits.upper_x + param.limits.lower_x)/2.0;
                                    }
                                    if(j==1){
                                         xvar2 = ( global_optimal_data[j] - sum ) * 0.1 - r + ((param.limits.upper_y - param.limits.lower_y) * r + param.limits.lower_y) * 0.1;
                                          stay_in_range(xvar2,param.limits.lower_y, param.limits.upper_y);
                                          newxvar2 = param.limits.upper_y + param.limits.lower_y - xvar2;
                                         dij = (param.limits.upper_y + param.limits.lower_y)/2.0;
                                    }
                                    if(j==2){
                                        xvar2 = ( global_optimal_data[j] - sum ) * 0.1 - r + ((param.limits.upper_z - param.limits.lower_z) * r + param.limits.lower_z) * 0.1;
                                        stay_in_range(xvar2,param.limits.lower_z, param.limits.upper_z);
                                        newxvar2 = param.limits.upper_z + param.limits.lower_z - xvar2;
                                        dij = (param.limits.upper_z + param.limits.lower_z)/2.0;
                                    }

                                    if(xvar2 < dij){
                                        newxvar2 = dij + (newxvar2 - dij) * r; 
                                    }else{
                                        newxvar2 = newxvar2 +( dij - newxvar2)* r;
                                    }

                                    if(j==0)
                                        stay_in_range(newxvar2,param.limits.lower_x, param.limits.upper_x);
                                    if(j==1)
                                        stay_in_range(newxvar2,param.limits.lower_y, param.limits.upper_y);
                                    if(j==2)
                                        stay_in_range(newxvar2,param.limits.lower_z, param.limits.upper_z);

                                    op_population->set_data(newxvar2, i, j);
                                }else{
                                    double q_f = pow(tvar, (double)(2*r-1)/ (double) pow(1-param.iter_max,2));
                                    double G1 = 2*r -1;
                                    double G2 = 2*(1-(double)tvar/(double)param.iter_max);

                                    xvar2 = q_f * global_optimal_data[j] - (G1 * op_population->get_data(i, j)* r) - G2*levy + r*G1;

                                    if(j==0){
                                        stay_in_range(xvar2,param.limits.lower_x, param.limits.upper_x);
                                        newxvar2 = param.limits.upper_x + param.limits.lower_x - xvar2;
                                        dij = (param.limits.upper_x + param.limits.lower_x)/2.0;
                                    }
                                    if(j==1){
                                        stay_in_range(xvar2,param.limits.lower_y, param.limits.upper_y);
                                        newxvar2 = param.limits.upper_y + param.limits.lower_y - xvar2;
                                        dij = (param.limits.upper_y + param.limits.lower_y)/2.0;
                                    }
                                    if(j==2){
                                        stay_in_range(xvar2,param.limits.lower_z, param.limits.upper_z);
                                        newxvar2 = param.limits.upper_z + param.limits.lower_z - xvar2;
                                        dij = (param.limits.upper_z + param.limits.lower_z)/2.0;
                                    }

                                    if(xvar2 < dij){
                                        newxvar2 = dij + (newxvar2 - dij) * r; 
                                    }else{
                                        newxvar2 = newxvar2 +( dij - newxvar2)* r;
                                    }

                                    if(j==0)
                                        stay_in_range(newxvar2,param.limits.lower_x, param.limits.upper_x);
                                    if(j==1)
                                        stay_in_range(newxvar2,param.limits.lower_y, param.limits.upper_y);
                                    if(j==2)
                                        stay_in_range(newxvar2,param.limits.lower_z, param.limits.upper_z);

                                    op_population->set_data(newxvar2, i, j);
                                }
                        }
                    }
                }

            // update the new optimal cost
            op_population->evaluate_costval(function_id);
            op_population->calculate_min_costval();
            // if necessary , update the global min
            if(op_population->get_min_costval() < global_optimal_cost){
                global_optimal_cost = op_population ->get_min_costval();
                global_optimal_set_data(op_population->get_min_costval_data());
            }
            r = mstw2(rng);

            debug2(save_result(global_optimal_cost,global_optimal_data, "sinecosine_f"+std::to_string(function_id)+"_iter.csv"));
        }
        
        if (std::isnan(global_optimal_cost))
            return global_optimal_cost;
        debug3(global_optimal_cost);
        return global_optimal_cost;
    }
};

#endif