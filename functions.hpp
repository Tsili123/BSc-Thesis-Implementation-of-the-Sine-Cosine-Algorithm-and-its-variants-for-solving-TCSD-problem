#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include "populations.hpp"
#include <math.h>  
#include <cstdlib>

// class with 7 standard criterion functions
// input of function get_function_type is the function id type
// returns the function pointer

class Functions
{
    private:
    // f1: Schwefel's function
    // dataset array parameter
    // dimension Size of the dataset array, (number of dimensions) parameter
    // returns numeric result of the function(double)  parameter
    double func1_Schwefels(double* dataset, int dimension);

    // f2: 1st De Jong function
    double func2_DeJongs1(double* dataset, int dimension);

    // f3: Rastrigin function
    double func3_Rastrigin(double* dataset, int dimension);

    // f4: Griewank function
    double func4_Griewank(double* dataset, int dimension);

    // f5: Ackleys function
    double func5_Ackleys(double* dataset, int dimension);

    // f6: f6 step
    double func6_f6step(double* dataset, int dimension);

    // f7: Goldstein_price X
    double func7_Goldstein_price(double* dataset, int dimension);

    // f8: Goldstein_price
    double func8_6H_Camel_Back(double* dataset, int dimension);

    // f9: Rosenborg -
    double func9_Rosenbrock(double* dataset, int dimension);

    // f10: Quartic
    double func10_Quartic(double* dataset, int dimension);

    // f11: Schwefels2_21 
    double func11_Schwefels2_21(double* dataset, int dimension);

    // f12: Schwefels2_22 X
    double func12_Schwefels2_22(double* dataset, int dimension);

    // f13: Schwefels2_12 X
    double func13_Schwefels2_12(double* dataset, int dimension);

    // f14: Generalized_Penalize2
    double func14_Generalized_Penalized2(double* dataset, int dimension);

    // f15: func11_Schwefels2_12
    double func15_Generalized_Penalized(double* dataset, int dimension);

    // f16: spring problem
    double func16_spring(double* dataset, int dimension);

    public:

    // Type defenition of a function pointer
    typedef double (Functions::*functionctr_pointer) (double*, int); // function pointer
    
    // Receive the function pointer by parameter identifier of the criterion function
    /// @param id code identifier of the criterion function
    /// @return Function pointer of the  criterion function

    functionctr_pointer get_function_type(int id);
};


#endif