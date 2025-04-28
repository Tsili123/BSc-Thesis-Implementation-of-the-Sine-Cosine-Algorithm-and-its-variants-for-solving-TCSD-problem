#include "functions.hpp"

typename Functions::functionctr_pointer Functions::get_function_type(int type_id)
{
// return the address of a function
functionctr_pointer funcp;

switch(type_id)
{
       case 1:
         funcp = &Functions::func1_Schwefels;
       break;
       case 2:
         funcp = &Functions::func2_DeJongs1;
       break;
       case 3:
         funcp = &Functions::func3_Rastrigin;
       break;
       case 4:
         funcp = &Functions::func4_Griewank;
       break;
       case 5:
         funcp = &Functions::func5_Ackleys;
       break;
       case 6:
         funcp = &Functions::func6_f6step;
       break;
       case 7:
         funcp = &Functions:: func7_Goldstein_price;
       break;
       case 8:
         funcp = &Functions:: func8_6H_Camel_Back;
       break;
       case 9:
         funcp = &Functions::func9_Rosenbrock;
       break;
       case 10:
         funcp = &Functions::func10_Quartic;
       break;
       case 11:
         funcp = &Functions::func11_Schwefels2_21;
       break;
       case 12:
         funcp = &Functions::func12_Schwefels2_22;
       break;
       case 13:
         funcp = &Functions::func13_Schwefels2_12;
       break;
       case 14:
         funcp = &Functions::func14_Generalized_Penalized2;
       break;
       case 15:
         funcp = &Functions::func15_Generalized_Penalized;
       break;
       case 16:
         funcp = &Functions::func16_spring;
       break;
}
return funcp;
}

//multimodal 
double Functions::func1_Schwefels(double* dataset, int dimension){//8

  double offset = 418.9829 * dimension;
  double result_sum;

    for (int i = 0; i < dimension; i++) {
       result_sum += (dataset[i] * sin(sqrt(abs(dataset[i]))));
   }

    return (offset - result_sum);
}

double Functions::func2_DeJongs1(double* dataset, int dimension){
   double result_sum = 0.0;

    for (int i = 0; i < dimension; i++)
       result_sum += dataset[i]*dataset[i];
    return result_sum;
}

double Functions::func3_Rastrigin(double* dataset, int dimension){
    double result_sum = 0.0;

      for (int i = 0; i < dimension; i++)
       result_sum += (dataset[i] * dataset[i] - 10.0 * cos(2.0 * M_PI * dataset[i]));

     result_sum = result_sum+  10.0 * dimension ;

    return result_sum;
}


double Functions::func4_Griewank(double* dataset, int dimension){
    double result_data = 0.0;
    double totsum = 0.0;
    double totproduct = 1.0;

    for (int i = 0; i < dimension; i++)
    {
         totsum += dataset[i]*dataset[i] / 4000.0;
    }

    for (int i = 0; i < dimension; i++)
    {
       totproduct *= cos(dataset[i] / sqrt(i + 1.0));
    }

    result_data = 1.0 + totsum - totproduct;
    return result_data;
}


double Functions::func5_Ackleys(double* dataset, int dimension)
{  
    double result_sum=0.0;
    double prd1=0.0;
    double prd2=0.0;

   for (int i = 0; i < (dimension); i++)
   {
      prd1 += dataset[i] * dataset[i];
      prd2 += cos(2.0 * M_PI * dataset[i]);
   }

   return (-1.0)*20.0 * exp((-1.0) * 0.2 * sqrt((1.0/dimension) * prd1)) - exp((1.0/dimension) * prd2) + exp(1) + 20.0;

   return result_sum;
}

//unimodal
double Functions::func6_f6step(double* dataset, int dimension)
{
    double result_sum = 0.0;
    for (int i = 0; i < dimension; i++)
    {
        result_sum += (floor(dataset[i] + 0.5)) * (floor(dataset[i] + 0.5));
    }

    return result_sum;
}

//multimodal
// -2 <= x(i) <= 2 X
double Functions::func7_Goldstein_price(double* dataset, int dimension ){

  return(( 1.0 + pow(dataset[0]+dataset[1]+1,2)*(19.0-14.0*dataset[0]+3.0*pow(dataset[0],2)-14.0*dataset[1]+6.0*dataset[0]*dataset[1]+3.0*pow(dataset[1],2)))*(30+pow(2.0*dataset[0]-3.0*dataset[1],2)*(18.0-32.0*dataset[0]+12.0*pow(dataset[0],2)+48.0*dataset[1]-36.0*dataset[0]*dataset[1]+27*pow(dataset[1],2)) ));

}

// X
double Functions::func8_6H_Camel_Back(double* dataset, int dimension ){

  return(4.0*pow(dataset[0],2)-2.1*pow(dataset[0],4)+pow(dataset[0],6)/3.0+dataset[0]*dataset[1]-4.0*pow(dataset[1],2)+4.0*pow(dataset[1],4));

}

double Functions::func9_Rosenbrock(double* dataset, int dimension)
{
    double result_sum = 0.0;
    for (int i = 0; i < dimension-1; i++)
    {
        result_sum += 100.0*pow(dataset[i+1]-pow(dataset[i],2),2)+pow(1.0-dataset[i],2);
    }

    return result_sum;
}

double Functions::func10_Quartic(double* dataset, int dimension)
{
    std::mt19937 rng2; // mersenne twister random generator

    //initialize seed of random generator
    rng2.seed(time(NULL));

    std::uniform_real_distribution<double> mst(0.0,1.0);

    double result_sum = 0.0;
    for (int i = 0; i < dimension; i++)
    {
        result_sum += (i+1.0)*pow(dataset[i],4);
    }

    return result_sum+ mst(rng2);
}

//multimodal
double Functions::func11_Schwefels2_21(double* dataset, int dimension){

  double max = -1.0;

    for (int i = 0; i < dimension; i++) {
       if(max< fabs(dataset[i]) )
            max = fabs(dataset[i]);
    }

    return (max);
}

//multimodal X
double Functions::func12_Schwefels2_22(double* dataset, int dimension){

  double result_sum=0.0,result_mul=1.0;

    for (int i = 0; i < dimension; i++) {
       result_sum+= fabs(dataset[i]);
       result_mul*= fabs(dataset[i]);
    }

    result_sum+=result_mul;
    return (result_sum);
}

//multimodal
double Functions::func13_Schwefels2_12(double* dataset, int dimension){

  double result_sum=0.0,result_sum2=0.0;

    for (int i = 0; i < dimension; i++) {
      result_sum2=0.0;
       for (int j = 0; j < i; j++) {
          result_sum2+= dataset[j];
        }
        result_sum+=pow(result_sum2,2);
    }

    return (result_sum);
}

double y(double dataset) {
        return 1.0 + (dataset + 1.0)/4.0;
}

double u(double dataset, double a, double k, double m) {
        if (dataset > a)
            return k * pow((dataset - a), m);
        else if (dataset < a)
            return k * pow((-dataset- a), m);
        else 
            return 0;
}

double Functions::func14_Generalized_Penalized2(double* dataset, int dimension){

      double result_sum = 0.0;
      result_sum += 10.0*pow(sin(3.0 * M_PI * dataset[0]), 2);
      
      for (int i = 0; i < dimension-1; i++)
        result_sum += pow(dataset[i]-1.0, 2) * (1 + pow(sin(3.0 * M_PI * dataset[i+1]), 2));

      result_sum += (dataset[dimension-1]-1.0) * (dataset[dimension-1]-1.0) ;
      result_sum *= 0.1;

      for (int i = 0; i < dimension; i++)
        result_sum += u(dataset[i], 5, 100, 4);

      return result_sum;
}

// Generalized penalized function 
double Functions::func15_Generalized_Penalized(double *dataset,int dimension) {

    double result_sum = 0.0;
    result_sum += 10.0 * pow(sin( M_PI * y(dataset[0]) ), 2);

    for (int i = 0; i < dimension -1; i++)
      result_sum += pow(y(dataset[i])-1.0, 2) * (1.0 + 10.0* pow( sin(M_PI * y(dataset[i+1])) , 2) );

    result_sum += pow( y(dataset[dimension]) - 1.0 , 2);
    result_sum *= M_PI / dimension ;

    for (int i = 0; i < dimension ; i++)
      result_sum += u(dataset[i], 10, 100, 4);

    return result_sum;
  }

double Functions::func16_spring(double* dataset,int dimension){
  //dataset[0] : x , dataset[1] : y , dataset[2] : z

    //regularization term 
    double R = 10000;
    double p1{0.0}, p2{0.0}, p3{0.0}, p4{0.0};
    double d1 = 1.0 - ((double)(std::pow(dataset[1], 3) * dataset[2]) / (double)( 71785.0 * std::pow(dataset[0], 4) ) );

    double d2 = ( (double)(4.0 * std::pow(dataset[1], 2) - dataset[0] * dataset[1]) / (double)(12556.0 * (dataset[1] * std::pow(dataset[0], 3) - std::pow(dataset[0], 4)))) + (1.0 / (double)(5108.0 * std::pow(dataset[0], 2))) -1.0;

    double d3 = 1.0 -((double)(140.45 * dataset[0]) / (double)(std::pow(dataset[1], 2) * dataset[2]));

    double d4 = ((double)(dataset[0] + dataset[1]) / (double)1.5) - 1.0;

    if (d1 > 0)
        p1 = std::pow(d1, 2);

    if (d2 > 0)
        p2 = std::pow(d2, 2);

    if (d3 > 0)
        p3 = std::pow(d3, 2);

    if (d4 > 0)
        p4 = std::pow(d4, 2);

    return (dataset[2] + 2) * dataset[1] * std::pow(dataset[0], 2) + R * p1 + R * p2 + R * p3 + R * p4;
}



