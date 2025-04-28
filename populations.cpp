#include "populations.hpp"

/// @brief populations constructor

/// save information about number of items and the dimension
/// memory allocation and random seed generation based on time
/// @param num_items number of elements
/// @param dimension dimension or size of each element

Populations::~Populations(){
    // free memory for dataset
    if(dataset){
        for (int i = 0; i<num_items; i++)
        {
            if(dataset[i])
            delete [] dataset[i]; //delete array
        }

        delete [] dataset; // delete array of array
    }
}

Populations::Populations(int num_items, int dimension){
    this->num_items = num_items;
    this->dimension = dimension;

    // random generator seed
    rng.seed(time(NULL));

    // allocate memory for the dataset
    dataset = new double*[num_items];

    for (int i = 0; i<num_items; i++)
        dataset[i] = new double[dimension];
    
}

void Populations::assign_random(double range_low, double range_high)
{       
    std::uniform_real_distribution<double> unipop(range_low,range_high);

    for (int i = 0; i < num_items; i++) {  
        for (int j = 0; j < dimension; j++) {
                dataset[i][j] = unipop(rng);
        }
    }
}

void Populations::assign_random2(double range_low, double range_high,int offset){       
    std::uniform_real_distribution<double> unipop(range_low,range_high);

    for (int i = 0; i < num_items; i++) {  
                dataset[i][offset] = unipop(rng);
    }
}

void Populations::print_population(){

    std::cout<< "---populations---" <<std::endl;
    for (int i = 0; i < num_items; i++) {  
        for (int j = 0; j < dimension; j++) {
            std::cout<< dataset[i][j] << " ";
        }
        std::cout<<std::endl;
    }
}

double* Populations::get_data(int i){
    return dataset[i];
}

double Populations::get_data(int i, int j){
    return dataset[i][j];
}

void Populations::set_data(double *item, int i){
    for(int j=0; j<dimension; j++)
        dataset[i][j] = item[j];
}


void Populations::set_data(double data_i, int i, int j){
    dataset[i][j] = data_i;
    
}
