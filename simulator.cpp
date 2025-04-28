#include "simulator.hpp"
#define DEBUG3
#define MACRO2
#define SPRING 1


Simulator::Simulator( int num_runs){
    this->num_runs = num_runs;
    // allocate memory for output array
    output = new double[num_runs];
}

void Simulator::printoutput(){
    for (int i = 0; i < num_runs; i++) 
        std::cout<< output[i] << std::endl;
}

Simulator::~Simulator(){
    delete [] output; // free memory
}

void Simulator::Optimizer( double min_range, double max_range,  int function_type, int dim , std::string conf_file)
{
    this->min_range_rand = min_range;
    this->max_range_rand = max_range;
    this->dimension = dim;
    this->function_type = function_type;

    ClockClass clock_var; 

    SCAParserFile(conf_file); // parse parameters from configuration file :iterations r1 , r3, population number
    #ifdef SPRING
    sca_params.limits.lower = min_range; 
    sca_params.limits.upper = max_range;
    sca_params.dimension = this->dimension; // pass them to sca parameter class
    #endif
    sca_algo = new SineCosineClass(sca_params);

    debug1(sca_algo->population->print_population());

    clock_var.tick();
    for (int i = 0; i < num_runs; i++){
        output[i] = sca_algo->execute(function_type);
        //std::cout<<"optimal result of sine cosine algorithm: " << output[i]<<std::endl;
    }

    compute_stats(clock_var.tock()); // compute the statistics

    #ifdef MACRO1 
        save_stats("sca_algorithm");
    #endif

    #ifdef MACRO2
        save_stats("modified_sca_algorithm");
    #endif

    #ifdef MACRO3 
        save_stats("chaotic_sca_algorithm");
    #endif

    #ifdef MACRO4
        save_stats("modified_exp_sca_algorithm");
    #endif

    #ifdef MACRO5
        save_stats("greedy_levy_sca_algorithm");
    #endif

    #ifdef MACRO6
        save_stats("centroid_opposition_sca_algorithm");
    #endif

    #ifdef MACRO7
        save_stats("opposition_sca_algorithm");
    #endif

    #ifdef MACRO8
        save_stats("oblaquila_sca_algorithm");
    #endif

    #ifdef MACRO9
        save_stats("springproblem_sca_algorithm");
    #endif

    #ifdef MACRO10
        save_stats("springproblem_chaotic_sca_algorithm");
    #endif

    #ifdef MACRO11
        save_stats("springproblem_levy_sca_algorithm");
    #endif

    #ifdef MACRO12
        save_stats("springproblem_aquila_sca_algorithm");
    #endif


    delete sca_algo;

    debug1(std::cout<<"output of sine cosine algorithm: " <<std::endl);
    debug1(printoutput());
}

//https://stackoverflow.com/questions/2114797/compute-median-of-values-stored-in-vector-c
void Simulator::compute_stats(double time_milsec)
{
    statistical_analysis.time_milsec = time_milsec;
    
    // sort the array of outputs
    std::sort(output, output + num_runs);
    
    // find the median_value
    // even case 
    if (num_runs % 2 != 0) 
        statistical_analysis.median_value = output[num_runs / 2]; 
    else
        statistical_analysis.median_value = (output[(num_runs - 1) / 2] + output[num_runs / 2]) / 2.0;
    
    // calculate range
    statistical_analysis.max_range = output[num_runs-1];
    statistical_analysis.min_range = output[0];
    statistical_analysis.range = statistical_analysis.max_range - statistical_analysis.min_range;
    
    // calculate mean_value
    double mean_value = 0.0;
    for(int i=0 ; i<num_runs; i++)
    {
        mean_value+=output[i];
    }
    mean_value = mean_value/num_runs;
    statistical_analysis.mean_value = mean_value;
    

    // calculate standard deviation
    double standardDeviation = 0.0;
    for(int i=0 ; i<num_runs; i++)
    {
        standardDeviation += (output[i] - mean_value)*(output[i] - mean_value);
    }
    standardDeviation /= num_runs;
    standardDeviation = sqrt(standardDeviation);
    statistical_analysis.standard_deviation = standardDeviation;
}

void Simulator::SCAParserFile(std::string conf_file)
{
    // pointer of file
    std::fstream instream; 

    // Open an SCA file which is already created 
    instream.open(conf_file, std::ios::in); 
    int counter_row = 0; 

    // Read the words from the file 
    // as strings (vectors)
    std::vector<std::string> row_vec; 
    std::string line_var, word_var; 

    while (instream) { 
        counter_row++;
        row_vec.clear();

        // read one row and store it in line variable
        std::getline(instream, line_var); 

        // used for splitting words 
        std::stringstream strvar(line_var); 
        
        while (std::getline(strvar, word_var, ',')) { 
            // add data to a vector
            row_vec.push_back(word_var); 
        } 
        
        // Compare the row number 
        if (counter_row==2) { 
            sca_params.iter_max = std::stoi(row_vec[0]); 
            sca_params.a = std::stod(row_vec[1]); 
            sca_params.population_size = std::stod(row_vec[2]);
            break; 
        } 
    } 

    if(counter_row!=2)
        std::cout <<conf_file<< " failed to find config file \n"; 

    instream.close();
}

void Simulator::save_stats(const std::string& input)
{
    // file pointer 
    std::fstream outstream; 
  
    // opens a created csv file or creates a new file. 
    std::string sca_file_names = "sca_algorithm_f"+std::to_string(function_type)+".csv";
  
    if(!std::ifstream(sca_file_names)) // if file does not exist, create new one
    {
        outstream.open(sca_file_names, std::fstream::in | std::fstream::out | std::fstream::app);//(append flag)
        outstream    
                << input << ","
                << "function_type" << ","
                << "num_runs" << ","
                << "dimension" << ","
                << "min_range_random" << ","
                << "max_range_random" << ","
                << "mean_value" << ", "
                << "standard_deviation" << ", "
                << "range" << ", "
                << "median_value" << ", "
                << "time_milsec" << ", "
                << "min_range" << ", " 
                << "max_range" << ", " 
                << "\n"; 
        outstream    
                << input << ","
                << function_type << ","
                << num_runs << ","
                << dimension << ","
                << min_range_rand << ","
                << max_range_rand << ","
                << statistical_analysis.mean_value << ", "
                << statistical_analysis.standard_deviation << ", "
                << statistical_analysis.range << ", "
                << statistical_analysis.median_value << ", "
                << statistical_analysis.time_milsec << ", "
                << statistical_analysis.min_range << ", " 
                << statistical_analysis.max_range << ", " 
                << "\n"; 
        outstream.close();
    }
    else
    {
        outstream.open(sca_file_names, std::fstream::in | std::fstream::out | std::fstream::app);
        outstream    
                << input << ","
                << function_type << ","
                << num_runs << ","
                << dimension << ","
                << min_range_rand << ","
                << max_range_rand << ","
                << statistical_analysis.mean_value << ", "
                << statistical_analysis.standard_deviation << ", "
                << statistical_analysis.range << ", "
                << statistical_analysis.median_value << ", "
                << statistical_analysis.time_milsec << ", "
                << statistical_analysis.min_range << ", "  
                << statistical_analysis.max_range << ", " 
                << "\n"; 
        outstream.close();
    }
}