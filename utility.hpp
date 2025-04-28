#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <fstream>
#include <ctime> //time
#include <iostream>
#include "debuglvl.hpp"


/// @brief this class  calculates and register  the CPU processing time 

class ClockClass{
    private:   
    std::clock_t start_clock, stop_clock; 

    public:

    /// @brief start clock before code block execution
    void tick(){
        start_clock = std::clock();
    }

    /// @brief register the clock when the code end, also return the elapsed cpu time 
    /// @return return the cpu elapsed time 
    double tock()
    {
        stop_clock = std::clock();
         // CLOCKS_PER_SEC=1000000 in linux enviroment
        return  ((double)stop_clock - (double)start_clock)/CLOCKS_PER_SEC*1000.0; // convert to milisecond
    }
};

#endif