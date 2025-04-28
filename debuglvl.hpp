#ifndef _DEBUGLVL_H_
#define _DEBUGLVL_H_

//comment out to use
//#define DEBUG1
#define DEBUG2
#define DEBUG3

    //std::cout << y << std::endl;
    #ifdef DEBUG1
        #define debug1(y) y
    #else
        #define debug1(y) 
    #endif 

    #ifdef DEBUG3
        #define debug3(y) std::cout << #y << " = " << y << std::endl; // var
    #else
        #define debug3(y)
    #endif 

    #ifdef DEBUG2 // file
        #define debug2(y) y
    #else
        #define debug2(y) 
    #endif 



#endif