#include "timer.hpp"
#include <iostream>


Timer::Timer(){
    //struct timeval t_start;
    //struct timeval t_now;
    t_cumul=0.;
    running=0;
}

void Timer::start(){
    if(running){
        std::cerr<<"timer already running."<<std::endl;
        return;
    }
    gettimeofday(&t_start,NULL);
    running=1;
}

void Timer::reset(){
    t_cumul=0.;
    running=0;
}

double Timer::stop(){
    gettimeofday(&t_now,NULL);
    if(running){
        t_cumul+=(t_now.tv_sec + 1e-6 * t_now.tv_usec) -
                    (t_start.tv_sec + 1e-6 * t_start.tv_usec);
        running=0;
        return t_cumul;
    }
    else{
        std::cerr<<"timer is not running."<<std::endl;
        return 0;
    }
}

double Timer::elapsed(){
    return t_cumul;
}


