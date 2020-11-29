#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>

class Timer{
    struct timeval  t_start;
    struct timeval  t_now;

    double  t_cumul;
    int running;

    public:
    Timer();
    void start();
    void reset();
    double stop();
    double elapsed();

};


#endif // TIMER_H
