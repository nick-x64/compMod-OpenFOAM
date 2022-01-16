#ifndef TIMER_H
#define TIMER_H

#include "stdafx.h"

class Timer
{
public:

    Timer()
    {
        start = std::chrono::high_resolution_clock::now();
    }

    ~Timer()
    {
        end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<float, std::ratio<1, 1>> duration = end - start;

        std::cout << "DURATION: " << duration.count() << " seconds\n";
    }
private:

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
};



#endif // TIMER_H
