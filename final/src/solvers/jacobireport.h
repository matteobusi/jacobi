#ifndef JACOBI_JACOBIREPORT_H
#define JACOBI_JACOBIREPORT_H


#include <ostream>

class JacobiReport
{
public:
    friend std::ostream &operator<<(std::ostream &os, const JacobiReport &report)
    {
        os << report.nWorkers << "," << report.nIterations << "," << report.updateTime << "," << report.convTime << "," << report.latency << "," << report.error;
        return os;
    }

    int nIterations;
    float error;
    int nWorkers;
    float updateTime;
    float convTime;
    float latency;
};


#endif //JACOBI_JACOBIREPORT_H
