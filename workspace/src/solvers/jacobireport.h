#ifndef JACOBI_JACOBIREPORT_H
#define JACOBI_JACOBIREPORT_H


#include <ostream>

class JacobiReport
{
public:
    friend std::ostream &operator<<(std::ostream &os, const JacobiReport &report)
    {
        float latency = report.compTime + report.updateTime + report.convTime;
        os << report.nWorkers << "," << report.nIterations << "," << report.compTime << "," \
                                        << report.updateTime << "," << report.convTime << "," << latency << "," << report.error;
        return os;
    }

    int nIterations;
    float error;
    int nWorkers;
    float compTime;
    float updateTime;
    float convTime;
};


#endif //JACOBI_JACOBIREPORT_H
