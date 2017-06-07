//
// Created by caos on 07/02/17.
//

#ifndef JACOBI_JACOBIREPORT_H
#define JACOBI_JACOBIREPORT_H


#include <ostream>

class JacobiReport
{
public:
    friend std::ostream &operator<<(std::ostream &os, const JacobiReport &report)
    {
        double latency = report.compTime + report.updateTime + report.convTime;
        os << report.nWorkers << "," << report.nIterations << "," << report.compTime << "," \
                                        << report.updateTime << "," << report.convTime << "," << latency << "," << report.error;
        return os;
    }

    int nIterations;
    double error;
    int nWorkers;
    double compTime;
    double updateTime;
    double convTime;
};


#endif //JACOBI_JACOBIREPORT_H
