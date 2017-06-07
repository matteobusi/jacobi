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
        os << report.nWorkers << ", " << report.nIterations << ", " << report.time << ", " << report.error;
        return os;
    }

public:
    int nIterations;
    double error;
    int nWorkers;
    double time;
};


#endif //JACOBI_JACOBIREPORT_H
