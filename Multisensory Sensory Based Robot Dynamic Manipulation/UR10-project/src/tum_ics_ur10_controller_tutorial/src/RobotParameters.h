#ifndef ROBOTPARAMETES_H
#define ROBOTPARAMETES_H

#include<Eigen/Core>

using namespace Eigen;

namespace param {
    // Lengths
    double L1 =      0.1280;
    double L2 =      0.1639;
    double L3 =      0.6127;
    double L4 =      0.0922;
    double L5 =      0.5716;
    double L6 =      0.1000;
    double L7 =      0.1500;
    double L8 =      0.3064;
    double L9 =      0.0819;
    double L10 =     0.3437;
    double L11 =     0.1157;
    double L12 =     0.0922;

    // Masses
    double m1 =     7.1000;
    double m2 =    12.7000;
    double m3 =     4.2700;
    double m4 =     2.0000;
    double m5 =     2.0000;
    double m6 =     0.3650;

    // Inertias

    // q1
    double I111 =     0.0400;
    double I112 =     0.0200;
    double I113 =     0.0200;
    double I122 =     0.0400;
    double I123 =     0.0200;
    double I133 =     0.0400;

    // q2
    double I211 =     0.0400;
    double I212 =     0.0200;
    double I213 =     0.0200;
    double I222 =     0.0400;
    double I223 =     0.0200;using namespace Eigen;
    double I233 =     0.0400;

    // q3
    double I311 =     0.0400;
    double I312 =     0.0200;
    double I313 =     0.0200;
    double I322 =     0.0400;
    double I323 =     0.0200;
    double I333 =     0.0400;

    // q4
    double I411 =     0.0400;
    double I412 =     0.0200;
    double I413 =     0.0200;
    double I422 =     0.0400;
    double I423 =     0.0200;
    double I433 =     0.0400;

    // q5
    double I511 =     0.0400;
    double I512 =     0.0200;
    double I513 =     0.0200;
    double I522 =     0.0400;
    double I523 =     0.0200;
    double I533 =     0.0400;

    // q6
    double I611 =     0.0400;
    double I612 =     0.0200;
    double I613 =     0.0200;
    double I622 =     0.0400;
    double I623 =     0.0200;
    double I633 =     0.0400;

    // Gravity
    double gx = 0;
    double gy = 0;
    double gz = -1;

    double g = 9.81;

    int count = 0;



}



#endif // ROBOTPARAMETES_H