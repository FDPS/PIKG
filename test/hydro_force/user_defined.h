#pragma once
#include <pikg_vector.hpp>

class Force_hydro {
public:
    PIKG::F32vec acc;
    PIKG::F32 eng_dot;
    PIKG::F32 dt;
};

class EP_hydro {
public:
    PIKG::F64vec pos;
    PIKG::F64    mass;
    PIKG::F32vec vel;
    PIKG::F32    eng;
    PIKG::F32    h;
    PIKG::F32    dens;
    PIKG::F32    pres;
    PIKG::F32    gradh;
    PIKG::F32    snds;
    PIKG::F32    BalSW;
    PIKG::F32    alpha;
};

