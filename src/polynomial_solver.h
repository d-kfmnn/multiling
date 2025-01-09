/*------------------------------------------------------------------------*/
/*! \file polynomial_solver.h
    \brief contains the polynomial solving routine

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#ifndef MULTILING_SRC_POLYNOMIAL_SOLVER_H_
#define MULTILING_SRC_POLYNOMIAL_SOLVER_H_
/*------------------------------------------------------------------------*/
#include "elimination.h"
#include "witness.h"
/*------------------------------------------------------------------------*/
// / If final remainder is not equal to zero a counter example is generated and
// / printed to file <input_name>.wit, default is true, can be turned of
//  using command line input
extern bool gen_witness;
/*------------------------------------------------------------------------*/

/**
    Calls the preprocessing, slicing, reduction routines
    If certify is true, files need to be provided to store the proof.

    @param inp_f name of input file

    @returns boolean whether circuit is correct (1) or incorrect (0)
*/
bool verify(const char * inp_f = 0, Polynomial * spec = 0);

Polynomial * mult_spec_poly();
Polynomial * miter_spec_poly();

#endif  // MULTILING_SRC_POLYNOMIAL_SOLVER_H_
