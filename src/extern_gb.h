/*------------------------------------------------------------------------*/
/*! \file extern_gb.h
    \brief contains functions used in the polynomial solver

  This file contains all functions used for preprocessing the Gröbner basis
  and for reducing the specification by the slices.

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#ifndef MULTILING_SRC_GB_H_
#define MULTILING_SRC_GB_H_
/*------------------------------------------------------------------------*/
#include <string.h>
#include <vector>
#include <set>
#include <sstream>

#include "polyparser.h"



/*------------------------------------------------------------------------*/

bool   linearize_via_gb(Gate *g, int depth, bool pre, mpz_t coeff, bool full);
Polynomial *flip_var_in_poly(const Polynomial *p1, Var *v);
Polynomial *unflip_poly(Polynomial *p);

#endif  // MULTILING_SRC_GB_H_
