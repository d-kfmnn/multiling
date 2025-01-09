/*------------------------------------------------------------------------*/
/*! \file elimination.h
    \brief contains functions used in the polynomial solver

  This file contains all functions used for preprocessing the Gröbner basis
  and for reducing the specification by the slices.

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#ifndef MULTILING_SRC_ELIMINATION_H_
#define MULTILING_SRC_ELIMINATION_H_
/*------------------------------------------------------------------------*/
#include <string.h>
#include <vector>

#include "dual_rewriting.h"
#include "extern_gb.h"
/*------------------------------------------------------------------------*/
// / 1 for pac
// / 2 for hybrid
// / 3 for nss
extern int proof;


/*------------------------------------------------------------------------*/

void preprocessing_elimination();

Polynomial * linearize_spec(Polynomial * spec);

Polynomial * reduce(Polynomial * spec);


#endif  // MULTILING_SRC_ELIMINATION_H_
