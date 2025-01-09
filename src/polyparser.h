/*------------------------------------------------------------------------*/
/*! \file parser.h
    \brief core functions for parsing

  Part of Pacheck 2.0 : PAC proof checker.
  Copyright(C) 2020 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#ifndef PACHECK2_SRC_PARSER_H_
#define PACHECK2_SRC_PARSER_H_
/*------------------------------------------------------------------------*/
#include <stdarg.h>
#include <vector>

#include "polynomial.h"
#include "gate.h"


/*------------------------------------------------------------------------*/




Polynomial * parse_specification_polynomial(const char * file_name);
//void parse_linear_gb(const char * file_name);
/*------------------------------------------------------------------------*/
#endif  // PACHECK2_SRC_PARSER_H_
