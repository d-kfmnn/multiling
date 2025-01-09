/*------------------------------------------------------------------------*/
/*! \file polynomial_solver.cpp
    \brief contains the polynomial solving routine

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#include "polynomial_solver.h"
/*------------------------------------------------------------------------*/
// Global variable
bool gen_witness = 1;
/*------------------------------------------------------------------------*/
// ERROR CODES:
static int err_rem_poly =42; // remainder poly no witness
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
Polynomial * mult_spec_poly(){
  mpz_t coeff;
  mpz_init(coeff);

  // outputs
  for (int i = NN-1; i >= 0; i--) {
    Var * v = gates[i+M-1]->get_var();

    mpz_pow_ui(coeff, base, i);
    mpz_neg(coeff, coeff);
    if (i == static_cast<int>(NN-1) && signed_mult) mpz_neg(coeff, coeff);

    Term * t1 = new_term(v, 0);
    Monomial * m1 = new Monomial(coeff, t1);
    push_mstack(m1);
  }

  // partial products
  for (int i = NN/2-1; i >= 0; i--) {
    Var * a = gates[a0+i*ainc]->get_var();

    for (int j = NN/2-1; j >= 0; j--) {
      Var * b = gates[b0+j*binc]->get_var();
      mpz_pow_ui(coeff, base, i+j);
      if (i == static_cast<int>(NN/2-1) && signed_mult) mpz_neg(coeff, coeff);
      if (j == static_cast<int>(NN/2-1) && signed_mult) mpz_neg(coeff, coeff);
      add_to_vstack(b);
      add_to_vstack(a);
      Term * t1 = build_term_from_stack();
      Monomial * m1 = new Monomial(coeff, t1);
      push_mstack(m1);
    }
  }
  mpz_clear(coeff);

  Polynomial * p = build_poly();
  return p;
}
/*------------------------------------------------------------------------*/
Polynomial * miter_spec_poly(){
  mpz_t coeff;
  mpz_init(coeff);

  // outputs
  for (int i = MM-1; i >= 0; i--) {
    Var * v = gates[i+M-1]->get_var();

    mpz_pow_ui(coeff, base, i);
    mpz_neg(coeff, coeff);
    if (i == static_cast<int>(MM-1) && signed_mult) mpz_neg(coeff, coeff);

    Term * t1 = new_term(v, 0);
    Monomial * m1 = new Monomial(coeff, t1);
    push_mstack(m1);
  }

  Monomial * m1 = new Monomial(one, 0);
  push_mstack(m1);
  mpz_clear(coeff);

  Polynomial * p = build_poly();
  return p;
}
/*------------------------------------------------------------------------*/

bool verify(const char * inp_f, Polynomial * spec) {




  //init_slices();


  Polynomial * rem = linearize_spec(spec);
  

  Polynomial * rest = reduce(rem);
  rem = rest;
 
 
 

  
  bool res;
  if (rem && !rem->is_constant_zero_poly())  {
    if (!check_inputs_only(rem)){
      msg("REMAINDER IS");
      fputs("[mltlng] ", stdout);
      rem->print(stdout);
      msg("");
      die(err_rem_poly, "sorting error - remainder polynomial contains non-inputs");

    }

    res = 0;
    msg("INCORRECT MULTIPLIER");
    msg("");

    if (inp_f && gen_witness) {
      msg("REMAINDER IS");
      fputs("[mltlng] ", stdout);
      rem->print(stdout);
      msg("");
      generate_witness(rem, inp_f);
    }
  } else {
    res = 1;
    msg("");
    msg("CORRECT MULTIPLIER");

  }
  delete(rem);



  return res;
}
