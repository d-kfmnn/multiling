/*------------------------------------------------------------------------*/
/*! \file monomial.cpp
    \brief contains the class Monomial and further functions to
    manipulate monomials

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#include "monomial.h"
/*------------------------------------------------------------------------*/

Monomial::Monomial(mpz_t _c, Term * _t): ref(1) {
  mpz_init_set(coeff, _c);
  if (mpz_sgn(_c) == 0) term = 0;
  else  term = _t;
}

/*------------------------------------------------------------------------*/

Monomial * Monomial::copy() {
  Monomial * m = this;
  if (!m) return 0;
  assert(ref > 0);

  ++ref;
  assert(ref);
  return this;
}

/*------------------------------------------------------------------------*/

Monomial::~Monomial() {
  assert(ref == 0);

  mpz_clear(coeff);
  deallocate_term(term);
}

/*------------------------------------------------------------------------*/

void Monomial::print(FILE * file, bool lm) const {
  int sign = mpz_sgn(coeff);
  if (!sign) return;
  else if (!lm && sign > 0) fputc_unlocked('+', file);

  if (term) {
    if (mpz_cmp_si(coeff, -1) == 0) { fputc_unlocked('-', file);
    } else if (mpz_cmp_si(coeff, 1) != 0) {
      mpz_out_str(file, 10, coeff);
      fputc_unlocked('*', file);
    } term->print(file);
  }
  else  mpz_out_str(file, 10, coeff);
}

/*------------------------------------------------------------------------*/
static std::string mpz_to_string(const mpz_t num, int base = 10) {
    // Convert mpz_t to a C-style string
    char* cstr = mpz_get_str(nullptr, base, num);

    // Create a std::string from the C-style string
    std::string result(cstr);

    // Free the allocated memory for the C-style string
    // mpz_get_str allocates memory that needs to be freed with free()
    free(cstr);

    return result;
}


std::string Monomial::print_to_string(bool lm) const{
  std::string res;
  int sign = mpz_sgn(coeff);
  if (!sign) return res;
  else if (!lm && sign > 0) res += "+";

  if (term) {
    if (mpz_cmp_si(coeff, -1) == 0) { res += "-";
    } else if (mpz_cmp_si(coeff, 1) != 0) {
      res += mpz_to_string(coeff);
      res += "*";
    } res+= term->print_to_string();
  }
  else  res+= mpz_to_string(coeff);

  return res;
}
/*------------------------------------------------------------------------*/

Monomial * multiply_monomial(const Monomial * m1, const Monomial *m2) {
  assert(m1);
  assert(m2);

  mpz_t coeff;
  mpz_init(coeff);
  mpz_mul(coeff, m1->coeff, m2->coeff);

  Term * t;
  if (m1->get_term() && m2->get_term())
    t = multiply_term(m1->get_term(), m2->get_term());
  else if (m2->get_term()) t = m2->get_term_copy();
  else if (m1->get_term()) t = m1->get_term_copy();
  else
    t = 0;

  Monomial * mon = new Monomial(coeff, t);

  mpz_clear(coeff);
  return mon;
}

/*------------------------------------------------------------------------*/

void deallocate_monomial(Monomial * m) {
  assert(m->get_ref() > 0);
  if (m->dec_ref() > 0) return;

  delete(m);
}
