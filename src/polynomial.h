/*------------------------------------------------------------------------*/
/*! \file polynomial.h
    \brief contains arithmetic operations for polynomials

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#ifndef MULTILING_SRC_POLYNOMIAL_H_
#define MULTILING_SRC_POLYNOMIAL_H_
/*------------------------------------------------------------------------*/
#include <list>
#include <deque>
#include <cstring>

#include "monomial.h"
/*------------------------------------------------------------------------*/

/** \class Polynomial
    This class is used to polynomials.
*/
class Polynomial {
 
  Monomial ** mon;  

  size_t num_mon = 0;

  int deg = 0;

 public:
  /** Constructor */
  Polynomial();

  Polynomial (Monomial ** m, size_t len, int d);


  Monomial * get_mon(size_t i) const {
     if(i < num_mon) return mon[i];
     else return 0;
   }

   size_t size() const { return num_mon;}

  size_t degree() const {return deg;}


   

  /** Returns the leading term

      @return Term*
  */
  Term * get_lt() const {return mon[0]->get_term();}

  Monomial * get_lm() const {return mon[0];}

  Term * get_tail_term() const {return mon[1]->get_term();}
  Term * get_largest_term() const;
  Monomial * get_largest_mon() const;

  Var* contains_dual_var() const;

  /**
      Copy routine

      @return A hard copy of the current polynomial
  */
  Polynomial * copy() const;

Polynomial * reorder_lex() const;
  Polynomial * reorder_dlex() const;


  /**
      Printing routine

      @param file Output file
      @param end if true we print trailing ";"
  */
  void print(FILE * file, bool end = 1) const;

   std::string print_to_string(bool end = 1) const;

  /** Destructor */
  ~Polynomial();


  /**
      Returns whether the polynomial is the constant zero polynomial

      @return bool
  */
  bool is_constant_zero_poly() const;

  /**
      Returns whether the polynomial is the constant one polynomial

      @return bool
  */
  bool is_constant_one_poly() const;

  /**
      Returns the size of the smallest term

      @return integer
  */
  int min_term_size() const;
};
/*------------------------------------------------------------------------*/
// Polynomials are generated using a sorted array "mstack"

/**
    Enlarges the allocated size of mstack
*/
void enlarge_mstack();

/**
    Sets the number of elements in the stack to 0
*/
void clear_mstack();

/**
    Deallocates mstack
*/
void deallocate_mstack();

/**
    Checks whether mstack is empty

    @return true if mstack is empty
*/
bool mstack_is_empty();



/**
    Pushes a monomial to the end of the stack

    @param t monomial to be added to the mstack
*/
void push_mstack_end(Monomial *t);

/**
    Pushes a monomial to the stack such that mstack remains sorted

    @param t monomial to be added to the mstack
*/
void push_mstack(Monomial *t);



void push_mstack_dlex(Monomial *m);


/**
    Generates a polynomial from mstack and clears mstack

    @return Polynomial*
*/
Polynomial * build_poly();


/*------------------------------------------------------------------------*/

/**
    Add two polynomials p1+p2

    @param p1 Polynomial*
    @param p2 Polynomial*

    @return Polynomial*, sum of p1+p2
*/
Polynomial * add_poly(const Polynomial *p1, const Polynomial *p2);


/**
    Add two polynomials p1+p2

    @param p1 Polynomial*
    @param p2 Polynomial*

    @return Polynomial*, sum of p1-p2
*/
Polynomial * sub_poly(const Polynomial *p1, const Polynomial *p2);

/**
    Multiplies two polynomials p1*p2

    @param p1 Polynomial*
    @param p2 Polynomial*

    @return Polynomial*, product of p1*p2
*/
Polynomial * multiply_poly(const Polynomial *p1, const Polynomial *p2);

/**
    Multiplies a polynomial p1 with a constant c

    @param p1: Polynomial*
    @param c:  mpz_t object

    @return Polynomial*, product of c*p1
*/
Polynomial * multiply_poly_with_constant(const Polynomial *p1, mpz_t c);

/**
    Multiplies a polynomial p1 with a constant c

    @param p1: Polynomial*
    @param c:  mpz_t object

    @return Polynomial*, product of c*p1
*/
Polynomial * multiply_poly_with_term(const Polynomial *p1, Term * t);
/**
    Returns the quotient of dividing a polynomial p1 by a term t

    @param p1 Polynomial*
    @param t  Term *

    @return Polynomial defining the quotient of p1/t
*/
Polynomial * divide_by_term(const Polynomial * p1, const Term * t);
Polynomial * build_one_minus_x_times_y_poly(Term * x, Term *y);
 
Polynomial * reduce_by_one_poly(Polynomial * p1, Polynomial * p2);
Polynomial * substitute_linear_poly(Polynomial * p1, Polynomial * p2);

bool equal_poly (Polynomial * p1, Polynomial * p2);

/*---------------------------------------------------------------------------*/
// / gmp for 1
extern mpz_t one;

// / gmp for -1
extern mpz_t minus_one;

// / gmp for 2
extern mpz_t base;

// / gmp for 2^NN
extern mpz_t mod_coeff;
/*---------------------------------------------------------------------------*/
/**
    Initializes all global gmp objects

    @param exp unsigned exponent for mod coeff
*/
void init_mpz(unsigned exp);

/**
    Clears all globally allocated gmp objects
*/
void clear_mpz();

#endif  // MULTILING_SRC_POLYNOMIAL_H_
