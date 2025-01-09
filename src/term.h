/*------------------------------------------------------------------------*/
/*! \file term.h
    \brief contains the class Term and further functions to
    manipulate terms

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#ifndef MULTILING_SRC_TERM_H_
#define MULTILING_SRC_TERM_H_
/*------------------------------------------------------------------------*/
#include <list>
#include <vector>
#include <algorithm>

#include "variable.h"
#include "signal_statistics.h"
/*------------------------------------------------------------------------*/

/** \class Term
    This class is used to represent terms in a polynomial.
    Terms are represented as ordered linked lists of variables.
*/

class Term {
  // / head variable
  Var * variable;

  // / tail in linked list
  Term * rest;

  // / reference counter
  uint64_t ref;

  // / hash value
  const uint64_t hash;

  // / hash collision chain link
  Term * next;

  size_t deg;

 public:
  /** Constructor

      @param _v Var*
      @param _r Term*
      @param _hash uint64_t
      @param _n Term*
  */
  Term(Var * _v, Term * _r, uint64_t _hash, Term * _n);

  /** Getter for member variable

      @return Var*
  */
  Var * get_var() const {return variable;}

  /** Getter for level of variable

      @return integer
  */
  int get_var_level() const {return variable->get_level();}

  /** Getter for num of variable

      @return integer
  */
  int get_var_num() const {return variable->get_num();}

  /** Getter for name of variable

      @return char*
  */
  const char * get_var_name() const {return variable->get_name();}

  /** Getter for member rest

      @return Term*
  */
  Term * get_rest() const {return rest;}

  /** Getter for member hash

      @return uint64_t
  */
  uint64_t get_hash() const {return hash;}

  /** Getter for member next

      @return Term*
  */
  Term * get_next() const {return next;}

  size_t degree() const {return deg;}

  /** Setter for member term

      @param t Term*
  */
  void set_next(Term * t) {next = t;}

  /** Getter for member ref

      @return uint64_t
  */
  uint64_t get_ref() const {return ref;}

  /** Increases ref

      @return uint64_t
  */
  uint64_t inc_ref() {return ++ref;}

  /** Decreases ref

      @return uint64_t
  */
  uint64_t dec_ref() {return --ref;}

  /**
      Copy routine

      @return A copy of the current term
  */
  Term * copy();


  /**
      Printing routine

      @param file Output file
  */
  void print(FILE * file) const;

  std::string print_to_string() const;

  /**
      Returns the number of variables in a term

      @return unsigned integer
  */
 // unsigned size() const;

  /**
      Compares this term to term t using the levels of the variables, based on lex ordering

      @param t Term*

      @return +1 if this > t
              -1 if this < t
              0  if this = t
  */
  int cmp(const Term *t) const;

  bool equal_up_to_duality(const Term *t) const;
 /**
      Compares this term to term t using the levels of the variables, based on dlex ordering

      @param t Term*

      @return +1 if this > t
              -1 if this < t
              0  if this = t
  */
  int cmp_dlex(const Term *t) const;

  /**
      Checks whether v is contained in Term

      @param v Var*

      @return true if v is contained in term
  */
  bool contains(Var * v) const;

  bool contains_input() const;

  Var * first_input() const;

  Var * extract_first_dual_var() const;

  Var * last_input() const;

  bool contains_t(const Term * t) const;
};

/*------------------------------------------------------------------------*/
// We organize terms in a hash table that is dynamically enlarged.
// Every time a new term is defined, we compute a hash value and insert
// the term. Terms are counted using a reference counter, which is incremented
// and decremented depending how often the term occurs in polynomials.

/**
    Compute hash_values

    @param variable Variable*
    @param rest Term*

    @return computed hash value for the term(variable, rest)
*/
uint64_t compute_hash_term(Var * variable, const Term * rest);

/**
    Enlarges the hash table
*/
void enlarge_terms();

/**
    Builds a term, where variable is added at the front of rest

    @param variable Variable*
    @param rest     Term*

    @return Term*
*/
Term * new_term(Var * variable, Term * rest);

Term * new_quadratic_term(Var * v1, Var* v2);



/**
    removes elements from vstack
*/
void clear_vstack();

Term * sort_and_build_term_from_vector(std::vector<Var*> v);
Term * gen_dual_term_from_vector(std::vector<Term*> t);


Var * first_dual_different_var(Term * t0, Term* t1);

/**
    Decrements the reference count of a term, and actually deletes a
    term if its reference count goes to zero.  In this case it also
    removes it from the hash table and applies the same procedure to the
    suffix 'rest'.

    @param t Term*
*/
void deallocate_term(Term * t);


/**
    Deallocates the hash table "term_table"
*/
void deallocate_terms();

/*------------------------------------------------------------------------*/
// Terms are generated using a stack "vstack"


void add_to_vstack(Var* v);

Term * divide_by_var(Term *t, Var*v);

Term * divide_by_term(Term *t, Term * t1);


/**
    Generates a term from the variable stack

    @return Term* generated from the variable stack
*/
Term * build_term_from_stack(bool sort = 0);


Term * multiply_term(Term * t1, Term * t2);
Term * multiply_term_by_var(Term * t1, Var * v);

  Term * reorder_term_lex(Term * t1);


Term * remainder(const Term * t, Var * v);

Term * remainder_t(const Term * t, const Term * t1);

#endif  // MULTILING_SRC_TERM_H_
