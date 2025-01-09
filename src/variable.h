/*------------------------------------------------------------------------*/
/*! \file variable.h
    \brief contains the class Var

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#ifndef MULTILING_SRC_VARIABLE_H_
#define MULTILING_SRC_VARIABLE_H_
/*------------------------------------------------------------------------*/
#include <string>
#include <cstdlib>

#include "hash_val.h"
/*------------------------------------------------------------------------*/

/** \class Var
    represents a variable is assigned to a gate(see <gate.h>) and is
    used to represent variables in terms, see <term.h>
*/

class Var {
    // / name of variable
    const std::string name;

    // / Hash value of variables, used for storing terms
    int hash;

    // / Increasing value that indicates the order of the variable
    int level;
    int distance = 0;

    // / corresponding value used to relate AIG gates to Gate class
    const int num;

    bool inp;

    // / links dual variable
    Var * dual;

    bool d = 0;

    


 public:
  /** Constructor

     @param name_ name
     @param level_ level
     @param num_ num, default is 0

  */
  Var(std::string name_, int level_, int num_ = 0, bool inp_ = 0, bool dual_ = 0):
     name(name_),  level(level_), num(num_), inp(inp_), d(dual_) {
       hash = hash_string(name_);
     }

   /** Getter for member name, and converts string to char*

       @return const char *
   */
  const char * get_name() const {return name.c_str();}

  /** Getter for member hash

      @return integer
  */
  int get_hash() const {return hash;}


  /** Getter for member d

      @return bool
  */
  int is_dual() const {return d;}

   int is_input() const {return inp;}

  /** Getter for member level

      @return integer
  */
  int get_level() const {return level;}

  void set_level(int l) {level = l;}

    int get_dist() const {return distance;}

  void set_dist(int l) {distance = l;}

  /** Getter for member num

      @return integer
  */
  int get_num() const {return num;}


  /** Getter for member dual

    @return integer
  */
  Var * get_dual() const {return dual;}

  /** Setter for dual variable
      @param d Variable *
  */
  void set_dual_var(Var *d) { dual = d;}

  

};


#endif  // MULTILING_SRC_VARIABLE_H_