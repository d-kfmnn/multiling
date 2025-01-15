/*------------------------------------------------------------------------*/
/*! \file gate.h
    \brief contains the class Gate and further functions to
    organize the gate structure, such as initializing the gate constraints

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#ifndef MULTILING_SRC_GATE_H_
#define MULTILING_SRC_GATE_H_
/*------------------------------------------------------------------------*/
#include <list>
#include <map>
#include <queue>
#include <string>

#include "aig.h"
#include "polynomial.h"
/*------------------------------------------------------------------------*/
// / set to true when a signed or unsigned multiplier is verified
extern bool signed_mult;
extern int add_var;
extern int max_dist;
/*------------------------------------------------------------------------*/

/** \class Gate
  Internal structure to represent the AIG graph.
*/
class Gate {
  // / Variable of the gate, as used in the polynomials
  Var * v;

  Var* replace_var = 0;

  // / True if gate is an input
  bool input;

  // / True if the gate is an output s_i
  bool output;

  // / True if the gate is an output in the aig
  bool aig_output = 0;

  // / True if gate is identified as a partial product
  bool partial_product = 0;


  // / True if gate is redinated during preprocessing
  bool red = 0;
  bool rew = 0;
  bool ext = 0;
  

  // / Polynomial implied by the aig gate
  Polynomial * gate_constraint = 0;
  Polynomial * expanded_constraint = 0;

   // / list of gates that create vanishing monomials
  std::list<Gate*> van_twins;
   std::list<Gate*> inverse_van_twins;

  std::list<Gate*> neutrals;

  std::list<Gate*> subset_of;


  // / list of gates that are parents
  std::list<Gate*> parents;

  // / list of gates that are children
  std::list<Gate*> children;

 public:
  /**
      Constructor
      Calls constructor of Var

      @param n_ value, corresponding to aiger value
      @param name_ string name of the variable
      @param level position in order of the variable
  */
  Gate(
    int n_, std::string name_, int level_, bool input_ = 0, bool output_ = 0);

  /**
      Getter for v

      @return member v
  */
  Var * get_var() const {return v;}
  Var * get_rep_var() const {return replace_var;}

  void set_rep_var(Var *v)  {replace_var = v;}

  /**
      Getter for number of v

      @return integer
  */
  int get_var_num() const {return v->get_num();}
   int get_var_level() const {return v->get_level();}

void set_var_level(int l){
    v->set_level(l);
    v->get_dual()->set_level(l+1);
}

  std::list<Gate*> get_children(){return children;}
  void set_children(std::list<Gate*> c){children =c;}
  void delete_children() {children = std::list<Gate*>();}

  Gate * children_front() const {return children.front();}
  Gate * children_back() const {return children.back();}


  std::list<Gate*> get_parents(){return parents;}
  Gate * parents_front() const {return parents.front();}
  Gate * parents_back() const {return parents.back();}

  /**
      Getter for name of v

      @return char*
  */
  const char * get_var_name()const {return v->get_name();}

  /**
      Getter for input

      @return member input
  */
  bool get_input()  const {return input;}

  void add_children_to_neutrals(Gate * g); 

  void update_gate_poly(Polynomial * p);


  /**
      Getter for output

      @return member output
  */
  bool get_output() const {return output;}


  /**
      Sets aig_output to true
  */
  void mark_aig_output() {aig_output = 1;}

  /**
      Getter for  aig_output

      @return member  aig_output
  */
  bool get_aig_output() const {return aig_output;}



  /**
      Getter for partial_product

      @return member partial_product
  */
  bool get_pp() const {return partial_product;}

  /**
      Sets partial_product to true
  */
  void mark_pp() {partial_product = 1;}

  bool some_children_are_inputs() const;

  /**
      Getter for red

      @return member red
  */
  bool get_red() const {return red;}
  bool get_ext() const {return ext;}
  void set_ext()  {ext = 1;}
  
  /**
      Sets red to true
  */
  void mark_red();
    /**
      Getter for red

      @return member red
  */
  bool get_rew() const {return rew;}

  /**
      Sets red to true
  */
  void mark_rew() {rew = 1;}

  /**
      Getter for gate_constraint

      @return member gate_constraint
  */
  Polynomial * get_gate_constraint() const;
 
 Polynomial * get_expanded_constraint() const {return expanded_constraint;}
  /**
      Setter for gate_constraint

      @param p Polynomial *
  */
  void set_gate_constraint(Polynomial * p) {gate_constraint = p;}
void set_expanded_constraint(Polynomial * p) {expanded_constraint = p;}
  /**
      Prints the gate constraint

      @param file output file
  */
  void print_gate_constraint(FILE * file) const { gate_constraint->print(file);}


   /**
      Determines whether n is contained in the parents of this gate

      @param n Gate*

      @return True whether n is parent of this gate
  */
  bool is_van_twin(const Gate *n) const;
  bool is_inverse_van_twin(const Gate *n) const;

  bool is_neutral(const Gate *n) const;

    bool is_in_subsets(const Gate *n) const;

  /**
      Getter for size of parents

      @return size_t
  */
  size_t van_twins_size() const {return van_twins.size();}

  size_t neutrals_size() const {return neutrals.size();}

  size_t subset_of_size() const {return subset_of.size();}



  /**
      Appends gate n to the parents

      @param n Gate*
  */
  void neutrals_push_back(Gate * n) {neutrals.push_back(n);}

    void subset_of_push_back(Gate * n) {subset_of.push_back(n);}

 void inverse_van_twins_push_back(Gate * n) {inverse_van_twins.push_back(n);}
 void van_twins_push_back(Gate * n) {van_twins.push_back(n);}
  /**
      Determines whether n is contained in the parents of this gate

      @param n Gate*

      @return True whether n is parent of this gate
  */
  bool is_in_parents(const Gate *n) const;

  void update_children_parents_reduce(Gate * n);
  void update_children_parents_elim(Gate * n);

  void update_children_parents_tail_sub(Gate * n);
  /**
      Getter for begin of parents

      @return std::list<Gate*>::const_iterator
  */
  std::list<Gate*>::const_iterator parents_begin() const {
    return parents.begin();
  }

  std::list<Gate*>::const_iterator van_twins_begin() const {
    return van_twins.begin();
  }

    std::list<Gate*>::const_iterator inverse_van_twins_begin() const {
    return inverse_van_twins.begin();
  }
 

    std::list<Gate*>::const_iterator subset_begin() const {
    return subset_of.begin();
  }
  /**
      Getter for end of parents

      @return std::list<Gate*>::const_iterator
  */
  std::list<Gate*>::const_iterator parents_end()   const {
    return parents.end();
  }

    std::list<Gate*>::const_iterator van_twins_end() const {
    return van_twins.end();
  }
     std::list<Gate*>::const_iterator inverse_van_twins_end() const {
    return inverse_van_twins.end();
  }
 
    std::list<Gate*>::const_iterator subset_end()   const {
    return subset_of.end();
  }
   std::list<Gate*>::const_iterator children_begin() const {
    return children.begin();
  }

  /**
      Getter for end of parents

      @return std::list<Gate*>::const_iterator
  */
  std::list<Gate*>::const_iterator children_end()   const {
    return children.end();
  }

  /**
      Getter for size of parents

      @return size_t
  */
  size_t parents_size() const {return parents.size();}

  /**
      Appends gate n to the parents

      @param n Gate*
  */
  void parents_push_back(Gate * n) {parents.push_back(n);}

  /**
      Removes gate n from the parents

      @param n Gate*
  */
  void parents_remove(Gate * n) {parents.remove(n);}


  /**
      Determines whether n is contained in the children of this gate

      @param n Gate*

      @return True whether n is child of this gate
  */
  bool is_child(const Gate *n) const;


  bool equal_children(const Gate *n) const;

  /**
      Getter for size of children

      @return size_t
  */
  size_t children_size() const {return children.size();}


  /**
      Appends gate n to the children

      @param n Gate*
  */
  void children_push_back(Gate * n) { if(n) children.push_back(n);}

  /**
      Removes gate n from the children

      @param n Gate*
  */
  void children_remove(Gate * n) {children.remove(n);}

  /**
      Destructor
  */
  ~Gate();
};

/*------------------------------------------------------------------------*/
// / Gate ** where all gates are stored
extern Gate ** gates;

// / Counts the number of gates
extern unsigned num_gates;


/*------------------------------------------------------------------------*/

/**
    Returns the gate with aiger value 'lit'

    @param lit unsigned integer

    @returns Gate*
*/
Gate * gate(int lit);

/**
    Generates the constraint -v-v_+1

    @param v Var*
*/

Polynomial * gen_dual_constraint(Var * v);
 Polynomial * gen_gate_constraint(unsigned i) ;

Gate * search_for_parent(Term * t, Gate * exclude = 0);

Gate * search_for_parent_dual(Term * t);

Gate * neutral_child(Gate*parent, Gate*g);


void init_gates();
void enlarge_gates(int added_size);
Term * extend_var_gates(Term *t);
void adjust_level_of_extended_gates();
/**
    Deletes Gate** gates by calling destructur of Gate
*/
void delete_gates();

std::list<Gate*>  get_var_of_poly(Gate * g);

// / vector-list Gate* matrix to store slices
extern std::vector<Gate*> slices;
/**
    Allocates the memory for the slices
*/
void init_slices();

/**
    Prints the polynomials in the slices to stdout
*/
void print_slices();

void resort_slices_lex();

#endif  // MULTILING_SRC_GATE_H_
