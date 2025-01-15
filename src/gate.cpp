/*------------------------------------------------------------------------*/
/*! \file gate.cpp
    \brief contains the class Gate and further functions to
    organize the gate structure, such as initializing the gate constraints

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#include <string>
#include <list>

#include "gate.h"
/*------------------------------------------------------------------------*/
// Global variables
bool signed_mult = 0;
int add_var = 0;
int max_dist = 0;
/*------------------------------------------------------------------------*/
// ERROR CODES:
static int err_allocate       = 91; // failed to allocate gates
std::vector<Gate*> slices;
/*------------------------------------------------------------------------*/
static Polynomial * positive_poly(Var * v) {
  Term * t = new_term(v, 0);
  Monomial * m = new Monomial(one, t);
  push_mstack_end(m);

  return build_poly();
}

/*------------------------------------------------------------------------*/
std::list<Gate*>  get_var_of_poly(Gate * g){
  std::list<Gate*> res;
  Polynomial * p = g->get_gate_constraint();

  for (size_t i = 1 ; i < p->size(); i++) {
    Monomial * m = p->get_mon(i);
    Term * t = m->get_term();
    while (t) {
      Gate * tmp = gate(t->get_var_num());
      auto it = std::find (res.begin(), res.end(), tmp);
      if (it == res.end()){
        res.push_back(tmp);
      }
      t = t->get_rest();
    }
  }
  

  return res;
}

/*------------------------------------------------------------------------*/

static Polynomial * negative_poly(Var * v) {
  Term * t1 = new_term(v, 0);
  Monomial * m1 = new Monomial(minus_one, t1);
  push_mstack(m1);

  Monomial * m2 = new Monomial(one, 0);
  push_mstack(m2);
  
  return build_poly();
}

/*------------------------------------------------------------------------*/

static Polynomial * get_node_constraint(Gate * g, unsigned sign){
  if(g) {
    Var * v1 = g->get_var();
    if (sign) return positive_poly(v1->get_dual());
    else return positive_poly(v1);
  } else {
    if(sign) {
      push_mstack_end(new Monomial(one, 0));
      return build_poly();
    } else return 0;
  }
}

/*------------------------------------------------------------------------*/
Polynomial * gen_gate_constraint(unsigned i) {
  assert(i >= NN && i < M + MM - 1);
  Polynomial * p = 0;
  // gate constraint
  if (i < M-1) {
    Gate * n = gates[i];
    assert(!n->get_input());

    aiger_and * and1 = is_model_and(n->get_var_num());
    assert(and1);

    unsigned l = and1->rhs0, r = and1->rhs1;
    Gate * l_gate = gate(l), *r_gate = gate(r);

    Var * v = n->get_var();    
    Term * t1 = new_term(v, 0);
    Monomial * m1 = new Monomial(minus_one, t1);
    
    push_mstack_end(m1);
    Polynomial * p_h = build_poly();
    Polynomial * p1 = get_node_constraint(l_gate, aiger_sign(and1->rhs0));
    Polynomial * p2 = get_node_constraint(r_gate, aiger_sign(and1->rhs1));
    Polynomial * p_tl = multiply_poly(p1, p2);


    p = add_poly(p_h, p_tl);
    delete(p_h);
    delete(p_tl);

    delete(p1);
    delete(p2);
  } else {  // output
    Gate * n = gates[i];
    assert(n->get_output());

    unsigned lit = slit(i-M+1);
    Var * v = n->get_var();
    Term * t1 = new_term(v, 0);
    Monomial * m1 = new Monomial(minus_one, t1);
    push_mstack_end(m1);
    Polynomial * p_h = build_poly();

    Polynomial * p_tl;
     if(lit == 1){
       push_mstack_end(new Monomial(one, 0));
       p_tl = build_poly();
      }
      else if (aiger_sign(slit(i-M+1))){
        p_tl =  negative_poly(gate(slit(i-M+1))->get_var());
      } else {
        p_tl = positive_poly(gate(slit(i-M+1))->get_var());
      }

      p = add_poly(p_h, p_tl);
      delete(p_h);
      delete(p_tl);


    
  }

  return p;
}

//-------------------------------------------------------------
static void init_gate_constraint(unsigned i) {
  assert(i >= NN && i < M + MM - 1);
  Polynomial * p = gen_gate_constraint(i);
  Gate * n = gates[i];
  n->set_gate_constraint(p);  
}
//-------------------------------------------------------------

Gate::Gate(int n_, std::string name_, int level_, bool input_, bool output_):
  v(new Var(name_, level_, n_,input_)), input(input_), output(output_)  {
    std::string negname = name_;
    //negname.insert(1,1,'_');
    negname.insert(0,"(1-");
    negname.insert(negname.size(),")");
    Var * d = new Var (negname, level_+1, n_, input_, 1);
    v->set_dual_var(d);
    d->set_dual_var(v);
}

void Gate::mark_red() {
  if(red) return;
  red = 1; 
  if(verbose > 3) msg("eliminated %s", get_var_name());
  if (gate_constraint) {
    delete(gate_constraint);
    gate_constraint = 0;
  }
  if (expanded_constraint) { 
    delete(expanded_constraint);
    expanded_constraint = 0;
  }
}

void Gate::update_gate_poly(Polynomial * p){
      if(gate_constraint) delete(this->gate_constraint);
    
      this->set_gate_constraint(p);
  

      auto g_orig_children = this->get_children();
      for(auto& gc : g_orig_children){
        gc->parents_remove(this);
      }

      auto g_update_children = get_var_of_poly(this);
      this->set_children(g_update_children);
     

      for(auto& gp : g_update_children){ 
        gp->parents_push_back(this);
      }

      for(auto& gc : g_orig_children){
        if(gc->parents_size() == 0) gc->mark_red();
      }
     

}


/*------------------------------------------------------------------------*/

Gate::~Gate() {
  delete(v->get_dual());
  delete(v);
  if (gate_constraint) delete(gate_constraint);
  if (expanded_constraint) delete(expanded_constraint);
}

/*------------------------------------------------------------------------*/
bool Gate::is_neutral(const Gate * n) const {
  for (auto it=neutrals.begin(); it != neutrals.end(); ++it) {
    Gate * g = *it;
    if (g == n) return 1;
  }
  return 0;
}
/*------------------------------------------------------------------------*/
bool Gate::is_in_subsets(const Gate * n) const {
  for (auto it=subset_of.begin(); it != subset_of.end(); ++it) {
    Gate * g = *it;
    if (g == n) return 1;
  }
  return 0;
}

/*------------------------------------------------------------------------*/
void Gate::add_children_to_neutrals(Gate * g) {
  Gate * tmp = this;
  for (auto it=g->children.begin(); it != g->children.end(); ++it) {
    tmp->neutrals_push_back(*it);
  }
}

/*------------------------------------------------------------------------*/
bool Gate::is_van_twin(const Gate * n) const {
  for (auto it=van_twins.begin(); it != van_twins.end(); ++it) {
    Gate * g = *it;
    if (g == n) return 1;
  }
  return 0;
}

/*------------------------------------------------------------------------*/
bool Gate::is_inverse_van_twin(const Gate * n) const {
  for (auto it=inverse_van_twins.begin(); it != inverse_van_twins.end(); ++it) {
    Gate * g = *it;
    if (g == n) return 1;
  }
  return 0;
}

/*------------------------------------------------------------------------*/

bool Gate::is_in_parents(const Gate * n) const {
  for (auto it=parents.begin(); it != parents.end(); ++it) {
    Gate * parents = *it;
    if (parents == n) return 1;
  }
  return 0;
}
/*------------------------------------------------------------------------*/

bool Gate::equal_children(const Gate * n) const {
  if(children_size() != n->children_size()) return 0;

  for (auto it = children.begin(); it != children.end(); ++it) {
    Gate * child = *it;
    if (!n->is_child(child)) return 0;
  }
  return 1;
}



bool Gate::is_child(const Gate * n) const {
  for (auto it = children.begin(); it != children.end(); ++it) {
    Gate * child = *it;
    if (child == n) return 1;
  }
  return 0;
}

bool Gate::some_children_are_inputs()const {
    for (auto it = children.begin(); it != children.end(); ++it) {
    Gate * child = *it;
    if (child->get_input()) return 1;
  }
  return 0;
}

/*------------------------------------------------------------------------*/
Polynomial * Gate::get_gate_constraint() const {
  if (!gate_constraint) {
    // output aig are 0, -1, ...-NN+2
    if (output) init_gate_constraint(-1*get_var_num()+M-1);
    // gates are numbered 2,4,6,8,..
    else
      init_gate_constraint(get_var_num()/2-1);
  }
  return gate_constraint;
}

/*------------------------------------------------------------------------*/
void Gate::update_children_parents_reduce(Gate * n ){
  
  Gate * g = this;
  msg("update children and parents of %s and %s", g->get_var_name(), n->get_var_name());
  g->children_remove(n);

  n->parents_remove(g);
  for(auto it =n->children.begin(); it!=n->children.end(); it++){
    Gate * tmp = *it;
    g->children_push_back(tmp);
  }
   for(auto it =g->parents.begin(); it!=g->parents.end(); it++){
    Gate * tmp = *it;
    n->parents_push_back(tmp);
  }
}

/*------------------------------------------------------------------------*/
void Gate::update_children_parents_elim(Gate * n ){
  
  Gate * g = this;
  msg("update elim children and parents of %s and %s", g->get_var_name(), n->get_var_name());
  g->children_remove(n);
  n->parents_remove(g);

  for(auto& n_c : n->get_children()){
    g->children_push_back(n_c);
  }
}

void Gate::update_children_parents_tail_sub(Gate * n ){
  Gate * g = this;
  msg("update children and parents tail of %s and %s", g->get_var_name(), n->get_var_name());
  g->children_push_back(n);
  
  n->parents_push_back(g);
  for(auto it =n->children.begin(); it!=n->children.end(); it++){
    Gate * tmp = *it;
    g->children_remove(tmp);
  }
}

Gate ** gates;
unsigned num_gates;
unsigned size_gates;
unsigned extended_gates = 0;

/*------------------------------------------------------------------------*/

Gate * gate(int lit) {
  if(lit <= 0) return gates[M-lit-1];
  if(lit < 2) return 0;
  return gates[lit/2-1];
}

/*------------------------------------------------------------------------*/
 void enlarge_gates(int added_size) {
  uint64_t new_size_gates = size_gates+ added_size;
  Gate ** new_gates_table = new Gate*[new_size_gates];
  memcpy( new_gates_table, gates, size_gates * sizeof(Gate*) );
  delete[] gates;
  gates = new_gates_table;
  size_gates = new_size_gates;
}
/*------------------------------------------------------------------------*/
static void mark_aig_outputs() {
  for (unsigned i = 0; i < MM; i++) {
    unsigned lit = slit(i);
    if(lit < 2) continue;
    Gate * n = gate(lit);
    n->mark_aig_output();
  }
}
/*------------------------------------------------------------------------*/
static void allocate_gates() {
  unsigned aiger;
  num_gates = M + MM - 1;

  msg("allocating %i gates", num_gates);
  gates = new Gate*[num_gates + MM];
  size_gates=num_gates+MM;

  if (!gates) die(err_allocate, "failed to allocate gates");
  int level = 0;
  if(mult_inp){
  // inputs a
  for (unsigned i = a0; i <= al; i+=ainc) {
    aiger = 2*(i+1);
    assert(is_model_input(aiger));

    std::string name = "a" + std::to_string((i-a0)/ainc);
    gates[i] = new Gate(aiger, name, level+=2, 1);
  }

  // inputs b
  for (unsigned i = b0; i <= bl; i+=binc) {
    aiger = 2*(i+1);
    assert(is_model_input(aiger));

    std::string name = "b" + std::to_string((i-b0)/binc);
    gates[i] = new Gate(aiger, name, level+=2, 1);
  }
  } else {
    for (unsigned i = 0; i < NN; i++) {
      aiger = 2*(i+1);
      assert(is_model_input(aiger));

      std::string name = "i" + std::to_string((i-a0)/ainc);
      gates[i] = new Gate(aiger, name, level+=2, 1);
      if(verbose > 3) msg("allocated inp %s", gates[i]->get_var_name());
    }
  }

  // internal gates
  for (unsigned i = NN; i < M-1; i++) {
    aiger = 2*(i+1);
    assert(is_model_and(aiger));

    std::string name = "l" + std::to_string(aiger);
    gates[i] = new Gate(aiger, name, 0);
  }




  for (unsigned i = NN; i < M-1; i++) {
    Gate * n = gates[i];
    aiger_and * and1 = is_model_and(n->get_var_num());

    if (!and1) continue;
    unsigned l = and1->rhs0, r = and1->rhs1;
    Gate * l_gate = gate(l), *r_gate = gate(r);
    int dist_l = l_gate->get_var()->get_dist();
    int dist_r = r_gate->get_var()->get_dist();
    int dist_n = dist_l > dist_r ? dist_l+1 : dist_r+1;

    n->get_var()->set_dist(dist_n);
    n->get_var()->get_dual()->set_dist(dist_n);
    if(dist_n > max_dist) max_dist++;
  }
  msg("max dist is %i", max_dist);

  mark_aig_outputs();


  for(int dist_it = 1; dist_it <= max_dist; dist_it++){
    for (unsigned i = NN; i < M-1; i++) {
      Gate * n = gates[i];
      if(n->get_var()->get_dist() == dist_it){
        level+=2;
        n->set_var_level(level);
        slices.push_back(n);
      }
    }
  }

    // output s
  for (unsigned i = M-1; i < M-1+MM; i++) {
    aiger = i-M+1;
    std::string name = "s" + std::to_string(aiger);
    gates[i] = new Gate(M-i-1, name, 2*(i+1), 0, 1);
     if(verbose > 3) msg("allocated outp %s", gates[i]->get_var_name());
  }


}

/*------------------------------------------------------------------------*/
static void set_parents_and_children_of_extension_var(Term *t, Gate* g){
  while(t){
    Gate * tmp = gate(t->get_var_num());
    g->children_push_back(tmp);
    tmp->parents_push_back(g);
    t = t->get_rest();
  }
}

/*------------------------------------------------------------------------*/
Term * extend_var_gates(Term *t){
  if(num_gates == size_gates) die(2, "gates too small");

  int level = -2-2*extended_gates;
  std::string name = "t" + std::to_string(extended_gates++);
  gates[num_gates] = new Gate(M-num_gates-1, name, level, 0,0);
  num_gates++;
  Gate * g = gates[num_gates-1];

  Monomial ** intern_mstack = new Monomial*[2];
  
  intern_mstack[0] = new Monomial(minus_one, new_term(g->get_var(), 0));
  intern_mstack[1] = new Monomial(one, t);
  Polynomial * p = new Polynomial(intern_mstack, 2, t->degree());
  g->set_gate_constraint(p);
  g->set_ext();
  set_parents_and_children_of_extension_var(t, g);

  if( verbose >= 2) {
    msg("added extension poly to gates");
    p->print(stdout);
  }

  return p->get_lt();


  

}
/*------------------------------------------------------------------------*/
void adjust_level_of_extended_gates(){
  for(unsigned i = M-1+MM; i < num_gates; i++){
    Gate * g = gates[i];
    g->set_var_level(g->get_var_level()+2*NN+2);
  }

  for (unsigned i = 0; i < NN; i++) {
    Gate * g = gates[i];
    g->set_var_level(g->get_var_level()-2*extended_gates-2);
  }
}



/*------------------------------------------------------------------------*/
Gate * search_for_parent(Term * t, Gate * exclude){
  assert(t);

  Gate * g = gate(t->get_var_num());
  for(auto it = g->parents_begin(); it != g->parents_end(); ++it){
    Gate * parent = *it;
    if(parent == exclude) continue;
    Polynomial * gc = parent->get_gate_constraint();
    if (gc->size() != 2) continue;
    
    if (t == gc->get_tail_term()) return parent;
  }
  return 0;
}
/*------------------------------------------------------------------------*/
Gate * search_for_parent_dual(Term * t){
  assert(t);
  Gate * g = gate(t->get_var_num());
 
  for(auto it = g->parents_begin(); it != g->parents_end(); ++it){
    Gate * parent = *it;
    Polynomial * gc = parent->get_gate_constraint();
  
    if (gc->size() != 2) continue;
    if (t->equal_up_to_duality(gc->get_tail_term())) return parent;
    
  }
  return 0;
}
/*------------------------------------------------------------------------*/
Gate * neutral_child(Gate*parent, Gate*g){
  for (auto it = parent->children_begin(); it!= parent->children_end(); it++){
    Gate * tmp = *it;
    msg("testing %s", tmp->get_var_name());
    if(g->is_neutral(tmp)) return tmp;
  }
  return 0;
}


Polynomial * gen_dual_constraint(Var * v) {
  Var * d = v->get_dual();

  push_mstack_end(new Monomial(minus_one, new_term(v, 0)));
  push_mstack_end(new Monomial(minus_one, new_term(d, 0)));
  push_mstack_end(new Monomial(one, 0));

  Polynomial * p = build_poly();
  return p;
}

/*------------------------------------------------------------------------*/


static void init_gate_constraints() {
  for (unsigned i = NN; i < M-1; i++) {
    init_gate_constraint(i);
  }

  for (unsigned i = 0; i < MM; i++) {
    init_gate_constraint(i+M-1);
  }

  
}


/*------------------------------------------------------------------------*/
static void set_parents_and_children() {
  unsigned pp = 0;

  for (unsigned i = NN; i < M; i++) {
    Gate * n = gates[i];
    assert(!n->get_input());

    aiger_and * and1 = is_model_and(n->get_var_num());

    if (!and1) continue;
    unsigned l = and1->rhs0, r = and1->rhs1;
    Gate * l_gate = gate(l), *r_gate = gate(r);
    n->children_push_back(l_gate);
    n->children_push_back(r_gate);
    if (verbose >= 4) msg("node %s has children %s, %s", n->get_var_name(), l_gate->get_var_name(), r_gate->get_var_name());
  
    if(l_gate && r_gate) {
       if (l_gate->get_input() && r_gate->get_input()
            && !aiger_sign(l) && !aiger_sign(r)) {
        n->mark_pp();
        pp++;
        if (verbose >= 4) msg("partial product %s", n->get_var_name());
      }
    }
    if(l_gate) l_gate->parents_push_back(n);
    if(r_gate) r_gate->parents_push_back(n);
    
  }

  // set children for extra outputs
  for (unsigned i = 0; i < MM; i++) {
    Gate * n = gates[i+M-1];
    assert(n->get_output());
    unsigned lit = slit(i);
    if(lit < 2) continue;
    Gate * model_output_gate = gate(lit);
    n->children_push_back(model_output_gate);
    if (verbose >= 4) msg("node %s has child %s", n->get_var_name(), model_output_gate->get_var_name());
    model_output_gate->parents_push_back(n);
  }

  if (verbose >= 1) msg("found %i partial products", pp);
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
void init_gates(){
  allocate_gates();
  
  init_gate_constraints();
  set_parents_and_children();
}

/*------------------------------------------------------------------------*/

void delete_gates() {
  for (unsigned i = 0; i < num_gates; i++) {
    delete(gates[i]);
  }
  delete[] gates;
}

/*------------------------------------------------------------------------*/
void init_slices() {

  for (unsigned j = M+MM-2; j >= NN; j--) {
    slices.push_back(gates[j]);
  }
}
/*------------------------------------------------------------------------*/
void resort_slices_lex(){
  for (size_t i = 0; i < slices.size(); i++) {
    Polynomial * p = slices[i]->get_gate_constraint();
    if(!p) continue;
    Polynomial * tmp = p->reorder_lex();
    delete(p);
    slices[i]->set_gate_constraint(tmp);
    }

}
/*------------------------------------------------------------------------*/

void print_slices(){
  msg("print slices");

  for (auto it=slices.begin(); it != slices.end(); ++it) {
    Gate * g = *it;
    if(g->get_red()) continue;
    else if(g->parents_size() == 0){
      g->mark_red();
      continue;
    } 
    g->get_gate_constraint()->print(stdout);
  }
  msg("");
}

/*------------------------------------------------------------------------*/
