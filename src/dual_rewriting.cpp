/*------------------------------------------------------------------------*/
/*! \file dual_rewriting.cpp
    \brief contains function to identify the final stage adder

  Part of TeluMA : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2021 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#include <algorithm>

#include "dual_rewriting.h"
/*------------------------------------------------------------------------*/

std::vector<Gate*> resolved;

std::vector<Polynomial*> resolved_gc;

std::vector<Polynomial*> init_gc;

std::vector<Polynomial*> proof_poly;
Polynomial * target_proof;

std::vector<Term*> mult1v;
std::vector<Term*> remv;
std::vector<Term*> outerv;
std::vector<Polynomial *> target_tmp;
/*------------------------------------------------------------------------*/



static bool is_unit(Gate * g){
  if(g->get_red()) return 0;
  Polynomial * gc = g->get_gate_constraint();

  if (gc->size() > 2) return 0;
  if (gc->size() == 1) return 1;

  Term * t = gc->get_tail_term();
  if(t->degree() == 1) return 1;
  return 0;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

static void remove_var_from_stored_terms(std::vector<Term*> vec, Var * v){
  for(auto it = vec.begin(); it != vec.end(); ++it){
    Term * t = *it;
    if(!t->contains(v) && !t->contains(v->get_dual())) continue;
    while(t){
      if(t->get_var() != v && t->get_var() != v->get_dual()){
        add_to_vstack(t->get_var());
      } else if (t->get_var() == v){
        clear_vstack();
        break;

      }
      t = t->get_rest();
    }
    Term * new_t = build_term_from_stack();
    *it = new_t;
  }
}

static void replace_var_from_stored_terms(std::vector<Term*> vec, Var * v, Var * w){
  for(auto it = vec.begin(); it != vec.end(); ++it){
    Term * t = *it;
    if(!t->contains(v) && !t->contains(v->get_dual())) continue;
    while(t){
      if(t->get_var() != v && t->get_var() != v->get_dual()){
        add_to_vstack(t->get_var());
      } else if (t->get_var() == v){
        add_to_vstack(w);
      } else {
        add_to_vstack(w->get_dual());
      }
      t = t->get_rest();
    }
    Term * new_t = build_term_from_stack();
    *it = new_t;
  }
}


static void update_stored_terms(Gate * g){
  Polynomial * gc = g->get_gate_constraint();
  if(gc->size() == 1){
    remove_var_from_stored_terms(mult1v, g->get_var());
    remove_var_from_stored_terms(remv, g->get_var());
  }
  else if (is_unit(g)){
    replace_var_from_stored_terms(mult1v, g->get_var(), gc->get_tail_term()->get_var());
    replace_var_from_stored_terms(remv, g->get_var(), gc->get_tail_term()->get_var());
  }
}


static void eliminate_by_one_gate(Gate * n1, Gate *n2) {
  Polynomial * tmp_p1 = n1->get_gate_constraint();
  Polynomial * flip = gen_dual_constraint(n2->get_var());
  Polynomial * p1 = reduce_by_one_poly(tmp_p1, flip);
  Polynomial * p2 = n2->get_gate_constraint();
  const Polynomial * negfactor = divide_by_term(p1, p2->get_lt());


  if (negfactor->is_constant_zero_poly()) {
    delete(p1);
    delete(negfactor);
    delete(flip);
    return;
  }


  Polynomial * mult   = multiply_poly(negfactor, p2);
  Polynomial * rem    = add_poly(p1, mult);

  n1->update_gate_poly(rem);

  delete(flip);
  delete(mult);
  delete(negfactor);
  delete(p1);
  delete(tmp_p1);
}


static void eliminate_unit_gate(Gate * n){

  update_stored_terms(n);

  for (auto it_c = n->children_begin(); it_c != n->children_end(); ++it_c) {
    Gate * n_child = *it_c;
    n_child->parents_remove(n);
  }

  for (auto it_p = n->parents_begin(); it_p != n->parents_end(); ++it_p) {
    Gate * n_parent = *it_p;
     msg("before her1");
    eliminate_by_one_gate(n_parent, n);
    n_parent->children_remove(n);

    for (auto it_c = n->children_begin(); it_c != n->children_end(); ++it_c) {
      Gate * n_child = *it_c;
      if (!n_parent->is_child(n_child)) n_parent->children_push_back(n_child);
      if (!n_child->is_in_parents(n_parent)) n_child->parents_push_back(n_parent);
    }

    if(is_unit(n_parent)){
      eliminate_unit_gate(n_parent);
    }
    else if(n_parent->children_size() == 1 && n_parent->get_gate_constraint()->size()==3){
      auto it = n_parent->children_begin();
      Gate * tmp = *it;
      Polynomial * flip = gen_dual_constraint(tmp->get_var());
      Polynomial * rem1 = reduce_by_one_poly(n_parent->get_gate_constraint(), flip);
      delete(flip);

      if(rem1->size()!=2){
        delete(rem1);
        flip = gen_dual_constraint(tmp->get_var()->get_dual());
        rem1 = reduce_by_one_poly(n_parent->get_gate_constraint(), flip);
        delete(flip);
      }

      n_parent->update_gate_poly(rem1);
      eliminate_unit_gate(n_parent);
    }

  }


  if (verbose>2) msg("removed unit %s", n->get_var_name());
}

/*----------------------------------------------------------------------------*/

//
static void remove_only_positives(size_t parent_limit = 0) { //checked
  int counter = 0;
 // for (unsigned i=M-1; i >= NN; i--) { //do not even consider changing direction
 //    Gate * n = gates[i];
  for (int i = slices.size()-1; i >= 0; i--) {
    Gate * n = slices[i];
    if(parent_limit > 0 && n->parents_size() > parent_limit) continue;
    if(!parent_limit && n->parents_size() == 1) continue;
  //  if(n->get_aig_output() && n == gate(slit(2))) break;
    if (n->get_pp()) break;
    if(n->get_input()) break; 
    
    if (n->get_red()) continue;
    if (n->get_output() || n->get_aig_output()) continue;
    
    if (n->get_gate_constraint()->size() > 2) continue;
    bool flag = 0;

    for(auto & n_parent : n->get_parents()){
      if(n_parent->get_gate_constraint()->size() > 2) {
        flag = 1;
        break;
      }

      Monomial * m = n_parent->get_gate_constraint()->get_mon(1);
      if(!m->get_term()->contains(n->get_var())) {
        flag = 1;
        break;
      }
    }
    if (flag) continue;
    

    for(auto & n_child : n->get_children()){
      n_child->parents_remove(n);
    }

    for(auto & n_parent : n->get_parents()){
      Polynomial * rem = reduce_by_one_poly(n_parent->get_gate_constraint(), n->get_gate_constraint());

      n_parent->set_gate_constraint(rem);
      for (auto it_c = n->children_begin(); it_c != n->children_end(); ++it_c) {
        Gate * n_child = *it_c;
        n_child->parents_push_back(n_parent);
        n_parent->children_push_back(n_child);
      }
      n_parent->children_remove(n);
      
    }


    n->mark_red();
    delete(n->get_gate_constraint());
    n->set_gate_constraint(0);
   
    counter++;
    if(verbose > 2) msg("eliminated %s", n->get_var_name());

  }
  if (verbose >= 1) msg("removed %i positive gates", counter);
  elim_pos_nodes+= counter;
}
/*============================================================================*/


static bool do_backward_substitution(Gate * outer){ //checked
  Gate * inner;
  Polynomial * inner_gc;
  Monomial * inner_m;
  Term * inner_t;

  Polynomial * outer_gc = outer->get_gate_constraint();
  Monomial * outer_m = outer_gc->get_mon(outer_gc->size()-1);
  Term * outer_t = outer_m->get_term();
  Term * outer_t_it = outer_t;
  Term * div_term = 0;

  

  while(outer_t_it){
    div_term = divide_by_var(outer_t, outer_t_it->get_var());
      
    if(div_term -> get_ref() > 1){
      Gate * tmp = gate(div_term->get_var_num());
  
      for(auto it = tmp->parents_begin(); it != tmp->parents_end(); ++it) {
        inner = *it;
            
        if (inner == outer) continue;
        if (inner->get_red()) continue;
      
        inner_gc = inner->get_gate_constraint();
     
        if (inner_gc->size() != 2) continue;
        if (inner_gc->degree() > outer_gc->degree()) continue;
        
        inner_m = inner_gc->get_mon(1);
        inner_t = inner_m->get_term();
        
        if (inner_t != div_term) continue;

        Term * t0 = new_term(outer_t_it->get_var(), 0);
        Term * t1 = new_term(inner_gc->get_lt()->get_var(), 0);
        Term * t2;
        if(t0) t2 = multiply_term(t0, t1);
        else t2 = t1;



        push_mstack_end(outer_gc->get_mon(0)->copy());
        Monomial * tmp = new Monomial(one, t2);
        push_mstack_end(tmp);
        Polynomial * rewr = build_poly(); 


   
        outer->update_gate_poly(rewr);
    
    

        if(verbose > 3)
         msg("substituted %s in %s", inner->get_var_name(), outer->get_var_name());

        return 1;

      }
    
    }
 
    deallocate_term(div_term);
    div_term = 0;
    outer_t_it = outer_t_it->get_rest();
  
  }
  return 0;
}

/*------------------------------------------------------------------------*/
std::vector<Gate*> sub;


/*------------------------------------------------------------------------*/

static void backward_substitution() { //checked
  msg("started backward substitution");

  int counter = 0;

  Gate * outer;
  Polynomial * outer_gc;
  Monomial * outer_m;
  Term * outer_t;

  for (unsigned i = M-2; i >= NN; i--) {
    outer = gates[i];

    if (outer->get_red()) continue;
    if (outer->get_pp()) continue;
    
    outer_gc = outer->get_gate_constraint();
    if (outer_gc->size() != 2) continue;

    outer_m = outer_gc->get_mon(1);
    outer_t = outer_m->get_term();
    if(!outer_t) continue;

    if (outer_t->degree() < 3) continue;

    bool rewr = do_backward_substitution(outer);
    while(rewr){
      rewr = do_backward_substitution(outer);
    }
    if(is_unit(outer)) eliminate_unit_gate(outer);
    sub.push_back(outer);
    counter++;
  }

  if (verbose >= 1) msg("backwards substitution done", counter);
  //linearize_backward_sub();
}
/*============================================================================*/

/*----------------------------------------------------------------------------*/
static bool compare_poly_with_tail_of_gc(Polynomial * gc, Polynomial * tail){
    if(gc->size() != 3) return 0;
    if(tail->size() != 2) return 0;

    Monomial * m1 = tail->get_mon(0);
    for (size_t i = 1; i < gc->size(); i++){
      Monomial * m2 = gc->get_mon(i);
      if (m1->get_term()->cmp(m2->get_term())) return 0;
    
    }
    return 1;
}

static void update_multv_remv_outerv(Term *m, Term *r, Term *o){ 
  mult1v.push_back(m);
  remv.push_back(r);
  outerv.push_back(o);
}

static void add_to_resolved(Gate *g){
  resolved.push_back(g);
  Polynomial * p = gen_gate_constraint(g->get_var_num()/2-1);
  resolved_gc.push_back(p);
}

static int resolved_contain_tail(Term * t){
  Gate * g = gate(t->get_var_num());
  auto it = std::find(resolved.begin(), resolved.end(), g);
  if(it != resolved.end()) return it- resolved.begin();
  return -1;
}

static Term * search_for_poly_in_resolved(Polynomial * rem){
  for (auto i = resolved.begin(); i != resolved.end(); ++i){
    Gate * g = *i;
    Polynomial * gc = g->get_gate_constraint();
    if(gc->size() != 3) continue;
    if(compare_poly_with_tail_of_gc(gc, rem)) return gc->get_lt();
  }
  return 0;
}



static Var * search_for_tail(Term * t){
  assert(t);
  Gate * g = gate(t->get_var_num());
  for(auto it = g->parents_begin(); it != g->parents_end(); ++it){
    Gate * parent = *it;
    Polynomial * gc = parent->get_gate_constraint();
    if (gc->size() != 2) continue;
    if (t == gc->get_tail_term()) return parent->get_var();
  }
  return 0;
}
/*-----------------------------------------------------------------------------*/
static Var * search_for_tail_in_init(Term * t){
  assert(t);
  for(auto it = init_gc.begin(); it != init_gc.end(); ++it){
    Polynomial * gc = *it;
    if (t == gc->get_tail_term()) return gc->get_lt()->get_var();
  }
  return 0;
}

static Var * search_for_tail_in_resolved(Term * t){
  for(auto it = resolved.begin(); it != resolved.end(); ++it){
    Gate * g = *it;
    int index = it - resolved.begin();
    Polynomial * gc = resolved_gc[index];
    if (gc->size() != 2) continue;

    if (t == gc->get_tail_term()) return g->get_var();
  }
  return 0;
}

static Term * backward_rewrite_term(Term *t){
  if(t->degree() == 1) return t->copy();
  Var * tail;
  std::vector<Var*> remainder;
  std::vector<Var*> remainder2;


  while(t){
    if(t->get_ref() > 1){
      tail = search_for_tail(t);
      if(tail) {
        remainder.push_back(tail);
        break;
      }
    }
    remainder.push_back(t->get_var());
    remainder2.push_back(t->get_var());
    t = t->get_rest();
  }
  Term * quo = sort_and_build_term_from_vector(remainder);
  return quo;
}

Term * backward_rewrite_term_until_completion(Term * t){
  if(t->degree() == 1) return t;
  Term * tmp = backward_rewrite_term(t);
  while(tmp != t) {
    deallocate_term(t);
    t  = tmp;
    tmp = backward_rewrite_term(t);
  }
  deallocate_term(t);
  return tmp;
}


static size_t get_dual_count(Term *t){
  size_t cnt = 0;
  while (t){
    if(t->get_var()->is_dual()) cnt++;
    t = t->get_rest();
  }
  return cnt;
}
/*-----------------------------------------------------------------------------*/
static void remove_positives_in_carry (Gate * n) { // check

    Polynomial * n_gc = n->get_gate_constraint();
    Term * n_t = n_gc->get_tail_term();
    if (n_t->degree() == get_dual_count(n_t)) return;
    const Var * n_v;
    Gate * child;

    bool flag = 0;

    while(n_t && !flag){
      n_v = n_t->get_var();
      child = gate(n_v->get_num());

      if(child->get_input()) break;

      if(!n_v->is_dual()  && !child->get_pp() &&
        child->get_gate_constraint()->size() == 2) {
        Polynomial * rem = reduce_by_one_poly(n_gc, child->get_gate_constraint());
        n->update_gate_poly(rem);
    
        flag = 1;
        if(verbose >= 3)
          msg("substituted %s in %s", child->get_var_name(), n->get_var_name());

      }
      n_t = n_t->get_rest();
    }

    if(flag) remove_positives_in_carry(n);
}
/*-----------------------------------------------------------------------------------*/
static Term * expand_term(Term *t){
  Term * tail;
  std::vector<Var*> remainder;

  while(t){
    if(!t->get_var()->is_dual()){
      Polynomial * p = gate(t->get_var_num())->get_gate_constraint();

      if(p->size() != 2) return t;
      tail = p->get_tail_term();
      while(tail) {
        remainder.push_back(tail->get_var());
        tail = tail->get_rest();
      }
    }
    else remainder.push_back(t->get_var());
    t = t->get_rest();
  }
  Term * quo = sort_and_build_term_from_vector(remainder);
  return quo;
}
/*-----------------------------------------------------------------------------------*/
static Term * expand_term_until_completion(Term * t){
  Term * tmp = expand_term(t);
  while(tmp != t) {
    deallocate_term(t);
    t  = tmp;
    tmp = expand_term(t);
  }
  deallocate_term(t);
  return tmp;
}
/*-----------------------------------------------------------------------------------*/
static Term * divide_tail_term_by_pivot(Gate *g, Term * t){
  Polynomial * g_gc = g->get_gate_constraint();
  if (g_gc->size() != 2) return 0;
  Term * t0 = g_gc->get_tail_term();

  Term * rem = remainder_t (t0, t);

  if(rem && rem == t0) return 0;

   if(rem){
     if(rem->degree() > 1)
      rem = backward_rewrite_term_until_completion(rem);
   } 
  return rem;
}
/*-----------------------------------------------------------------------------------*/
static bool term_vector_contains_only_singles(std::vector<Term*> vec){
  for(auto i = vec.begin(); i != vec.end(); ++i){
    Term * tmp = *i;
    if(tmp->degree()!= 1) return 0;
  }
  return 1;
}

/*-----------------------------------------------------------------------------*/
static void expand_tail_term_and_update_gc (Gate *g){ //check
  Polynomial * gc = g->get_gate_constraint();

  if(gc->size() != 2) return;
  Term * t = gc->get_tail_term();
  Term * t_exp  = expand_term(t);

    push_mstack(new Monomial(minus_one, t->copy()));
    push_mstack(new Monomial(one, t_exp->copy()));

    Polynomial * p = build_poly();
    Polynomial * add = add_poly(gc, p);
   
    delete(p);
    delete(gc);
    g->update_gate_poly(add);
  
}
/*----------------------------------------------------------------------------*/
static Term * calc_pivot(Term * t0, Term *t1){ //checked
  while(t0 && t1) {
    if (t0->get_var() == t1->get_var()){
      add_to_vstack(t0->get_var());
      t0 = t0->get_rest();
      t1 = t1->get_rest();
    } else if (t0->get_var_level() > t1->get_var_level()){
      t0 = t0->get_rest();
    } else {
      t1 = t1->get_rest();
    }
  }
  Term * res = build_term_from_stack();
  return res;
}




static Term * carry_skip_rewriting(Term * head, Term * tail, Term * pivot, Term * rem){

  auto it = find(outerv.begin(), outerv.end(), tail);
  if(it == outerv.end()) return 0;
  int index = it - outerv.begin();
  Term * unrolled_tail = remv[index];
  Term * unrolled_head = mult1v[index];


  if(unrolled_tail->degree()==1 && unrolled_head->degree() == 2 &&
     head->contains(unrolled_tail->get_var()->get_dual())){

     Polynomial * p = build_one_minus_x_times_y_poly(unrolled_head, unrolled_tail);
     Term * found = search_for_poly_in_resolved(p);
     delete(p);
     if(found) return found;
  }

  Term * unrolled_pivot = calc_pivot(head, unrolled_head);
  if(!unrolled_pivot) return 0;


  Term * t0 = remainder_t(unrolled_head, unrolled_pivot);

  if(t0 && (t0->degree() != 1 || !t0->get_var()->is_dual())) return 0;

  if(t0) {
    Term * t0_flipped = new_term(t0->get_var()->get_dual(),0);
    t0 = t0_flipped; // TODO flip?

 } else {
    if (unrolled_head->degree() == 2) {
      Polynomial * p = build_one_minus_x_times_y_poly(unrolled_head, unrolled_tail);

      Term * found = search_for_poly_in_resolved(p);
      delete(p);
      if(found) return found;

    }
  }

  Term * t1 = remainder_t(head, unrolled_pivot);
  Term * quo = 0;



  if(t1->degree() > 1){
    Term * u_p = multiply_term(pivot, unrolled_pivot);
    Term * u_r = multiply_term(rem, unrolled_tail);
     Term * found = carry_skip_rewriting(t1, t0, u_p, u_r);
     if(!found) return 0;
     quo = new_term(found->get_var()->get_dual(), 0);

  } else {

    std::vector<Var*> tmp;
    tmp.push_back(t0->get_var());
    tmp.push_back(t1->get_var()->get_dual());
    quo = sort_and_build_term_from_vector(tmp);

    int idx = resolved_contain_tail(t0);
    Gate * parent = resolved[idx];
    Polynomial * parent_gc = resolved_gc[idx];
    Var * child = quo->get_rest()->get_var();
    if(parent_gc->size() == 2 && parent_gc->get_tail_term()->contains(child)){
        deallocate_term(quo);
        quo = new_term(parent->get_var()->get_dual(), 0);
    }
  }

  Term * new_head = multiply_term(unrolled_pivot, quo);



  if (new_head->degree() == 2) {
    Polynomial * p = build_one_minus_x_times_y_poly(new_head, unrolled_tail);
    Term * found = search_for_poly_in_resolved(p);
    delete(p);
    if(found) return found;

  }
  return 0;

}

static Term * get_all_but_last_var_of_term(Term *t){
  while(t->get_rest()){
    add_to_vstack(t->get_var());
    t = t->get_rest();
  }
  return build_term_from_stack(0);
}

static Term * get_last_var_of_term(Term *t){
  while(t->get_rest()){
    t = t->get_rest();
  }
  return new_term(t->get_var(), 0);
}


static Term * expand_head_to_find(Term * head, Term * tail){
  if (head->degree() == 1) return 0;
  if (tail->degree() != 1) return 0;
  if (tail->get_var()->is_dual()) return 0;

  Term * tail_exp = expand_term_until_completion(tail);
  Term * head_h = get_all_but_last_var_of_term(head);
  Term * head_t = get_last_var_of_term(head);

  int idx = resolved_contain_tail(head_t);
  if(idx == -1) return 0;

  Polynomial * gc = resolved_gc[idx];
  Term * gc_t = gc->get_tail_term();

  Term * t0 = multiply_term(head_h, new_term(gc_t->get_var()->get_dual(), 0));
  Term * t1 = multiply_term(head_h, new_term(gc_t->get_rest()->get_var()->get_dual(), 0));

  t0 = backward_rewrite_term_until_completion(t0);
  if(t0->degree() > 1) return 0;
  t0 = new_term(t0->get_var()->get_dual(),0);

  t1 = backward_rewrite_term_until_completion(t1);
  if(t1->degree() > 1) return 0;
  t1 = new_term(t1->get_var()->get_dual(),0);



  Term * t2 = multiply_term(t0, t1);
//  deallocate_term(t0);
//  deallocate_term(t1);

  Term * final = multiply_term(tail_exp, t2);
//  deallocate_term(t2);

  Var * v = search_for_tail_in_init(final);
//  deallocate_term(final);

  if(v)  return new_term(v,0);
  return 0;
} 


static Term * second_level_rewriting(std::vector<Term*> vec, Term * pivot, Term * rem){
  Term * head = vec.front()->degree() > vec.back()->degree() ? vec.front() : vec.back();
  Term * tail = vec.front()->degree() < vec.back()->degree() ? vec.front() : vec.back();

  if (tail->degree() != 1) return 0;

  Term * tail_flip = new_term(tail->get_var()->get_dual(), 0);
  deallocate_term(tail);
  tail = tail_flip;

  if(!tail->get_var()->is_dual() && resolved_contain_tail(tail) != -1){
     Term * found = carry_skip_rewriting(head, tail, pivot, rem);
     if(found) return found;
  }

  if (head->degree() == 2) {
    Polynomial * p = build_one_minus_x_times_y_poly(head, tail);
    Term * found = search_for_poly_in_resolved(p);
    delete(p);
    if(found) return found;
  }

  auto it = find(mult1v.begin(), mult1v.end(), head);
  if(it == mult1v.end()) {
    head = expand_term(head);
    it = find(mult1v.begin(), mult1v.end(), head);
    if(it == mult1v.end()) {
      return expand_head_to_find(head, tail);
    }
  }
  int index = it - mult1v.begin();

  Polynomial * p = gate(tail->get_var_num())->get_gate_constraint();
  if(p->size() != 2) return 0;
  if (remv[index] == p->get_tail_term()) {
    return outerv[index];
  }

  Term * tail_exp = expand_term_until_completion(tail);
  if (remv[index] == tail_exp) {
    return outerv[index];
  }
  return 0;

}

/*-----------------------------------------------------------------------------*/
static bool size_three_vec_first_has_size_two(std::vector<Term*> vec){
  if(vec.size() != 3) return 0;
  if(vec.front()->degree() != 2) return 0;
  if(vec[1]->degree()!= 1 || vec[2]->degree()!= 1) return 0;
  return 1;
}

static bool size_three_vec_first_has_size_three(std::vector<Term*> vec){
  if(vec.size() != 3) return 0;
  if(vec.front()->degree() != 3) return 0;
  if(vec[1]->degree()!= 1 || vec[2]->degree()!= 1) return 0;
  return 1;
}
/*----------------------------------------------------------------------------*/


static void expand_first_two(Term * t){
  if (t->degree() <= 2) return;

  Gate * g0 = gate(t->get_var_num());
  Gate * g1 = gate(t->get_rest()->get_var_num());
  expand_tail_term_and_update_gc(g0);
  expand_tail_term_and_update_gc(g1);
}
/*-----------------------------------------------------------------------------*/

static std::vector<Term*> internal_pivot_calc (Term * outer, Term * pivot, Term * carry = 0){
  std::vector<Var*> remainder;
  std::vector<Term*> quotient;
  std::vector<Term*> res;
  Term * outer_t = outer;



  if(carry) quotient.push_back(carry);

  while(outer_t){
    Term * t0 = divide_tail_term_by_pivot(gate(outer_t->get_var_num()), pivot);
    
    if(t0) {
     
      quotient.push_back(t0);
      
    } else {
      remainder.push_back(outer_t->get_var());
    }
    outer_t = outer_t->get_rest();
  }
  Term * rem = sort_and_build_term_from_vector(remainder);
  


 
  Term * quo;
  if(term_vector_contains_only_singles(quotient)){
    quo = gen_dual_term_from_vector(quotient);
  } else {
 
    Term * tmp_t = 0;
    if(quotient.size() == 2) {
  
      tmp_t = second_level_rewriting(quotient, pivot, rem);
      if(!tmp_t) return res;
      quo = new_term(tmp_t->get_var()->get_dual(), 0);


    } else if (size_three_vec_first_has_size_two(quotient) && !carry) {
      Gate * g = gate(outer->get_rest()->get_var_num());
      Polynomial * g_gc = g->get_gate_constraint();
      Gate * g1 = gate(g_gc->get_tail_term()->get_rest()->get_var_num());
    
      eliminate_by_one_gate(g, g1); // was guarded by xor
      return res;


    } else return res;
  }
 
 if(quo && quo->degree() == 2 && !quo->get_var()->is_dual()){
   Gate * parent = gate(quo->get_var_num());
   Polynomial * parent_gc = parent->get_gate_constraint();
   Var * child = quo->get_rest()->get_var();
   if(parent_gc->size() == 2 && parent_gc->get_tail_term()->contains(child)){
     deallocate_term(quo);
     quo = new_term(parent->get_var()->get_dual(), 0);
   }
 }

 if(quo && quo->degree() > 1) {

   Var * v = search_for_tail(quo);
   
  if(!v){
    v = search_for_tail_in_resolved(quo);
    
  }

   if(!v){

     Term * tmp = backward_rewrite_term(quo);

     v = search_for_tail(tmp);

     if(!v) {
       v = search_for_tail_in_resolved(tmp);
     }


     if(!v) {
       v = search_for_tail_in_init(tmp);

     }
   }

   if(v) {
     deallocate_term(quo);
     quo = new_term(v->get_dual(), 0);
   }
 }

 if(quo && quo->degree() == 1){
   Gate * g = gate(quo->get_var_num());
   if(g->get_gate_constraint()->size() == 1){
     deallocate_term(quo);
     quo = 0;
   }
 }
 res.push_back(quo);
 res.push_back(rem);
 return res;

}

/*-----------------------------------------------------------------------------*/


static std::vector<Term *> larger_pivot_rewriting(Term * mul, Term * rem){
  mul = backward_rewrite_term_until_completion(mul);
  std::vector<Term*> rewritten;

  Term * pivot =  calc_pivot(mul, gate(rem->get_var_num())->get_gate_constraint()->get_tail_term());

  if(!pivot) {
    rewritten.push_back(mul->copy());
    rewritten.push_back(rem->copy());
    return rewritten;
  }

  Term * carry = remainder_t(mul, pivot);

  std::vector<Term*> res = internal_pivot_calc(rem, pivot, carry);
  if(res.empty()) return rewritten;

  Term * new_sub;

  if(res.front()) new_sub = multiply_term(pivot, res.front());
  else new_sub = pivot;
  new_sub = backward_rewrite_term_until_completion(new_sub);

  rewritten.push_back(new_sub);
  rewritten.push_back(res.back());
  return rewritten;
}

/*----------------------------------------------------------------------------*/
static int pivot_rewriting(Gate * outer, Term * pivot){
  
  
  Polynomial * outer_gc = outer->get_gate_constraint();

  Term * outer_t = outer_gc->get_tail_term();
  

  // first step where all come from outer_t
  std::vector<Term*> res = internal_pivot_calc(outer_t, pivot);
 
  if(res.empty()) return 0;
  Term * sub = res.front();
  Term * rem = res.back();

  Term * mult1 = 0, * mult2, * mult3 = 0;
  bool larger_sub = 0;

  if(sub && sub->degree() == 1 && is_unit(gate(sub->get_var_num()))){
    Polynomial * gc = gate(sub->get_var_num())->get_gate_constraint();
    Term * tail = gc->get_tail_term();
    Term * sub_tmp = 0;
    if (sub->get_var()->is_dual())
      sub_tmp = new_term(tail->get_var()->get_dual(), 0);
    else sub_tmp = tail;
    deallocate_term(sub);
    sub = sub_tmp;
  }

 

  if(sub && sub->degree() == 1){
     mult1 = multiply_term(pivot, sub);
  } else if (sub && pivot->degree() == 1) {
     mult1 = multiply_term(pivot, sub);
     mult3 = multiply_term(pivot,rem);
     larger_sub = 1;
  } else if (sub) {
    deallocate_term(sub);
    return -1;
  } else mult1 = pivot->copy();

  

  if(pivot->degree() > 1){
    Term * old_mul = 0;
    while(old_mul != mult1){
      update_multv_remv_outerv(mult1->copy(), rem->copy(), outer_gc->get_lt());
      std::vector<Term*> rewritten = larger_pivot_rewriting(mult1, rem);
      if(rewritten.empty()) return 0;

      deallocate_term(rem);
      old_mul = mult1;
      mult1 = rewritten.front();
      rem = rewritten.back();
    }
    deallocate_term(old_mul);
  }


  if(rem) {
    mult2 = multiply_term(mult1, rem);
    deallocate_term(mult1);
  }  else mult2 = mult1;

  push_mstack(outer_gc->get_lm()->copy());
  if(rem != pivot){
    if(larger_sub){
      push_mstack(new Monomial(one, mult2));
      push_mstack(new Monomial(minus_one, mult3));
    } else {
      push_mstack(new Monomial(minus_one, mult2));
    }
    push_mstack(new Monomial(one, rem));
  }

  Polynomial * new_gc = build_poly();

  outer->delete_children();

  if(rem != pivot){
    while(mult2){
      Gate * n_child = gate(mult2->get_var_num());
      outer->children_push_back(n_child);
      n_child->parents_push_back(outer);
      mult2 = mult2->get_rest();
    }
  }



  //delete(outer_gc);
  outer->update_gate_poly(new_gc);


  if(is_unit(outer)){
    outer->delete_children();
    eliminate_unit_gate(outer);
  }

  return 1;
}
/*-----------------------------------------------------------------------------*/
static void removing_f_chains() {
  msg("removing negative chains");

  bool change = 1, first = 1;
  while(change) {
    change = 0;
    for (auto i = sub.begin(); i != sub.end(); ++i) {

      Gate * outer = *i;
      if (outer->get_red()) continue;
      if (outer->get_pp()) continue;

      Polynomial * outer_gc = outer->get_gate_constraint();
      if (outer_gc->size() != 2) continue;
      if (outer_gc->degree() < 4) continue;

      if (verbose > 1) {
        fputs("\n----- REWRITING -----\n\n", stdout);
        outer->get_gate_constraint()->print(stdout);
        msg("parentsize %i", outer->parents_size());
        fputs("\n",stdout);
      }


     
      Term * outer_t = outer_gc->get_tail_term();


    
      if(first){
        
        init_gc.push_back(outer_gc->copy());

        if(outer_t->degree() != get_dual_count(outer_t)){
          remove_positives_in_carry(outer);   // clean
          outer_gc = outer->get_gate_constraint();
          init_gc.push_back(outer_gc->copy());
          outer_t = outer_gc->get_tail_term();
        }
      }
   
      outer_t = outer_gc->get_tail_term();
      int flag = 0;
      target_proof = outer->get_gate_constraint()->copy();

      if(outer_t->degree() < 2) { msg("cannot calculate pivot of unit");
      } else if(outer_t->degree() == 2) {

        flag = 1;
      } else {

        Gate * g0 = gate(outer_t->get_var_num());
        Gate * g1 = gate(outer_t->get_rest()->get_var_num());

        Polynomial * g0_gc = g0->get_gate_constraint();
        if (g0_gc->size() != 2) continue;
     

        Polynomial * g1_gc = g1->get_gate_constraint();
        if (g1_gc->size() != 2) continue;


        Term * t0 = g0_gc->get_tail_term();
        Term * t1 = g1_gc->get_tail_term();

        Term * pivot =  calc_pivot(t0, t1);
    
        if(!pivot || pivot->degree()>1) continue;
 
        flag = pivot_rewriting(outer, pivot);
        
        if (flag == -1) expand_first_two(outer_t);

        deallocate_term(pivot);
      }

      if (flag == 1) {
        add_to_resolved(outer);
        sub.erase(i--);
        change = 1;
      }

      if(verbose > 1 && !outer->get_red())
        outer->get_gate_constraint()->print(stdout);

    }

    if(first && !sub.empty()) {
      std::reverse(sub.begin(), sub.end());
      first = 0;
    }
  }

  for(auto it = init_gc.begin(); it!= init_gc.end(); ++it){
    Polynomial * p = *it;
    delete(p);
  }
}



/*============================================================================*/

void dual_preprocessing(){
    if(mult_inp) {
      remove_only_positives(1);
      remove_only_positives(0);
    }
    backward_substitution();
    removing_f_chains();
    if(verbose > 3) print_slices();
}

