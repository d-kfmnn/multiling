/*------------------------------------------------------------------------*/
/*! \file term.cpp
    \brief contains the class Term and further functions to
    manipulate terms

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#include "term.h"
/*------------------------------------------------------------------------*/

Term::Term(Var * _v,  Term * _r, uint64_t _hash, Term * _n):
  variable(_v), ref(1), hash(_hash), next(_n) {
  if (_r) {
    rest = _r->copy();
    deg = (_r->degree())+1;
  } else {
    rest = 0;
    deg = 1;
  }
}

/*------------------------------------------------------------------------*/

Term * Term::copy() {
  assert(ref > 0);
  ++ref;
  assert(ref);
  return this;
}



/*------------------------------------------------------------------------*/

void Term::print(FILE * file) const {
  const Term * tmp = this;
  if (!tmp) fputc_unlocked('0', file);
  while (tmp) {
    fputs_unlocked(tmp->get_var_name(), file);
    tmp = tmp->get_rest();
    if (tmp) fputc_unlocked('*', file);
  }
}


std::string Term::print_to_string() const{
  std::string res;
  const Term * tmp = this;
  if (!tmp) res = "0";
  while (tmp) {
    for (const char* ptr = tmp->get_var_name(); *ptr != '\0'; ++ptr) {
        res += *ptr; // Append each character to the string
    }
   
    tmp = tmp->get_rest();
    if (tmp) res += "*";
  }
  return res;
}


/*------------------------------------------------------------------------*/
bool Term::equal_up_to_duality(const Term *t) const {
  const Term * tmp1 = this;
  const Term * tmp2 = t;
  if (tmp1 == tmp2) return 1;

  while (tmp1 && tmp2) {

      if (tmp1->get_var() != tmp2->get_var() && tmp1->get_var() != tmp2->get_var()->get_dual()) return 0;
      tmp1 = tmp1->get_rest();
      tmp2 = tmp2->get_rest();
  }  
  if (tmp1) return 0;
  else if (tmp2) return 0;
  

  return 1;
}

/*------------------------------------------------------------------------*/

int Term::cmp(const Term *t) const {
  const Term * tmp1 = this;
  if (!tmp1) return -1;
  const Term * tmp2 = t;
  if (!tmp2) return 1;

  if (tmp1 != tmp2) {
    while (tmp1 && tmp2) {
      if (tmp1->get_var_level() > tmp2->get_var_level()) return 1;
      else if (tmp1->get_var_level() < tmp2->get_var_level()) return -1;
      tmp1 = tmp1->get_rest();
      tmp2 = tmp2->get_rest();
    }
    if (tmp1) return 1;
    else if (tmp2) return -1;
  }

  return 0;
}

/*------------------------------------------------------------------------*/

int Term::cmp_dlex(const Term *t) const {
  const Term * tmp1 = this;
  if (!tmp1) return -1;
  const Term * tmp2 = t;
  if (!tmp2) return 1;

  if (tmp1 != tmp2) {
    const Term * tmp1l = tmp1;
    const Term * tmp2l = tmp2;

    while (tmp1l && tmp2l) {
      tmp1l = tmp1l->get_rest();
      tmp2l = tmp2l->get_rest();
    }

    if (tmp1l) return 1;
    else if (tmp2l) return -1;

   while (tmp1 && tmp2) {
      if (tmp1->get_var_level() > tmp2->get_var_level()) return 1;
      else if (tmp1->get_var_level() < tmp2->get_var_level()) return -1;
      tmp1 = tmp1->get_rest();
      tmp2 = tmp2->get_rest();
    }   
  }

  return 0;
}


/*------------------------------------------------------------------------*/

bool Term::contains(Var *v) const {
  assert(v);
  const Term * t = this;
  while (t) {
    if (t->get_var() == v) return 1;
   // else if (t->get_var_level() < v->get_level()) return 0;
    t = t->get_rest();
  }
  return 0;
}

/*------------------------------------------------------------------------*/

bool Term::contains_t(const Term *t) const {
  assert(t);
  const Term * tmp = t;
  while (tmp) {
    if (!this->contains(tmp->get_var())) return 0;
    tmp = tmp->get_rest();
  }
  return 1;
}

/*------------------------------------------------------------------------*/
 bool Term::contains_input() const{
  const Term * t = this;
  while (t) {
    if (t->get_var()->is_input()) return 1;
    t = t->get_rest();
  }
  return 0;
 }
/*------------------------------------------------------------------------*/
 Var * Term::extract_first_dual_var() const {
  const Term * t = this;
  while (t) {
      if (t->get_var()->is_dual()) return t->get_var();
    t = t->get_rest();
  }
  return 0;
 }
/*------------------------------------------------------------------------*/

  Var * Term::first_input() const {
    const Term * t = this;
    while (t) {
    if (t->get_var()->is_input()) return t->get_var();
    t = t->get_rest();
  }
  return 0;

  }

  Var * Term::last_input() const {
    const Term * t = this;
    while(t->get_rest()){
      t = t->get_rest();
    }
    if (t->get_var()->is_input()) return t->get_var();
    return 0;

  }

uint64_t size_terms;
uint64_t current_terms;
Term ** term_table;

/*------------------------------------------------------------------------*/

uint64_t compute_hash_term(Var * variable, const Term * rest) {
  assert(variable);
  uint64_t res = rest ? rest->get_hash() : 0;
  res *= get_nonces_entry(0);
  res += variable->get_hash();
  res *= get_nonces_entry(1);
  return res;
}

/*------------------------------------------------------------------------*/

void enlarge_terms() {
  uint64_t new_size_terms = size_terms ? 2*size_terms : 1;
  Term ** new_term_table = new Term*[new_size_terms]();
  for (uint64_t i = 0; i < size_terms; i++) {
    for (Term * m = term_table[i], * n; m; m = n) {
      uint64_t h = m->get_hash() &(new_size_terms - 1);
      n = m->get_next();
      m->set_next(new_term_table[h]);
      new_term_table[h] = m;
    }
  }
  delete[] term_table;
  term_table = new_term_table;
  size_terms = new_size_terms;
}

/*------------------------------------------------------------------------*/

Term * new_term(Var * variable, Term * rest) {
  if (current_terms == size_terms) enlarge_terms();
  const uint64_t hash = compute_hash_term(variable, rest);
  const uint64_t h = hash &(size_terms - 1);

  Term * res;
  for (res = term_table[h];
       res &&(res->get_var() != variable || res->get_rest() != rest);
       res = res->get_next())
       {}

  if (res) {
    res->inc_ref();  // here we extend that we found term once more
  } else {
    res = new Term(variable, rest, hash, term_table[h]);
    term_table[h] = res;
    current_terms++;
  }
  return res;
}

Term * new_quadratic_term(Var * v1, Var * v2) {
  Term * t1 = new_term(v1, 0);
  Term * t2 = new_term(v2, 0);
  return multiply_term(t1,t2);
}

/*------------------------------------------------------------------------*/

void deallocate_term(Term * t) {
  while (t) {
    assert(t->get_ref() > 0);
    if (t->dec_ref() > 0) break;  // t is still used
    Term * rest = t->get_rest();
    const uint64_t h = t->get_hash() &(size_terms - 1);
    Term * p = term_table[h];
    if (p == t) { term_table[h] = t->get_next();
    } else {
      Term * p2;
      while ((p2 = p->get_next()) != t) p = p2;
      p->set_next(t->get_next());
    }

    assert(current_terms);
    current_terms--;
    delete(t);
    t = rest;
  }
}

/*------------------------------------------------------------------------*/

void deallocate_terms() {
  for (uint64_t i = 0; i < size_terms; i++) {
    for (Term * m = term_table[i], *n; m; m = n) {
      n = m->get_next();
      assert(current_terms);
      current_terms--;

      delete(m);
    }
  }
  delete[] term_table;
}

/*------------------------------------------------------------------------*/
static std::vector<Var*> vstack;  // /< used to build a term

struct {
    bool operator()(const Var * a, const Var * b) const {
      return a->get_level() > b->get_level();
     }
} cmpVarLvl;

/*------------------------------------------------------------------------*/

void add_to_vstack(Var* v) {
  assert(v);
  vstack.push_back(v);
}
/*------------------------------------------------------------------------*/
void clear_vstack(){
  vstack = std::vector<Var*>();
}


/*------------------------------------------------------------------------*/
Term * build_term_from_stack(bool sort) {
  if (sort) std::sort(vstack.begin(), vstack.end(), cmpVarLvl);
  Term * res = 0;
  while (!vstack.empty()) {
    Term * t = new_term(vstack.back(), res);
    assert(t);

    vstack.pop_back();
    deallocate_term(res);
    res = t;
  }
  return res;
}

/*------------------------------------------------------------------------*/

Term * sort_and_build_term_from_vector(std::vector< Var*> v) {
  std::sort(v.begin(), v.end(), cmpVarLvl);

  for (unsigned i = 0; i < v.size(); i++) {
    if(i < v.size()-1 && v[i] == v[i+1]->get_dual()) i+=2;
    else add_to_vstack(v[i]);
  }

  Term * t = build_term_from_stack(0);
  return t;

}
/*------------------------------------------------------------------------*/
Term * gen_dual_term_from_vector(std::vector<Term*> t){
  std::vector< Var*> tmp;
  for(auto i = t.begin(); i != t.end(); ++i){
    Term * tmp_t = *i;
    tmp.push_back(tmp_t->get_var()->get_dual());
  }
  return sort_and_build_term_from_vector(tmp);
}
/*------------------------------------------------------------------------*/
Term * divide_by_var(Term *t, Var*v) {
  while(t){
    if(t->get_var() != v) add_to_vstack(t->get_var()); 
    t = t->get_rest();
  }
  Term * res = build_term_from_stack();
  return res;
}
/*------------------------------------------------------------------------*/
Term * divide_by_term(Term *t, Term *t1) {
  while(t && t1){
    if(t->get_var() != t1->get_var()) add_to_vstack(t->get_var()); 
    else t1 = t1->get_rest();
    t = t->get_rest();
  }
  while(t){
    add_to_vstack(t->get_var());
    t = t->get_rest();
  }
  Term * res = build_term_from_stack();
  return res;
}
/*------------------------------------------------------------------------*/

Term * reorder_term_lex(Term * t1) {
  while(t1){
    vstack.push_back(t1->get_var());
    t1 = t1->get_rest();
  }
  std::sort(vstack.begin(), vstack.end(), cmpVarLvl);
  return build_term_from_stack();
}

/*------------------------------------------------------------------------*/

Term * multiply_term(Term * t1, Term * t2) {
  if (!t1 || !t2) return 0;
  if (t1 == t2) return t1->copy();

  const Term * tmp1 = t1;
  const Term * tmp2 = t2;


  while (tmp1 && tmp2) {
    if (tmp1->get_var_level() > tmp2->get_var_level()) {
      vstack.push_back(tmp1->get_var());
      tmp1 = tmp1->get_rest();
    } else if (tmp1->get_var_level() < tmp2->get_var_level()) {
      vstack.push_back(tmp2->get_var());
      tmp2 = tmp2->get_rest();
    } else {
      vstack.push_back(tmp1->get_var());
      tmp1 = tmp1->get_rest();
      tmp2 = tmp2->get_rest();
    }
  }

  while (tmp1) {
    vstack.push_back(tmp1->get_var());
    tmp1 = tmp1->get_rest();
  }

  while (tmp2) {
    vstack.push_back(tmp2->get_var());
    tmp2 = tmp2->get_rest();
  }
  Term * t = build_term_from_stack();

  return t;
}
/*------------------------------------------------------------------------*/
Term * multiply_term_by_var(Term * t1, Var * v){
  if(t1->contains(v->get_dual())) return 0;
  Term * tmp = new_term(v,0);
  Term * res = multiply_term(t1, tmp);
  //delete(tmp);
  return res;
}

/*------------------------------------------------------------------------*/
Var * first_dual_different_var(Term * t0, Term* t1){

  while (t0 && t1) {

      if (t0->get_var() == t1->get_var()){
        t0 = t0->get_rest();
        t1 = t1->get_rest();
      }
      else if (t0->get_var() == t1->get_var()->get_dual()) return t0->get_var();
      else return 0;
  }  
  return 0;
}

Term * remainder(const Term * t, Var * v) {
  assert(v);


  while (t) {
    Var * v1 = t->get_var();

    if (v1 == v) {
      t = t->get_rest();
    } else {
      vstack.push_back(v1);
      t = t->get_rest();
    }
  }

  Term * res = build_term_from_stack();
  return res;
}

/*--------------------------------------------------------*/
Term * remainder_t(const Term * t, const Term * t1) {
assert(t1);

  while (t) {
    if(!t1 || t->get_var() != t1->get_var()) {
      vstack.push_back(t->get_var());
      t = t->get_rest();
    } else {
      t = t->get_rest();
      t1 = t1->get_rest();
    }
  }

  Term * res = build_term_from_stack();

  return res;
}