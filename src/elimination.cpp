/*------------------------------------------------------------------------*/
/*! \file elimination.cpp
    \brief contains functions used in the polynomial solver

  This file contains all functions used for preprocessing the Gröbner basis
  and for reducing the specification by the slices.

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#include <algorithm>
#include <list>

#include "variable.h"
#include "elimination.h"
/*------------------------------------------------------------------------*/
// Global variables
int proof = 0;

/*------------------------------------------------------------------------*/
static Polynomial *unfold_linear_poly(Polynomial *p)
{
  assert(p->degree() == 1);
  Polynomial *tmp_p = p->copy();

  for (unsigned j = 0; j < p->size(); j++)
  {
    Term *t = p->get_mon(j)->get_term();

    if (t && t->get_var()->is_dual())
    {
      Polynomial *tmp = flip_var_in_poly(tmp_p, t->get_var());
      delete (tmp_p);
      tmp_p = tmp;
    }
  }
  return tmp_p;
}

/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
static bool rewrite_l_f_pairs(Gate *g, Gate *g_c_p)
{

  if (!g->equal_children(g_c_p))
    return 0;

  Polynomial *g_poly = g->get_gate_constraint();
  if (g_poly->size() != 2)
    return 0;

  Polynomial *g_c_p_poly = g_c_p->get_gate_constraint();
  if (g_c_p_poly->degree() == 1)
    return 0;
  if (g_c_p_poly->size() != 2)
    return 0;

  Polynomial *g_c_p_poly_tmp = g_c_p_poly->copy();

  Term *m_tail = g_poly->get_tail_term();
  Term *n_tail = g_c_p_poly->get_tail_term();
  assert(m_tail->degree() == n_tail->degree());
  bool flipped = 0;
  while (m_tail)
  {
    Var *m_var = m_tail->get_var();
    Var *n_var = n_tail->get_var();

    if (m_var != n_var)
    {
      Polynomial *tmp = flip_var_in_poly(g_c_p_poly_tmp, n_var);
      delete (g_c_p_poly_tmp);
      g_c_p_poly_tmp = tmp;
      flipped = !flipped;
    }

    m_tail = m_tail->get_rest();
    n_tail = n_tail->get_rest();
  }

  Polynomial *res = flipped ? add_poly(g_c_p_poly_tmp, g_poly) : sub_poly(g_c_p_poly_tmp, g_poly);

  delete (g_c_p_poly_tmp);

  Polynomial *res_unfold = unfold_linear_poly(res);
  delete (res);

  Gate *big = g_c_p->get_var_level() > g->get_var_level() ? g_c_p : g;
  Gate *small = (big == g_c_p) ? g : g_c_p;

  big->update_gate_poly(res_unfold);

  big->van_twins_push_back(small);
  if (verbose >= 4)
    msg("found van pair %s %s", big->get_var_name(), small->get_var_name());
  return 1;
}

/*------------------------------------------------------------------------*/

static Gate *rewrite_l_f(Gate *g)
{
  if (verbose >= 3)
    msg("rewrite l f %s", g->get_var_name());


  auto g_update_children = get_var_of_poly(g);
  g->set_children(g_update_children);
  Gate *g_c = *(g->children_begin());

  for (auto it = g_c->parents_begin(); it != g_c->parents_end(); it++)
  {
    Gate *g_c_p = *it;
    if (g_c_p == g)
      continue;

    if (rewrite_l_f_pairs(g, g_c_p))
      return g_c_p;
  }
  return 0;
}


/*------------------------------------------------------------------------*/
// tailored for -p+x*y
static void remove_vanishing_mon_single_gate(Gate *g)
{
  Polynomial *p = g->get_gate_constraint();
  if (p->degree() != 2)
    return;
  Term *p_tail = p->get_tail_term();

  Gate *g0 = gate(p_tail->get_var_num());
  Gate *g1 = gate(p_tail->get_rest()->get_var_num());

  if (!g0->is_van_twin(g1))
    return;

  p = p->copy();

  if (p_tail->get_var()->is_dual())
  {
    Polynomial *tmp = flip_var_in_poly(p, p_tail->get_var());
    delete (p);
    p = tmp;
  }

  if (p_tail->get_rest()->get_var()->is_dual())
  {
    Polynomial *tmp = flip_var_in_poly(p, p_tail->get_rest()->get_var());
    delete (p);
    p = tmp;
  }

  Term *t0 = new_term(g1->get_var(), 0);
  Term *t1 = new_term(g0->get_var(), t0);
  Monomial *red = new Monomial(p->get_mon(1)->coeff, t1);
  push_mstack_end(red);
  Polynomial *sub = build_poly();
  Polynomial *sub_neg = multiply_poly_with_constant(sub, minus_one);
  Polynomial *res = add_poly(p, sub_neg);
  delete (sub);
  delete (sub_neg);
  delete (p);

  g->update_gate_poly(res);
}
/*------------------------------------------------------------------------*/
// CHILDREN ARE L_F_PAIRS
/*------------------------------------------------------------------------*/
static void children_l_f_pairs(Gate *g)
{
  if (verbose >= 3)
    msg("child l f pairs %s", g->get_var_name());
  Polynomial *p = g->get_gate_constraint();
  if (p->degree() != 2 || p->size() != 2)
    return;
  Gate *big_child = gate(p->get_tail_term()->get_var_num());
  if (big_child->get_input())
    return;
  Gate *little_child = gate(p->get_tail_term()->get_rest()->get_var_num());
  if (little_child->get_input())
    return;
  if (rewrite_l_f_pairs(big_child, little_child))
    remove_vanishing_mon_single_gate(g);
}
/*------------------------------------------------------------------------*/
// CHILDREN share L_F
/*------------------------------------------------------------------------*/
static bool share_l_f_child(Gate *g, Gate *q)
{
  if (g->get_input() || q->get_input())
    return 0;
  if (g->get_gate_constraint()->degree() == 1)
    return 0;
  if (q->get_gate_constraint()->degree() == 1)
    return 0;
  if (g->get_var_level() < q->get_var_level())
    std::swap(g, q);
  Term *gt = g->get_gate_constraint()->get_tail_term();
  Term *qt = q->get_gate_constraint()->get_tail_term();

  while (gt && qt)
  {
    Var *vg = gt->get_var();
    Var *vq = qt->get_var();
    if (vg->get_dual() == vq)
    {
      g->van_twins_push_back(q);
      return 1;
    }
    if (!vg->is_dual() && share_l_f_child(q, gate(vg->get_num())))
    {
      g->van_twins_push_back(q);
      return 1;
    }
    if (!vq->is_dual() && share_l_f_child(g, gate(vq->get_num())))
    {
      g->van_twins_push_back(q);
      return 1;
    }

    if (vg->get_num() > vq->get_num())
      gt = gt->get_rest();
    else
      qt = qt->get_rest();
  }
  return 0;
}
/*------------------------------------------------------------------------*/
static void children_share_l_f(Gate *g)
{
  if (verbose >= 3)
    msg("child share l_f %s", g->get_var_name());
  Polynomial *p = g->get_gate_constraint();
  if (p->degree() != 2 || p->size() != 2)
    return;
  Gate *big_child = gate(p->get_tail_term()->get_var_num());
  if (big_child->get_input())
    return;
  Gate *little_child = gate(p->get_tail_term()->get_rest()->get_var_num());
  if (little_child->get_input())
    return;
  if (share_l_f_child(big_child, little_child))
    remove_vanishing_mon_single_gate(g);
}

/*------------------------------------------------------------------------*/
// VAN MON DEPTH
/*------------------------------------------------------------------------*/
static void filter_van_mon_pairs(Gate *g, unsigned l)
{
  if (aiger_sign(l))
  {
    Gate *lg = gate(l - 1);
    g->van_twins_push_back(lg);
    if (verbose >= 4)
      msg("van twin %s %s", g->get_var_name(), lg->get_var_name());
    return;
  }

  Gate *lg = gate(l);
  if (!lg->get_input())
  {
    aiger_and *and1 = is_model_and(lg->get_var_num());
    unsigned ll = and1->rhs0, rr = and1->rhs1;

    filter_van_mon_pairs(g, ll);
    filter_van_mon_pairs(g, rr);
  }
}

/*------------------------------------------------------------------------*/
static void filter_van_mon_depth(Gate *g)
{
  if (verbose >= 3)
    msg("filter van mon depth");
  Polynomial *p = g->get_gate_constraint();
  Gate *big_child = gate(p->get_tail_term()->get_var_num());
  if (big_child->get_input())
    return;
  aiger_and *and1 = is_model_and(big_child->get_var_num());
  unsigned l = and1->rhs0, r = and1->rhs1;

  filter_van_mon_pairs(big_child, l);
  filter_van_mon_pairs(big_child, r);
  remove_vanishing_mon_single_gate(g);
}

/*------------------------------------------------------------------------*/
// GRANDCHILD ARE EQUAL
/*------------------------------------------------------------------------*/
static bool is_among_dual_children(Gate *g, Gate *q)
{
  if (verbose >= 3)
    msg("searching for %s", q->get_var_name());
  Polynomial *p = g->get_gate_constraint();
  Term *t = p->get_tail_term();
  while (t)
  {
    Var *v = t->get_var();
    if (v->is_dual() && v->get_dual() == q->get_var())
      return 1;
    else if (!v->is_dual())
      return is_among_dual_children(gate(v->get_num()), q);
    t = t->get_rest();
  }

  return 0;
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
// Equivalent gate

/*------------------------------------------------------------------------*/
static void local_rewriting(Gate *g, mpz_t coeff)
{
  
  if(verbose > 1) msg("locally rewrite %s", g->get_var_name());
  rewrite_l_f(g);
  if (g->get_gate_constraint()->degree() == 1)
  {
    l_f_count++;
    totalcount++;
    if (verbose >= 3)
      g->get_gate_constraint()->print(stdout);
    return;
  }

  children_l_f_pairs(g);
  if (g->get_gate_constraint()->degree() == 1)
  {
    child_l_f_count++;
    totalcount++;
    if (verbose >= 3)
      g->get_gate_constraint()->print(stdout);
    return;
  }

  children_share_l_f(g);
  if (g->get_gate_constraint()->degree() == 1)
  {
    children_share_l_f_count++;
    totalcount++;
    if (verbose >= 3)
      g->get_gate_constraint()->print(stdout);
    return;
  }



  filter_van_mon_depth(g);
  if (g->get_gate_constraint()->degree() == 1)
  {
    van_mon_depth_count++;
    totalcount++;
    if (verbose >= 3)
      g->get_gate_constraint()->print(stdout);
    return;
  }

  linearize_via_gb(g, 3, 0, coeff,0);

  if (!g->get_gate_constraint())
    die(2, "g lost gate constraint");
  if (g->get_gate_constraint()->degree() == 1)
  {
    lin_gb_count++;
    totalcount++;
    return;
  }
  if (verbose >= 3)
    g->get_gate_constraint()->print(stdout);
}
/*------------------------------------------------------------------------*/
static void replace_p1_by_p2(Gate *g1, Gate *g2)
{
  if (g1->get_var_level() < g2->get_var_level())
    std::swap(g1, g2);
  if (verbose > 1)
    msg("need to implement replacing of gate %s instead of %s", g2->get_var_name(), g1->get_var_name());
  g2->mark_rew();

  g1->set_rep_var(g2->get_var());
  Var *v = 0, *vd = 0;

  mpz_t coeff;
  mpz_init(coeff);

  for (auto &g1_p : g1->get_parents())
  {
    Polynomial *g1_gc = g1_p->get_gate_constraint();

    assert(g1_gc->size() == 2);
    Term *t = g1_gc->get_tail_term();
    Monomial *m = g1_gc->get_mon(1);
    v = t->contains(g1->get_var()) ? g1->get_var() : g1->get_var()->get_dual();
    vd = v->is_dual() ? g2->get_var()->get_dual() : g2->get_var();

    Term *t1 = divide_by_var(t, v);

    Term *t2 = multiply_term_by_var(t1->copy(), vd);

    mpz_neg(coeff, m->coeff);

    Monomial *m_neg = new Monomial(coeff, m->get_term_copy());
    Monomial *m_new = new Monomial(m->coeff, t2);
    push_mstack(m_neg);
    push_mstack(m_new);
    Polynomial *rep = build_poly();

    Polynomial *res = add_poly(g1_gc, rep);
    delete (rep);
    delete (g1_gc);
    g1_p->set_gate_constraint(res);
  }

  for (auto &g1_p : g1->get_parents())
  {
    g1_p->children_remove(g1);
    g1_p->children_push_back(g2);
    g1->parents_remove(g1_p);
    g2->parents_push_back(g1_p);
  }

  std::list<Gate *> red;
  mpz_clear(coeff);
  for (auto it = g2->parents_begin(); it != g2->parents_end(); it++)
  {

    Gate *outer = *it;

    if (!outer)
      continue;
    if (outer->get_gate_constraint()->size() > 2)
      continue;
    for (auto it2 = it; it2 != g2->parents_end(); it2++)
    {

      Gate *inner = *it2;
      if (outer == inner)
        continue;
      if (inner->get_gate_constraint()->size() > 2)
        continue;

      if (outer->get_gate_constraint()->get_tail_term() == inner->get_gate_constraint()->get_tail_term())
      {
        replace_p1_by_p2(outer, inner);
        if (outer->get_var_level() < inner->get_var_level())
          red.push_back(inner);
        else
          red.push_back(outer);
      }
    }
  }

  while (!red.empty())
  {
    Gate *redg = red.front();
    red.pop_front();
    for (auto &g1_c : redg->get_children())
    {
      g1_c->parents_remove(redg);
    }
  }
}


/*------------------------------------------------------------------------*/
static Var *get_sibling(Gate *gp, Gate *g)
{
  if (gp->children_size() != 2)
    return 0;

  Gate *gc = (g == gp->children_front()) ? gp->children_back() : gp->children_front();
  if (gp->get_gate_constraint()->get_tail_term()->contains(gc->get_var()))
    return gc->get_var();
  else
    return gc->get_var()->get_dual();
}

/*------------------------------------------------------------------------*/

static int try_rewrite_single_pos(Gate *g)
{
  Gate *gp = g->parents_front();
  Var *g_sibling = get_sibling(gp, g);

  Polynomial *p = gp->get_gate_constraint();
  Polynomial *g_poly = g->get_gate_constraint();
  Polynomial *tmp_poly = reduce_by_one_poly(p, g->get_gate_constraint());

  Var *v0 = g_poly->get_tail_term()->get_var();
  Term *t0 = new_quadratic_term(v0, g_sibling);

  Var *v1 = g_poly->get_tail_term()->get_rest()->get_var();
  Term *t1 = new_quadratic_term(v1, g_sibling);

  Gate *gp_q = search_for_parent(t0);
  Term *rem = divide_by_term(tmp_poly->get_tail_term(), t0);
  Var *rep_var = 0;
  if (gp_q && gp_q->get_var_level() > g->get_var_level() && gp_q->get_rep_var())
  {
    rep_var = gp_q->get_rep_var();
  }

  if (!gp_q || (gp_q->get_var_level() > g->get_var_level() && !rep_var))
  {
    gp_q = search_for_parent(t1);
    deallocate_term(rem);
    rem = divide_by_term(tmp_poly->get_tail_term(), t1);
  }
  if (gp_q && gp_q->get_var_level() > g->get_var_level() && gp_q->get_rep_var())
  {
    rep_var = gp_q->get_rep_var();
  }

  if (!gp_q || (gp_q->get_var_level() > g->get_var_level() && !rep_var))
    return 0;
  assert(rem->degree() == 1);

  //  msg("found %s", gp_q->get_var_name());

  Polynomial *multi = multiply_poly_with_term(gp_q->get_gate_constraint(), rem);
  Polynomial *res = sub_poly(tmp_poly, multi);

  delete (tmp_poly);
  tmp_poly = res;

  if (rep_var)
  {
    push_mstack(gp_q->get_gate_constraint()->get_lm()->copy());
    Term *t = new_term(rep_var, 0);
    push_mstack(new Monomial(one, t));
    Polynomial *rep_v = build_poly();

    Polynomial *multi2 = multiply_poly_with_term(rep_v, rem);

    Polynomial *res2 = add_poly(tmp_poly, multi2);
    delete (tmp_poly);
    delete (multi2);
    tmp_poly = res2;
  }

  delete (multi);

  delete (p);
  gp->set_gate_constraint(tmp_poly);
  if (verbose >= 2)
  {
    msg("single pos rewriting changed %s %s", gp->get_var_name(), g->get_var_name());
    tmp_poly->print(stdout);
  }

  for (auto &gc : gp->get_children())
  {
    gc->parents_remove(gp);
    gp->children_remove(gc);
  }
  
  for (auto &gc : g->get_children())
  {
    gc->parents_remove(g);
    g->children_remove(gc);
  }

  if (rep_var)
    gp_q = gate(rep_var->get_num());

  gp_q->parents_push_back(gp);
  gp->children_push_back(gp_q);
  Gate * gp_equiv = 0;
   Gate *gp_q_sib = 0;
  if(rem) {
    gp_q_sib = gate(rem->get_var_num());
    gp_q_sib->parents_push_back(gp);
    gp->children_push_back(gp_q_sib);

    Term *gp_tail = gp->get_gate_constraint()->get_tail_term();
    gp_equiv = search_for_parent(gp_tail, gp);
  }
  
  if (gp_equiv)
  {
    replace_p1_by_p2(gp, gp_equiv);
  }

  if (g->parents_size() == 0)
  {
    g->mark_red();
    if(verbose > 2) msg("eliminated %s", g->get_var_name());
  }

  if (gp_q_sib) return std::min(gp_q->get_var_level(), gp_q_sib->get_var_level());
  else return gp_q->get_var_level();
}
/*------------------------------------------------------------------------*/
static int rewrite_single_pos(unsigned lowerbound)
{
  msg("rewrite single pos");
  int min = 0;
  int count = 0;
  for (unsigned i = lowerbound; i < num_gates; i++)
  {
    Gate *g = gates[i];
    if (g->parents_size() != 1)
      continue;
    if (g->get_input())
      continue;
    if (g->get_pp())
      continue;
    if (g->get_aig_output())
      continue;

    Gate *gp = g->parents_front();
    Term *t = gp->get_gate_constraint()->get_tail_term();

    if (!t->contains(g->get_var()))
      continue; // This means t contains dual of g

    int flag = try_rewrite_single_pos(g);
    if (flag) count++;
     if (flag && !min)
      min = i;
  }

  msg("removed %i single pos", count);
  return min;
}

/*------------------------------------------------------------------------*/
static bool all_children_are_singles(Gate *g)
{
  for (auto &gc : g->get_children())
  {
    if (gc->parents_size() > 1)
      return 0;
  }
  return 1;
}
/*------------------------------------------------------------------------*/
static void replace_p1_by_p2_dual(Gate *g1, Gate *g2)
{
  bool flag = 0;
  if (all_children_are_singles(g2) && !all_children_are_singles(g1))
  {
    std::swap(g1, g2);
    flag = 1;
  }
  if (g1->parents_size() < 1)
    return;
  g2->mark_rew();

  if (verbose > 1)
    msg("need to implement replacing of dual gate %s instead of %s", g2->get_var_name(), g1->get_var_name());
  g1->set_rep_var(g2->get_var()->get_dual());
  Var *v = 0, *vd = 0;

  mpz_t coeff;
  mpz_init(coeff);

  for (auto &g1_p : g1->get_parents())
  {
    Polynomial *g1_gc = g1_p->get_gate_constraint();
    g1_gc->print(stdout);
    assert(g1_gc->size() == 2);
    Term *t = g1_gc->get_tail_term();
    Monomial *m = g1_gc->get_mon(1);
    v = t->contains(g1->get_var()) ? g1->get_var() : g1->get_var()->get_dual();
    vd = v->is_dual() ? g2->get_var() : g2->get_var()->get_dual();

    Term *t1 = divide_by_var(t, v);
    Term *t2 = multiply_term_by_var(t1->copy(), vd);

    mpz_neg(coeff, m->coeff);

    Monomial *m_neg = new Monomial(coeff, m->get_term_copy());
    Monomial *m_new = new Monomial(m->coeff, t2);
    push_mstack(m_neg);
    push_mstack(m_new);
    Polynomial *rep = build_poly();
    Polynomial *res = add_poly(g1_gc, rep);
    res->print(stdout);

    delete (rep);
    // delete(g1_gc);
    g1_p->update_gate_poly(res);
  }

  for (auto &g1_p : g1->get_parents())
  {
    g1_p->children_remove(g1);
    g1_p->children_push_back(g2);

    g1->parents_remove(g1_p);
    g2->parents_push_back(g1_p);
  }

  if (flag)
  {
    int level = g1->get_var_level();
    g1->set_var_level(g2->get_var_level());
    g2->set_var_level(level);
    msg("swapped levels of %s and %s", g2->get_var_name(), g1->get_var_name());

    Gate *g1_c1 = g1->children_front();
    Gate *g1_c2 = g1->children_back();

    Gate *g2_c1 = g2->children_front();
    Gate *g2_c2 = g2->children_back();

    level = g1_c1->get_var_level();
    g1_c1->set_var_level(g2_c1->get_var_level());
    g2_c1->set_var_level(level);
    msg("swapped levels of %s and %s", g2_c1->get_var_name(), g1_c1->get_var_name());

    level = g1_c2->get_var_level();
    g1_c2->set_var_level(g2_c2->get_var_level());
    g2_c2->set_var_level(level);
    msg("swapped levels of %s and %s", g2_c2->get_var_name(), g1_c2->get_var_name());

    resort_slices_lex();
  }

  mpz_clear(coeff);
  std::list<Gate *> red;
  for (auto it = g2->parents_begin(); it != --(g2->parents_end()); it++)
  {
    Gate *outer = *it;
    if (!outer)
      continue;
    for (auto it2 = it; it2 != g2->parents_end(); it2++)
    {
      Gate *inner = *it2;
      if (outer == inner)
        continue;
      if (outer->get_gate_constraint()->get_tail_term() == inner->get_gate_constraint()->get_tail_term())
      {
        replace_p1_by_p2(outer, inner);
        if (outer->get_var_level() < inner->get_var_level())
          red.push_back(inner);
        else
          red.push_back(outer);
      }
    }
  }

  while (!red.empty())
  {
    Gate *redg = red.front();
    red.pop_front();
    for (auto &g1_c : redg->get_children())
    {
      g1_c->parents_remove(redg);
    }
  }
}


/*------------------------------------------------------------------------*/
static void clean_unused_poly()
{
  msg("call clean unused");
  for (unsigned i = M - 2; i >= NN; i--)
  {
    Gate *g = gates[i];
    if (g->parents_size() == 0 && !g->get_red())
    {
      g->mark_red();
      for (auto &gc : g->get_children())
      {
        gc->parents_remove(g);
      }
    }
    else
    {
      for (auto &gp : g->get_parents())
      {
        if (gp->parents_size() == 0)
          gp->mark_red();
        // else msg("parent of %s is %s", g->get_var_name(), gp->get_var_name());
      }
    }
  }
}

static Gate *common_parent(Gate *g1, Gate *g2)
{
  for (auto &gp : g1->get_parents())
  {
    if (gp->is_child(g2))
      return gp;
  }
  return 0;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
static Polynomial *check_for_duality(Gate *g1, Gate *g2)
{
  Polynomial *p1 = g1->get_gate_constraint();
  Polynomial *p1_pos = unflip_poly(p1);

  for (auto &g1_c : g1->get_children())
  {
    Polynomial *tmp = reduce_by_one_poly(p1_pos, g1_c->get_gate_constraint());
    delete (p1_pos);
    p1_pos = tmp;
  }

  Polynomial *p1_unfolded_pos = unflip_poly(p1_pos);
  delete (p1_pos);

  Polynomial *p2 = g2->get_gate_constraint();
  Polynomial *p2_pos = unflip_poly(p2);

  for (auto &g2_c : g2->get_children())
  {
    Polynomial *tmp = reduce_by_one_poly(p2_pos, g2_c->get_gate_constraint());
    delete (p2_pos);
    p2_pos = tmp;
  }

  Polynomial *p2_unfolded_pos = unflip_poly(p2_pos);
  delete (p2_pos);

  Polynomial *add = add_poly(p1_unfolded_pos, p2_unfolded_pos);
  delete (p1_unfolded_pos);
  delete (p2_unfolded_pos);

  if (add->size() != 3 || add->degree() != 1)
    return 0;
  if (add->get_mon(0)->get_term()->get_var() != g1->get_var())
    return 0;
  if (mpz_cmp_si(add->get_mon(0)->coeff, -1) != 0)
    return 0;
  if (add->get_mon(1)->get_term()->get_var() != g2->get_var())
    return 0;
  if (mpz_cmp_si(add->get_mon(1)->coeff, -1) != 0)
    return 0;
  if (add->get_mon(2)->get_term())
    return 0;
  if (mpz_cmp_si(add->get_mon(2)->coeff, 1) != 0)
    return 0;

  return add;
}
/*------------------------------------------------------------------------*/

static void try_regrouping_of_children(Gate *g)
{

  Polynomial *p = g->get_gate_constraint();
  assert(p->size() == 2);
  assert(g->children_size() == 2);
  Gate *g0 = g->children_front();
  Gate *g1 = g->children_back();
  if (g0->get_var_level() < g1->get_var_level())
    std::swap(g0, g1);

  for (auto &g0_p : g0->get_parents())
  {
    if (g0_p->get_var_level() < g->get_var_level())
      continue;
    for (auto &g1_p : g1->get_parents())
    {
      if (g0_p->is_child(g1_p))
      {

        Polynomial *g0_poly = g0_p->get_gate_constraint();
        Polynomial *g1_poly = g1_p->get_gate_constraint();
        assert(g0_poly->size() == 2 && g0_poly->degree() == 2);
        assert(g1_poly->size() == 2 && g1_poly->degree() == 2);
        if (!g0_poly->get_tail_term()->contains(g1_p->get_var()))
          return;
        if (!g0_poly->get_tail_term()->contains(p->get_tail_term()->get_var()))
          return;
        if (!g1_poly->get_tail_term()->contains(p->get_tail_term()->get_rest()->get_var()))
          return;

        Term *tmp = divide_by_var(g1_poly->get_tail_term(), p->get_tail_term()->get_rest()->get_var());

        assert(tmp->degree() == 1);
        Gate *rest = gate(tmp->get_var_num());

        Term *tmp2 = multiply_term_by_var(tmp, g->get_var());
        Monomial *m = new Monomial(minus_one, tmp2);
        push_mstack(m);
        push_mstack(g0_poly->get_mon(1)->copy());
        Polynomial *res = build_poly();

        Polynomial *add = sub_poly(g0_poly, res);

        delete (g0_poly);
        g0_p->set_gate_constraint(add);

        if (verbose > 2)
        {
          msg("new gate constraint for %s", g0_p->get_var_name());
          g0_p->get_gate_constraint()->print(stdout);
        }

        g0_p->children_remove(g0);
        g0_p->children_remove(g1);
        g0_p->children_push_back(g);
        g0_p->children_push_back(rest);
        g0->parents_remove(g0_p);
        g1->parents_remove(g1_p);
        g->parents_push_back(g0_p);
        rest->parents_push_back(g0_p);
      }
    }
  }
}
/*------------------------------------------------------------------------*/

static bool check_for_double_xor(Gate *g)
{
  if (verbose > 1)
    msg("checking for double xor %s", g->get_var_name());
  bool flag = 0;
  std::list<Gate *> nodes;
  if (g->parents_size() == 4)
  {
    nodes = g->get_parents();
    Gate *g_p0 = nodes.front();
    for (auto &g_p : nodes)
    {

      if (!g_p->equal_children(g_p0))
        return 0;
    }
  }
  else
  {

    for (auto it = g->parents_begin(); it != --(g->parents_end()); it++)
    {
      nodes.clear();
      Gate *outer = *it;

      for (auto it2 = it; it2 != g->parents_end(); it2++)
      {
        Gate *inner = *it2;

        if (outer->equal_children(inner))
        {

          nodes.push_back(inner);
        }
      }
      if (nodes.size() == 4)
        break;
    }
  }
  if (nodes.size() != 4)
    return 0;

  // find grandparents
  Gate *p1 = 0, *p2 = 0;

  Gate *g0 = nodes.front();
  nodes.pop_front();
  Gate *g1 = nodes.front();
  nodes.pop_front();
  Gate *g2 = nodes.front();
  nodes.pop_front();
  Gate *g3 = nodes.front();
  nodes.pop_front();
  assert(nodes.empty());

  for (auto &g0_p : g0->get_parents())
  {
    if (g0_p->is_child(g1))
    {
      p1 = g0_p;
      p2 = common_parent(g2, g3);
      break;
    }

    if (g0_p->is_child(g2))
    {
      p1 = g0_p;
      p2 = common_parent(g1, g3);
      break;
    }

    if (g0_p->is_child(g3))
    {
      p1 = g0_p;
      p2 = common_parent(g1, g2);
      break;
    }
  }
  if (!p1 || !p2)
    return 0;
  if (p1 == p2)
    return 0;
  if (p2->parents_size() == 0)
    return 0;
  if (p1->get_var_level() < p2->get_var_level())
    std::swap(p1, p2);

  try_regrouping_of_children(p1);

  Polynomial *duality_poly = check_for_duality(p1, p2);
  if (duality_poly)
  {
    if (p2->get_rep_var())
      return 0;
    replace_p1_by_p2_dual(p1, p2);
    flag = 1;
  }

  return flag;
}

/*------------------------------------------------------------------------*/


void preprocessing_elimination()
{

  // rewrite_non_pp();
  rewrite_single_pos(0);

  // merge_gates_with_equal_children(); // buggy and causes non integer coeffs
  int singleposmin = 1;
  for (unsigned i = singleposmin - 1; i < num_gates; i++)
  {
    Gate *g = gates[i];
    if (g->get_red())
      continue;
    if (g->get_input())
      continue; // TODO maybe not necessary
    if (g->parents_size() >= 4)
      check_for_double_xor(g);
    // if (g->parents_size() == 3)  check_for_double_triple(g);
  }

  dual_preprocessing();
  clean_unused_poly();

  if (verbose > 3) print_slices();
}



/*------------------------------------------------------------------------*/
Polynomial *linearize_spec(Polynomial *spec)
{
  msg("started reducing non linear terms in spec");
  if (spec->degree() == 1)
    return spec;

  Polynomial *rem = spec;
  std::vector<size_t> indices;
  bool enlarged = 0;
  for (size_t i = 0; i < rem->size(); i++)
  {
    Monomial *m = rem->get_mon(i);
    Term *t = m->get_term();
    if (t->degree() == 1)
    {
      push_mstack(m->copy());
    }
    else if (t->get_ref() > 1)
    {
      Gate *g = gate(t->get_var_num());

      for (auto &sub : g->get_parents())
      {

        if (sub->get_red())
          continue;
        Polynomial *sub_gc = sub->get_gate_constraint();
        if (sub_gc->size() != 2)
          continue;

        Monomial *inner_m = sub_gc->get_mon(1);
        Term *inner_t = inner_m->get_term();
        if (t != inner_t)
          continue;
        Term *lt = sub_gc->get_lt()->copy();
        Monomial *tmp = new Monomial(m->coeff, lt);
        push_mstack(tmp);
        break;
      }
    }
    else
    {
      
      if (!enlarged)
      {
        enlarge_gates(rem->size());
        enlarged = 1;
      }
      Term *rep_t = extend_var_gates(t);
      Monomial *tmp = new Monomial(m->coeff, rep_t->copy());
      push_mstack(tmp);
    }
  }
  Polynomial *tmp = build_poly();
  if(verbose >= 2) tmp->print(stdout);
  adjust_level_of_extended_gates();


  return tmp;
}
/*------------------------------------------------------------------------*/
static Polynomial *mod_poly(const Polynomial *p1)
{

  int exp = NN;

  mpz_t coeff;
  mpz_init(coeff);

  for (size_t i = 0; i < p1->size(); i++)
  {
    Monomial *m = p1->get_mon(i);
    mpz_tdiv_r_2exp(coeff, m->coeff, exp);
    if (mpz_sgn(coeff) != 0)
    {
      Monomial *tmp;
      if (m->get_term())
        tmp = new Monomial(coeff, m->get_term_copy());
      else
        tmp = new Monomial(coeff, 0);
      push_mstack_end(tmp);
    }
  }
  mpz_clear(coeff);
  Polynomial *out = build_poly();

  delete (p1);
  return out;
}
/*------------------------------------------------------------------------*/

Polynomial *reduce(Polynomial *spec)
{
  preprocessing_elimination();
  msg("");
  msg("");
  msg("started reducing");

  Polynomial *rem = spec, *tmp = spec;
  if (verbose > 2)
    rem->print(stdout);

  while (!gate(rem->get_lt()->get_var_num())->get_input())
  {
    Gate *g = gate(rem->get_lt()->get_var_num());
    Polynomial *n = g->get_gate_constraint();

    if (g->get_gate_constraint()->degree() > 1)
    {
      local_rewriting(g, rem->get_mon(0)->coeff);
      n = g->get_gate_constraint();
    }
    if (n->degree() > 1 && g->get_var()->get_dist() <= 6)
    {
      msg("distance close to inputs, switching to non-linear rewriting");
      break;
    }

    if (n->degree() > 1 && g->get_var()->get_dist() > 6)
    {
      msg("failed to linearize gate poly");
      g->get_gate_constraint()->print(stdout);
      msg("remainder is");
      rem->print(stdout);

      die(32, "lin gb failure dist %i", g->get_var()->get_dist());
    }

    if (verbose >= 2 && n)
    {
      fputs("[mltlng] reducing by ", stdout);
      n->print(stdout);
    }

    tmp = substitute_linear_poly(rem, n);
    linear_count++;

    if (tmp)
      tmp = mod_poly(tmp);
    g->mark_red();

    if (tmp != rem)
    {
      delete (rem);
      rem = tmp;
    }

    if (!rem)
    {
      msg("remainder is 0");
      return rem;
    }

    if (verbose >= 3)
    {
      fputs("[mltlng] remainder is ", stdout);
      rem->print(stdout);
      msg(" ");
    }
  }

  while (!gate(rem->get_lt()->get_var_num())->get_input())
  {
    Gate *g = gate(rem->get_lt()->get_var_num());
    Polynomial *n = g->get_gate_constraint();

    Polynomial *n_unflip = unflip_poly(n);
    n = n_unflip;

    if (verbose >= 2 && n)
    {
      fputs("[mltlng] reducing by ", stdout);
      n->print(stdout);
    }

    tmp = reduce_by_one_poly(rem, n);
    non_linear_count++;

    if (tmp)
      tmp = mod_poly(tmp);
    g->mark_red();

    if (tmp != rem)
    {
      delete (rem);
      rem = tmp;
    }

    if (!rem)
    {
      msg("remainder is 0");
      return rem;
    }

    if (verbose >= 3)
    {
      fputs("[mltlng] remainder is ", stdout);
      rem->print(stdout);
      msg(" ");
    }
  }

  return rem;
}

/*------------------------------------------------------------------------*/
