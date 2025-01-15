/*------------------------------------------------------------------------*/
/*! \file extern_gb.cpp
    \brief contains functions used in the polynomial solver

  This file contains all functions used for preprocessing the Gröbner basis
  and for reducing the specification by the slices.

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#include "extern_gb.h"



/*------------------------------------------------------------------------*/
struct LargerGate {
  bool operator()(const Gate* a, const Gate* b) const {
    return a->get_var_level() > b->get_var_level();
  }
};

 std::set<Gate*, LargerGate> var;
 std::set<Gate*, LargerGate> gate_poly;
/*------------------------------------------------------------------------*/
static void  add_children(Gate *g, int depth, bool pre){
  if(!depth) return;
  if(g->get_input()) return;
  
  if(pre && g->get_gate_constraint()->degree() == 1) return;
 // if(verbose >= 4) msg("linearize via gb - push_back %s", g->get_var_name());
  gate_poly.insert(g);
  var.insert(g);
  
  for(auto &gc : g->get_children()){
 //   if(verbose >= 4) msg("linearize via gb - push_back child %s", gc->get_var_name());
    if(gc->get_red()) continue;
    var.insert(gc);
  }
  
  for(auto& gc : g->get_children()){
    add_children(gc, depth-1, pre);
  }
}
/*------------------------------------------------------------------------*/
static void  add_all_remaining(Gate *g){

   for (int i = num_gates-1; i>=0; i--)
   {
    Gate *gc = gates[i];
    if(gc->get_input()) {
      var.insert(gc);
      continue;
    }
    if (gc->get_var_level() > g->get_var_level()) continue;

    gate_poly.insert(gc);
    var.insert(gc);
   }
}
/*------------------------------------------------------------------------*/
static void add_spouses(Gate *g, bool pre ){
  
  for(auto& gc : g->get_children()){
    for(auto & g_sib : gc->get_parents()){
      if(g == g_sib) continue;
      if(!pre && (g->get_var_level() < g_sib->get_var_level())) continue;
      if(g_sib->get_red()) continue;
      bool flag = 0;
      for(auto& gsib_c : g_sib->get_children()){
        
        if(!var.contains(gsib_c)) {
          flag = 1;
          break;
        }
      }
      if(flag) continue;

//       if(verbose >= 4)msg("linearize via gb add spouse- push_back %s", g_sib->get_var_name());
      gate_poly.insert(g_sib);
      var.insert(g_sib);

      for(auto& gsib_c : g_sib->get_children()){
 //       if(verbose >= 4)msg("linearize via gb add spuse - push_back child %s", gsib_c->get_var_name());
        var.insert(gsib_c);
      }
   
    }
  }

}

/*------------------------------------------------------------------------*/
static void add_parents(Gate * node, Gate *g ){
  

  for(auto& node_p : node->get_parents()){
    if(node_p->get_var_level() > g->get_var_level()) continue;
    if(node_p->get_output()) continue;
    if(node_p->get_red()) continue;

    bool flag = 0;
    for(auto& node_p_c : node_p->get_children()){
      if(node == node_p_c) continue;
      if(!var.contains(node_p_c)){
          flag = 1;
          break;
        }
      }
      if(flag) continue;
   //    if(verbose >= 4)msg("linearize via gb add parents- push_back %s", node_p->get_var_name());
      gate_poly.insert(node_p);
      var.insert(node_p);
      add_parents(node_p, g);
  }
  
}

/*------------------------------------------------------------------------*/
static void add_common_ancestors(Gate *g, bool pre ){
  
  for(auto& node : var){
    for(auto& node_p : node->get_parents()){
      if(node_p == g) continue;
     // msg("investigating %s %s", node->get_var_name(), node_p->get_var_name());
      if(!pre && (node_p->get_var_level() > g->get_var_level())) continue;
      if(node_p->get_output()) continue;
      if(node_p->get_red()) continue;
      bool flag = 0;
      for(auto& node_p_c : node_p->get_children()){
        if(node == node_p_c) continue;
        if(!var.contains(node_p_c)){
          flag = 1;
          break;
        }
      }
      if(flag) continue;
       
   //   if(verbose >= 4)msg("linearize via gb add common anc- push_back %s", node_p->get_var_name());
      //node_p->get_gate_constraint()->print(stdout);
      gate_poly.insert(node_p);
      var.insert(node_p);
      add_parents(node_p, g);
    }
  }
}


/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
Polynomial *flip_var_in_poly(const Polynomial *p1, Var *v)
{

  Polynomial *flip = gen_dual_constraint(v);
  const Polynomial *negfactor = divide_by_term(p1, flip->get_lt());
  Polynomial *rem;

  if (negfactor->is_constant_zero_poly())
  {
    rem = p1->copy();
  }
  else
  {
    Polynomial *mult = multiply_poly(negfactor, flip);
    rem = add_poly(p1, mult);
    delete (mult);
  }

  delete (flip);
  delete (negfactor);
  return rem;
}
/*------------------------------------------------------------------------*/

Polynomial *unflip_poly(Polynomial *p)
{

  Polynomial *res = p->copy();
  Var *v = res->contains_dual_var();
  while (v)
  {
    Polynomial *tmp = flip_var_in_poly(res, v);
    delete (res);
    res = tmp;
    v = res->contains_dual_var();
  }

  return res;
}



/*------------------------------------------------------------------------*/
static bool print_to_msolve(Gate *g){

  std::string vars;
  std::string poly;
  std::string red;
  std::vector<Gate *> gate_vec;

  int randomNumber = std::rand();

    // Construct the string
  std::string output = "/tmp/" + std::to_string(randomNumber) + ".ms";
  std::string respath = "/tmp/" + std::to_string(randomNumber) + ".out";
  std::string polypath = "/tmp/" + std::to_string(randomNumber) + ".line";
   

  FILE * f = fopen(output.c_str(), "w");
  if(!f) die(2, "cannot open file %s", output.c_str());
    
  for (const auto& tag : var) {
    fprintf(f, "%s", tag->get_var_name());
    if(tag != *var.rbegin()) fprintf(f,",");
  }

  fprintf(f, "\n");
  fprintf(f, "1073741827\n");
  
  for (const auto& gatep : gate_poly) {
    Polynomial * gpol = gatep->get_gate_constraint();
    Polynomial * tmp = unflip_poly(gpol);
    tmp->print(f, 0);
    fprintf(f,",\n");
    delete(tmp);
  }

  for (const auto& gatep : var) {
    fprintf(f,"-%s^2+%s", gatep->get_var_name(), gatep->get_var_name());
    if(gatep != *var.rbegin()) fprintf(f,",\n");
  }
  fclose(f);
  
  
  std::string msolvecall = "msolve -f " + output + " -g 2 "; // this is hardcoded for the artifact
  std::string add_grep = msolvecall + "| grep -m2 " + g->get_var_name() + " | tail -n1 | ";
  std::string remove_ones = add_grep + "sed 's/\\(\\^\\)1\\b//g; s/+1073741826/-1/g; ; s/+1073741825/-2/g' > " + respath;
  std::string call = remove_ones;
 
  system(call.c_str()); 
  Polynomial * target = parse_specification_polynomial(respath.c_str());
  std::string delcall = "rm " + output + " " + respath;
  system(delcall.c_str());
  if(target->degree() > 1) return 0;

  g->update_gate_poly(target);

  return 1;

}

/*------------------------------------------------------------------------*/


bool linearize_via_gb(Gate *g, int depth, bool pre, mpz_t coeff, bool full){
       var.clear();
  gate_poly.clear();
  if(verbose >= 3) msg("calling GB in msolve for %s", g->get_var_name());

  int max_depth = g->get_var()->get_dist();
  bool prep = 0; // REMOVE THIS IF PRE IS FIXED TO BE FAST AND CHANGE EVERY PREP TO PRE
  
  if(!full){
    add_children(g, depth, prep);
    if(!pre) add_spouses(g, prep);
    add_common_ancestors(g, prep );
  } else {
    msg("max depth reached for %s - add all of the remaining circuit", g->get_var_name());
    add_all_remaining(g);
  }
  
  bool res = 0;
  msg("calling msolve for %s, d= %i, size= %i, var= %i", g->get_var_name(), depth, gate_poly.size(),  var.size());

  count_gb_calls++;
  if(depth == 3) unique_gb_calls++;
  count_arr[depth-3]++;
  size_arr[depth-3]+=gate_poly.size();
  var_size_arr[depth-3]+=var.size();

  call_init_time = process_time();
  res = print_to_msolve(g);
  call_end_time = process_time();
  double call_time = call_end_time - call_init_time;
  msg("used time for %s, d=%i: %.5f seconds", g->get_var_name(), depth, call_time);
  gb_time+=call_time;
  time_arr[depth-3]+=call_time;
  msg("");

  if(!res && depth < max_depth) {

      return linearize_via_gb(g, depth+1, pre, coeff,0);
    }
  else if(!res && depth == max_depth){
     return linearize_via_gb(g, depth+1, pre, coeff,1);
  } else if (!res) {
      return 0;
    }   

 
  var.clear();
  gate_poly.clear();
  return res;
  
}

