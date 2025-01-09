/*------------------------------------------------------------------------*/
/*! \file multiling.cpp
    \brief main file of our tool MultiLing

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#define VERSION "1.0"
/*------------------------------------------------------------------------*/
// / Manual of MultiLing, will be printed with command line '-h'
static const char * USAGE =
"\n"
"### USAGE ###\n"
"usage : multiling <input.aig> <mode | spec_file> [<option> ...] \n"
"\n"
"Either a <mode> with an automatically generated specification has to \n"
"be selected or a user-defined spec has to be provided in a file <spec_file>\n"
"\n"
"\n"
"<mode> = -mult:\n"
"    assumes that the AIG represents an unsigned n*n = 2n multiplier\n"
"    automatically generated spec: 2^(2n-1)s_(2n-1) + ...+s_0 - \n"
"      (2^(n-1)a_(n-1)+...+a_0)*(2^(n-1)b_(n-1)+...+b_0) \n"
"\n"
"<mode> = -miter:\n"
"    assumes that the AIG represents an unsigned n*n = 2n multiplier\n"
"    automatically generated spec: s-1 \n"
"\n"
"<spec_file> =    file containing a single polynomial used as the specification\n"
"\n"
"\n"
"<option> = the following options are available \n"
"   -h | --help           print this command line summary \n"
"   -v<1,2,3,4>           different levels of verbosity(default -v1) \n"
"   -no-counter-examples  do not generate and write counter examples\n"
"     \n"
"     \n"
;
/*------------------------------------------------------------------------*/
#include "parser.h"
#include "gate.h"
#include "polynomial_solver.h"
#include "polyparser.h"

#include <ctime>
/*------------------------------------------------------------------------*/
// / Name of the input file
static const char * input_name = 0;
static const char * spec_name = 0;

static int mode;
/*------------------------------------------------------------------------*/
// ERROR CODES:

static int err_no_file    = 10; // no input file given
static int err_mode_sel   = 11; // mode has already been selected/not selected
static int err_wrong_arg  = 12; // wrong number of arguments given
static int err_proof_form = 13; // too many proof formats selected


/*------------------------------------------------------------------------*/
/**
    Calls the deallocaters of the involved data types
    @see reset_all_signal_handlers()
    @see delete_gates()
    @see deallocate_terms()
    @see deallocate_mstack()
    @see clear_mpz()
*/
static void reset_all() {
  reset_all_signal_handlers();
  delete_gates();
  deallocate_terms();
  deallocate_mstack();
  clear_mpz();
}
/*------------------------------------------------------------------------*/

//----------------------------------------------------------------------


/**
    Main Function of MultiLing.
    Reads the given AIG and depending on the selected mode, either
    calls the substution engine or the polynomial solver.

    Prints statistics to stdout after finishing.
*/
int main(int argc, char ** argv) {
  std::srand(std::time(0));
  msg("MultiLinG Version " VERSION);
  msg("AIG Verification using Linear Extractions from Groebner Bases");
  msg("Copyright(C) 2024, 2025, Daniela Kaufmann, TU Wien");

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-h") ||
    !strcmp(argv[i], "--help")) {
      fputs(USAGE, stdout);
      fflush(stdout);
      exit(0);
    } else if (!strcmp(argv[i], "-v0")) { verbose = 0;
    } else if (!strcmp(argv[i], "-v1")) { verbose = 1;
    } else if (!strcmp(argv[i], "-v2")) { verbose = 2;
    } else if (!strcmp(argv[i], "-v3")) { verbose = 3;
    } else if (!strcmp(argv[i], "-v4")) { verbose = 4;
    } else if (!strcmp(argv[i], "-miter")) {
      if (mult_inp) die(err_proof_form, "too many modes selected(try '-h')");
      miter_inp = 1;
    } else if (!strcmp(argv[i], "-mult")) {
      if (miter_inp) die(err_proof_form, "too many modes selected(try '-h')");
      mult_inp = 1;      
    } else if (!strcmp(argv[i], "-no-counter-examples")) {
      gen_witness = 0;
    } else if (spec_name) {
      die(err_wrong_arg, "too many arguments '%s', '%s', and '%s'(try '-h')",
        input_name, spec_name, argv[i]);
    } else if (input_name) { spec_name = argv[i];
    } else {
      input_name = argv[i];
      bench = input_name;
    }
  }
  mode = 2;
   
  if (!input_name)  die(err_no_file, "no input file given(try '-h')");
  if (!miter_inp && !mult_inp && !spec_name)  die(err_no_file, "no spec file given(try '-h')");
  if(!miter_inp && !mult_inp) die(err_mode_sel, "select mode(try -h for more information)");

  msg("using msolve for calculating the Gröbner basis");
  
  init_all_signal_handers();
  init_nonces();

  parse_aig(input_name);
  bool res;
  init_mpz(NN);
  init_gates();

  Polynomial * spec;
  if(mult_inp) spec = mult_spec_poly();
  else if (miter_inp) spec = miter_spec_poly();
  else spec = parse_specification_polynomial(spec_name);
  if(verbose >= 2) spec->print(stdout);
  res = verify(input_name, spec);
  reset_aig_parsing();
  reset_all();

  reset_time = process_time();
  print_statistics();

  return res;
}
