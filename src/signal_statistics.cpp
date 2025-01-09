/*------------------------------------------------------------------------*/
/*! \file signal_statistics.cpp
    \brief used to handle signals, messages and statistics

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/

#include "signal_statistics.h"

/*------------------------------------------------------------------------*/
// Global variable
void(*original_SIGINT_handler)(int);
void(*original_SIGSEGV_handler)(int);
void(*original_SIGABRT_handler)(int);
void(*original_SIGTERM_handler)(int);
/*------------------------------------------------------------------------*/
const char * bench = 0;
int van_mon_depth_count = 0;
int child_superset_count= 0;
int grandchild_is_subset_of_child_count= 0;
int grandchild_is_superset_child_count= 0;
int two_subsets_count= 0;
int grand_children_are_equal_count=0;
int l_f_count=0;
int child_l_f_count=0;
int children_share_l_f_count = 0;
int totalcount = 0;
int equiv_gate_count = 0;
int elim_pos_nodes = 0;
int lin_gb_count = 0;
int count_gb_calls = 0;
bool miter_inp = 0;
bool mult_inp = 0;


const char * signal_name(int sig) {
  switch (sig) {
    case SIGINT: return "SIGINT";
    case SIGSEGV: return "SIGSEGV";
    case SIGABRT: return "SIGABRT";
    case SIGTERM: return "SIGTERM";
    default: return "SIGUNKNOWN";
  }
}

/*------------------------------------------------------------------------*/

void init_all_signal_handers() {
  original_SIGINT_handler  = signal(SIGINT,  catch_signal);
  original_SIGSEGV_handler = signal(SIGSEGV, catch_signal);
  original_SIGABRT_handler = signal(SIGABRT, catch_signal);
  original_SIGTERM_handler = signal(SIGTERM, catch_signal);
}

/*------------------------------------------------------------------------*/

void catch_signal(int sig) {
  printf("c\nc caught signal '%s'(%d)\nc\n", signal_name(sig), sig);
  printf("c\nc raising signal '%s'(%d) again\n", signal_name(sig), sig);
  reset_all_signal_handlers();
  fflush(stdout);
  raise(sig);
}

/*------------------------------------------------------------------------*/

void reset_all_signal_handlers() {
  (void) signal(SIGINT, original_SIGINT_handler);
  (void) signal(SIGSEGV, original_SIGSEGV_handler);
  (void) signal(SIGABRT, original_SIGABRT_handler);
  (void) signal(SIGTERM, original_SIGTERM_handler);
}

/*------------------------------------------------------------------------*/
// Global variable
int verbose = 1;

/*------------------------------------------------------------------------*/

void msg(const char *fmt, ...) {
  va_list ap;
  fputs_unlocked("[mltlng] ", stdout);
  va_start(ap, fmt);
  vfprintf(stdout, fmt, ap);
  va_end(ap);
  fputc_unlocked('\n', stdout);
  fflush(stdout);
}

/*------------------------------------------------------------------------*/

void die(int error_code, const char *fmt, ...) {
  fflush(stdout);
  va_list ap;
  fprintf(stderr, "*** [mltlng] error code %i \n", error_code);
  fputs_unlocked("*** [mltlng] ", stderr);
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fputc('\n', stderr);
  fflush(stderr);
  exit(error_code);
}

/*------------------------------------------------------------------------*/
// Global variables

double init_time, slicing_elim_time, reduction_time, reset_time;
double substitution_time;

/*------------------------------------------------------------------------*/

size_t maximum_resident_set_size() {
  struct rusage u;
  if (getrusage(RUSAGE_SELF, &u)) return 0;
  return((size_t) u.ru_maxrss) << 10;
}

/*------------------------------------------------------------------------*/

double process_time() {
  struct rusage u;
  if (getrusage(RUSAGE_SELF, &u)) return 0;
  double res = u.ru_utime.tv_sec + 1e-6 * u.ru_utime.tv_usec;
  res += u.ru_stime.tv_sec + 1e-6 * u.ru_stime.tv_usec;
  return res;
}

static double percent(unsigned a, unsigned b) { return b ? 100.0*a/b : 0; }
/*------------------------------------------------------------------------*/

void print_statistics() {
  msg("");
  msg("maximum resident set size:     %18.2f MB",
  maximum_resident_set_size() / static_cast<double>((1<<20)));
 // msg("van_mon_depth_count: %32i (%.2f %)",van_mon_depth_count, percent(van_mon_depth_count,totalcount)); 
 // msg("l_f_count: %42i (%.2f %)",l_f_count, percent(l_f_count,totalcount));
  msg("merged_nodes: %35i ",child_l_f_count, percent(child_l_f_count,totalcount));
 // msg("children_share_l_f_count: %27i (%.2f %)",children_share_l_f_count, percent(children_share_l_f_count,totalcount));
 // msg("lin_gb_count: %39i (%.2f %)",lin_gb_count, percent(lin_gb_count,totalcount));
  msg("eliminated_pos_nodes: %27i", elim_pos_nodes);
  msg("gb_calls_count: %33i",count_gb_calls);
  msg("");
  msg("total process time:            %18.2f seconds",  process_time());
}
