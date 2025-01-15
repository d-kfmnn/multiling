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
int unique_gb_calls = 0;
bool miter_inp = 0;
bool mult_inp = 0;
int count_arr[6];
int size_arr[6];
int var_size_arr[6];
double time_arr[6];
double gb_time=0;
double call_init_time =0;
double call_end_time=0;
int non_linear_count=0;
int linear_count=0;


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
  
  fprintf(stdout,"\n\n\nprinting statistics until interruption:\n");
  double interruption = process_time();
  double last_gb_time = interruption - call_init_time;
  gb_time+=last_gb_time;
  print_statistics();

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



size_t maximum_resident_set_size() {
  struct rusage u;
  if (getrusage(RUSAGE_SELF, &u)) return 0;
  return((size_t) u.ru_maxrss) << 10;
}
/*------------------------------------------------------------------------*/
struct timeval start_tv;
void init_time(){
  gettimeofday(&start_tv, NULL);
}

/*------------------------------------------------------------------------*/

double process_time() {
  struct timeval tv;
  double elapsed = 0.0;

  gettimeofday(&tv, NULL);
  elapsed = (tv.tv_sec - start_tv.tv_sec) +   (tv.tv_usec - start_tv.tv_usec) / 1000000.0;
  return elapsed;

 /* struct rusage u;
  if (getrusage(RUSAGE_SELF, &u)) return 0;
  double res = u.ru_utime.tv_sec + 1e-6 * u.ru_utime.tv_usec;
  res += u.ru_stime.tv_sec + 1e-6 * u.ru_stime.tv_usec;
  return res;*/
}

static double percent(double a, double b) { return b ? 100.0*a/b : 0; }
static double avg_time(double a, double b) { return b ? a/b : 0; }
/*------------------------------------------------------------------------*/

void print_statistics() {

  int red_count = linear_count + non_linear_count;
  msg("");
  msg("--------------------------------------------------");
  msg("STATISTICS:");
  msg("");
  msg("PREPROCESSING: ");
  msg("merged_nodes: %35i ",child_l_f_count);
  msg("eliminated_pos_nodes: %27i", elim_pos_nodes);
  msg("");
  msg("MSOLVE:");
  msg("total msolve calls: %29i",count_gb_calls);
  msg("unique msolve calls: %28i", unique_gb_calls);
  msg("calls with depth 3: %29i (avg: %3.f gate poly; %3.f var; %5.3f seconds       total: %5.3f seconds)", count_arr[0], avg_time(size_arr[0],count_arr[0]), avg_time(var_size_arr[0],count_arr[0]),avg_time(time_arr[0],count_arr[0]),time_arr[0]);
  msg("calls with depth 4: %29i (avg: %3.f gate poly; %3.f var;  %5.3f seconds      total: %5.3f seconds)", count_arr[1], avg_time(size_arr[1],count_arr[1]), avg_time(var_size_arr[1],count_arr[1]),avg_time(time_arr[1],count_arr[1]),time_arr[1]);
  msg("calls with depth 5: %29i (avg: %3.f gate poly; %3.f var;  %5.3f seconds      total: %5.3f seconds)", count_arr[2], avg_time(size_arr[2],count_arr[2]), avg_time(var_size_arr[2],count_arr[2]),avg_time(time_arr[2],count_arr[2]),time_arr[2]);
  msg("calls with depth 6: %29i (avg: %3.f gate poly; %3.f var;  %5.3f seconds      total: %5.3f seconds)", count_arr[3], avg_time(size_arr[3],count_arr[3]), avg_time(var_size_arr[3],count_arr[3]),avg_time(time_arr[3],count_arr[3]),time_arr[3]);
  msg("calls with depth 7: %29i (avg: %3.f gate poly; %3.f var;  %5.3f seconds      total: %5.3f seconds)", count_arr[4], avg_time(size_arr[4],count_arr[4]), avg_time(var_size_arr[4],count_arr[4]),avg_time(time_arr[4],count_arr[4]),time_arr[4]);
  msg("calls with depth 8: %29i (avg: %3.f gate poly; %3.f var;  %5.3f seconds      total: %5.3f seconds)", count_arr[5], avg_time(size_arr[5],count_arr[5]), avg_time(var_size_arr[5],count_arr[5]),avg_time(time_arr[5],count_arr[5]),time_arr[5]);
  msg("");
  msg("REDUCTIONS: ");
  msg("total reductions: %31i", red_count);
  msg("linear reductions: %30i (%2.2f %)", linear_count, percent(linear_count, red_count));
  msg("non-linear reductions: %26i (%2.2f %)", non_linear_count, percent(non_linear_count, red_count));
  msg("");
  msg("TIME AND MEMORY: ");
  msg("maximum resident set size:     %18.2f MB",
  maximum_resident_set_size() / static_cast<double>((1<<20)));
  msg("total process time:            %18.3f seconds",  process_time());
  msg("including msolve time:         %18.3f seconds (%.2f %)",  gb_time, percent(gb_time, process_time()));
  msg("--------------------------------------------------");

}
