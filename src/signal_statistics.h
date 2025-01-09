/*------------------------------------------------------------------------*/
/*! \file signal_statistics.h
    \brief used to handle signals, messages and statistics

  Part of MultiLing : AIG Verification using Linear Groebner Bases.
  Copyright(C) 2024, 2025 Daniela Kaufmann, TU Wien
*/
/*------------------------------------------------------------------------*/
#ifndef MULTILING_SRC_SIGNAL_STATISTICS_H_
#define MULTILING_SRC_SIGNAL_STATISTICS_H_
/*------------------------------------------------------------------------*/
#include <stdarg.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <iostream>
/*------------------------------------------------------------------------*/
extern void(*original_SIGINT_handler)(int);
extern void(*original_SIGSEGV_handler)(int);
extern void(*original_SIGABRT_handler)(int);
extern void(*original_SIGTERM_handler)(int);

extern int van_mon_depth_count;
extern int child_superset_count;
extern int grandchild_is_subset_of_child_count;
extern int grandchild_is_superset_child_count;
extern int grand_children_are_equal_count;
extern int two_subsets_count;
extern int l_f_count;
extern int lin_gb_count;
extern int child_l_f_count;
extern int children_share_l_f_count;
extern int elim_pos_nodes;
extern int totalcount;
extern int equiv_gate_count;
extern bool miter_inp;
extern bool mult_inp;
extern const char * bench;




/**
    Returns name of signal

    @param sig integer

    @return char *
*/
const char * signal_name(int sig);

/**
    Initialize all signals
*/
void init_all_signal_handers();


/**
    Catch signal and prints corresponding message

    @param sig integer
*/
void catch_signal(int sig);

/**
    Resets all signal handlers
*/
void reset_all_signal_handlers();

/*------------------------------------------------------------------------*/
// / Level of output verbosity, ranges from 0 to 4
extern int verbose;


/**
    Prints an error message to stderr and exits the program

    @param char* fmt message
    @param int error_code
*/
void die(int error_code, const char *fmt,...);

/**
    Prints a message to stdout

    @param char* fmt message
*/
void msg(const char *fmt, ...);
/*------------------------------------------------------------------------*/

// / Time measures used for verify/certify modus
extern double init_time;          // /< measure initializing time
extern double slicing_elim_time;  // /< measure time used to eliminate & slice
extern double reduction_time;     // /< measure time used to reduce
extern double reset_time;         // /< measure resetting time
extern double substitution_time;  // /< measure time used in substitution
extern int count_gb_calls;
/**
    Determines max used memory
*/
size_t maximum_resident_set_size();

/**
    Determines the used process time
*/
double process_time();

/**
    Print statistics of maximum memory and used process time depending on
    selected modus

    @param modus integer
*/
void print_statistics();

#endif  // MULTILING_SRC_SIGNAL_STATISTICS_H_
