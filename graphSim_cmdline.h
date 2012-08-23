/** @file graphSim_cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.3
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef GRAPHSIM_CMDLINE_H
#define GRAPHSIM_CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE "graphSim"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "graphSim"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "0.1"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  int verbose_lvl_arg;	/**< @brief Level of verbosity (0 = nothing, 3 = full) (default='0').  */
  char * verbose_lvl_orig;	/**< @brief Level of verbosity (0 = nothing, 3 = full) original value given at command line.  */
  const char *verbose_lvl_help; /**< @brief Level of verbosity (0 = nothing, 3 = full) help description.  */
  int rates_flag;	/**< @brief Create rates for treekin (default=off).  */
  const char *rates_help; /**< @brief Create rates for treekin help description.  */
  char * rates_file_arg;	/**< @brief File where to write rates (default='rates.out').  */
  char * rates_file_orig;	/**< @brief File where to write rates original value given at command line.  */
  const char *rates_file_help; /**< @brief File where to write rates help description.  */
  double temp_arg;	/**< @brief Temperature in Celsius (default='37.0').  */
  char * temp_orig;	/**< @brief Temperature in Celsius original value given at command line.  */
  const char *temp_help; /**< @brief Temperature in Celsius help description.  */
  int print_graph_flag;	/**< @brief Generate xmgrace graph from best paths? (default=off).  */
  const char *print_graph_help; /**< @brief Generate xmgrace graph from best paths? help description.  */
  int RNAstructs_flag;	/**< @brief Assume RNA structs input from RNAsubopt (default=on).  */
  const char *RNAstructs_help; /**< @brief Assume RNA structs input from RNAsubopt help description.  */
  int dot_print_flag;	/**< @brief Generate dot file? (default=off).  */
  const char *dot_print_help; /**< @brief Generate dot file? help description.  */
  int dot_graph_flag;	/**< @brief Generate also neato file from dot output? (forces dot-print) (default=off).  */
  const char *dot_graph_help; /**< @brief Generate also neato file from dot output? (forces dot-print) help description.  */
  int dot_threshold_arg;	/**< @brief Maximal number of sequences in graph (0=infty) (default='0').  */
  char * dot_threshold_orig;	/**< @brief Maximal number of sequences in graph (0=infty) original value given at command line.  */
  const char *dot_threshold_help; /**< @brief Maximal number of sequences in graph (0=infty) help description.  */
  int dot_struct_arg;	/**< @brief Print neato graph only around one structure (graph will be saved to graph<num>.pdf).  */
  char * dot_struct_orig;	/**< @brief Print neato graph only around one structure (graph will be saved to graph<num>.pdf) original value given at command line.  */
  const char *dot_struct_help; /**< @brief Print neato graph only around one structure (graph will be saved to graph<num>.pdf) help description.  */
  int dot_radius_arg;	/**< @brief Radius of --dot-struct print (default='8').  */
  char * dot_radius_orig;	/**< @brief Radius of --dot-struct print original value given at command line.  */
  const char *dot_radius_help; /**< @brief Radius of --dot-struct print help description.  */
  int all_from_arg;	/**< @brief Construct all paths from node number\n(all-to must be specified too) (default='0').  */
  char * all_from_orig;	/**< @brief Construct all paths from node number\n(all-to must be specified too) original value given at command line.  */
  const char *all_from_help; /**< @brief Construct all paths from node number\n(all-to must be specified too) help description.  */
  int all_to_arg;	/**< @brief Construct all paths to node number.  */
  char * all_to_orig;	/**< @brief Construct all paths to node number original value given at command line.  */
  const char *all_to_help; /**< @brief Construct all paths to node number help description.  */
  int all_len_arg;	/**< @brief Maximal length of paths (default='8').  */
  char * all_len_orig;	/**< @brief Maximal length of paths original value given at command line.  */
  const char *all_len_help; /**< @brief Maximal length of paths help description.  */
  int simulate_flag;	/**< @brief Simulate? (simulates from all-from to all-to) (default=off).  */
  const char *simulate_help; /**< @brief Simulate? (simulates from all-from to all-to) help description.  */
  int sim_len_arg;	/**< @brief Maximal length of simulation (default='10000').  */
  char * sim_len_orig;	/**< @brief Maximal length of simulation original value given at command line.  */
  const char *sim_len_help; /**< @brief Maximal length of simulation help description.  */
  int sim_cnt_arg;	/**< @brief Number of simulated trajectories (default='1000').  */
  char * sim_cnt_orig;	/**< @brief Number of simulated trajectories original value given at command line.  */
  const char *sim_cnt_help; /**< @brief Number of simulated trajectories help description.  */
  int sim_print_len_arg;	/**< @brief Number of best paths to print (-1 means all of them) (default='-1').  */
  char * sim_print_len_orig;	/**< @brief Number of best paths to print (-1 means all of them) original value given at command line.  */
  const char *sim_print_len_help; /**< @brief Number of best paths to print (-1 means all of them) help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int verbose_lvl_given ;	/**< @brief Whether verbose-lvl was given.  */
  unsigned int rates_given ;	/**< @brief Whether rates was given.  */
  unsigned int rates_file_given ;	/**< @brief Whether rates-file was given.  */
  unsigned int temp_given ;	/**< @brief Whether temp was given.  */
  unsigned int print_graph_given ;	/**< @brief Whether print-graph was given.  */
  unsigned int RNAstructs_given ;	/**< @brief Whether RNAstructs was given.  */
  unsigned int dot_print_given ;	/**< @brief Whether dot-print was given.  */
  unsigned int dot_graph_given ;	/**< @brief Whether dot-graph was given.  */
  unsigned int dot_threshold_given ;	/**< @brief Whether dot-threshold was given.  */
  unsigned int dot_struct_given ;	/**< @brief Whether dot-struct was given.  */
  unsigned int dot_radius_given ;	/**< @brief Whether dot-radius was given.  */
  unsigned int all_from_given ;	/**< @brief Whether all-from was given.  */
  unsigned int all_to_given ;	/**< @brief Whether all-to was given.  */
  unsigned int all_len_given ;	/**< @brief Whether all-len was given.  */
  unsigned int simulate_given ;	/**< @brief Whether simulate was given.  */
  unsigned int sim_len_given ;	/**< @brief Whether sim-len was given.  */
  unsigned int sim_cnt_given ;	/**< @brief Whether sim-cnt was given.  */
  unsigned int sim_print_len_given ;	/**< @brief Whether sim-print-len was given.  */

  char **inputs ; /**< @brief unamed options (options without names) */
  unsigned inputs_num ; /**< @brief unamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* GRAPHSIM_CMDLINE_H */
