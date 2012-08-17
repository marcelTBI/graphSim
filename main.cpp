#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <algorithm>

#include "graph.h"

extern "C" {
  #include "graphSim_cmdline.h"
}

using namespace std;

bool sort_gilles (const SimplePath &first, const SimplePath &second) { return first.gilles>second.gilles; }
bool sort_pathprob (const SimplePath &first, const SimplePath &second) { return first.path_prob>second.path_prob; }
bool sort_simpleprob (const SimplePath &first, const SimplePath &second) { return first.simple_prob>second.simple_prob; }

int main(int argc, char **argv)
{
  // parse arguments
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info) != 0) {
    fprintf(stderr, "Argument parsing problem.\n");
    exit(EXIT_FAILURE);
  }

  Graph graph;
  graph.LoadStructsFromFile(stdin);
  //graph.LoadFromFile(stdin);
  if (args_info.dot_print_flag || args_info.dot_graph_flag) graph.PrintToDot("graph.dot", true, args_info.dot_graph_flag, args_info.dot_threshold_arg);

  //graph.PrintRates("rates.txt", 200.0);

  /*graph.AssignWeights(37.0);
  vector<int> path = graph.Dijkstra(2, 4);

  for (int i=0; i<path.size(); i++) {
    printf("%d ", path[i]);
  }
  printf("\n");*/


  // all paths generation
  if (args_info.all_to_given) {

    float temp = args_info.temp_arg;

    vector<SimplePath> paths;
    paths = graph.ConstructAllPaths(args_info.all_from_arg, args_info.all_to_arg, args_info.all_len_arg, temp);

    /*// erase first some and keep only some
    paths.erase(paths.begin(), paths.begin()+10);
    paths.erase(paths.begin()+10, paths.end());
    */

    int minprint = 10;

    for (unsigned int i=0; i<min(0, (int)paths.size()); i++) paths[i].Print(true);
    if (args_info.print_graph_flag) PrintGraph("path", graph, paths, true, true, 10);

    printf("\nBest paths due to energy barriers:\n");
    for (int i=0; i<min(minprint, (int)paths.size()); i++) paths[i].Print(true);
    sort(paths.begin(), paths.end(), sort_gilles);
    printf("Best paths due to gillespie:\n");
    for (int i=0; i<min(minprint, (int)paths.size()); i++) paths[i].Print(true);
    sort(paths.begin(), paths.end(), sort_pathprob);
    printf("Best paths due to path probability:\n");
    for (int i=0; i<min(minprint, (int)paths.size()); i++) paths[i].Print(true);
    /*sort(paths.begin(), paths.end(), sort_simpleprob);
    printf("Best paths due to simple probability:\n");
    for (int i=0; i<min(minprint, (int)paths.size()); i++) paths[i].Print(true);*/
  }

  // Simulation
  if (args_info.simulate_flag) {
    vector<SimplePath> paths = graph.Simulate(args_info.all_from_arg, args_info.all_to_arg, args_info.sim_cnt_arg, args_info.temp_arg, args_info.sim_len_arg);
    //root.Print();
    printf("Simulation:\n");
    int len = (args_info.sim_print_len_arg!=-1 ? min(args_info.sim_print_len_arg, (int)paths.size()) : paths.size());
    for (unsigned int i=0; i<len; i++) paths[i].Print(true, true, false);
  }

  //graph.Print(true);

  cmdline_parser_free(&args_info);
  return 0;
}
