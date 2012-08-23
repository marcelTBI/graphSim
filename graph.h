#ifndef __GRAPH_H
#define __GRAPH_H

#include <math.h>

#include <vector>
#include <string>
#include <assert.h>
#include <map>

class Node;

#include "tree.h"

using namespace std;

class Edge;
class SimplePath;

class Node
{
public:
  string name;
  int num;
  float energy;

  bool sink; // if this node is a sink

  // number of the minima where this one belongs (basin)
  int belongs;

  // maybe structure
  short *structure;

  // coincidence (list + matrix)
  vector<Edge*> edges;

  vector<Edge*> coincidence;

  // probabilities for simulations
  vector<double> prob;
  vector<double> logprob;

  vector<double> co_prob;
  vector<double> co_logprob;

  double sumprob;
  double logsumprob;

public:
  Node(int num, float energy, char *name, short *structure = NULL);
  Node() {
    num =-1;
  };
  ~Node();
  void CreateProb(double temp, int max_edges, bool norm = false);

  void PrintDeb(bool probab = false);

  bool operator==(Node &right) {
    return num==right.num;
  }
};

class Edge
{
public:
  Node *src;
  Node *dest;
  float saddle;

  // some path measurements stuff
  double weight;

public:
  Edge(Node *src, Node *dest, float saddle);

  bool coincide(int num);
  bool coincide(int num, int num2);
  void PrintDeb();
  int goesTo(int from);
  Node *goesTo(Node *from);

private:
  Edge() {};
};

static bool compEdge(Edge *frst, Edge *scnd) {
  if (frst->src->num == scnd->src->num) return frst->dest->num < scnd->dest->num;
  if (frst->src->num == scnd->dest->num) return frst->dest->num < scnd->src->num;
  if (frst->dest->num == scnd->dest->num) return frst->src->num < scnd->src->num;
  if (frst->dest->num == scnd->src->num) return frst->src->num < scnd->dest->num;

  assert(false);
  return true;
}

class Graph
{
public:
  // nodes
  vector<Node*> nodes;
  map<short*, int> nodes_map;

  // edges
  vector<Edge*> edges;

  // constant for Simulation
  int max_edges;

  // sequence
  char *sequence;

public:
  void BuildAdjacency();
  void ResolveDegeneracy();
  bool Paint(vector<bool> nodes_paint, int node_num);

  void CreateProb(double temp, bool norm = false);

  int LoadFromFile(const char *filename);
  int LoadStructsFromFile(const char *filename);
  int LoadFromFile(FILE *file);
  int LoadStructsFromFile(FILE *file);
  int PrintToDot(const char *filename, bool cluster = true, bool neato = false, int threshold = -1);  // cluster = nice clustering + colouring
  int PrintDotStruct(int structure, int radius, bool cluster = true, int threshold = -1);  // print only around one structure
  int PrintRates(const char *filename, float temp);

  // functions to assign weights
  void AssignWeights(double temp = 37.0);

  vector<int> Dijkstra(int src, int target);

  vector<SimplePath> Simulate(int source, int dest, int how_many = 1000, double temp = 37.0, int threshold = 100);

  vector<SimplePath> ConstructAllPaths(int source, int dest, int max_length, double temp = 37.0, double threshold = 0.0);
  void ConstructPath(vector<SimplePath> &paths, SimplePath &path, int dest, int max_length);

  // score a single path
  void Score(SimplePath &path);

  // print graph
  void Print(bool prob = false);

private:
  void addEdge(int src, int dest, float saddle = INFINITY);
public:
  Graph();
  ~Graph();
};


class SimplePath {
public:
  vector<int> points;
  double gilles;
  double max_energy;
  double path_prob;
  double simple_prob;

  // how many of them
  int number;

  // closed?
  bool closed;

  bool operator<(const SimplePath &left) const {
    if (max_energy == left.max_energy) return points.size() < left.points.size();
    else return max_energy < left.max_energy;
  }

public:
  SimplePath();
  SimplePath(const SimplePath &path);

  void AddPoint(int num);
  void Close();

  bool ContainsNode(int num);
  int FindNode(int num);
  void Print(bool whole_path, bool with_count = false, bool force_print = true, FILE *out = stdout);
};

// print xmgrace file of "how_many" paths into "filename"
void PrintGraph(char *filename, Graph &gr, vector<SimplePath> &paths, bool stretch = true, bool xmgrace = true, int how_many = 20);

#endif
