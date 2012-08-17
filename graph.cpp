#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "math.h"

#include <algorithm>
#include <queue>
#include <set>
#include <utility>

extern "C" {
  #include "pair_mat.h"
  //#include "fold.h"
}

#include "graph.h"

inline double absD(double a) {return (a>0.0?a:-a);}

bool neighbours(short *first, short *second);

/* reads a line no matter how long*/
char* my_getline(FILE *fp)
{
  char s[512], *line, *cp;
  line = NULL;
  do {
    if(fgets(s, 512, fp) == NULL) break;
    cp = strchr(s, '\n');
    if(cp != NULL) *cp = '\0';
    if(line == NULL) line = (char *) calloc(strlen(s) + 1, sizeof(char));
    else line = (char *) realloc(line, strlen(s) + strlen(line) + 1);
    strcat (line, s);
  } while (cp == NULL);
  return (line);
}

Node::Node(int num, float energy, char *name, short *structure)
{
  this->num = num;
  this->energy = energy;
  this->name = name;
  this->structure = structure;
  this->sink = true;
}

Node::~Node()
{
  if (structure) free(structure);
}

void Node::PrintDeb(bool probab)
{
  fprintf(stderr, "%d %s %f\n", num, name.c_str(), energy);
  if (probab) {
    for (unsigned int i=0; i<prob.size(); i++) {
      fprintf(stderr, "%d:%.4f (%.4f) ", edges[i]->goesTo(num), prob[i], logprob[i]);
    }
    fprintf(stderr, "\n%.4f (%.4f)\n", sumprob, logsumprob);
  }
  for (unsigned int i=0; i<edges.size(); i++) {
    edges[i]->PrintDeb();
  }
  fprintf(stderr, "\n");
}

void Node::CreateProb(double temp, int max_edges, bool norm)
{
  double _kT = 0.00198717*(273.15 + temp);
  prob.resize(edges.size());
  logprob.resize(edges.size());
  sumprob = 0.0;
  for (unsigned int i=0; i<edges.size(); i++) {
    prob[i]=min(1.0, exp((energy - edges[i]->goesTo(this)->energy)/_kT))/(double)max_edges;
    sumprob += prob[i];
  }
  for (unsigned int i=0; i<edges.size(); i++) {
    if (norm) prob[i]=prob[i]/sumprob;

    logprob[i]=log(prob[i]);

    // and coincidence prob
    co_prob[edges[i]->goesTo(num)] = prob[i];
    co_logprob[edges[i]->goesTo(num)] = logprob[i];
  }
  if (norm) sumprob = 1.0;
  logsumprob = log(sumprob);

  // and coincidence prob

}

//########################################################################

Edge::Edge(Node *src, Node *dest, float saddle)
{
  this->src = src;
  this->dest = dest;
  this->saddle = saddle;
}

void Edge::PrintDeb()
{
  fprintf(stderr, "\t%d %d", src->num, dest->num);
  if (!isinf(saddle)) fprintf(stderr, " %f", saddle);
  fprintf(stderr, "\n");
}

int Edge::goesTo(int from)
{
  if (src->num == from) return dest->num;
  if (dest->num == from) return src->num;
  return -1;
}

Node* Edge::goesTo(Node *from)
{
  if (src == from) return dest;
  if (dest == from) return src;
  return NULL;
}

bool Edge::coincide(int num)
{
  return (src->num == num || dest->num == num);
}

bool Edge::coincide(int num, int num2)
{
  return ((src->num == num && dest->num == num2) || (src->num == num2 && dest->num == num));
}

//########################################################################

Graph::Graph()
{
  srand( time(NULL) );
  //srand( 12 );

  max_edges = 0;

  sequence = NULL;
}

Graph::~Graph()
{
  for (unsigned int i=0; i<nodes.size(); i++) delete nodes[i];
  for (unsigned int i=0; i<edges.size(); i++) delete edges[i];
  if (sequence) free(sequence);
}

void Graph::addEdge(int src, int dest, float saddle)
{
  if (dest>=(int)nodes.size() || !nodes[dest]) {
    fprintf(stderr, "Non-existent node (%d)!!!\n", dest);
    return;
  }
  if (src>=(int)nodes.size() || !nodes[src]) {
    fprintf(stderr, "Non-existent node (%d)!!!\n", src);
    return;
  }
  Edge *e = new Edge(nodes[src], nodes[dest], saddle);
  edges.push_back(e);
  nodes[src]->edges.push_back(e);
  nodes[dest]->edges.push_back(e);

  // not sink anymore
  if (nodes[src]->energy < nodes[dest]->energy) nodes[dest]->sink = false;
  if (nodes[src]->energy > nodes[dest]->energy) nodes[src]->sink = false;
  // maybe we created some LM - !!!!  must resolve degeneracy - call ResolveDegeneracy

  // fill coincidence matrix
  if (nodes[src]->coincidence.size() != nodes.size()) nodes[src]->coincidence.resize(nodes.size());
  if (nodes[dest]->coincidence.size() != nodes.size()) nodes[dest]->coincidence.resize(nodes.size());
  nodes[src]->coincidence[dest] = e;
  nodes[dest]->coincidence[src] = e;
}

void Graph::ResolveDegeneracy()
{
  vector<bool> nodes_paint (nodes.size(), false);
  for (int i=0; i<nodes.size(); i++) {
    if (nodes[i] && nodes[i]->sink) {
      bool lower_found = Paint(nodes_paint, i);
      nodes[i]->sink = !lower_found;
    }
  }
}

// find lower?
bool Graph::Paint(vector<bool> nodes_paint, int node_num)
{
  nodes_paint[node_num] = true;
  bool lower_found = false;
  for (int i=0; i<nodes[node_num]->edges.size(); i++) {
    Edge *e = nodes[node_num]->edges[i];
    int new_node = e->goesTo(node_num);
    // degeneracy:
    if (nodes[new_node]->energy == nodes[node_num]->energy && !nodes_paint[new_node]) {
      nodes[new_node]->sink = false; // definitely false...
      lower_found = Paint(nodes_paint, new_node) || lower_found;
    }
    // found lower = not LM
    if (nodes[new_node]->energy < nodes[node_num]->energy) lower_found = true;
  }
  return lower_found;
}

void Graph::CreateProb(double temp, bool norm)
{
  for (unsigned int i=0; i<nodes.size(); i++) if (nodes[i]) {
    nodes[i]->co_prob.resize(nodes.size(), 0.0);
    nodes[i]->co_logprob.resize(nodes.size(), 0.0);
    nodes[i]->CreateProb(temp, max_edges, norm);
  }
}

int Graph::LoadFromFile(FILE *file)
{
  int num;
  char name[500];
  char line[600];
  float energy;
  while (gets(line) && sscanf(line, "%d %s %f\n", &num, name, &energy)==3) {
    if ((int)nodes.size()<num+1) nodes.resize(num+1, NULL);
    nodes[num] = new Node(num, energy, name);
  }

  int dest;
  int src;
  energy = 1.0/0.0;
  while (gets(line) && sscanf(line, "%d %d %f\n", &src, &dest, &energy)>=2) {
    addEdge(src, dest, energy);
  }

  // sort edges in nodes
  for (unsigned int i=0; i<nodes.size(); i++) {
    if (!nodes[i]) continue;
    sort(nodes[i]->edges.begin(), nodes[i]->edges.end(), compEdge);
    max_edges = max(max_edges, (int)nodes[i]->edges.size());
  }
  return 0;
}

int Graph::LoadStructsFromFile(FILE *file)
{
  char structure[500];
  char line[600];
  float energy;
  int num = 0;

  int belongs = 0;

  sequence = my_getline(file);

  while (gets(line) && sscanf(line, "%s %f %d\n", structure, &energy, &belongs)>=2) {
    if ((int)nodes.size()<num+1) nodes.resize(num+1, NULL);
    // make name + structure + store it
    char name[100];
    sprintf(name, "%d", num);
    short *str = make_pair_table(structure);
    nodes[num] = new Node(num, energy, name, str);
    nodes[num]->belongs = belongs;
    nodes_map.insert(make_pair(str, num));
    num++;

    belongs = 0;
  }

  BuildAdjacency();

  return 0;
}

int Graph::LoadStructsFromFile(const char *filename)
{
  FILE *file;
  file = fopen(filename, "r");
  if (file==NULL) {
    fprintf(stderr, "File \"%s\" not found!!!\n", filename);
    return -1;
  }

  LoadStructsFromFile(file);

  fclose(file);
  return 0;
}

int Graph::LoadFromFile(const char *filename)
{
  FILE *file;
  file = fopen(filename, "r");
  if (file==NULL) {
    fprintf(stderr, "File \"%s\" not found!!!\n", filename);
    return -1;
  }

  LoadFromFile(file);

  fclose(file);
  return 0;
}

void Graph::BuildAdjacency()
{
  for (unsigned int i=0; i<nodes.size(); i++) {
    int edges_added = 0;
    for (unsigned int j=i+1; j<nodes.size(); j++) {
      if (neighbours(nodes[i]->structure, nodes[j]->structure)) {  // very naive
        // create and add edge
        addEdge(i, j);
        edges_added++;
      }
    }
    if (edges_added > max_edges) max_edges = edges_added;
  }

  ResolveDegeneracy();
}

int Graph::PrintToDot(const char *filename, bool cluster, bool neato, int threshold)
{
  FILE *file;
  file = fopen(filename, "w");
  if (file==NULL) {
    fprintf(stderr, "Couldn't open file \"%s\"!!!\n", filename);
    return -1;
  }

  fprintf(file, "Graph G {\n\tnode [width=0.1, height=0.1, shape=circle];\n");

  int cnt = 0;

  // nodes
  double min_energy = 1e10;
  double max_energy = -1e10;
  for (unsigned int i=0; i<nodes.size(); i++) {
    if (nodes[i]) {
      if (nodes[i]->energy<min_energy) min_energy = nodes[i]->energy;
      if (nodes[i]->energy>max_energy) max_energy = nodes[i]->energy;

      cnt++;
      if (threshold>0 && cnt>=threshold) {
        threshold = nodes[i]->num;
        break;
      }
    }
  }

  for (unsigned int i=0; i<nodes.size(); i++) {
    if (!nodes[i]) continue;
    if (threshold>0 && nodes[i]->num>threshold) {
        break;
      }

    if (cluster)  {
      double value = (nodes[i]->sink? 0.000 : 0.750);
      fprintf(file, "\"%d\"\t[color = \"%.3f %.3f 1.000\", fontcolor=\"1.000 0.000 %.3f\", fontsize=8.0];\n", nodes[i]->num, 0.666*(nodes[i]->energy-min_energy)/(max_energy-min_energy), 1.0-value, value);
    } else          fprintf(file, "\t\"%d\"\t[pos=\"%d,%.2f!\"];\n", nodes[i]->num, i, nodes[i]->energy);
  }

  // edges
  double max_difference = 0.0;
  if (cluster) {
    for (unsigned int i=0; i<edges.size(); i++) {
      if (threshold>0 && (edges[i]->src->num>threshold || edges[i]->dest->num>threshold)) continue;
      double src = nodes[edges[i]->src->num]->energy;
      double out = nodes[edges[i]->dest->num]->energy;
      if (max_difference<absD(out-src)) max_difference = absD(out-src);
    }
  }

  for (unsigned int i=0; i<edges.size(); i++) {
    if (threshold>0 && (edges[i]->src->num>threshold || edges[i]->dest->num>threshold)) continue;
    if (cluster) {
      double src = nodes[edges[i]->src->num]->energy;
      double out = nodes[edges[i]->dest->num]->energy;
      double difference = 0.666*(1.0-absD(out-src)/max_difference);
      string direction = (src>out? "forward":(src==out?"both":"back"));
      fprintf(file, "\t\"%d\" -- \"%d\"\t[color=\"%.3f 1.000 1.000\", dir=%s];\n", edges[i]->src->num, edges[i]->dest->num, difference, direction.c_str());
    } else fprintf(file, "\t\"%d\" -- \"%d\";\n", edges[i]->src->num, edges[i]->dest->num);
  }

  fprintf(file, "}\n");

  fclose(file);

  // start neato
  if (neato) {
    char syst[200];
    char f3[100];
    strcpy(f3, filename);
    char *f2 = strtok(f3, " .");
    sprintf(syst, "neato -Tps < %s > %s.eps", filename, (f2==NULL? "graph2":f2));
    system(syst);
  }

  return 0;
}

int Graph::PrintRates(const char *filename, float temp)
{
  double _kT = 0.00198717*(273.15 + temp);
  FILE *rates;
  rates = fopen(filename, "w");
  if (rates==NULL) {
    fprintf(stderr, "Couldn't open file \"%s\"!!!\n", filename);
    return -1;
  }

/*
  for (unsigned int i=0; i<edges.size(); i++) {
    edges[i].PrintDeb();
  }
*/

  for (unsigned int i=0; i<nodes.size(); i++) {
    if (!nodes[i]) continue;
    //nodes[i]->PrintDeb();
    int e_point=0;
    for (unsigned int j=0; j<nodes.size(); j++) {
      if (!nodes[j]) continue;
      double rate = 0.0;
      if (e_point<(int)nodes[i]->edges.size() && nodes[i]->edges[e_point]->coincide(i, j)) {
        float saddle = nodes[i]->edges[e_point]->saddle;
        if (!isinf(saddle)) { // saddle heights specified
          float e_diff = saddle - nodes[i]->energy;
          rate = 1.0*exp(-e_diff/_kT);
        } else {
          rate = min(1.0, exp((nodes[i]->energy - nodes[j]->energy)/_kT));
        }
        e_point++;
      }
      fprintf(rates, "%10.4g ", rate);
    }
    fprintf(rates, "\n");
  }

  fclose(rates);
  return 0;
}

struct classcomp {
  bool operator() (const SimplePath& lhs, const SimplePath& rhs) const {
    int i=0;
    if (lhs.points.size()!=rhs.points.size()) return lhs.points.size()<rhs.points.size();
    while (i<lhs.points.size() && lhs.points[i]==rhs.points[i]) i++;
    if (i==lhs.points.size()) return false;
    return lhs.points[i]<rhs.points[i];
  }
};
bool numsort (const SimplePath &i, const SimplePath &j) {
  if (i.number == j.number) return i.points.size()<j.points.size();
  return (i.number>j.number);
}

vector<SimplePath> Graph::Simulate(int source, int dest, int how_many, double temp, int threshold)
{
  // create probabilities
  CreateProb(temp, true);

  // create map
  map<SimplePath, int, classcomp> paths;

  int all_time = 0;

  // simulate
  for (int i=0; i<how_many; i++) {
    Node *state = nodes[source];
    SimplePath path;
    path.AddPoint(source);
    int count = 0;
    while (count<threshold && state->num != dest) {
      assert(state);
      //printf("%4d", state->num);
      double rnd = rand()/(double)RAND_MAX;
      unsigned int edge = 0;
      while (edge<state->prob.size()) {
        rnd -= state->prob[edge];
        if (rnd<=0.0) break;
        edge++;
      }
      if (edge!=state->prob.size()) {
        state = state->edges[edge]->goesTo(state);
        path.AddPoint(state->num);
      }
      count++;
    }
    if (state->num==dest) path.Close();
    //printf("\n%d\n", count);
    all_time += count;
    if (paths.count(path)) paths[path]++;
    else paths.insert(make_pair(path, 1));
  }

  // FPT
  //printf ("fpt = %f\n", all_time/(double)how_many);

  vector<SimplePath> vec;
  map<SimplePath, int, classcomp>::iterator it;
  for (it=paths.begin(); it!=paths.end(); it++) {
    vec.push_back(it->first);
    vec[vec.size()-1].number = it->second;
  }

  for (int i=0; i<vec.size(); i++) {
    Score(vec[i]);
  }

  sort(vec.begin(), vec.end(), numsort);

  return vec;
}

vector<SimplePath> Graph::ConstructAllPaths(int source, int dest, int max_length, double temp, double threshold)
{
  // create probabilities
  CreateProb(temp);

  vector<SimplePath> paths;

  SimplePath path;
  path.points.push_back(source);

  // construct all path recursively
  ConstructPath(paths, path, dest, max_length-1);

  // measure them
  for (int i=0; i<paths.size(); i++) {
    Score(paths[i]);
  }

  // sort
  sort(paths.begin(), paths.end());

  return paths;
}

void Graph::ConstructPath(vector<SimplePath> &paths, SimplePath &path, int dest, int max_length)
{
  int num = path.points[path.points.size()-1];

  if (num == dest) {
    paths.push_back(path);
    //if (paths.size()%100==1) printf("found paths: %d\n", (int)paths.size());
    return ;
  }

  if (max_length == 0) return;

  // all edges
  for (int i=0; i<nodes[num]->edges.size(); i++) {
    Edge *edge = nodes[num]->edges[i];
    int goesTo = edge->goesTo(num);
    if (!path.ContainsNode(goesTo)) {
      path.points.push_back(goesTo);
      ConstructPath(paths, path, dest, max_length-1);
      path.points.pop_back();
    }
  }
}

void Graph::Score(SimplePath &path)
{
  double sumprob = 0.0;

  path.max_energy = -1e10;
  path.path_prob = 0.0;
  path.simple_prob = 0.0;
  path.gilles = 0.0;

  for (int i=1; i<path.points.size(); i++){
    int old = path.points[i-1];
    int curr = path.points[i];

    // find prob:
    double prob = nodes[old]->co_prob[curr];
    double logprob = nodes[old]->co_logprob[curr];
    if (prob == 0.0) {fprintf(stderr, "smth is wrong"); return;}

    // max energy
    path.max_energy = max(path.max_energy, (double)nodes[curr]->energy);

    // gillespie
    sumprob += max_edges*nodes[old]->sumprob;
    path.gilles += (logprob - nodes[old]->logsumprob);

    // path_prob
    path.path_prob += logprob;

    // simple_prob
    path.simple_prob += logprob/nodes[old]->sumprob;
  }
  path.gilles -= log(sumprob);
  path.path_prob += log(path.points.size()-1);
}

// functions to assign weights (log path prob.)
void Graph::AssignWeights(double temp)
{
  CreateProb(temp);
  // create probabilities
  for (unsigned int i=0; i<nodes.size(); i++) if (nodes[i]) {
    //nodes[i]->PrintDeb(true);
    for (unsigned int j=0; j<nodes[i]->edges.size(); j++) {
      double logprob = nodes[i]->logprob[j];
      nodes[i]->edges[j]->weight = logprob;
    }
  }
}

void Graph::Print(bool prob)
{
  for (unsigned int i=0; i<nodes.size(); i++) if (nodes[i]) {
    nodes[i]->PrintDeb(prob);
  }
}

inline int minimum_index(vector<bool> done, vector<double> arr) {
  int min_ind = -1;
  double min = 1.0/0.0;
  for (unsigned int i=1; i<arr.size(); i++) {
    if (!done[i] && arr[i]<=min) {
      min_ind = i;
      min = arr[i];
    }
  }
  return min_ind;
}

vector<int> Graph::Dijkstra(int src, int dst)
{
  vector<int> path;
  if ((int)nodes.size()+1 < max(src, dst)) return path;

  vector<bool>  done;
  vector<double> dist;
  vector<int> prev;
  dist.resize(nodes.size(), 1.0/0.0);
  prev.resize(nodes.size(), -1);
  done.resize(nodes.size(), false);

  dist[src] = 0;
  int count = 0;
  for (unsigned int i=0; i<nodes.size(); i++) if (nodes[i]) count++;

  while (count--) {
    int u = minimum_index(done, dist);
    if (u == -1 || dist[u] > 1e100) break;

    done[u] = true;

    if (u==dst) break;

    for (unsigned int i=0; i<nodes[u]->edges.size(); i++) {
      int v = nodes[u]->edges[i]->goesTo(u);
      if (done[v]) continue;
      Edge *edge = nodes[u]->edges[i];

      double alt = dist[u] + edge->weight;
      if (alt < dist[v]) {
        dist[v] = alt;
        prev[v] = u;
      }
    }
  }

  // retreive path:
  int u = dst;
  while (prev[u] != -1) {
    path.push_back(u);
    u = prev[u];
  }
  path.push_back(u);

  reverse(path.begin(), path.end());
  return path;
}

// only non-shift neighbours
bool neighbours(short *first, short *second)
{
  if (first[0]!=second[0]) {
    fprintf(stderr, "Lengths are not equal!!! (%d, %d)", first[0], second[0]);
    return false;
  }
  int diffs = 0;
  for (int i=1; i<=first[0]; i++) {
    if (first[i]!=second[i]) {
      diffs++;
      if (diffs>2) break;
    }
  }

  if (diffs != 2) return false;
  else return true;
}

SimplePath::SimplePath()
{
  gilles = 0.0;
  path_prob = 0.0;
  simple_prob = 0.0;
  max_energy = -1e10;
  closed = false;
}

void SimplePath::Close()
{
  closed = true;
}

SimplePath::SimplePath(const SimplePath &path)
{
  points.assign(path.points.begin(), path.points.end());
  gilles = path.gilles;
  path_prob = path.path_prob;
  simple_prob = path.simple_prob;
  max_energy = path.max_energy;
  number = path.number;
  closed = path.closed;
}

void SimplePath::AddPoint(int num)
{
  int h;
  if ((h = FindNode(num)) != -1) {
    points.erase(points.begin()+h+1, points.end());
  } else {
    points.push_back(num);
  }
}

int SimplePath::FindNode(int num)
{
  for (int i=0; i<points.size(); i++) {
    if (points[i]==num) return i;
  }
  return -1;
}

bool SimplePath::ContainsNode(int num)
{
  for (int i=0; i<points.size(); i++) {
    if (points[i]==num) return true;
  }
  return false;
}

void SimplePath::Print(bool whole_path, bool with_count, bool force_print, FILE *out)
{
  if (!force_print && !closed) return;
  if (with_count) fprintf(out, "%4d: ", number);
  fprintf(out, "e(%6.2f) ", max_energy);
  if (!with_count) fprintf(out, "p(%7.3f) ", path_prob);
  if (!with_count) fprintf(out, "g(%7.3f) ", gilles);
  if (with_count) fprintf(out, "                ");
  //fprintf(out, "(%8.4f) ", simple_prob);
  fprintf(out, "%4d:", (int)points.size());

  if (whole_path) {
    for (unsigned int i=0; i<points.size(); i++) {
      fprintf(out, "%4d ", points[i]);
    }
  }

  fprintf(out, "\n");
}

void PrintGraph(char *filename, Graph &gr, vector<SimplePath> &paths, bool stretch, bool xmgrace, int how_many)
{
  int max_size = 0;
  for (int i=0; i<paths.size(); i++) {
    if (max_size<paths[i].points.size()) max_size = paths[i].points.size();
  }

  char syst[10000] = "xmgrace -nxy ";

  for (int i=0; i<min(how_many, (int)paths.size()); i++) {

    FILE *graph;
    char f[100];
    sprintf(f, "%s.%d", filename, i);
    strcat(syst, f);
    strcat(syst, " ");
    graph = fopen(f, "w");
    if (!graph) {
      fprintf(stderr, "WARNING: Cannot open file \"%s\" for writing!\n", filename);
      return;
    }

    for (int j=0; j<paths[i].points.size(); j++) {
      double value = j;
      if (stretch) value *= (max_size-1)/((double)paths[i].points.size()-1);
      if (!stretch && j == paths[i].points.size()-1) fprintf(graph, "%7d %6.2f \n", max_size-1, gr.nodes[paths[i].points[j]]->energy);
      else fprintf(graph, "%7g %6.2f \n", value, gr.nodes[paths[i].points[j]]->energy-0.01*i);
    }

    fclose(graph);
  }

  // start xmgrace
  if (xmgrace) {
    strcat(syst, "&");
    system(syst);
  }
}

