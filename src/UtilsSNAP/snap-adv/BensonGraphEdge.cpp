#include "BensonGraphEdge.h"

// SEQUENTIAL //////////////////////////////////////////////////////////////////
/******************************************************************************/
/*                         Some Utilities Functions                           */

BensonGraphEdge::~BensonGraphEdge(){
}

// Dummy Constructor
BensonGraphEdge::BensonGraphEdge() : BensonGraph(), motif() {
}

// Constructor: If decomp is true, the MAM is built
BensonGraphEdge::BensonGraphEdge(PNGraph& one_graph,MotifTypeEdge& one_motif,bool decomp, bool verb) : BensonGraph(one_graph,decomp,verb), motif(one_motif) {
  if (decomp){
    int tt = motif.getType();
    switch (tt){
      case DIR:
        DirAdjacency(); // iof the directed Graph is asked
        break;
      case SYM:
        SymAdjacency(); // if the symmetrised graph is asked
        break;
      default:
        TExcept::Throw("Unknown edge type.");
    }
  }
}

/******************************************************************************/
/*                             Symmetrised Graph                              */

// For incrementing weight in the symmetrised graph
void BensonGraphEdge::IncrementWeight(int const i, int const j){
  {THash<TInt, TInt>& edge_weights = weight[i];
  edge_weights(j) += 1;}
  {THash<TInt, TInt>& edge_weights = weight[j];
  edge_weights(i) += 1;}
}

void BensonGraphEdge::SymProcOneNode(const TNGraph::TNodeI& NI,const int& src){
  for (int i = 0; i < NI.GetOutDeg(); i++) { // For each out neighbour of a node
    int nbr = NI.GetOutNId(i);
    IncrementWeight(src,  nbr); // symmetrically update the MAM
  }
}

void BensonGraphEdge::SymAdjacency(){
  if (verbose) {printf("Building the symmetrised graph sequentially...\n");}
  int cpt = 0;
  for (TNGraph::TNodeI NI = graph->BegNI(); NI < graph->EndNI(); NI++) { // Process each node
    if (verbose){
      if (cpt%10000 == 0){
        printf("Processing node %d\n",cpt);
      }
      cpt++;
    }
    int src = NI.GetId();
    SymProcOneNode(NI, src);
  }
}

/******************************************************************************/
/*                               Directed Graph                               */

// For incrementing weight in the directed graph
void BensonGraphEdge::IncrementDirWeight(int const src, int const tgt){
  THash<TInt, TInt>& edge_weights = weight[src];
  edge_weights(tgt) += 1;
}

void BensonGraphEdge::DirProcOneNode(const TNGraph::TNodeI& NI,const int& src){
  for (int i = 0; i < NI.GetOutDeg(); i++) { // For each out neighbour of a node
    int nbr = NI.GetOutNId(i);
    IncrementDirWeight(src,nbr); // Non symmetrically update the MAM
  }
}

void BensonGraphEdge::DirAdjacency(){
  if (verbose) {printf("Building the directed graph sequentially...\n");}
  int cpt = 0;
  for (TNGraph::TNodeI NI = graph->BegNI(); NI < graph->EndNI(); NI++) { // Process each node
    if (verbose){
      if (cpt%10000 == 0){
        printf("Processing node %d\n",cpt);
      }
      cpt++;
    }
    int src = NI.GetId();
    DirProcOneNode(NI, src);
  }
}

// MULTITHREADED ///////////////////////////////////////////////////////////////
/******************************************************************************/
/*                         Some Utilities Functions                           */

// Destructor :  Destroy all the openMP locks
BensonGraphEdgeMP::~BensonGraphEdgeMP(){
  for (TNGraph::TNodeI VI = graph->BegNI(); VI < graph->EndNI(); VI++) {
    int vi = VI.GetId();
    omp_destroy_lock(&TLocks[vi]);
  }
}

BensonGraphEdgeMP::BensonGraphEdgeMP() : BensonGraphEdge() {
  TLocks = TVec<omp_lock_t>(0);
}

BensonGraphEdgeMP::BensonGraphEdgeMP(PNGraph& one_graph,MotifTypeEdge& one_motif,const TInt& nThreads,bool verb) : BensonGraphEdge(one_graph,one_motif,false,verb), nbProc(nThreads) {

  weight = WeightVHM(graph->GetMxNId());
  TLocks = TVec<omp_lock_t>(graph->GetMxNId()); // Create and init one lock per node
  for (TNGraph::TNodeI VI = graph->BegNI(); VI < graph->EndNI(); VI++) {
    int vi = VI.GetId();
    weight[vi] = THash<TInt,TInt>();
    omp_lock_t one_lock;
    omp_init_lock(&one_lock);
    TLocks[vi] = one_lock;
  }

  int tt = motif.getType();
  switch (tt){
    case DIR:
      DirAdjacency();
      break;
    case SYM:
      SymAdjacency();
      break;
    default:
      TExcept::Throw("Unknown edge type.");
  }
}

/******************************************************************************/
/*                             Symmetrised Graph                              */

// For incrementing weight in the symmetrised graph, avoiding concurrential accesses
void BensonGraphEdgeMP::IncrementWeight(int i, int j) {
  { omp_lock_t& lock_i = TLocks[i];
    omp_set_lock(&lock_i);
    THash<TInt, TInt>& edge_weights = weight[i];
    edge_weights(j) += 1;
    omp_unset_lock(&lock_i);
  }
  {
    omp_lock_t& lock_j = TLocks[j];
    omp_set_lock(&lock_j);
    THash<TInt, TInt>& edge_weights = weight[j];
    edge_weights(i) += 1;
    omp_unset_lock(&lock_j);
  }
}

void BensonGraphEdgeMP::SymAdjacency(){
  if (verbose) {printf("Building the symmetrised graph using %d threads...\n",nbProc);}
  int nbNodes = graph->GetMxNId()+1;

  omp_set_num_threads(nbProc);
  int tid(-1),cpt(0);

  #pragma omp parallel for firstprivate(cpt) private(tid) schedule(dynamic,1)
  for (int src=0;src<nbNodes;src++){ // One proc process one node, in parallel
    if (not graph->IsNode(src)){
      continue;
    }
    if (verbose){
      tid = omp_get_thread_num();
      if (cpt%10000 == 0){
        printf("Thread %d processes its %dth Node...\n",tid,cpt);
      }
      cpt ++;
    }
    TNGraph::TNodeI NI = graph->GetNI(src);
    SymProcOneNode(NI,src);
  }
}


/******************************************************************************/
/*                               Directed Graph                               */

// For incrementing weight in the directed graph, avoiding concurrential accesses
void BensonGraphEdgeMP::IncrementDirWeight(int i, int j) {
  { omp_lock_t& lock_i = TLocks[i];
    omp_set_lock(&lock_i);
    THash<TInt, TInt>& edge_weights = weight[i];
    edge_weights(j) += 1;
    omp_unset_lock(&lock_i);
  }
}

void BensonGraphEdgeMP::DirAdjacency(){
  if (verbose) {printf("Building the symmetrised graph using %d threads...\n",nbProc);}
  int nbNodes = graph->GetMxNId()+1;

  omp_set_num_threads(nbProc);
  int tid(-1),cpt(0);

  #pragma omp parallel for firstprivate(cpt) private(tid) schedule(dynamic,1)
  for (int src=0;src<nbNodes;src++){ // One proc process one node, in parallel
    if (not graph->IsNode(src)){
      continue;
    }
    if (verbose){
      tid = omp_get_thread_num();
      if (cpt%10000 == 0){
        printf("Thread %d processes its %dth Node...\n",tid,cpt);
      }
      cpt ++;
    }
    TNGraph::TNodeI NI = graph->GetNI(src);
    DirProcOneNode(NI,src);
  }
}
