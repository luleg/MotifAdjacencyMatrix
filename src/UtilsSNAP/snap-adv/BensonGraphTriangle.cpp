#include "BensonGraphTriangle.h"

// SEQUENTIAL //////////////////////////////////////////////////////////////////
/******************************************************************************/
/*                         Some Utilities Functions                           */

BensonGraphTriangle::~BensonGraphTriangle(){
  motif = MotifTypeTriangle(); // To unhandle the motif, which is linked to a functor
}

// Dummy Constructor
BensonGraphTriangle::BensonGraphTriangle() : BensonGraph(), motif() {
}


BensonGraphTriangle::BensonGraphTriangle(PNGraph& one_graph,MotifTypeTriangle& one_motif,bool decomp,bool verb) : BensonGraph(one_graph,decomp,verb), motif(one_motif) {

  motif.SetGraph(graph);

  if (decomp){// if decomp is true, the MAM is computed in the constructor
    int tt = motif.getType();
    switch (tt){
      case TRIANGLE:
        TriangleMotifAdjacency();
        break;
      case WEDGE:
        WedgeMotifAdjacency();
        break;
      default:
        TExcept::Throw("Unknown triangle type.");
    }
  }
}

// For incrementing weight in the MAM
void BensonGraphTriangle::IncrementWeight(int const i, int const j, const int k){
  {
    THash<TInt, TInt>& edge_weights = weight[i];
    edge_weights(j) += 1;
    edge_weights(k) += 1;
  }
  {
    THash<TInt, TInt>& edge_weights = weight[j];
    edge_weights(i) += 1;
    edge_weights(k) += 1;
  }
  {
    THash<TInt, TInt>& edge_weights = weight[k];
    edge_weights(i) += 1;
    edge_weights(j) += 1;
  }
}
/******************************************************************************/
/*                        Closed Triangle Weighting                           */

void BensonGraphTriangle::TriangleProcOneNode(const TNGraph::TNodeI& NI, const int& src, const int& src_pos){
  // Algorithm from Chiba and NISHIZEKI for listing triangles
  TIntV neighbors_higher;
  for (int i = 0; i < NI.GetOutDeg(); i++) {
    int nbr = NI.GetOutNId(i);
    if (order[nbr] > src_pos) {
      neighbors_higher.Add(nbr);
    }
  }
  for (int i = 0; i < NI.GetInDeg(); i++) {
    int nbr = NI.GetInNId(i);
    if (!NI.IsOutNId(nbr) && order[nbr] > src_pos) {
      neighbors_higher.Add(nbr);
    }
  }

  for (int ind1 = 0; ind1 < neighbors_higher.Len(); ind1++) {
    for (int ind2 = ind1 + 1; ind2 < neighbors_higher.Len(); ind2++) {
      int dst1 = neighbors_higher[ind1];
      int dst2 = neighbors_higher[ind2];
      // Check for triangle formation
      if (graph->IsEdge(dst1, dst2) || graph->IsEdge(dst2, dst1)) {
        bool motif_occurs = motif.IsInstance(src, dst1, dst2);
        if (motif_occurs) { // we update MAM if the ibnduced subgraph is an instance of the motif
          IncrementWeight(src,  dst1, dst2);
        }
      }
    }
  }
}


void BensonGraphTriangle::TriangleMotifAdjacency() {
  // Algorithm from Chiba and NISHIZEKI for listing triangles
  if (verbose) {printf("Building the MAM of a closed triangle sequentially...\n");}
  int cpt = 0;
  for (TNGraph::TNodeI NI = graph->BegNI(); NI < graph->EndNI(); NI++) {
    if (verbose){
      if (cpt%10000 == 0){
        printf("Processing node %d\n",cpt);
      }
      cpt++;
    }
    int src = NI.GetId();
    int src_pos = order[src];
    TriangleProcOneNode(NI, src,src_pos);
  }
}

/******************************************************************************/
/*                            Wedges Weighting                               */

void BensonGraphTriangle::WedgeProcOneNode(const TNGraph::TNodeI& NI, const int& center){
  TIntV neighbors;
  for (int i = 0; i < NI.GetOutDeg(); i++) {
    int nbr = NI.GetOutNId(i);
    neighbors.Add(nbr);
  }
  for (int i = 0; i < NI.GetInDeg(); i++) {
    int nbr = NI.GetInNId(i);
    if (!NI.IsOutNId(nbr)) {
      neighbors.Add(nbr);
    }
  }

  for (int ind1 = 0; ind1 < neighbors.Len(); ind1++) {
    for (int ind2 = ind1 + 1; ind2 < neighbors.Len(); ind2++) {
      int dst1 = neighbors[ind1];
      int dst2 = neighbors[ind2];
      bool motif_occurs = motif.IsInstance(center, dst1, dst2);
      // Increment weights of (center, dst1, dst2) if it occurs.
      if (motif_occurs) {
        IncrementWeight(center,  dst1, dst2);
      }
    }
  }
}


void BensonGraphTriangle::WedgeMotifAdjacency() {
  if (verbose) {printf("Building the MAM of a wedge sequentially...\n");}
  int cpt = 0;
  for (TNGraph::TNodeI NI = graph->BegNI(); NI < graph->EndNI(); NI++) {
    if (verbose){
      if (cpt%10000 == 0){
        printf("Processing node %d\n",cpt);
      }
      cpt++;
    }

    int center = NI.GetId();
    WedgeProcOneNode(NI, center);
  }
}

// MULTITHREADED ///////////////////////////////////////////////////////////////
/******************************************************************************/
/*                         Some Utilities Functions                           */

 // Destructor :  Destroy all the openMP locks
BensonGraphTriangleMP::~BensonGraphTriangleMP(){
  for (TNGraph::TNodeI VI = graph->BegNI(); VI < graph->EndNI(); VI++) {
    int vi = VI.GetId();
    omp_destroy_lock(&TLocks[vi]);
  }
}

BensonGraphTriangleMP::BensonGraphTriangleMP() : BensonGraphTriangle() {
  TLocks = TVec<omp_lock_t>(0);
}

BensonGraphTriangleMP::BensonGraphTriangleMP(PNGraph& one_graph,MotifTypeTriangle& one_motif,const TInt& nThreads, bool verb) : BensonGraphTriangle(one_graph,one_motif,false,verb), nbProc(nThreads) {

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
    case TRIANGLE:
      TriangleMotifAdjacency();
      break;
    case WEDGE:
      WedgeMotifAdjacency();
      break;
    default:
      TExcept::Throw("Unknown triangle type.");
  }

}

// For incrementing weight in the MAM avoiding concurrential accesses :
void BensonGraphTriangleMP::IncrementWeight(int const i, int const j, const int k){
  {
    omp_lock_t& lock_i = TLocks[i];
    omp_set_lock(&lock_i);
    THash<TInt, TInt>& edge_weights = weight[i];
    edge_weights(j) += 1;
    edge_weights(k) += 1;
    omp_unset_lock(&lock_i);
  }
  {
    omp_lock_t& lock_j = TLocks[j];
    omp_set_lock(&lock_j);
    THash<TInt, TInt>& edge_weights = weight[j];
    edge_weights(i) += 1;
    edge_weights(k) += 1;
    omp_unset_lock(&lock_j);
  }
  {
    omp_lock_t& lock_k = TLocks[k];
    omp_set_lock(&lock_k);
    THash<TInt, TInt>& edge_weights = weight[k];
    edge_weights(i) += 1;
    edge_weights(j) += 1;
    omp_unset_lock(&lock_k);
  }
}

/******************************************************************************/
/*                        Closed Trinagles Weighting                          */
void BensonGraphTriangleMP::TriangleMotifAdjacency() {
  if (verbose) {printf("Building MAM of closed triangle using %d threads...\n",nbProc);}
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
    int src_pos = order[src];
    TriangleProcOneNode(NI, src,src_pos);
  }
}


/******************************************************************************/
/*                            Wedges Weighting                               */

void BensonGraphTriangleMP::WedgeMotifAdjacency() {
  if (verbose) {printf("Building MAM of a wedge using %d threads...\n",nbProc);}

  int nbNodes = graph->GetMxNId()+1;

  omp_set_num_threads(nbProc);
  int tid(-1),cpt(0);

  #pragma omp parallel for private(tid) schedule(dynamic,1)
  for (int center=0;center<nbNodes;center++){ // One proc process one node, in parallel

    if (not graph->IsNode(center)){
      continue;
    }
    if (verbose){
      tid = omp_get_thread_num();
      if (cpt%10000 == 0){
        printf("Thread %d processes its %dth Node...\n",tid,cpt);
      }
      cpt ++;
    }

    TNGraph::TNodeI NI = graph->GetNI(center);
    WedgeProcOneNode(NI, center);
  }

}
