#include "BensonGraphQuad.h"


// SEQUENTIAL //////////////////////////////////////////////////////////////////
/******************************************************************************/
/*                         Some Utilities Functions                           */

// A destructor
BensonGraphQuad::~BensonGraphQuad(){
  // printf("~BensonGraphQuad destructor\n");
}
// A constructor
BensonGraphQuad::BensonGraphQuad() : BensonGraph(), motif() {
}

// the real constructor
BensonGraphQuad::BensonGraphQuad(PNGraph& one_graph,MotifTypeQuad& one_motif,bool decomp, bool verb) : BensonGraph(one_graph,decomp,verb), motif(one_motif), ungraph(TSnap::ConvertGraph<PUNGraph>(one_graph)) {

  if (decomp){
    int qq = motif.getType();
    switch (qq){
      case 4:
        Q4MotifAdjacency();
        break;
      case 5:
        Q5MotifAdjacency();
        break;
      case 6:
        Q6MotifAdjacency();
        break;
      default:
        TExcept::Throw("Unknown Quadrangle.");
    }
  }
}

// For incrementing weight in the MAM
void BensonGraphQuad::IncrementAllWeight(int const v,int const w,int const n1,int const n2) {
  {
    THash<TInt, TInt>& edge_weights = weight[v];
    edge_weights(w) += 1;
    edge_weights(n1) += 1;
    edge_weights(n2) += 1;
  }{
    THash<TInt, TInt>& edge_weights = weight[w];
    edge_weights(v) += 1;
    edge_weights(n1) += 1;
    edge_weights(n2) += 1;
  }{
    THash<TInt, TInt>& edge_weights = weight[n1];
    edge_weights(v) += 1;
    edge_weights(w) += 1;
    edge_weights(n2) += 1;
  }{
    THash<TInt, TInt>& edge_weights = weight[n2];
    edge_weights(v) += 1;
    edge_weights(w) += 1;
    edge_weights(n1) += 1;
  }
}
/* Checks if the kind of edge (uni/bi) between a node and anotehr one is known.
Computes and stores it if not, and returns it */
const void BensonGraphQuad::UpdateOneEdge(TIntV& TwW, const int& w_id, const int& v_id, int& RelWV){
  if (TwW.Last() > -2){
    RelWV = 1 * ( graph->GetNI(w_id).IsOutNId(v_id) ? 1 : 0) +
            2 * ( graph->GetNI(w_id).IsInNId(v_id) ? 1 : 0);

    TwW.Add(RelWV);
    TwW.Add(0);
    TwW.Add(-(w_id+2));
  }
  else {
    RelWV = TwW[TwW.Len()-3];
  }
}
/* Checks if the kind of edges (uni/bi) between a node and two others is known.
Computes and stores them, if not, and returns them */
const void BensonGraphQuad::UpdateTwoEdges(TIntV& TwU, const int& u_id, const int& v_id, const int& w_id, int& RelUV, int& RelUW){
  if (TwU.Empty() || TwU.Last() > -2){
    RelUV = 1 * ( graph->GetNI(u_id).IsOutNId(v_id) ? 1 : 0) +
            2 * ( graph->GetNI(u_id).IsInNId(v_id) ? 1 : 0);

    if (order[u_id]<order[w_id]){
      RelUW += 1 * ( graph->GetNI(u_id).IsOutNId(w_id) ? 1 : 0) +
               2 * ( graph->GetNI(u_id).IsInNId(w_id) ? 1 : 0);
    }
    else {
      RelUW += 1 * ( graph->GetNI(w_id).IsInNId(u_id) ? 1 : 0) +
               2 * ( graph->GetNI(w_id).IsOutNId(u_id) ? 1 : 0);
    }
    TwU.Add(RelUV);
    TwU.Add(RelUW);
    TwU.Add(-(w_id+2));
  }
  else if (TwU.Last() != -(w_id+2)){
    int NbElts = TwU.Len();
    RelUV = TwU[NbElts-3];
    if (order[u_id]<order[w_id]){
      RelUW += 1 * ( graph->GetNI(u_id).IsOutNId(w_id) ? 1 : 0) +
               2 * ( graph->GetNI(u_id).IsInNId(w_id) ? 1 : 0);
    }
    else {
      RelUW += 1 * ( graph->GetNI(w_id).IsInNId(u_id) ? 1 : 0) +
               2 * ( graph->GetNI(w_id).IsOutNId(u_id) ? 1 : 0);
    }
    TwU[NbElts-2] = RelUW;
    TwU[NbElts-1] = -(w_id+2);
  }
  else {
    int NbElts = TwU.Len();
    RelUV = TwU[NbElts-3];
    RelUW = TwU[NbElts-2];
  }
}

/******************************************************************************/
/*                       Q4 quadrangle MAM Weighting                          */

void BensonGraphQuad::Q4ProcOneNode(const int& vi, const int& lab_vi){
  // Algorithm from Niba and Nishizeki for listing every quadrangles.
  // Adapted to 4-node quadrangles of type Q4 (with 4 (undirected) edges)

  TVec< TIntV > Tw = TVec<TIntV>(graph->GetMxNId() + 1); // 2-hop neighbours of VI: listing their paths to VI
  TIntV w_match; // To store those with more than 2 paths, that  are elligible nodes

  TUNGraph::TNodeI VI = ungraph->GetNI(vi);
  // Get all neighbors of VI who come before in the ordering
  for (int s = 0; s < VI.GetDeg(); s++) {  // u : u ~ vi

    int u = VI.GetNbrNId(s);

    if ( order[u] >= lab_vi) {
      continue;
    }
    TUNGraph::TNodeI UI = ungraph->GetNI(u);

    // Get all the neighbors of UI (2-hop neighbours of VI)
    for (int tmp = 0; tmp< UI.GetDeg();tmp++){ // w : w ~ u
      int w = UI.GetNbrNId(tmp);
      bool first = false;

      // Preprocessing :
      if (Tw[w].Empty()){ // If this is the 1st time we meet this node :
        if (order[w] >= lab_vi) { // is its degree above or equal to VI ?
          Tw[w].Add(-1); // If true, w will never be an admissible node: this is
          // indicated by the negative element in the vector first coordinate
          continue;
        }
        if (ungraph->GetNI(w).IsNbrNId(vi)){ // Otherwise, is there any edge between VI and w ?
          Tw[w].Add(-1); // If true, w will never be an admissible node: this is
          // indicated by the negative element in the vector first coordinate
          continue;
        }
        // otherwise, u is the first 2-path from w to vi
        Tw[w].Add(1);
        Tw[w].Add(u);
        first = true;
      }
      int ffirst = Tw[w][0];
      if (ffirst<0){
        continue; // w is not an admissible point, we leave
      }
      if (not first){ //this means u is at least the second 2-path between VI and w
        // So we have to add w to marked nodes, if it has not been done alreadyy.
        Tw[w].Add(u);
        Tw[w][0] += 1;

        if (Tw[w][0]==2){
          w_match.Add(w);
        }
      }
    }
  }

  // Check all marked nodes
  for (int cpt=0; cpt<w_match.Len();cpt++){
    int w = w_match[cpt];
    int nbElts = Tw[w][0];

    for (int s=1;s<nbElts+1;s++){
      int u_s = Tw[w][s];
      int id_s = 0;

      for (int t=s+1;t<nbElts+1;t++){ //(u_s,u_t,vi,w) may be a graphlet:
        int u_t = Tw[w][t];

        if (order[u_s] < order[u_t]){ // if u_s ~ u_t, we leave
          if (graph->GetNI(u_s).IsNbrNId(u_t)){
            continue;
          }
        }
        else{
          if (graph->GetNI(u_t).IsNbrNId(u_s)){
            continue;
          }
        }
        // Otherwise, we check the graphlet.
        if (id_s == 0){
          int RelUsV(0),RelUsW(0);
          TIntV& TwUs = Tw[u_s];
          UpdateTwoEdges(TwUs, u_s,vi,w,RelUsV,RelUsW); //
          id_s = 8 * (RelUsV == 1 ? 1 : 0) + 4096 * (RelUsV == 2 ? 1 : 0) + 4104 * (RelUsV == 3 ? 1 : 0) +
                 4 * (RelUsW == 1 ? 1 : 0) + 256  * (RelUsW == 2 ? 1 : 0) + 260  * (RelUsW == 3 ? 1 : 0);
        }

        int RelUtV(0),RelUtW(0);
        TIntV& TwUt = Tw[u_t];
        UpdateTwoEdges(TwUt, u_t,vi,w,RelUtV,RelUtW);

        int id_st = id_s +
                    128 * (RelUtV == 1 ? 1 : 0) + 8192 * (RelUtV == 2 ? 1 : 0) + 8320 * (RelUtV == 3 ? 1 : 0) +
                     64 * (RelUtW == 1 ? 1 : 0) + 512  * (RelUtW == 2 ? 1 : 0) + 576  * (RelUtW == 3 ? 1 : 0);
        if (motif.IsInstance(id_st)){
            IncrementAllWeight(vi,w,u_s,u_t);
        }
      }
    }
  }
}

void BensonGraphQuad::Q4MotifAdjacency() {
  if (verbose) {printf("Building a Q4 MAM  sequentially...\n");}
  TIntV NIdV;
  graph->GetNIdV(NIdV);
  for (int k = 0; k < NIdV.Len(); k ++) {
    if (verbose){
      if (k%1000 == 0){
        printf("Processing node %d\n",k);
      }
    }
    int vi = NIdV[k];
    int lab_vi = order[vi];
    Q4ProcOneNode(vi,lab_vi);
  }
}

/******************************************************************************/
/*                       Q5 quadrangle MAM Weighting                          */

void BensonGraphQuad::Q5ProcOneNode(const int& vi, const int& lab_vi){
  TVec< TIntV > Tw = TVec<TIntV>(graph->GetMxNId() + 1);
  TIntV w_match;
  TUNGraph::TNodeI VI = ungraph->GetNI(vi);

  // Get all neighbors who come before in the ordering
  for (int s = 0; s < VI.GetDeg(); s++) {  // u : u ~ vi

    int u = VI.GetNbrNId(s);

    if ( order[u] >= lab_vi) {
      continue;
    }
    TUNGraph::TNodeI UI = ungraph->GetNI(u);
    for (int tmp = 0; tmp< UI.GetDeg();tmp++){ // w : w ~ u

      int w = UI.GetNbrNId(tmp);
      if (order[w] >= lab_vi) {
        continue;
      }
      Tw[w].Add(u);
      if (Tw[w].Len()==2){
        w_match.Add(w);
      }
    }
  }

  // Check all marked nodes and update weight if needed
  for (int cpt=0; cpt<w_match.Len();cpt++){

    int w = w_match[cpt];
    TIntV& wTw = Tw[w];

    int nbElts = wTw.Len();

    int RelWV(0);
    UpdateOneEdge(wTw, w,vi,RelWV);
    int id_w =  2048 * (RelWV == 1 ? 1 : 0) + 16384 * (RelWV == 2 ? 1 : 0) + 18432 * (RelWV == 3 ? 1 : 0);

    if (nbElts == wTw.Len()){
      nbElts = nbElts - 3;
    }

    for (int s=0;s<nbElts;s++){
      int u_s = Tw[w][s];
      int id_s = 0;

      for (int t=s+1;t<nbElts;t++){
        int u_t = Tw[w][t];
        int id_st = 0;

        if (id_w == 0) { // w and v are not linked, u_s and u_t thus must be
          if (order[u_s] < order[u_t]){
            id_st  = 2 * (graph->GetNI(u_s).IsOutNId(u_t) ? 1 : 0) + 16 * (graph->GetNI(u_s).IsInNId(u_t) ? 1 : 0);
          }
          else{
            id_st  = 2 * (graph->GetNI(u_t).IsInNId(u_s) ? 1 : 0) + 16* (graph->GetNI(u_t).IsOutNId(u_s) ? 1 : 0);
          }
          if (id_st == 0){ // If not we leave
            continue;
          }
        }
        else{// w and v are linked, u_s and u_t thus must not be
          if (order[u_s] < order[u_t]){
            if (graph->GetNI(u_s).IsNbrNId(u_t)){
              continue;
            }
          }
          else if (graph->GetNI(u_t).IsNbrNId(u_s)){
            continue;
          }
        }

        if (id_s == 0){
          int RelUsV(0),RelUsW(0);
          TIntV& TwUs = Tw[u_s];
          UpdateTwoEdges(TwUs, u_s,vi,w,RelUsV,RelUsW); //
          id_s = 8 * (RelUsV == 1 ? 1 : 0) + 4096 * (RelUsV == 2 ? 1 : 0) + 4104 * (RelUsV == 3 ? 1 : 0) +
                 4 * (RelUsW == 1 ? 1 : 0) + 256  * (RelUsW == 2 ? 1 : 0) + 260  * (RelUsW == 3 ? 1 : 0);
        }

        int RelUtV(0),RelUtW(0);
        TIntV& TwUt = Tw[u_t];
        UpdateTwoEdges(TwUt, u_t,vi,w,RelUtV,RelUtW);

        id_st += id_w + id_s +
                  128 * (RelUtV == 1 ? 1 : 0) + 8192 * (RelUtV == 2 ? 1 : 0) + 8320 * (RelUtV == 3 ? 1 : 0) +
                   64 * (RelUtW == 1 ? 1 : 0) + 512  * (RelUtW == 2 ? 1 : 0) + 576  * (RelUtW == 3 ? 1 : 0);

        if (motif.IsInstance(id_st)){
          IncrementAllWeight(vi,  w,  u_s, u_t);
        }
      }
    }
  }
}
void BensonGraphQuad::Q5MotifAdjacency() {
  if (verbose) {printf("Building a Q5 MAM  sequentially...\n");}
  TIntV NIdV;
  graph->GetNIdV(NIdV);
  for (int k = 0; k < NIdV.Len(); k ++) {
    if (verbose){
      if (k%1000 == 0){
        printf("Processing node %d\n",k);
      }
    }
    int vi = NIdV[k];
    int lab_vi = order[vi];
    Q5ProcOneNode(vi,lab_vi);
  }
}

/******************************************************************************/
/*                       Q6 quadrangle MAM Weighting                          */


/* Checks if the kind of edges (uni/bi) between a node and two others is known.
Computes and stores them, if not, and returns them */
// This is different than the other as the ordering between u_id, v_id and w_id
// is known in Q6. This allows some simplifications
const void BensonGraphQuad::UpdateTwoEdgesQ6(TIntV& TwU, const int& u_id, const int& v_id, const int& w_id, int& RelUV, int& RelUW){
  if (TwU.Empty() || TwU.Last() > -2){
    RelUV = 1 * ( graph->GetNI(u_id).IsOutNId(v_id) ? 1 : 0) +
            2 * ( graph->GetNI(u_id).IsInNId(v_id) ? 1 : 0);

    RelUW += 1 * ( graph->GetNI(u_id).IsOutNId(w_id) ? 1 : 0) +
               2 * ( graph->GetNI(u_id).IsInNId(w_id) ? 1 : 0);

    TwU.Add(RelUV);
    TwU.Add(RelUW);
    TwU.Add(-(w_id+2));
  }
  else if (TwU.Last() != -(w_id+2)){
    int NbElts = TwU.Len();
    RelUV = TwU[NbElts-3];
    RelUW += 1 * ( graph->GetNI(u_id).IsOutNId(w_id) ? 1 : 0) +
             2 * ( graph->GetNI(u_id).IsInNId(w_id) ? 1 : 0);

    TwU[NbElts-2] = RelUW;
    TwU[NbElts-1] = -(w_id+2);
  }
  else {
    int NbElts = TwU.Len();
    RelUV = TwU[NbElts-3];
    RelUW = TwU[NbElts-2];
  }
}

void BensonGraphQuad::Q6ProcOneNode(const int& vi, const int& lab_vi){
  TVec< TIntV > Tw = TVec<TIntV>(graph->GetMxNId() + 1);
  TIntV w_match;


  TUNGraph::TNodeI VI = ungraph->GetNI(vi);

  // Get all neighbors who come before in the ordering
  for (int s = 0; s < VI.GetDeg(); s++) {  // u : u ~ vi

    int u = VI.GetNbrNId(s);

    int lab_u = order[u];

    if ( lab_u >= lab_vi) {
      continue;
    }
    TUNGraph::TNodeI UI = ungraph->GetNI(u);

    // Get all the neighbors of UI (2-hop neighbours of VI)
    for (int tmp = 0; tmp< UI.GetDeg();tmp++){ // w : w ~ u
      int w = UI.GetNbrNId(tmp);

      if (order[w] <= lab_u){ // We need this condition, otherwise we will coutn each motif several times.
        continue;
      }

      TIntV& wTw = Tw[w];
      bool first=false;

      // Preprocessing :
      if (wTw.Empty()){ // If this is the 1st time we meet this node :
        if (order[w] >= lab_vi) { // is its degree above or equal to U ?
          wTw.Add(-1); // If true, w will never be an admissible node: this is
          // indicated by the negative element in the vector first coordinate
          continue;
        }
        if (ungraph->GetNI(w).IsNbrNId(vi)){ // Otherwise, is there an edge between VI and w ?
          // If true, u is the first 2-path from w to vi
          wTw.Add(1);
          wTw.Add(u);
          first=true;
        }else{
          wTw.Add(-1); // otherwise, w will never be an admissible node: this is
          // indicated by the negative element in the vector first coordinate
          continue;
        }
      }
      int ffirst = wTw[0];
      if (ffirst<0){
        continue; // w is not an admissible point, we leave
      }
      else if (not first){ // if this is true, this means u is at least the second 2-path between VI and w
        // So we have to add w to marked nodes, if it has not been done alreadyy.
        wTw[0] +=1;
        wTw.Add(u);
        if (wTw[0]==2){
          w_match.Add(w);
        }
      }
    }
  }

  // Check all marked nodes and update weight if needed

  for (int cpt=0; cpt<w_match.Len();cpt++){
    int w = w_match[cpt];
    int nbElts = Tw[w][0];
    int id_w = 0;

    for (int s=1;s<nbElts+1;s++){
      int u_s = Tw[w][s];
      int id_s = 0;

      for (int t=s+1;t<nbElts+1;t++){
        int id_st = 0;
        int u_t = Tw[w][t];

        if (order[u_s] < order[u_t]){
          id_st = 2 * (graph->GetNI(u_s).IsOutNId(u_t) ? 1 : 0) + 16 * (graph->GetNI(u_s).IsInNId(u_t) ? 1 : 0);
        }
        else{
          id_st = 2 * (graph->GetNI(u_t).IsInNId(u_s) ? 1 : 0) + 16 * (graph->GetNI(u_t).IsOutNId(u_s) ? 1 : 0);
        }
        if (id_st == 0){// if u_s ~ u_t, we leave
          continue;
        }

        if (id_w==0){
          int RelWV(0);
          TIntV& TwW = Tw[w];
          UpdateOneEdge(TwW, w,vi,RelWV);
          id_w = 2048 * (RelWV == 1 ? 1 : 0) + 16384 * (RelWV == 2 ? 1 : 0) + 18432 * (RelWV == 3 ? 1 : 0);
        }

        if (id_s == 0){
          int RelUsV(0),RelUsW(0);
          TIntV& TwUs = Tw[u_s];
          UpdateTwoEdgesQ6(TwUs, u_s,vi,w,RelUsV,RelUsW); //
          id_s  = 8 * (RelUsV == 1 ? 1 : 0) + 4096 * (RelUsV == 2 ? 1 : 0) + 4104 * (RelUsV == 3 ? 1 : 0) +
                  4 * (RelUsW == 1 ? 1 : 0) + 256  * (RelUsW == 2 ? 1 : 0) + 260  * (RelUsW == 3 ? 1 : 0);
        }

        int RelUtV(0),RelUtW(0);
        TIntV& TwUt = Tw[u_t];
        UpdateTwoEdgesQ6(TwUt, u_t,vi,w,RelUtV,RelUtW);

        id_st += id_w + id_s +
                  128 * (RelUtV == 1 ? 1 : 0) + 8192 * (RelUtV == 2 ? 1 : 0) + 8320 * (RelUtV == 3 ? 1 : 0) +
                   64 * (RelUtW == 1 ? 1 : 0) + 512  * (RelUtW == 2 ? 1 : 0) + 576  * (RelUtW == 3 ? 1 : 0);
        if (motif.IsInstance(id_st)){
            IncrementAllWeight(vi,w,u_s,u_t);
        }
      }
    }
  }
}
void BensonGraphQuad::Q6MotifAdjacency() {
  if (verbose) {printf("Building a Q6 MAM  sequentially...\n");}
  TIntV NIdV;
  graph->GetNIdV(NIdV);
  for (int k = 0; k < NIdV.Len(); k ++) {
    if (verbose){
      if (k%1000 == 0){
        printf("Processing node %d\n",k);
      }
    }
    int vi = NIdV[k];
    int lab_vi = order[vi];
    Q6ProcOneNode(vi,lab_vi);
  }
}

// MULTITHREADED ///////////////////////////////////////////////////////////////
/******************************************************************************/
/*                         Some Utilities Functions                           */

 // Destructor :  Destroy all the openMP locks
BensonGraphQuadMP::~BensonGraphQuadMP(){
  // printf("~BensonGraphQuadMP destructor\n");
  for (TNGraph::TNodeI VI = graph->BegNI(); VI < graph->EndNI(); VI++) {
    int vi = VI.GetId();
    omp_destroy_lock(&TLocks[vi]);
  }
}

BensonGraphQuadMP::BensonGraphQuadMP() : BensonGraphQuad() {
  TLocks = TVec<omp_lock_t>(0);
}

BensonGraphQuadMP::BensonGraphQuadMP(PNGraph& one_graph,MotifTypeQuad& one_motif,const TInt& nThreads, bool verb) : BensonGraphQuad(one_graph,one_motif,false,verb), nbProc(nThreads) {

  weight = WeightVHM(graph->GetMxNId());
  TLocks = TVec<omp_lock_t>(graph->GetMxNId());
  for (TNGraph::TNodeI VI = graph->BegNI(); VI < graph->EndNI(); VI++) {
    int vi = VI.GetId();
    weight[vi] = THash<TInt,TInt>();
    omp_lock_t one_lock;
    omp_init_lock(&one_lock);
    TLocks[vi] = one_lock;
  }

  int qq = motif.getType();
  switch (qq){
    case 4:
      Q4MotifAdjacency();
      break;
    case 5:
      Q5MotifAdjacency();
      break;
    case 6:
      Q6MotifAdjacency();
      break;
    default:
      TExcept::Throw("Unknown Quadrangle.");
  }
}

// For incrementing weight in the MAM avoiding concurrential accesses :
void BensonGraphQuadMP::IncrementAllWeight(int const v,int const w,int const n1,int const n2) {
  {
    omp_lock_t& lock_v = TLocks[v];
    omp_set_lock(&lock_v);
    THash<TInt, TInt>& edge_weights = weight[v];
    edge_weights(w) += 1;
    edge_weights(n1) += 1;
    edge_weights(n2) += 1;
    omp_unset_lock(&lock_v);
  }{
    omp_lock_t& lock_w = TLocks[w];
    omp_set_lock(&lock_w);
    THash<TInt, TInt>& edge_weights = weight[w];
    edge_weights(v) += 1;
    edge_weights(n1) += 1;
    edge_weights(n2) += 1;
    omp_unset_lock(&lock_w);
  }{
    omp_lock_t& lock_n1 = TLocks[n1];
    omp_set_lock(&lock_n1);
    THash<TInt, TInt>& edge_weights = weight[n1];
    edge_weights(v) += 1;
    edge_weights(w) += 1;
    edge_weights(n2) += 1;
    omp_unset_lock(&lock_n1);
  }{
    omp_lock_t& lock_n2 = TLocks[n2];
    omp_set_lock(&lock_n2);
    THash<TInt, TInt>& edge_weights = weight[n2];
    edge_weights(v) += 1;
    edge_weights(w) += 1;
    edge_weights(n1) += 1;
    omp_unset_lock(&lock_n2);
  }
}

/******************************************************************************/
/*                       Q4 quadrangle MAM Weighting                          */
void BensonGraphQuadMP::Q4MotifAdjacency() {
  // For parallelism
  if (verbose) {printf("Building a Q4 MAM using %d threads...\n",nbProc);}

  TIntV NIdV;
  graph->GetNIdV(NIdV);

  omp_set_num_threads(nbProc);
  int tid(-1),cptNode(0);

  #pragma omp parallel for firstprivate(cptNode) private(tid) schedule(dynamic,1)
  for (int k=0;k<NIdV.Len();k++){
    if (verbose){
      if (cptNode%1000 == 0){
        tid = omp_get_thread_num();
        printf("Thread %d Processing its %dth node...\n",tid,cptNode);
      }
      cptNode++;
    }
    int vi = NIdV[k];

    int lab_vi = order[vi];


    Q4ProcOneNode(vi,lab_vi);

  }
}

/******************************************************************************/
/*                       Q5 quadrangle MAM Weighting                          */
void BensonGraphQuadMP::Q5MotifAdjacency() {
  // For parallelism
  if (verbose) {printf("Building a Q5 MAM using %d threads...\n",nbProc);}

  TIntV NIdV;
  graph->GetNIdV(NIdV);

  omp_set_num_threads(nbProc);
  int tid(-1),cptNode(0);

  #pragma omp parallel for firstprivate(cptNode) private(tid) schedule(dynamic,1)
  for (int k=0;k<NIdV.Len();k++){
    if (verbose){
      if (cptNode%1000 == 0){
        tid = omp_get_thread_num();
        printf("Thread %d Processing its %dth node...\n",tid,cptNode);
      }
      cptNode++;
    }
    int vi = NIdV[k];
    int lab_vi = order[vi];

    Q5ProcOneNode(vi,lab_vi);
  }
}

/******************************************************************************/
/*                       Q6 quadrangle MAM Weighting                          */
void BensonGraphQuadMP::Q6MotifAdjacency() {
  // For parallelism
  if (verbose) {printf("Building a Q6 MAM using %d threads...\n",nbProc);}

  TIntV NIdV;
  graph->GetNIdV(NIdV);

  omp_set_num_threads(nbProc);
  int tid(-1),cptNode(0);

  #pragma omp parallel for firstprivate(cptNode) private(tid) schedule(dynamic,1)
  for (int k=0;k<NIdV.Len();k++){
    if (verbose){
      if (cptNode%1000 == 0){
        tid = omp_get_thread_num();
        printf("Thread %d Processing its %dth node...\n",tid,cptNode);
      }
      cptNode++;
    }
    int vi = NIdV[k];
    int lab_vi = order[vi];

    Q6ProcOneNode(vi,lab_vi);
  }
}
