//// CONTAINS THE SPEC OF THE MOTHER CLASS OF BENSONGRAPHS (MAM),
// **************BensonGraph
// This is a template class to enable handle undirected graphs in the future
// (Snap directed and non directed graphs does not have same types)

#ifndef bensongraph_h
#define bensongraph_h
#include "motifType.h"

// Typedef to formally contain the MAM
typedef TVec< THash<TInt, TInt> > WeightVHM;


template<class TypeGraph>
class BensonGraph{
  protected:
    const TPt<TypeGraph> graph; // the initial graph
    WeightVHM weight; // the MAM
    TIntV order; // order of nodes, sorted by degree
    TBool verbose; // Displaying info when program is running, if true

    void DegreeOrdering(); // function to sort nodes by degree

    virtual void IncrementWeight(int const src,int const tgt); // to increment weight when an instance of a graphlet is found

  public:
    BensonGraph(TPt<TypeGraph>& one_graph,bool& decomp,bool verb=false);
    BensonGraph();
    virtual ~BensonGraph();
    void SaveAsTxt(const TStr& file) const; // To save the MAM in a text file
    void SaveAsLv(const TStr& file) const; //  to save the MAM in two binary files (structure and wieghts)
    WeightVHM& GetWeight();
    WeightVHM& GetWeight(int& singles);

};


template<class TypeGraph>
BensonGraph<TypeGraph>::~BensonGraph(){

}

template<class TypeGraph>
BensonGraph<TypeGraph>::BensonGraph() : graph(NULL)
{
  weight = WeightVHM(0);
  order = TIntV(0);
}

template<class TypeGraph>
BensonGraph<TypeGraph>::BensonGraph(TPt<TypeGraph>& one_graph,bool& decomp, bool verb) : graph(one_graph), verbose(verb)
{
  if (decomp){
    weight = WeightVHM(graph->GetMxNId());
    for (typename TypeGraph::TNodeI VI = graph->BegNI(); VI < graph->EndNI(); VI++) {
      int vi = VI.GetId();
      weight[vi] = THash<TInt,TInt>();
    }
  }
  DegreeOrdering();
}

template<class TypeGraph>
void BensonGraph<TypeGraph>::IncrementWeight(int i, int j) {
  // int minval = MIN(i, j);
  // int maxval = MAX(i, j);
  {THash<TInt, TInt>& edge_weights = weight[i];
  edge_weights(j) += 1;}
  {THash<TInt, TInt>& edge_weights = weight[j];
  edge_weights(i) += 1;}
}


// Extracted from Benson Code
template<class TypeGraph>
void BensonGraph<TypeGraph>::DegreeOrdering() {
  // Note: This is the most efficient when the nodes are numbered
  // 0, 1, ..., num_nodes - 1.
  int max_nodes = graph->GetMxNId() + 1;
  TVec< TKeyDat<TInt, TInt> > degrees(max_nodes);
  degrees.PutAll(TKeyDat<TInt, TInt>(0, 0));
  // Set the degree of a node to be the number of nodes adjacent to the node in
  // the undirected graph.
  for (typename TypeGraph::TNodeI NI = graph->BegNI(); NI < graph->EndNI(); NI++) {
    int src = NI.GetId();
    int num_nbrs = NI.GetOutDeg();
    // For each incoming edge that is not an outgoing edge, include it.
    for (int i = 0; i < NI.GetInDeg(); ++i) {
      int dst = NI.GetInNId(i);
      if (!NI.IsOutNId(dst)) {
        ++num_nbrs;
      }
    }
    degrees[src] = TKeyDat<TInt, TInt>(num_nbrs, src);
  }

  degrees.Sort();
  order = TIntV(graph->GetMxNId());
  for (int i = 0; i < degrees.Len(); ++i) {
      order[degrees[i].Dat] = i;
  }
}

template<class TypeGraph>
void BensonGraph<TypeGraph>::SaveAsTxt(const TStr& txt_output) const{
  int nnbNodes(0);
  int nbEdges(0);
  FILE *TxtOut = fopen(txt_output.CStr(),"w");

  for (typename TypeGraph::TNodeI NI = graph->BegNI(); NI < graph->EndNI();NI++){
    uint node = NI.GetId();
    const THash<TInt, TInt>& edge_list = weight[node];
    uint deg = edge_list.Len();
    if (deg>0){
      nnbNodes++;
      nbEdges+=deg;
    }

    for (THash<TInt, TInt>::TIter it = edge_list.BegI(); it < edge_list.EndI();it++) {
      fprintf(TxtOut,"%d %d %d\n",node,it->Key,it->Dat);
    }
  }
  fclose(TxtOut);
  if (verbose){
    printf("MAM ::: Number of nodes : %d; Number of edges : %d\n",nnbNodes,nbEdges);
  }
}

template<class TypeGraph>
void BensonGraph<TypeGraph>::SaveAsLv(const TStr& lv_output) const{
  int nnbNodes(0);
  int nbEdges(0);
  FILE *GraphOut = fopen((lv_output+".bin").CStr(),"wb");
  FILE *WeightOut = fopen((lv_output+".weight").CStr(),"wb");
  uint nbNodes = graph->GetNodes();
  fwrite((char*)(&nbNodes),sizeof(nbNodes),1,GraphOut);
  for (uint node = 0; node < weight.Len(); node ++){
    const THash<TInt, TInt>& edge_list = weight[node];
    uint deg = edge_list.Len();
    if (deg>0){
      nnbNodes++;
      nbEdges+=deg;
    }

    fwrite((char*)(&deg),sizeof(deg),1,GraphOut);

    for (THash<TInt, TInt>::TIter it = edge_list.BegI(); it < edge_list.EndI();it++) {
      uint nei = (uint)it->Key;
      uint weight = (uint)it->Dat;
      fwrite((char*)(&nei),sizeof(nei),1,GraphOut);
      fwrite((char*)(&weight),sizeof(weight),1,WeightOut);
    }
  }
  fclose(GraphOut);
  fclose(WeightOut);
  if (verbose){
    printf("MAM ::: Number of nodes : %d; Number of edges : %d\n",nnbNodes,nbEdges);
  }
}

template<class TypeGraph>
WeightVHM& BensonGraph<TypeGraph>::GetWeight(){
  return weight;
}

template<class TypeGraph>
WeightVHM& BensonGraph<TypeGraph>::GetWeight(int& singles){
  singles = 0;
  for (typename TypeGraph::TNodeI NI = graph->BegNI(); NI < graph->EndNI();NI++){
    int src = NI.GetId();
    if (weight[src].Len() == 0){
      singles++;
    }
  }
  return weight;
}

#endif // bensongraph_h
