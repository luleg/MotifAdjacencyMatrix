//// CONTAINS THE SPEC OF THE DAUGHTER CLASS OF BENSONGRAPH (MAM) WHEN "MOTIF"
// IS AN EDGE
// **************BensonGraphEdge
#ifndef bensongraphedge_h
#define bensongraphedge_h

#include "BensonGraph.h"


class BensonGraphEdge : public BensonGraph<TNGraph>{
private:
  void SymAdjacency(); // to build the MAM if the symmetrised version of the graph is asked
  void DirAdjacency(); // to build the MAM if the directed version of the graph is asked

protected:
  MotifTypeEdge motif;
  void SymProcOneNode(const TNGraph::TNodeI& NI,const int& src);
  void DirProcOneNode(const TNGraph::TNodeI& NI,const int& src);
  virtual void IncrementWeight(int const src, int const tgt);
  virtual void IncrementDirWeight(int const src, int const tgt);

public:
  BensonGraphEdge(PNGraph& one_graph, MotifTypeEdge& one_motif,bool decomp, bool verb = false);
  BensonGraphEdge();
  virtual ~BensonGraphEdge();
};

/// For MultiProcessing computation --- Only need to be carefull when updating the MAM
class BensonGraphEdgeMP : public BensonGraphEdge{
private:
  int nbProc;
  TVec<omp_lock_t> TLocks;

  void SymAdjacency();
  void DirAdjacency();
  virtual void IncrementWeight(int const src, int const tgt);
  virtual void IncrementDirWeight(int const src, int const tgt);

public:
  BensonGraphEdgeMP(PNGraph& one_graph, MotifTypeEdge& one_motif, const TInt& nThreads, bool verb = false);
  BensonGraphEdgeMP();
  virtual ~BensonGraphEdgeMP();
};

#endif  //  bensongraphedge_h
