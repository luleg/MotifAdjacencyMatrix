//// CONTAINS THE SPEC OF THE DAUGHTER CLASS OF BENSONGRAPH (MAM) WHEN THE MOTIF
// HAS 3 NODES
// **************BensonGraphTriangle
#ifndef bensongraphtriangle_h
#define bensongraphtriangle_h

#include "BensonGraph.h"


class BensonGraphTriangle : public BensonGraph<TNGraph>{
private:
  void TriangleMotifAdjacency(); // To build the MAM when the motif is closed (triangle)
  void WedgeMotifAdjacency(); // To build the MAM when the motif is open (wedge)

protected:
  MotifTypeTriangle motif;
  void TriangleProcOneNode(const TNGraph::TNodeI& NI, const int& src, const int& src_pos);
  void WedgeProcOneNode(const TNGraph::TNodeI& NI, const int& center);
  virtual void IncrementWeight(int const i, int const j, int const k);

public:
  BensonGraphTriangle(PNGraph& one_graph, MotifTypeTriangle& one_motif,bool decomp, bool verb =false);
  BensonGraphTriangle();
  virtual ~BensonGraphTriangle();
};

/// For MultiProcessing computation --- Only need to be carefull when updating the MAM

class BensonGraphTriangleMP : public BensonGraphTriangle{
private:
  int nbProc;
  TVec<omp_lock_t> TLocks;

  void TriangleMotifAdjacency();
  void WedgeMotifAdjacency();
  virtual void IncrementWeight(int const i, int const j, int const k);

public:
  BensonGraphTriangleMP(PNGraph& one_graph, MotifTypeTriangle& one_motif, const TInt& nThreads, bool verb = false);
  BensonGraphTriangleMP();
  virtual ~BensonGraphTriangleMP();
};

#endif  //  bensongraphtriangle_h
