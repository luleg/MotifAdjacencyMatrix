//// CONTAINS THE SPEC OF THE DAUGHTER CLASS OF BENSONGRAPH (MAM) WHEN THE MOTIF
// IS A QUADRANGLE (4-NODE MOTIF THAT CONTAINS A 4-NODE LOOP)
// **************BensonGraphQuad
#ifndef bensongraphquad_h
#define bensongraphquad_h

#include "BensonGraph.h"


class BensonGraphQuad : public BensonGraph<TNGraph>{
  private :
    PUNGraph ungraph;
    void Q4MotifAdjacency();
    void Q5MotifAdjacency();
    void Q6MotifAdjacency();

  protected:
    MotifTypeQuad motif;

    void Q4ProcOneNode(const int& vi, const int& lab_vi);
    void Q5ProcOneNode(const int& vi, const int& lab_vi);
    void Q6ProcOneNode(const int& vi, const int& lab_vi);

/* Utlise Functions */
    const void UpdateTwoEdges(TIntV& TwU, const int& u_id, const int& v_id, const int& w_id,
                            int& RelUV, int& RelUW);

    const void UpdateTwoEdgesQ6(TIntV& TwU, const int& u_id, const int& v_id, const int& w_id,
                            int& RelUV, int& RelUW);

    const void UpdateOneEdge(TIntV& TwW, const int& w_id, const int& v_id, int& RelWV);

    virtual void IncrementAllWeight(int const v,int const w,int const n1,int const n2);

  public:
    BensonGraphQuad(PNGraph& one_graph, MotifTypeQuad& one_motif,bool decomp, bool verb = false);
    BensonGraphQuad();
    virtual ~BensonGraphQuad();
  };

/// For MultiProcessing computation --- Only need to be carefull when updating the MAM
class BensonGraphQuadMP : public BensonGraphQuad{
  private:
    int nbProc;
    TVec<omp_lock_t> TLocks;

    void Q4MotifAdjacency();
    void Q5MotifAdjacency();
    void Q6MotifAdjacency();
    virtual void IncrementAllWeight(int const v,int const w,int const n1,int const n2);

  public:
    BensonGraphQuadMP(PNGraph& one_graph, MotifTypeQuad& one_motif,const TInt& nThreads, bool verb = false);

    BensonGraphQuadMP();
    virtual ~BensonGraphQuadMP();
};


#endif  //bensongraphquad_h
