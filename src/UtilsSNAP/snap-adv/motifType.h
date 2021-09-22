//// CONTAINS THE SPEC OF THE MOTHER CLASS OF MOTIFS,
// **************MotifType
// AND DAUGHTER CLASSES  FOR THE DIFFERENT KINDS OF GRAPHLETS, NAMELY:
// *****MotifTypeTriangle                  FOR 3-NODE,
// *****MotifTypeQuad                     FOR QUADRANGLES,
// *****MotifTypeEdge                     FOR SYMMETRISED OR DIRECTED EDGES.


#ifndef snap_motiftype_h
#define snap_motiftype_h
#include "my_functors.h" // Used for 3-node graphlets


// Some MACROS to characterise the motifs
#define NON 0
// MACROS for edge motifs
#define DIR 1
#define SYM 2
// MACROS for 3-node motifs
#define WEDGE 2
#define TRIANGLE 3
// MACROS for quadrangles
#define Q4 4
#define Q5 5
#define Q6 6


// Mother class of motifs
// type to specify the type:
//      for 3-node triangle or wedge;
//      for quadrangles Q4, Q5 or Q6;
//      for edges symmetrised or directed
// id to specify specifically the graphlets e.g.:
//      204 for bifan
//      12 for 3-node path

class MotifType {
  protected:
    int id;
    int type;
  public:
    MotifType() : id(0), type(NON) {
    };
    virtual ~MotifType() {};
    int getType() const {return type;};
};

/******************************
      QUADRANGLE
*******************************/
class MotifTypeQuad : public MotifType{
private:
  TIntV cpt2id; //to know if an induced subgraph is the required motif

public:
  MotifTypeQuad();
  MotifTypeQuad(const TStr& motif);

  virtual ~MotifTypeQuad();
  bool IsInstance(const int& cpt_id) const; // return true is cpt_id is a valid id of this quadrangle
};

/******************************
      TRIANGLE
*******************************/
class MotifTypeTriangle : public MotifType{

private:
  Func_IsInstance* FInstance; //to know if a induced subgraph is the required motif

public:
  MotifTypeTriangle(const TStr& motif);
  MotifTypeTriangle();
  virtual ~MotifTypeTriangle();
  bool IsInstance(const int& u, const int& v, const int& w) const; // return true if the subgraph induced by nodes u,v,w is isomorphic to the motif
  void SetGraph(const PNGraph& one_graph); // Links FInstance to a graph (required to explore induced subgraphs)
};



/******************************
      EDGES
*******************************/
class MotifTypeEdge : public MotifType{

public:
  MotifTypeEdge(const TStr& motif);
  MotifTypeEdge();
  virtual ~MotifTypeEdge();
};


#endif  // snap_motiftype_h
