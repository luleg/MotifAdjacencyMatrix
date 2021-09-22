#ifndef snap_fctor_triangle_h
#define snap_fctor_triangle_h
#include "Snap.h"

// Each kind of 3-node graphlet has its own FInstance, which is an functor
// containing a graph G and a function F. F return true if the induced subrgaph
// formed by (u,v,w) in G is an instance of the motif.

class Func_IsInstance {
  protected:
    PNGraph graph;
    bool IsNoEdge(int u, int v){
      return !graph->IsEdge(u, v) && !graph->IsEdge(v, u);
    };
    bool IsUnidirEdge(int u, int v){
      return graph->IsEdge(u, v) && !graph->IsEdge(v, u);
    };
    bool IsBidirEdge(int u, int v){
      return graph->IsEdge(u, v) && graph->IsEdge(v, u);
    };
  public:
    virtual ~Func_IsInstance() {
      // printf("~Func_IsInstance destructor\n");
    };
    Func_IsInstance() : graph(NULL){};
    void SetGraph(const PNGraph& one_graph){
      graph = one_graph;
    };
    virtual bool operator()(int u, int v, int w) = 0;
};

class Func_IsInstance_M1 : public Func_IsInstance {
  public:

    virtual ~Func_IsInstance_M1(){
      // printf("~Func_IsInstance_M1 destructor\n");
    };
    Func_IsInstance_M1() : Func_IsInstance() {};

    virtual bool operator()(int u, int v, int w) {
      return ((IsUnidirEdge(u, v) && IsUnidirEdge(v, w) &&
    	   IsUnidirEdge(w, u)) ||
              (IsUnidirEdge(u, w) && IsUnidirEdge(w, v) &&
    	   IsUnidirEdge(v, u)));
    };
};

class Func_IsInstance_M2 : public Func_IsInstance {
  public:

  virtual ~Func_IsInstance_M2(){
    // printf("~Func_IsInstance_M2 destructor\n");
  };
  Func_IsInstance_M2() : Func_IsInstance() {};

  virtual bool operator()(int u, int v, int w) {
    return ((IsBidirEdge(u, v) && IsUnidirEdge( u, w) &&
       IsUnidirEdge(w, v)) ||
            (IsBidirEdge( u, v) && IsUnidirEdge( w, u) &&
       IsUnidirEdge( v, w)) ||
            (IsBidirEdge( u, w) && IsUnidirEdge( u, v) &&
       IsUnidirEdge( v, w)) ||
            (IsBidirEdge( u, w) && IsUnidirEdge( v, u) &&
       IsUnidirEdge( w, v)) ||
            (IsBidirEdge( v, w) && IsUnidirEdge( v, u) &&
       IsUnidirEdge( u, w)) ||
            (IsBidirEdge( v, w) && IsUnidirEdge( u, v) &&
       IsUnidirEdge( w, u)));
  };
};

class Func_IsInstance_M3 : public Func_IsInstance {
public:

  virtual ~Func_IsInstance_M3(){
    // printf("~Func_IsInstance_M3 destructor\n");
  };
  Func_IsInstance_M3() : Func_IsInstance() {};

  virtual bool operator()(int u, int v, int w) {
    if ((IsBidirEdge( u, v) && IsBidirEdge(  v, w)) &&
        (IsUnidirEdge( u, w) || IsUnidirEdge( w, u))) { return true; }
    if ((IsBidirEdge(  u, w) && IsBidirEdge(  w, v)) &&
        (IsUnidirEdge( u, v) || IsUnidirEdge( v, u))) { return true; }
    if ((IsBidirEdge( w, u) && IsBidirEdge( u, v)) &&
        (IsUnidirEdge( w, v) || IsUnidirEdge(v, w))) { return true; }
    return false;
  };
};
class Func_IsInstance_M4 : public Func_IsInstance {
public:
  virtual ~Func_IsInstance_M4(){
    // printf("~Func_IsInstance_M4 destructor\n");
  };
  Func_IsInstance_M4() : Func_IsInstance() {};

  virtual bool operator()(int u, int v, int w) {
    return IsBidirEdge(u, v) && IsBidirEdge( u, w) && IsBidirEdge(v, w);
  };
};
class Func_IsInstance_M5 : public Func_IsInstance {
public:

  virtual ~Func_IsInstance_M5(){
    // printf("~Func_IsInstance_M5 destructor\n");
  };
  Func_IsInstance_M5() : Func_IsInstance() {};

  virtual bool operator()(int u, int v, int w) {
    if ((IsUnidirEdge(u, v) && IsUnidirEdge(u, w)) &&
        (IsUnidirEdge(v, w) || IsUnidirEdge(w, v))) { return true; }
    if ((IsUnidirEdge(v, u) && IsUnidirEdge(v, w)) &&
        (IsUnidirEdge(u, w) || IsUnidirEdge(w, u))) { return true; }
    if ((IsUnidirEdge(w, v) && IsUnidirEdge(w, u)) &&
        (IsUnidirEdge(v, u) || IsUnidirEdge(u, v))) { return true; }
    return false;
  };
};
class Func_IsInstance_M6 : public Func_IsInstance {
public:

  virtual ~Func_IsInstance_M6(){
    // printf("~Func_IsInstance_M6 destructor\n");
  };

  Func_IsInstance_M6() : Func_IsInstance() {};

  virtual bool operator()(int u, int v, int w) {
    return ((IsUnidirEdge(u, v) && IsUnidirEdge(u, w) &&
       IsBidirEdge(v, w)) ||
            (IsUnidirEdge(v, u) && IsUnidirEdge(v, w) &&
       IsBidirEdge(u, w)) ||
            (IsUnidirEdge(w, u) && IsUnidirEdge(w, v) &&
       IsBidirEdge(u, v)));
  };
};
class Func_IsInstance_M7 : public Func_IsInstance {
public:

  virtual ~Func_IsInstance_M7(){
    // printf("~Func_IsInstance_M7 destructor\n");
  };

  Func_IsInstance_M7() : Func_IsInstance() {};

  virtual bool operator()(int u, int v, int w) {
    return ((IsUnidirEdge(v, u) && IsUnidirEdge(w, u) &&
       IsBidirEdge(v, w)) ||
            (IsUnidirEdge(u, v) && IsUnidirEdge(w, v) &&
       IsBidirEdge(u, w)) ||
            (IsUnidirEdge(u, w) && IsUnidirEdge(v, w) &&
       IsBidirEdge(u, v)));
  };
};
class Func_IsInstance_M8 : public Func_IsInstance {
public:

  virtual ~Func_IsInstance_M8(){
    // printf("~Func_IsInstance_M8 destructor\n");
  };
  Func_IsInstance_M8() : Func_IsInstance() {};

  virtual bool operator()(int center, int v, int w) {
    return IsNoEdge(v, w) && IsUnidirEdge(center, v) &&
      IsUnidirEdge(center, w);
  };
};
class Func_IsInstance_M9 : public Func_IsInstance {
public:

  virtual ~Func_IsInstance_M9(){
    // printf("~Func_IsInstance_M9 destructor\n");
  };
  Func_IsInstance_M9() : Func_IsInstance() {};

  virtual bool operator()(int center, int v, int w) {
    return IsNoEdge(v, w) &&
      ((IsUnidirEdge(center, v) && IsUnidirEdge(w, center)) ||
       (IsUnidirEdge(center, w) && IsUnidirEdge(v, center)));
  };
};
class Func_IsInstance_M10 : public Func_IsInstance {
public:

  virtual ~Func_IsInstance_M10(){
    // printf("~Func_IsInstance_M10 destructor\n");
  };
  Func_IsInstance_M10() : Func_IsInstance() {};

  virtual bool operator()(int center, int v, int w) {
    return IsNoEdge(v, w) && IsUnidirEdge(v, center) &&
      IsUnidirEdge(w, center);
  };
};
class Func_IsInstance_M11 : public Func_IsInstance {
public:
  virtual ~Func_IsInstance_M11(){
    // printf("~Func_IsInstance_M11 destructor\n");
  };
  Func_IsInstance_M11() : Func_IsInstance() {};

  virtual bool operator()(int center, int v, int w) {
    return IsNoEdge(v, w) &&
      ((IsBidirEdge(center, v) && IsUnidirEdge(center, w)) ||
       (IsBidirEdge(center, w) && IsUnidirEdge(center, v)));
  };
};
class Func_IsInstance_M12 : public Func_IsInstance {
public:
  virtual ~Func_IsInstance_M12(){
    // printf("~Func_IsInstance_M12 destructor\n");
  };
  Func_IsInstance_M12() : Func_IsInstance() {};

  virtual bool operator()(int center, int v, int w) {
    return IsNoEdge(v, w) &&
      ((IsBidirEdge(center, v) && IsUnidirEdge(w, center)) ||
       (IsBidirEdge(center, w) && IsUnidirEdge(v, center)));
  };
};
class Func_IsInstance_M13 : public Func_IsInstance {
public:
  virtual ~Func_IsInstance_M13(){
    // printf("~Func_IsInstance_M13 destructor\n");
  };
  Func_IsInstance_M13() : Func_IsInstance() {};

  virtual bool operator()(int center, int v, int w) {
    return IsNoEdge(v, w) && IsBidirEdge(center, v)
      && IsBidirEdge(center, w);
  };
};
#endif //snap_fctor_triangle_h
