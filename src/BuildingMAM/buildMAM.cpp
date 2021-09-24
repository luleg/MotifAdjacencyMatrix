#include "stdafx.h"
#include "BensonGraphQuad.h"
#include "BensonGraphTriangle.h"
#include "BensonGraphEdge.h"
#include <cstring>
// #include <stdio>
#include <stdio.h>

TStr FileGraph = "";
TStr Motif = "E2";
TStr TextFileOutput = "";
TInt nbProc = -1;
TBool help = false;
TBool verbose = false;


TStr FileNGraph = "";
TStr BinaryFileOutput = "";

int usage(bool verbose){
  printf("PURPOSE ::\n\tBuilding the Motif Adjacency Matrix (MAM) of a graph, given a motif (graphlet).\n");
  printf("==========================================================================================\n");
  printf("USAGE   ::\n\tbuildMAM -i InputFile [-m Motif -nthreads NbTh -otxt OutputFile");
  if (verbose) printf(" -igraph InputSNAP -obin OutputBin ");
  printf(" -h -hv -v]\n");

  printf("\t\t-i InputFile :: a text file representing the graph as an edgelist, with integer nodes.\n");
  if (verbose) printf("\t\t\tNodes should be labelled with successive integers, starting with 0.\n");

  printf("\t\t-m Motif :: Motif to build the MAM. Listed in Motifs.pdf. Default is E2 (A+A^T).\n");
  if (verbose) {
    printf("\t\t\tMotif format \"<char><int>\": char for its size (E=2, T=3, Q=4), int for its identifier, e.g.:\n");
    printf("\t\t\t\tE1 returns the directed graph. E2 symmetrises the graph as A+A^T.\n");
    printf("\t\t\t\tT38 returns the MAM of the Feed Forward Loop. T12 gives the MAM of the 2-path (i->j->k).\n");
    printf("\t\t\t\tQ204 returns the MAM of the Bifan.\n");
  }
  printf("\t\t-nthreads NbTh :: Number of threads used to build the MAM. Default is sequential.\n");
  if (verbose) printf("\t\t\tWith only this program running on a laptop, best choice is the CPU number of threads (often 8).\n");

  printf("\t\t-otxt OutputFile :: the file in which the MAM is written as a weighted edgelist.\n");
  if (verbose) {
    printf("\t\t\t\tEach line corresponds to an edge : \"source target weight\",\n");
    printf("\t\t\t\teach edge appears twice (\"source target weight\" and \"target source weight\").\n");
  }
  if (verbose){
    printf("\t\t-igraph InputSNAP :: If the input graph is in the SNAP format (extension .ngraph).\n");
    printf("\t\t\tFor debugging purpose.\n");
    printf("\t\t-obin OutputBin :: Returns the MAM in two binary files (\"OutputBin\".bin and \"OutputBin\".weight).\n");
    printf("\t\t\tThese files are to be read by downstream partitioning methods.\n");
  }

  printf("\t\t-h  :: Plots the usage message.\n");
  printf("\t\t-hv :: Plots the extended usage message.\n");
  printf("\t\t-v  :: Runs the method in its verbose version.\n");
  printf("==========================================================================================\n");
  printf("OTHER   ::buildMAM -f ParamsFile\n");
  if (verbose){
    printf("\t\t-f ParamsFile :: A file that contains the above mentioned arguments, one by line, ended by a \"?\" char.\n");
    printf("\t\t\tSee file Data/inputs.dat for an example\.n");
  }
  return 0;
}

int readFromFile(char *file){
  FILE * Finput = fopen(file,"r");
  if (Finput == NULL){
    printf("ERROR :: FILE %S NOT FOUND\n",file);
    help = true;
    return 1;
  }
  long lSize;
  char * buffer;
  size_t result;
  fseek (Finput , 0 , SEEK_END);
  lSize = ftell (Finput);
  rewind (Finput);
  buffer = (char*) malloc (sizeof(char)*lSize);
  result = fread (buffer,1,lSize,Finput);
  if (result != lSize) {
    printf ("ERROR :: PROBLEM WHILE READING FILE %S\n",file);
    help = true;
    return 1;
  }
  fclose(Finput);

  char * oneArg;
  oneArg = strtok(buffer,"- ?");
  while (oneArg != NULL)
  {
    if (strcmp(oneArg,"i")==0){
      oneArg = strtok (NULL, "- ?");
      FileGraph = TStr(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg,"m")==0){
      oneArg = strtok (NULL, "- ?");
      Motif = TStr(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg,"nthreads")==0){
      oneArg = strtok (NULL, "- ?");
      nbProc = atoi(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg,"otxt")==0){
      oneArg = strtok (NULL, "- ?");
      TextFileOutput = TStr(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg,"igraph")==0){
      oneArg = strtok (NULL, "- ?");
      FileNGraph = TStr(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg,"obin")==0){
      oneArg = strtok (NULL, "- ?");
      BinaryFileOutput = TStr(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg,"h")==0){
      help =true;
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg,"v")==0){
      verbose = true;
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg,"hv")==0){
      help =true;
      verbose = true;
      oneArg = strtok (NULL, "- ?");
    }
    oneArg = strtok (NULL, "- ?");
  }
  free (buffer);

  return 0;

}


int read_args(int argc, char **argv){
  int num_arg = 1;
  while( num_arg < argc){
    if (strcmp(argv[num_arg], "-i") == 0){ // input file : an edgelist
      num_arg ++;
      FileGraph = argv[num_arg];
    }
    else if (strcmp(argv[num_arg], "-m") == 0){ // motif for building the MAM (Default: the graph is symmetrised like A+A^T).
      num_arg ++;
      Motif = argv[num_arg];
    }
    else if (strcmp(argv[num_arg], "-otxt") == 0){ // output file where the MAM edgelist will be stored
      num_arg ++;
      TextFileOutput = argv[num_arg];
    }
    else if (strcmp(argv[num_arg], "-nthreads") == 0){ // number of trheads to build the MAM (default: sequential)
      num_arg++;
      nbProc = atoi(argv[num_arg]);

    }
    else if (strcmp(argv[num_arg], "-h") == 0){ // if raised, displays usage info and exit
      help = true;
      return 0;
    }
    else if (strcmp(argv[num_arg], "-v") == 0){
      verbose = true;
    }
    else if (strcmp(argv[num_arg], "-hv") == 0){
      help = true;
      verbose = true;
      return 0;
    }
    else if (strcmp(argv[num_arg], "-igraph") == 0){ // if the graph is in a binary format (snap)
      num_arg++;
      FileNGraph = argv[num_arg];
    }
    else if (strcmp(argv[num_arg], "-obin") == 0){ // output files to store the MAM (binary format)
      num_arg++;
      BinaryFileOutput = argv[num_arg];
    }
    else if (strcmp(argv[num_arg], "-f") == 0){ // output files to store the MAM (binary format)
      num_arg++;
      readFromFile(argv[num_arg]);
      return 0;
    }
    else{
      printf("Unknown Parameter : %s\n",argv[num_arg]);
    }
    num_arg++;
  }
  return 0;
}




int main(int argc, char **argv){

  read_args(argc,argv); // Read the arguments -- input files, motif, output files, number of threads, flags.
  if (help){ // If flag -h, plot info about usage and returns
    usage(verbose);
    return 0;
  }

////////////////////////////////////////////////////////////////////////////////
// READ GRAPH:: FROM TEXT TO SNAP FORMAT
////////////////////////////////////////////////////////////////////////////////
  // If No input file is provided, then return.
  if (FileGraph.Empty() && FileNGraph.Empty()){
    printf("No input file. Nothing done.\n");
    return 0;
  }

  TSecTm TStart = TSecTm::GetCurTm();

  // Read the grah, and write the time taken for this operation
  TExeTm ExeTmReadGr;
  PNGraph Graph ;
  if (FileNGraph.Empty()){ // If the graph is in a textual edgelist file
    Graph = TSnap::LoadEdgeList<PNGraph>(FileGraph);
  }
  else{ // If the graph is in a snap format.
    TFIn FIn(FileNGraph);
    Graph = TNGraph::Load(FIn);
  }
  TStart = TSecTm::GetCurTm()-TStart;
  printf("%%%%%%%% CPU Time to load the Graph : \t\t%s\n%%%%%%%% Ellapsed Time : \t\t\t%s.\n%%%%%%%%\n",
      ExeTmReadGr.GetTmStr(),TStart.GetTmStr().CStr());


////////////////////////////////////////////////////////////////////////////////
// BUILDING THE MAM
////////////////////////////////////////////////////////////////////////////////
  TExeTm ExeTmBenson;
  BensonGraph<TNGraph> *Benson(0);
  // Decide which operation should be performed, provided the specified graphlet/motif
  char type_mot = Motif.GetLc()[0]; // Kind of motif is given by the first char (the "Letter")


  ////////// QUADRANGLE
  if (type_mot == 'q'){// The Motif is a quadrangle::
    MotifTypeQuad QuadMot(Motif);

    /// MULTIPROC
    if (nbProc>0){ // If Several processors
      TStart = TSecTm::GetCurTm();
      Benson = new BensonGraphQuadMP(Graph,QuadMot,nbProc, verbose);
      TStart = TSecTm::GetCurTm()-TStart;
    }
    /// SEQUENTIAL
    else{ // Otherwise
      TStart = TSecTm::GetCurTm();
      Benson = new BensonGraphQuad(Graph,QuadMot,true, verbose);
      TStart = TSecTm::GetCurTm()-TStart;
    }
  }

  ////////// TRIANGLE
  else if (type_mot == 't'){// The Motif is a triangle
    MotifTypeTriangle TriMot(Motif);

    /// MULTIPROC
    if (nbProc>0){
      TStart = TSecTm::GetCurTm();
      Benson = new BensonGraphTriangleMP(Graph,TriMot,nbProc, verbose);
      TStart = TSecTm::GetCurTm()-TStart;
    }
    /// SEQUENTIAL
    else{
      TStart = TSecTm::GetCurTm();
      Benson = new BensonGraphTriangle(Graph,TriMot,true, verbose);
      TStart = TSecTm::GetCurTm()-TStart;
    }
  }

  ////////// EDGE (Directed or Symmetrised)
  else if (type_mot == 'e'){
    MotifTypeEdge MEdge(Motif);

    /// MULTIPROC
    if (nbProc>0){ // If Several processors
      TStart = TSecTm::GetCurTm();
      Benson = new BensonGraphEdgeMP(Graph,MEdge,nbProc, verbose);
      TStart = TSecTm::GetCurTm()-TStart;
    }
    /// SEQUENTIAL
    else{ // Otherwise
      TStart = TSecTm::GetCurTm();
      Benson = new BensonGraphEdge(Graph,MEdge,true, verbose);
      TStart = TSecTm::GetCurTm()-TStart;
    }
  }
  ////////// UNKNOWN GRAPHLET
  else{
    TExcept::Throw("Unhandled kind of Motif.");
  }

  //// Print the time taken to build the Benson graph
  printf("%%%%%%%% CPU Time to compute the MAM for Motif %s : \t\t%s\n%%%%%%%% Ellapsed Time : \t\t\t%s.\n%%%%%%%%\n",
      Motif.CStr(),ExeTmBenson.GetTmStr(),TStart.GetTmStr().CStr());



////////////////////////////////////////////////////////////////////////////////
// WRITING OUTPUTS
////////////////////////////////////////////////////////////////////////////////

  if (not BinaryFileOutput.Empty()){  // If a Path to a binary file has been provided,...
    TExeTm ExeTmWrite;
    TStart = TSecTm::GetCurTm();
    Benson->SaveAsLv(BinaryFileOutput); //...save the binaries...
    TStart = TSecTm::GetCurTm()-TStart;
    printf("%%%%%%%% CPU Time to write the Benson graph in a Louvain compatible file : \t\t%s\n%%%%%%%% Ellapsed Time : \t\t\t%s.\n%%%%%%%%\n",
        ExeTmWrite.GetTmStr(),TStart.GetTmStr().CStr()); //...and displays the time takn by this operation.
  }
  if (not TextFileOutput.Empty()){ // If a text file has been provided,...
    TExeTm ExeTmWrite;
    TStart = TSecTm::GetCurTm();
    Benson->SaveAsTxt(TextFileOutput); //...save the textual edgelist...
    TStart = TSecTm::GetCurTm()-TStart;
    printf("%%%%%%%% CPU Time to write the Benson graph in a text file : \t\t%s\n%%%%%%%% Ellapsed Time : \t\t\t%s.\n%%%%%%%%\n",
        ExeTmWrite.GetTmStr(),TStart.GetTmStr().CStr());//...and displays the time takn by this operation.
  }

  // CLEAN POINTERS

  delete Benson;
  return 0;

}
