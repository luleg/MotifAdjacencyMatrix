#include "motifType.h"

/******************************
      QUADRANGLE
*******************************/

// Utility function that allows to link id of a 4-node subgraph to valid quadrangle ids.
// Required to know if an induced subgraph is an instance of a quadrangle
void loadIds(TStr Qfile,TIntV& cpt2id){
  // Given an 4-node induced subgraph of Adjacency amtrix A, and idA equal to
  // e'*(A.*B)*e) with B = reshape(2.^[0:15],[4,4])'; we have
  // cpt2id[idA] = idQ with idQ the minimal id of the quadrangle represented by A
  // #cpt2id is equal to 2^15-1=32767

  // Binary have been built that contains the correspondances.
  FILE* FQ = fopen(Qfile.CStr(),"rb");
  if (FQ == NULL ) TExcept::Throw("Unable to open the binary file containing QUADRANGLEs.");
  int idMax,nbIds;
  fread(&idMax,sizeof(idMax),1,FQ);
  cpt2id.Reserve(idMax,idMax);
  fread(&nbIds,sizeof(nbIds),1,FQ);
  for (int cpt_id = 0; cpt_id<nbIds;cpt_id++){
    int id,nb_others;
    fread(&id,sizeof(id),1,FQ);
    fread(&nb_others,sizeof(nb_others),1,FQ);
    for (int one = 0;one<nb_others;one++){
      int another;
      fread(&another,sizeof(another),1,FQ);
      cpt2id[another] = id;
    }
  }
}

// Dummy constructor
MotifTypeQuad::MotifTypeQuad() : MotifType(){
  cpt2id =TIntV(0);
}

//Constructor:
// given the id of the quadrangle, specify "type", "id", and cpt2id
// Since induced subgraphs to be tested as isomorphic to the quadrangle have the
// same type (ie non directed structure (from the algorithm)), cpt2id contains
//the ids of quandrangles of this type.
// There are 129 directed quadrangles, thus this unction contains 129 if.
MotifTypeQuad::MotifTypeQuad(const TStr& motif) {

  TStr motif_lc = motif.GetLc();
  if (motif_lc == "q204"){
      id=204;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q456"){
      id=456;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q460"){
      id=460;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q904"){
      id=904;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q908"){
      id=908;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q972"){
      id=972;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q4548"){
      id=4548;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q4556"){
      id=4556;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q4740"){
      id=4740;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q4748"){
      id=4748;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q4812"){
      id=4812;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q5004"){
      id=5004;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q5064"){
      id=5064;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q5068"){
      id=5068;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q13260"){
      id=13260;
      loadIds("src/UtilsSNAP/snap-adv/Q4.bin",cpt2id);
      type = 4;
  }
  else if (motif_lc == "q206"){
      id=206;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q222"){
      id=222;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q458"){
      id=458;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q462"){
      id=462;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q472"){
      id=472;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q474"){
      id=474;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q476"){
      id=476;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q478"){
      id=478;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q906"){
      id=906;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q910"){
      id=910;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q922"){
      id=922;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q924"){
      id=924;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q926"){
      id=926;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q974"){
      id=974;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q990"){
      id=990;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q2190"){
      id=2190;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q2204"){
      id=2204;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q2206"){
      id=2206;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q2252"){
      id=2252;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q2458"){
      id=2458;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q2462"){
      id=2462;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4546"){
      id=4546;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4550"){
      id=4550;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4558"){
      id=4558;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4562"){
      id=4562;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4564"){
      id=4564;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4566"){
      id=4566;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4572"){
      id=4572;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4574"){
      id=4574;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4742"){
      id=4742;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4750"){
      id=4750;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4758"){
      id=4758;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4764"){
      id=4764;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4766"){
      id=4766;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4814"){
      id=4814;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4830"){
      id=4830;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4994"){
      id=4994;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q4998"){
      id=4998;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5002"){
      id=5002;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5006"){
      id=5006;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5010"){
      id=5010;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5012"){
      id=5012;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5014"){
      id=5014;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5016"){
      id=5016;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5018"){
      id=5018;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5020"){
      id=5020;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5022"){
      id=5022;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5058"){
      id=5058;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5062"){
      id=5062;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5066"){
      id=5066;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5070"){
      id=5070;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5074"){
      id=5074;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5076"){
      id=5076;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5078"){
      id=5078;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5080"){
      id=5080;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5082"){
      id=5082;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5084"){
      id=5084;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q5086"){
      id=5086;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q6348"){
      id=6348;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q6550"){
      id=6550;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q6552"){
      id=6552;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q6554"){
      id=6554;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q6558"){
      id=6558;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q6604"){
      id=6604;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q6858"){
      id=6858;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q6874"){
      id=6874;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q13142"){
      id=13142;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q13146"){
      id=13146;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q13148"){
      id=13148;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q13150"){
      id=13150;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q13262"){
      id=13262;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q13278"){
      id=13278;
      loadIds("src/UtilsSNAP/snap-adv/Q5.bin",cpt2id);
      type = 5;
  }
  else if (motif_lc == "q2254"){
      id=2254;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q2270"){
      id=2270;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q2506"){
      id=2506;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q2510"){
      id=2510;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q2524"){
      id=2524;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q2526"){
      id=2526;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q3038"){
      id=3038;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6342"){
      id=6342;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6350"){
      id=6350;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6356"){
      id=6356;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6358"){
      id=6358;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6364"){
      id=6364;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6366"){
      id=6366;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6598"){
      id=6598;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6602"){
      id=6602;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6606"){
      id=6606;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6614"){
      id=6614;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6616"){
      id=6616;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6618"){
      id=6618;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6620"){
      id=6620;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6622"){
      id=6622;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6854"){
      id=6854;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6862"){
      id=6862;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6870"){
      id=6870;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6876"){
      id=6876;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q6878"){
      id=6878;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q7126"){
      id=7126;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q7128"){
      id=7128;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q7130"){
      id=7130;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q7134"){
      id=7134;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q14678"){
      id=14678;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q14686"){
      id=14686;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q14790"){
      id=14790;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q14798"){
      id=14798;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q14810"){
      id=14810;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q14812"){
      id=14812;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q14814"){
      id=14814;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q15258"){
      id=15258;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q15262"){
      id=15262;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q15310"){
      id=15310;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q15326"){
      id=15326;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else if (motif_lc == "q31710"){
      id=31710;
      loadIds("src/UtilsSNAP/snap-adv/Q6.bin",cpt2id);
      type = 6;
  }
  else {

     TExcept::Throw("Unknown motif");
   }

}


bool MotifTypeQuad::IsInstance(const int& cpt_id) const{
  return (cpt2id[cpt_id] == id);
}
//


MotifTypeQuad::~MotifTypeQuad(){
}

/******************************
      TRIANGLE
*******************************/

MotifTypeTriangle::MotifTypeTriangle(): MotifType(){
  FInstance = NULL;
}

void MotifTypeTriangle::SetGraph(const PNGraph& one_graph){
  FInstance->SetGraph(one_graph); // Links FInstance to a graph (required to explore induced subgraphs)
}


MotifTypeTriangle::~MotifTypeTriangle(){
  delete FInstance; // This is a pointer, needs to be delete
}

// Constructor. Understands motif as in Benson's Thesis (M1,M2...) and also
// minimal graph id (T98, T12, etc)
MotifTypeTriangle::MotifTypeTriangle(const TStr& motif) {
  TStr motif_lc = motif.GetLc();
  if ((motif_lc == "m1")||(motif_lc == "t98"))          {
    id=1;
    type = TRIANGLE;
    FInstance = new Func_IsInstance_M1();
  }
  else if ((motif_lc == "m2")||(motif_lc == "t102"))         {
    id=2;
    type = TRIANGLE;
    FInstance = new Func_IsInstance_M2();
  }
  else if ((motif_lc == "m3")||(motif_lc == "t110"))          {
    id=3;
    type = TRIANGLE;
    FInstance = new Func_IsInstance_M3();
  }
  else if ((motif_lc == "m4")||(motif_lc == "t238"))          {
    id=4;
    type = TRIANGLE;
    FInstance = new Func_IsInstance_M4();
  }
  else if ((motif_lc == "m5")||(motif_lc == "t38"))          {
    id=5;
    type = TRIANGLE;
    FInstance = new Func_IsInstance_M5();
  }
  else if ((motif_lc == "m6")||(motif_lc == "t108"))          {
    id=6;
    type = TRIANGLE;
    FInstance = new Func_IsInstance_M6();
  }
  else if ((motif_lc == "m7")||(motif_lc == "t46"))          {
    id=7;
    type = TRIANGLE;
    FInstance = new Func_IsInstance_M7();
  }
  else if ((motif_lc == "m8")||(motif_lc == "t6"))          {
    id=8;
    type = WEDGE;
    FInstance = new Func_IsInstance_M8();
  }
  else if ((motif_lc == "m9")||(motif_lc == "t12"))          {
    id=9;
    type = WEDGE;
    FInstance = new Func_IsInstance_M9();
  }
  else if ((motif_lc == "m10")||(motif_lc == "t36"))         {
    id=10;
    type = WEDGE;
    FInstance = new Func_IsInstance_M10();
  }
  else if ((motif_lc == "m11")||(motif_lc == "t14"))         {
    id=11;
    type = WEDGE;
    FInstance = new Func_IsInstance_M11();
  }
  else if ((motif_lc == "m12")||(motif_lc == "t74"))         {
    id=12;
    type = WEDGE;
    FInstance = new Func_IsInstance_M12();
  }
  else if ((motif_lc == "m13")||(motif_lc == "t78"))         {
    id=13;
    type = WEDGE;
    FInstance = new Func_IsInstance_M13();
  }
  else { TExcept::Throw("Unknown motif"); }
}

bool MotifTypeTriangle::IsInstance(const int& u, const int& v, const int& w) const{
  return (*FInstance)(u,v,w);
  // Each kind of 3-node graphlet has its own FInstance, which is an functor
  // containing a graph G and a function F. F return true if the induced subrgaph
  // formed by (u,v,w) in G is an instance of the motif
}

/******************************
      EDGES
*******************************/

MotifTypeEdge::MotifTypeEdge() : MotifType(){
}

MotifTypeEdge::MotifTypeEdge(const TStr& motif){
  if (motif[1] == '1'){
    type = DIR;
  }
  else if(motif[1] == '2'){
    type=SYM;
  }
}

MotifTypeEdge::~MotifTypeEdge(){
}
