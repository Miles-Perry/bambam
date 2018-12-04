#include <utils/bamtools_pileup_engine.h>

#include <string>
#include <sstream>
#include <vector>
#include <list>


enum dNTP { A, C, G, T };
std::string dNTPstr = "ACGT";

class Insert {
 public:
  std::string seq;
  int count;

  Insert (std::string _seq) : seq(_seq), count(1) {}
};


class ConsensusVisitor : public BamTools::PileupVisitor {
  float maf;
  float mif;
  int min;
  int start;
  int end;
  std::vector<std::string> results;
  std::vector<std::string> insertions;
  std::vector<int> coverage;
  
public:
 ConsensusVisitor() : BamTools::PileupVisitor() {
    start = 0;
    end = 1;
    maf = 1.0;
    min = 0;
    mif = 0.0;
    init(start, end);
  }

  ConsensusVisitor (int _start, int _end, float _maf, int _min, float _mif = 0.0) : BamTools::PileupVisitor(), maf(_maf), min(_min), mif(_mif) {
    init(_start, _end);
  }
    
  ~ConsensusVisitor() {
    clear();
  }

  // setup or reset
  void init(int _start, int _end) {
    start = _start;
    end = _end;
    results.resize(end-start);
    if (mif) insertions.resize(end-start);
    coverage.resize(end-start);

    clear();
  }

  // remove results
  void clear() {
    for (int i = 0; i < end-start; i++) {
      results[i] = "";
      if (mif) insertions[i] = "";
      coverage[i] = 0;
    }
  }
  


  void Visit (const BamTools::PileupPosition& pileup) {
    if (pileup.Position < start || pileup.Position >= end) return;

    int vars[] = {0,0,0,0};
    coverage[pileup.Position-start] = pileup.PileupAlignments.size();;
    int dels = 0;
    std::list<Insert> inserts;
    std::string varStr = "";


    // total base frequencies
    for (std::vector<BamTools::PileupAlignment>::const_iterator itr = pileup.PileupAlignments.begin(); itr != pileup.PileupAlignments.end(); itr++) {
      BamTools::PileupAlignment pup = *itr;

      // check for deletions
      if (pup.IsCurrentDeletion) {
	dels++;
	continue;
      }

      // check for insertions
      if (mif && pup.InsertionLength) {

	// build insertion string
	std::string seq = "";
	for (int i = 0; i < pup.InsertionLength; i++) seq += toupper(pup.Alignment.QueryBases[pup.PositionInAlignment+i+1]);
	bool found = false;

	// increment vote if already present
	for (std::list<Insert>::iterator itr = inserts.begin(); itr != inserts.end(); itr++) {
	  if (seq == itr->seq) {
	    itr->count++;
	    found = true;
	    break;
	  }
	}

	if (!found) {
	  Insert ins(seq);
	  inserts.push_back (ins);
	}
      }


      // total bases at point
      char base = toupper(pup.Alignment.QueryBases[pup.PositionInAlignment]);

      for (int i = 0; i < 4; i++) {
	if (base == dNTPstr[i]) {
	  vars[i]++;
	  break;
	}
      }
    }


    // find most common alleles
    float mafThreshold = maf * coverage[pileup.Position-start];
    if (mafThreshold < min) mafThreshold = min;

    for (int i = 0; i < 4; i++) {
      if (vars[i] >= mafThreshold) {
	varStr += dNTPstr[i];
      }
    }

    
    // show indels
    if (mif) {
      float mifThreshold = mif * (dels+coverage[pileup.Position-start]);
      if (mifThreshold < min) mifThreshold = min;
      
      // deletions
      if (dels >= mifThreshold) {
	varStr = "-";
      }
    
      else { // insertion
	for (std::list<Insert>::iterator itr = inserts.begin(); itr != inserts.end(); itr++) {
	  if (itr->count >= mifThreshold) {
	    insertions[pileup.Position-start] = itr->seq;
	    break;
	  }
	}
      }
    }

    results[pileup.Position-start] = varStr; // record results (note offset)
  }



  // access
  int getStart() {
    return start;
  }

  std::string& get (int i) {
    if (!results[i-start].length()) results[i-start] = "N";
    return results[i-start];
  }

  std::string& getInsert (int i) {
    return insertions[i-start];
  }

  int getCoverage (int i) {
    return coverage[i-start];
  }

  void set (int i, const char* val) {
    results[i-start] = val;
    coverage[i-start] = 0;
  }
};
