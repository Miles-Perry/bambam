#include <utils/bamtools_pileup_engine.h>

#include <string>
#include <vector>

class CoverageVisitor : public BamTools::PileupVisitor {
 private:
  int start;
  int end;
  
 public:
  std::vector<int> results;

 CoverageVisitor() : BamTools::PileupVisitor() {
    start = 0;
    end = 1;
    init(start, end);
  }
  
  CoverageVisitor (int _start, int _end) : BamTools::PileupVisitor() {
    init(_start, _end);
  }

  ~CoverageVisitor() {
    clear();
  }

  // setup or reset
  void init (int _start, int _end) {
    start = _start;
    end = _end;
    results.resize(end-start);
    clear();
  }

  // remove results
  void clear() {
    for (int i = 0; i < results.size(); i++) results[i] = 0;
  }

  void Visit (const BamTools::PileupPosition& pileup) {
    if (pileup.Position < start || pileup.Position >= end) return;

    results[pileup.Position-start] = 0;
    for (std::vector<BamTools::PileupAlignment>::const_iterator itr = pileup.PileupAlignments.begin(); itr != pileup.PileupAlignments.end(); itr++) {
      BamTools::PileupAlignment pup = *itr;
      if (!pup.IsCurrentDeletion) results[pileup.Position-start]++;
    }
  }

  int get (int i) {
    return results[i-start];
  }
};
