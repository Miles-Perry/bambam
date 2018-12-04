#include <utils/bamtools_pileup_engine.h>

#include <string>
#include <vector>


class MpileupVisitor : public BamTools::PileupVisitor {
 private:
  int start;
  int end;
  
 public:
  std::vector<std::string> results;
  
  MpileupVisitor (int _start, int _end) : BamTools::PileupVisitor() {
    init(_start, _end);
  }

  ~MpileupVisitor() {
    clear();
  }

  void init (int _start, int _end) {
    start = _start;
    end = _end;
    results.resize(end-start);
    clear();
  }

  void clear() {
    for (int i = 0; i < end-start; i++) results[i] = "";
  }

  void Visit (const BamTools::PileupPosition& pileup) {
    if (pileup.Position < start || pileup.Position >= end) return;

    for (std::vector<BamTools::PileupAlignment>::const_iterator itr = pileup.PileupAlignments.begin(); itr != pileup.PileupAlignments.end(); itr++) {
      BamTools::PileupAlignment pup = *itr;
      if (pup.IsCurrentDeletion) continue;

      BamTools::BamAlignment& aln = pup.Alignment;
      char base = aln.QueryBases[pup.PositionInAlignment];
      if (aln.IsReverseStrand()) base = tolower(base);

      results[pileup.Position-start] += base;
    }
  }

  std::string& get (int i) {
    return results[i-start];
  }
};
