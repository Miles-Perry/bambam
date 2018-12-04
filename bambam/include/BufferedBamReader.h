#include "BamUtils.h"

#include <api/BamReader.h>
#include <api/BamAlignment.h>


class BufferedBamReader {
 private:
  BamTools::BamReader* reader;
  BamTools::BamAlignment nextAln;


 public:
  bool hasNext;

  BufferedBamReader(std::string filename) {
    reader = openBam (filename, true);
    hasNext = true;
    next();
  }

  ~BufferedBamReader() {
    reader->Close();
    delete reader;
  }


  BamTools::BamReader* getReader() {
    return reader;
  }


  void next() {
    hasNext = (reader->GetNextAlignment (nextAln));
  }

  BamTools::BamAlignment& get() {
    return nextAln;
  }
};
