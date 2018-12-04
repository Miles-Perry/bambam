// HapHunt
// Justin Page
// jtpage68@gmail.com

#include "BamUtils.h"

#include <vector>
#include <string>


class Snp {
 public:
  int pos;
  std::string vars;

  Snp () : pos(-1), vars("") {}
  Snp (int _pos, std::string& _vars) : pos(_pos), vars(_vars) {}
};

class Cluster {
  friend class Solution;

 private:
  std::string consensus;
  std::vector<std::string> reads;

 public:
 Cluster(std::string& _consensus) : consensus(_consensus) {}
 Cluster() : consensus("") {}

  void addRead (std::string& read) {
    reads.push_back (read);
  }

  void update() {
    std::string newCons = "";
    
    if (reads.size()) {
      for (int i = 0; i < consensus.length(); i++) {

        int vars[4] = {0,0,0,0};

        // total up base counts
        for (std::vector<std::string>::iterator itr = reads.begin(); itr != reads.end(); itr++) {
          char base = itr->at(i);

          for (int j = 0; j < 4; j++) {
            if (base == dNTPstr[j]) {
              vars[j]++;
              break;
            }
          }
        }

	// find most represented base
        int max = -1;
        for (int j = 0; j < 4; j++) {
          if ((vars[j] > 0 && (max < 0 || vars[j] > vars[max])) ||
	      (vars[j] == vars[max] && rand()%2)) { // randomly assign if unclear
            max = j;
          }
        }

	// stick with current base unless clear majority
	char cons = 'N';
	if (max >= 0) cons = dNTPstr[max];
	newCons += cons;
      }

      // swap in updated consensus
      consensus = newCons;
    }
  }

  float distance () {
    float dist = 0.0;

    for (std::vector<std::string>::iterator itr = reads.begin(); itr != reads.end(); itr++) {
      dist += seqDistance (consensus, *itr, 0);
    }

    if (reads.size()) dist /= reads.size();
    return dist;
  }

  float distance (Cluster& other) {
    return seqDistance (consensus, other.consensus, 0);
  }
};
