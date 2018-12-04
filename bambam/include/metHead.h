// methHead
// Justin Page
// jtpage68@gmail.com

#include <api/BamAlignment.h>

#include "BamUtils.h"
#include "MultiIndex.h"


void determineStrand (BamTools::BamAlignment& alignment, int* strandVotes, MultiIndex* cyts, MultiIndex* snps = NULL) {
  for (Snp* cyt = cyts->get (alignment.RefID, alignment.Position); cyt && cyt->pos <= alignment.GetEndPosition(); cyt = cyt->next) {
    char base = offsetBase (alignment, cyt->pos);
    if (!base) continue;

    // don't count strand vote at SNP locus
    if (snps) {
      Snp* nextSnp = snps->get (alignment.RefID, cyt->pos);
      if (nextSnp && nextSnp->pos == cyt->pos) continue;
    }

    if (cyt->vars[0][0] == 'C' && base == 'T') strandVotes[0]++;
    else if (cyt->vars[0][0] == 'G' && base == 'A') strandVotes[1]++;
  }
}
