#pragma once

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <sys/time.h>
#include <string>
#include <vector>

#include <api/BamReader.h>



// return current time in seconds
double time() {
  struct timeval tp;
  gettimeofday (&tp, NULL);
  return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}


// open BAM file
BamTools::BamReader* openBam (std::string bamFile, bool exception=false) {

  // extract base of name
  std::string bamBase = bamFile;
  bamBase.erase (bamBase.length()-4);

  // open BAM and index
  BamTools::BamReader* newBam = new BamTools::BamReader();
  if (!newBam->Open (bamFile)) throw std::invalid_argument ("Failed to open BAM file " + bamFile);
  if (!exception && !newBam->OpenIndex (bamFile + ".bai")) {
    if (!newBam->OpenIndex (bamBase + ".bai")) throw std::invalid_argument ("Failed to open BAM index " + bamFile + ".bai");
  }

  return newBam;
}


// generate a cigar string from a BAM alignment
char* getCigar (BamTools::BamAlignment& aln) {
  char* result = (char*) malloc (aln.CigarData.size()*3);
  sprintf (result, "");
  
  for (int i = 0; i < aln.CigarData.size(); i++) {
    sprintf(result+strlen(result), "%d%c", aln.CigarData[i].Length, aln.CigarData[i].Type);
  }
  
  return result;
}


// make sure all BAMs point to the same reference and (unless exception) have an index
bool checkBams (std::vector<std::string>& bamFiles, bool exception=false) {
  
  // open BAMs
  std::vector<BamTools::BamReader*> bams;
  for (int i = 0; i < bamFiles.size(); i++) {
    bams.push_back (openBam(bamFiles[i], exception));
  }


  // compare BAMs
  bool result = true;

  // check number
  int count = bams[0]->GetReferenceCount();
  for (int i = 1; result && i < bams.size(); i++) {
    if (bams[i]->GetReferenceCount() != count) result = false;
  }

  // check identity
  for (int i = 0; result && i < bams[0]->GetReferenceCount() && result; i++) {
    std::string name = bams[0]->GetReferenceData()[i].RefName;
    int length = bams[0]->GetReferenceData()[i].RefLength;
    
    for (int j = 1; result && j < bams.size(); j++) {
      if (name != bams[j]->GetReferenceData()[i].RefName || length != bams[j]->GetReferenceData()[i].RefLength) result = false;
    }
  }

  // cleanup
  for (int i = 0; i < bams.size(); i++) {
    bams[i]->Close();
    delete bams[i];
  }

  return result;
}


// return base at given position in alignment
char offsetBase (BamTools::BamAlignment& alignment, int pos) {
  int refOffset = alignment.Position;
  int seqOffset = 0;
  char result = 0;
  if (refOffset > pos) return result;


  // parse CIGAR string
  std::vector<BamTools::CigarOp>::const_iterator iter = alignment.CigarData.begin();
  std::vector<BamTools::CigarOp>::const_iterator end = alignment.CigarData.end();
  char lastOp = 'N';

  while (iter != end) {
    BamTools::CigarOp op = *iter;
    int length = op.Length;

    if (op.Type == BamTools::Constants::BAM_CIGAR_MATCH_CHAR ||
        op.Type == BamTools::Constants::BAM_CIGAR_DEL_CHAR ||
        op.Type == BamTools::Constants::BAM_CIGAR_REFSKIP_CHAR) {

      if (refOffset + length > pos) length = pos - refOffset; // don't pass target
      refOffset += length;
    }

    if (op.Type == BamTools::Constants::BAM_CIGAR_MATCH_CHAR ||
        op.Type == BamTools::Constants::BAM_CIGAR_INS_CHAR ||
        op.Type == BamTools::Constants::BAM_CIGAR_SOFTCLIP_CHAR ||
        op.Type == BamTools::Constants::BAM_CIGAR_HARDCLIP_CHAR) {

      seqOffset += length;
    }

    if (refOffset == pos) {
      lastOp = op.Type;
      break;
    }

    ++iter;
  }


  // extract actual base
  std::string seq = alignment.QueryBases;
  if (lastOp == BamTools::Constants::BAM_CIGAR_MATCH_CHAR && seqOffset < seq.length()) result = toupper(seq.at (seqOffset));

  return result;
}


// find nucleotide distance between two sequences
float seqDistance (std::string seq1, std::string seq2, float punish) {
  assert (seq1.length() == seq2.length());

  float diff = 0.0;
  float count = 0;

  for (int i = 0; i < seq1.length(); i++) {

    if (seq1[i] != 'N' && seq2[i] != 'N') {
      if (seq1[i] != seq2[i]) diff++;
      count++;
    }
    else diff += punish;
  }

  if (count) diff /= count;
  return diff;
}
