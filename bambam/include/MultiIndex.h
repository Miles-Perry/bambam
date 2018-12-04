#pragma once

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include <string>
#include <vector>

#include "api/BamConstants.h"
#include "api/BamReader.h"


// Snp contains the info for a single SNP
// also acts as a node in a singly-linked list
class Snp {
 public:
  int pos;
  int endPos;
  char** vars;
  Snp* next;

  Snp () {
    pos = endPos = -1;
    vars = NULL;
    next = NULL;
  }
  
  Snp (int _pos, std::vector<std::string> _vars, int _endPos = -1) {
    pos = _pos;
    endPos = _endPos;
    next = NULL;

    vars = (char**) malloc (sizeof(char*) * _vars.size());
    for (int i = 0; i < _vars.size(); i++) {
      vars[i] = (char*) malloc (sizeof(char) * 5);
      strcpy (vars[i], _vars.at(i).c_str());
    }
  }


  ~Snp() {
    free (vars);
  }
};


// ChrNode is a list of SnpNode* that also has head and tail pointers to manage growth of the list
// each SnpNode* indexes into a linked list of SnpNode, at the first SNP after a certain position
class ChrNode {
  friend class MultiIndex;

 public:
  ChrNode() {
    head = tail = NULL;
  }

  ~ChrNode () {
    Snp* current = head;

    while (current) {
      Snp* next = current->next;
      delete current;
      current = next;
    }

    for (int i = 0; i < bins.size(); i++) delete bins[i];
  }


  Snp* set (int i, Snp* newNode) {
    return bins[i] = newNode;
  }

  Snp* get (int i) {
    if (!head) return NULL;
    return bins[i];
  }


 private:
  Snp* head;
  Snp* tail;
  std::vector<Snp*> bins;
};



// SnpIndex is a map of chromosome IDs (from the BamReader) to ChrNodes
// also contains loadSnps, which builds whole index
class MultiIndex {
 public:
  MultiIndex(std::string file, BamTools::BamReader* bam, int _step) {
    index.resize(bam->GetReferenceCount());
    step = _step;
    ploidy = -1;

    loadSnps (file, bam);
  }


  // returns first Snp after the queried position
  Snp* get (int chr, int pos) {
    Snp* result = index[chr].get(pos/step);
    return result;
  }

  int getPloidy() {
    return ploidy;
  }

  std::vector<std::string>& getNames() {
    return names;
  }
  
  
 private:
  std::vector<ChrNode> index;
  int ploidy;
  std::vector<std::string> names;
  int step;


  // builds index from snpFile
  void loadSnps (std::string snpFile, BamTools::BamReader* bam) {
    int MAX_LENGTH = 512;
    char delim[5] = " \t\r\n";
    char buffer[MAX_LENGTH];

    int nextStep = -1;
    Snp* lastIndel = NULL;

    FILE* fh = fopen (snpFile.c_str(), "r");
    if (!fh) throw std::runtime_error ("Cannot open file " + snpFile);

    while (fgets(buffer, MAX_LENGTH, fh)) {
      
      char* chrStr = strtok (buffer, delim);
      char* posStr = strtok (NULL, delim);
      if (!chrStr || !posStr) throw std::runtime_error ("Unrecognized file format for " + snpFile);

      // build vars masks
      std::vector<std::string> varList;
      while (char* base = strtok (NULL, delim)) varList.push_back (base);
      if (varList.size() < 1) throw std::runtime_error ("No alleles found in " + snpFile);


      // get genome names
      if (buffer[0] == '#') {
	if (ploidy < 0) {
	  ploidy = varList.size();
	  for (int i = 0; i < varList.size(); i++) names.push_back (varList[i]);
	}
      }


      else {
	int chr = bam->GetReferenceID(chrStr);
	if (chr < 0) throw std::invalid_argument ("Unrecognized chromosome name " + std::string(chrStr));

	char* posStr1 = strtok (posStr, ",");
	int pos = atoi(posStr1) - 1; // shift to 0-based

	Snp* newSnp = new Snp(pos, varList);

	if (char* posStr2 = strtok (NULL, ",")) {
	  newSnp->endPos = atoi(posStr2) - 1;
	  lastIndel = newSnp;
	}

	
	if (!index[chr].head) {
	  nextStep = 0;
	  index[chr].bins.resize(bam->GetReferenceData()[chr].RefLength/step+1);
	  index[chr].head = index[chr].tail = newSnp;
	}
	else {
	  index[chr].tail->next = newSnp;
	  index[chr].tail = newSnp;
	}
	
	int actualStep = pos / step;
	while (nextStep <= actualStep) {
	  if (lastIndel && lastIndel->endPos >= nextStep*step && lastIndel->pos <= (nextStep+1)*step) index[chr].set(nextStep, lastIndel);
	  else {
	    index[chr].set(nextStep, newSnp);
	    lastIndel = NULL;
	  }
	  ++nextStep;
	}
      }
    }
    fclose (fh);
  }
};
