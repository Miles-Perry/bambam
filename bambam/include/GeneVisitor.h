// GeneVisitor
// Justin Page
// jtpage68@gmail.com

#include <iostream>
#include <string>
#include <api/BamReader.h>

#define MAX_LINE 2048


class Gene {
 public:
  int chr;
  int start;
  int end;
  std::string id;

  Gene (int _chr, int _start, int _end, std::string& _id) : chr(_chr), start(_start), end(_end), id(_id) {}
};


// visit each region denoted in annotFile (gff/bed) and perform function fn
void visitGenes (BamTools::BamReader* bam, void (*fn)(Gene g), const std::string annotFile = "0") {

  if (annotFile != "0") {
    // determine file type                                               
    int pos = annotFile.find_last_of(".");
    if (pos == std::string::npos) throw std::invalid_argument ("Cannot identify annotation file format for " + annotFile);
    
    std::string type = annotFile.substr(pos+1,3);
    for (int i = 0; i < type.length(); i++) type[i] = tolower(type[i]);
    if (type != "gff" && type != "bed") throw std::invalid_argument ("Unrecognized annotation file format " + type);
    
    
    std::ifstream in;
    in.open(annotFile.c_str());
    char buffer[MAX_LINE];
    while (in.getline (buffer, MAX_LINE)) {
      
      // parse values
      std::string chrStr = strtok (buffer, "\t");
      if (chrStr.substr(0,5) == "track" || chrStr.substr(0,7) == "browser") continue; // header types
      if (chrStr.substr(0,1) == "#") continue; // GFF comment

      int chr = bam->GetReferenceID (chrStr);
      if (chr < 0) throw std::invalid_argument ("Unrecognized chromosome name " + chrStr);
      
      if (type == "gff") {
	strtok (NULL, "\t");
	strtok (NULL, "\t");
      }

      std::string startStr = strtok (NULL, "\t");
      std::string endStr = strtok (NULL, "\t");

      if (type == "gff") {
	strtok (NULL, "\t");
	strtok (NULL, "\t");
	strtok (NULL, "\t");
      }

      char* idStr = strtok (NULL, "\t");
      std::string id = "";
      if (idStr) id = idStr;
      else id = chrStr + ":" + startStr + "-" + endStr;

      int start = atoi(startStr.c_str())-1;
      int end = atoi(endStr.c_str());
      
      // call function
      Gene g (chr, start, end, id);
      fn(g);
    }
    in.close();
  }

  // use BAM file for reference sequences
  else {
    for (int chr = 0; chr < bam->GetReferenceCount(); chr++) {
      
      // get chromosome information
      std::string id = bam->GetReferenceData()[chr].RefName;
      int start = 0;
      int end = bam->GetReferenceData()[chr].RefLength;

      // call function
      Gene g (chr, start, end, id);
      fn(g);
    }
  }
}
