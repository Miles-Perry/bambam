#pragma once

#include <cstdlib>
#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <cassert>


class Arg {
 public:
  std::string id;
  std::string def;
  std::string val;
  std::string help;

  Arg (const char* _id, const char* _def, const char* _help) : id(_id), def(_def), val(_def), help(_help) {
    if (def.length() == 0) def = "required";
  }
};


class ArgList {
 public:
  std::string program;
  std::string description;
  std::vector<Arg> params;
  std::vector<Arg> opts;
  int cur;
  bool good;
  
  void printHelp() {
    // print version
    printf ("\nBamBam Version 1.3\nContact Justin (jtpage68@gmail.com) with any problems\n");

    // print description
    if (description.length()) printf ("\nSUMMARY: %s\n", description.c_str());

    // print usage statement
    printf ("\nUSAGE: %s <options>", program.c_str());
    for (int i = 0; i < params.size(); i++) {
      printf (" %s", params[i].id.c_str());
    }
    printf ("\n");
    
    // print parameter information
    for (int i = 0; i < params.size(); i++) {
      printf ("\t%s\t%s\n", params[i].id.c_str(), params[i].help.c_str());
    }
    printf ("\nOPTIONS:\n");
    for (int i = 0; i < opts.size(); i++) {
      printf ("\t%s\t(%s) %s\n", opts[i].id.c_str(), opts[i].def.c_str(), opts[i].help.c_str());
    }
    printf ("\n");
  }
  
 ArgList(int argc, char* argv[], int optc, const char* optv[], const std::string& descrip = "") : good(true), cur(0), description(descrip) {
    program = argv[0];
    
    for (int i = 0; i < 3*optc; i+=3) {
      Arg a (optv[i], optv[i+1], optv[i+2]);
      if (optv[i][0] == '-') {
	assert (strlen(optv[i]) == 2);
	opts.push_back (a);
      }
      else params.push_back (a);
    }
    
    // parse args
    for (int i = 1; good && i < argc; i++) {
      
      if (argv[i][0] == '-') { // optional parameter

	// if parameter doesn't have a space between the -x and the value
	if (strlen(argv[i]) > 2) {
	  std::string temp = argv[i];
	  std::string id = temp.substr(0,2);
	  std::string val = temp.substr(2,temp.length()-2);

	  for (int i = 0; i < opts.size(); i++) {
	    if (opts[i].id == id) opts[i].val = val;
	  }
	}

	// if parameter does have a space
	else if (i+1 >= argc) good = false; // make sure there's a value
	else {
	  std::string id = argv[i];
	  std::string val = argv[++i];
	  
	  for (int i = 0; i < opts.size(); i++) {
	    if (opts[i].id == id) opts[i].val = val;
	  }
	}
      }

      else {
	if (cur >= params.size()) params.push_back (params[params.size()-1]); // allow multiple files
	params[cur++].val = argv[i]; // core parameter
      }
    }
    
    // check completeness
    if (cur < params.size()) {
      good = false;
    }
    for (int i = 0; good && i < opts.size(); i++) {
      if (opts[i].val.length() < 1) {
	good = false;
      }
    }
    

    if (!good) printHelp();
  }
    
  const char* get (int i) {
    if (i < params.size()) return params[i].val.c_str();
    else return opts[i-params.size()].val.c_str();
  }

  const char* get (const char* query) {
    for (int i = 0; i < opts.size(); i++) {
      if (!strcmp (query, opts[i].id.c_str())) return opts[i].val.c_str();
    }
    throw std::invalid_argument ("Unrecognized argument query " + std::string(query));
  }

  std::vector<std::string> getAll() {
    std::vector<std::string> result;
    for (int i = 0; i < params.size(); i++) result.push_back (params[i].val);
    return result;
  }

  int size() {
    return params.size();
  }
};
