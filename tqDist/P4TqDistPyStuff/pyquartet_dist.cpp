#include <Python.h>

#include <iostream>
#include <cstring>

#include "int_stuff.h"
#include "QuartetDistanceCalculator.h"
#include "newick_parser.h"


#ifndef _MSC_VER
#define _stricmp strcasecmp
#endif



int qdistFromStrings(char* treestring1, char* treestring2) {

  UnrootedTree *ut1 = NULL;
  UnrootedTree *ut2 = NULL;
  NewickParser parser;

  // std::cout << "Got treestring 1 " << treestring1 << std::endl;
  // std::cout << "Got treestring 2 " << treestring2 << std::endl;



  ut1 = parser.parseStr(treestring1);
  if (ut1 == NULL || parser.isError()) {
    std::cerr << "Error: Parsing of treestring 1\"" << treestring1 << "\" failed." << endl;
    std::cerr << "Aborting!" << endl;
    return -1;
  }

  ut2 = parser.parseStr(treestring2);
  if(ut2 == NULL || parser.isError()) {
    cerr << "Parsing of tree string 2 \"" << treestring2 << "\" failed." << endl;
    cerr << "Aborting!" << endl;
    return -1;
  }

  QuartetDistanceCalculator quartetCalc;
  INTTYPE_N4 dist = quartetCalc.calculateQuartetDistance(ut1, ut2);

  if(dist == -1)
    exit(-1);

  // std::cout << "qdistFromStrings() " << dist << std::endl;
  return dist;
}

int qdistFromStringsVerbose(char* treestring1, char* treestring2) {

  UnrootedTree *ut1 = NULL;
  UnrootedTree *ut2 = NULL;
  NewickParser parser;

  // std::cout << "Got treestring 1 " << treestring1 << std::endl;
  // std::cout << "Got treestring 2 " << treestring2 << std::endl;


  // File reading is bypassed, but so is checking for spaces.
  ut1 = parser.parseStr(treestring1);
  if (ut1 == NULL || parser.isError()) {
    std::cerr << "Error: Parsing of treestring 1\"" << treestring1 << "\" failed." << endl;
    std::cerr << "Aborting!" << endl;
    return -1;
  }

  ut2 = parser.parseStr(treestring2);
  if(ut2 == NULL || parser.isError()) {
    cerr << "Parsing of tree string 2 \"" << treestring2 << "\" failed." << endl;
    cerr << "Aborting!" << endl;
    return -1;
  }

  QuartetDistanceCalculator quartetCalc;
  INTTYPE_N4 dist = quartetCalc.calculateQuartetDistance(ut1, ut2);

  if(dist == -1)
    exit(-1);

    INTTYPE_N4 resolvedQuartetsAgree = quartetCalc.get_resolvedQuartetsAgree();
    INTTYPE_N4 resolvedQuartetsAgreeDiag = quartetCalc.get_resolvedQuartetsAgreeDiag();
    INTTYPE_N4 resolvedQuartetsDisagree = quartetCalc.get_resolvedQuartetsDisagree();
    INTTYPE_N4 resolvedQuartetsDisagreeDiag = quartetCalc.get_resolvedQuartetsDisagreeDiag();
    INTTYPE_N4 resolvedQuartetsAgreeUpper = quartetCalc.get_resolvedQuartetsAgreeUpper();
    INTTYPE_N4 resolvedQuartetsDisagreeUpper = quartetCalc.get_resolvedQuartetsDisagreeUpper();

    INTTYPE_N4 n = quartetCalc.get_n();
    INTTYPE_N4 totalNoQuartets = quartetCalc.get_totalNoQuartets();
    double dist_norm = double(dist) / double(totalNoQuartets);
    INTTYPE_N4 resAgree = resolvedQuartetsAgree + resolvedQuartetsAgreeDiag + resolvedQuartetsAgreeUpper;
    double resAgree_norm = double(resAgree) / double(totalNoQuartets);
    INTTYPE_N4 unresolvedQuartetsAgree = quartetCalc.get_unresolvedQuartets();
    double unresolvedQuartetsAgree_norm = double(unresolvedQuartetsAgree) / double(totalNoQuartets);
    
    std::cout << n                            << "\t"
	      << totalNoQuartets              << "\t"
	      << dist                         << "\t"
	      << dist_norm                    << "\t"
	      << resAgree                     << "\t"
	      << resAgree_norm                << "\t"
	      << unresolvedQuartetsAgree      << "\t"
	      << unresolvedQuartetsAgree_norm << std::endl;


  // std::cout << "qdistFromStrings() " << dist << std::endl;
  return dist;
}



int qdistFromFiles(char* filename1, char* filename2) {


  // Check, then turn off
  // #ifdef quartetsToo
  //   std::cout << "quartetsToo is defined" << std::endl;
  // #endif

  QuartetDistanceCalculator quartetCalc;
  INTTYPE_N4 dist = quartetCalc.calculateQuartetDistance(filename1, filename2);

  if(dist == -1)
    exit(-1);

  return dist;
}

int qdistFromFilesVerbose(char* filename1, char* filename2) {

  QuartetDistanceCalculator quartetCalc;
  INTTYPE_N4 dist = quartetCalc.calculateQuartetDistance(filename1, filename2);

  if(dist == -1)
    exit(-1);

    INTTYPE_N4 resolvedQuartetsAgree = quartetCalc.get_resolvedQuartetsAgree();
    INTTYPE_N4 resolvedQuartetsAgreeDiag = quartetCalc.get_resolvedQuartetsAgreeDiag();
    INTTYPE_N4 resolvedQuartetsDisagree = quartetCalc.get_resolvedQuartetsDisagree();
    INTTYPE_N4 resolvedQuartetsDisagreeDiag = quartetCalc.get_resolvedQuartetsDisagreeDiag();
    INTTYPE_N4 resolvedQuartetsAgreeUpper = quartetCalc.get_resolvedQuartetsAgreeUpper();
    INTTYPE_N4 resolvedQuartetsDisagreeUpper = quartetCalc.get_resolvedQuartetsDisagreeUpper();

    INTTYPE_N4 n = quartetCalc.get_n();
    INTTYPE_N4 totalNoQuartets = quartetCalc.get_totalNoQuartets();
    double dist_norm = double(dist) / double(totalNoQuartets);
    INTTYPE_N4 resAgree = resolvedQuartetsAgree + resolvedQuartetsAgreeDiag + resolvedQuartetsAgreeUpper;
    double resAgree_norm = double(resAgree) / double(totalNoQuartets);
    INTTYPE_N4 unresolvedQuartetsAgree = quartetCalc.get_unresolvedQuartets();
    double unresolvedQuartetsAgree_norm = double(unresolvedQuartetsAgree) / double(totalNoQuartets);
    
    std::cout << n                            << "\t"
	      << totalNoQuartets              << "\t"
	      << dist                         << "\t"
	      << dist_norm                    << "\t"
	      << resAgree                     << "\t"
	      << resAgree_norm                << "\t"
	      << unresolvedQuartetsAgree      << "\t"
	      << unresolvedQuartetsAgree_norm << std::endl;

  return dist;
}

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
 
BOOST_PYTHON_MODULE(pytqdist)
{
    using namespace boost::python;
    scope().attr("__doc__") = "A wrapper on tqDist quartet_dist from BIRC.  Includes an interface via strings.  If interfacing via strings, then spaces in the strings are not allowed.";
    def("qdistFromFiles", qdistFromFiles);
    def("qdistFromStrings", qdistFromStrings);
    def("qdistFromFilesVerbose", qdistFromFilesVerbose);
    def("qdistFromStringsVerbose", qdistFromStringsVerbose);
}
