#include <Python.h>

// This stuff is modified from quartet_dist.cpp in the tqdist code
#include <iostream>
#include <cstring>

#include "int_stuff.h"
#include "QuartetDistanceCalculator.h"
#include "newick_parser.h"

#ifndef _MSC_VER
#define _stricmp strcasecmp
#endif

int qdist(char* treestring1, char* treestring2) {

  UnrootedTree *ut1 = NULL;
  UnrootedTree *ut2 = NULL;
  NewickParser parser;


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

  return dist;
}


#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
using namespace boost::python;
 
BOOST_PYTHON_MODULE(tqdist)
{
    def("qdist", qdist);
}

