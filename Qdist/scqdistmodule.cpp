#include <Python.h>

// This stuff is mostly from main.cpp in the qdist code
#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <set>

// #include "Util.hpp"
#include "NewickParser.hpp"
#include "Tree.hpp"
#include "TreeUtil.hpp"
#include "QDist.hpp"

static const std::string &ExtractLeafLabel(const LeafNode *n) {
    return n->GetLabel();
}

int qdist(const char* input1, const char* input2) {


    Tree* tree1;
    Tree* tree2;

    //load newick strings from files
    //std::string filename1 = argv[1];
    //std::string input1 = Util::LoadFileToString(filename1);

    //std::string filename2 = argv[2];
    //std::string input2 = Util::LoadFileToString(filename2);

    //parse strings
    NewickParser* parser = new NewickParser();
    tree1 = parser->Parse(input1);
    tree2 = parser->Parse(input2);

    // Check that the leaf lists match
    std::set<std::string> leaves1, leaves2;
    
    std::transform(tree1->GetLeafNodes().begin(), tree1->GetLeafNodes().end(), 
                   std::inserter(leaves1, leaves1.begin()), ExtractLeafLabel);
    std::transform(tree2->GetLeafNodes().begin(), tree2->GetLeafNodes().end(), 
                   std::inserter(leaves2, leaves2.begin()), ExtractLeafLabel);
    
    if (leaves1 != leaves2) {
        std::cout << "The two trees do not have the same leaf sets!" << std::endl;
        return 1;
    }
    
    
    //make identical numbering of the leaves
    TreeUtil::RenumberTreeAccordingToOther(tree2, tree1);    

    TreeUtil::CheckTree(tree1);
    TreeUtil::CheckTree(tree2);

    //long n = leaves1.size();
    //long max_qdist = Util::Choose(n, 4);
    
    long myqdist, b1, b2, shared_b, diff_b;
    myqdist = SubCubicQDist(tree1, tree2, b1, b2, shared_b, diff_b);

    return myqdist;
    
    //long min_b = std::min(b1, b2);
    
    //std::cout << "N\tB1\tB2\tS\tD\tNorm B\tQ\tNorm Q" << std::endl;
    //std::cout << n << '\t' << b1 << '\t' << b2 << '\t' << shared_b << '\t' << diff_b << '\t' << (double(shared_b) / min_b)  << '\t' << qdist << '\t' << (double(qdist) / max_qdist) << std::endl;

}
 
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
using namespace boost::python;
 
BOOST_PYTHON_MODULE(scqdist)
{
    def("qdist", qdist);
}

