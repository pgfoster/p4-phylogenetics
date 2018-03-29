#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <vector>
#include <cmath>

namespace p = boost::python;
namespace np = boost::python::numpy;

class InTrNo       // inTree Node
{
public:
    boost::dynamic_bitset<> splitKey;
    double support;

    InTrNo(std::string splKeyStr, double support): 
        splitKey(splKeyStr), support(support)
        {
        }
};

class InTr
{
public:
    int tNum;
    int nTax;
    int nAllTax;
    boost::dynamic_bitset<> taxBits;
    int firstTaxPos;
    std::vector<InTrNo*> internals;

    InTr(int tNum, int nTax, int nAllTax, std::string taxBitsStr, int firstTaxPos): 
        tNum(tNum), nTax(nTax), nAllTax(nAllTax), taxBits(taxBitsStr), firstTaxPos(firstTaxPos)
        {
        } 

    void setInTrNo(std::string splKeyStr, double support) {
        InTrNo* myInTrNo = new InTrNo(splKeyStr, support);
        //std::cout << "InTr setInTrNo  splitKey "  << myInTrNo->splitKey << std::endl;  doesn't work!
        internals.push_back(myInTrNo);
    }
};

class BigTNo
{
public:
    int nodeNum;
    boost::dynamic_bitset<> spl; 
    boost::dynamic_bitset<> spl2;
    boost::dynamic_bitset<> theSpl;
    boost::dynamic_bitset<> maskedSplitWithFirstTaxOne;


    BigTNo(int nodeNum):
        nodeNum(nodeNum)
        {
        }

    void setSpl(std::string splStr, std::string splStr2) {
        boost::dynamic_bitset<> splx(splStr);       // Why is this so awkward? Why can't I set spl directly?
        spl = splx;
        boost::dynamic_bitset<> splxx(splStr2);
        spl2 = splxx;
        //boost::dynamic_bitset<>spl(splStr);      // This does not work.
        //boost::dynamic_bitset<>spl2(splStr2);
        //std::cout << "BigTNo setSpl(), setting spl to " << spl  << std::endl;
        //std::cout << "BigTNo setSpl(), setting spl2 to " << spl2  << std::endl;
    }

};


class BigT
{
public:
    int nNo;
    int nTax;
    //int nAllTax;
    //BigTNo* root;
    np::ndarray postOrder;
    np::ndarray spaQ;
    std::vector<BigTNo*> nodes;
    int btNum;

    BigT(int nNo, int nTax, np::ndarray postOrder, np::ndarray spaQ, int btNum): 
        nNo(nNo), nTax(nTax), postOrder(postOrder), spaQ(spaQ), btNum(btNum)
        {

            //std::cout << "BigT init: " << p::extract<char const *>(p::str(postOrder)) << std::endl;
            for(int i=0; i < nNo; i++) {
                BigTNo* newNo = new BigTNo(i);
                nodes.push_back(newNo);
            }
        } 
        
};  // end of BigT class

class FastSpa
{
public:
    std::vector<BigT*> bigTT;
    std::vector<InTr*> inTrees;
    int useSplitSupport;
    double logLike;

    FastSpa(int useSplitSupport):
        useSplitSupport(useSplitSupport)
        {
        }

    void summarizeInTrs() {
        InTr* inTr;

        for(size_t i = 0; i < inTrees.size(); i++){
            inTr = inTrees[i];
            std::cout << i << "      =inTr " << inTr->tNum << " inTr->taxBits " << inTr->taxBits << std::endl;
            for(size_t i2 = 0; i2 < inTr->internals.size(); i2++){
                std::cout << "               =splitKey " << inTr->internals[i2]->splitKey << std::endl;
            }
        }
    }

    void setInTr(int tNum, int nTax, int nAllTax, std::string taxBitsStr, int firstTaxPos) {
        InTr* inTr;
        inTr = new InTr(tNum, nTax, nAllTax, taxBitsStr, firstTaxPos);
        inTrees.push_back(inTr);
    }

    void setInTrNo(int tNum, std::string splitKeyStr, double support) {
        InTr* inTr;

        inTr = inTrees[tNum];
        inTr->setInTrNo(splitKeyStr, support);
    }

    void setBigT(int nNo, int nTax, np::ndarray postOrder, np::ndarray spaQ) {
        BigT* bigT;
        int btNum;

        //std::cout << "FastSpa setBigT: " << p::extract<char const *>(p::str(postOrder)) << std::endl;
        btNum = bigTT.size();
        bigT = new BigT(nNo, nTax, postOrder, spaQ, btNum);
        bigTT.push_back(bigT);
    }

    void setBigTNoSpl(int btNum, int nNum, std::string splStr, std::string spl2Str) {
        BigT* bigT;
        BigTNo* bigTNo;

        bigT = bigTT[btNum];
        bigTNo = bigT->nodes[nNum];
        bigTNo->setSpl(splStr, spl2Str);
    }

    double calcLogLike(int btNum) {
        //std::cout << std::endl << "------Likelihood calculator here -----"  << std::endl;
        BigT* bigT = bigTT[btNum];
        std::vector<BigTNo*> relevantSplits;
        std::map<boost::dynamic_bitset<>, BigTNo*> nodeForSplitDict;
        InTr* inTr;
        InTrNo* inTrNo;
        BigTNo* bigTNo;
        int* myPostOrder;
        double* mySpaQ;
        int nNum;
        size_t i;
        int j, onesCount, upperGood;
        double S, S_st, q, r, R, logq, logr;

        myPostOrder = (int *)bigT->postOrder.get_data(); 
        mySpaQ = (double *)bigT->spaQ.get_data();
        // std::cout << "myPostOrder: ";
        // for(j=0; j<bigT->nNo; j++) {
        //     std::cout << " " << myPostOrder[j];
        // }
        // std::cout << std::endl;

        //summarizeInTrs();
        //std::cout << std::endl;

        logLike = 0.0;
        for(i = 0; i < inTrees.size(); i++){
            inTr = inTrees[i];
            //std::cout << "InTr " << i << std::endl;
            upperGood = inTr->nTax - 2;
            relevantSplits.clear();
            nodeForSplitDict.clear();
            for(j = 0; j < bigT->nNo; j++) {
                nNum = myPostOrder[j];
                if(nNum == -10000) {
                    break;
                }
                bigTNo = bigT->nodes[nNum];
                if(bigTNo->spl.size() > 1) {
                    //for(int k=0; k < nTax; k++) {
                    //    std::cout << "bigTNo->spl[" << k << "] is " << bigTNo->spl[k] << std::endl;
                    //}
                    //std::cout << "    bigTNo spl " << bigTNo->spl << ", spl2 " << bigTNo->spl2;
                    //std::cout << ", inTr->firstTaxPos " << inTr->firstTaxPos << std::endl;
                    // dbitset[0] is the least significant bit, ie on the right.
                    // That is opposite to Python bitset
                    if(bigTNo->spl[(bigT->nTax-1) - (inTr->firstTaxPos)]) {    
                        bigTNo->theSpl = bigTNo->spl;
                    }
                    else {
                        bigTNo->theSpl = bigTNo->spl2;
                    }
                    //std::cout << "    theSpl " << bigTNo->theSpl << " taxBits " << inTr->taxBits << std::endl;
                    bigTNo->maskedSplitWithFirstTaxOne = bigTNo->theSpl;
                    bigTNo->maskedSplitWithFirstTaxOne &= inTr->taxBits;
                    //std::cout << "    maskedSplitWithFirstTaxOne " << bigTNo->maskedSplitWithFirstTaxOne << std::endl;
                    onesCount =  bigTNo->maskedSplitWithFirstTaxOne.count();
                    if(onesCount >= 2 && onesCount <= upperGood) {
                        //std::cout << "          -> relevant" << std::endl;
                        relevantSplits.push_back(bigTNo);
                        nodeForSplitDict[bigTNo->maskedSplitWithFirstTaxOne] = bigTNo;
                    }
                    //std::cout << "got " << relevantSplits.size() << "relevantSplits" << std::endl;
                }
            }
            //std::cout << std::endl;

            // S is the number of possible splits in an it-sized tree
            //  S_st is the number of splits in the reduced supertree
            // S_st nodeForSplitDict.size()
            //S = 2 ** (inTr->nTax - 1) - (inTr->nTax + 1);
            S = pow(2.0, inTr->nTax - 1) - ((double)inTr->nTax + 1.0);
            S_st = (double)nodeForSplitDict.size(); 
            //std::cout << "S= " << S << ", S_st= " << S_st << ", S - S_st " << S-S_st << std::endl;
            if(S_st) {
                //q = bigT->spaQ / S_st;
                //R = 1. - bigT->spaQ;
                q = mySpaQ[0] / S_st;
                R = 1. - mySpaQ[0];
                r = R / (S - S_st);
                //std::cout << "q " <<  q << std::endl;
                logq = log(q);
            }
            else {
                R = 1.;
                r = R / S;
            }
            //std::cout << "r " <<  r << std::endl;
            logr = log(r);
            // if(S_st) {
            //     std::cout << "logq " << logq << " logr " << logr << std::endl;
            // }
            // else {
            //     std::cout << "logq is undefined, "  << " logr " << logr << std::endl;
            // }
            for(j = 0; j < (int)inTr->internals.size(); j++) {
                inTrNo = inTr->internals[j];
                if(nodeForSplitDict[inTrNo->splitKey]) {
                    if(useSplitSupport && inTrNo->support >= 0.0) {
                        logLike += log(r + (inTrNo->support * (q - r)));
                    }
                    else {
                        logLike += logq;
                    }
                }
                else {
                    logLike += logr;
                }
            }
        }
        //std::cout << "logLike " << logLike << std::endl;
        return(logLike);
    }

};  // end of FastSpa class


using namespace boost::python;
 
BOOST_PYTHON_MODULE(fastspa)
{
    Py_Initialize();
    np::initialize();

    class_<InTr>("InTr", init<int, int, int, std::string, int>())
        .def("setInTrNo", &InTr::setInTrNo)
        ;

    class_<BigT>("BigT", init<int, int, np::ndarray, np::ndarray, int>())
        ;

    class_<FastSpa>("FastSpa", init<int>())
        .def("summarizeInTrs", &FastSpa::summarizeInTrs)
        .def("setInTr", &FastSpa::setInTr)
        .def("setInTrNo", &FastSpa::setInTrNo)
        .def("setBigT", &FastSpa::setBigT)
        .def("setBigTNoSpl", &FastSpa::setBigTNoSpl)
        .def("calcLogLike", &FastSpa::calcLogLike)
        ;

}

