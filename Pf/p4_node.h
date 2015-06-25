p4_node *p4_newNode(int nodeNum, p4_tree *aTree, int seqNum, int isLeaf, int inTree);
void p4_freeNode(p4_node *aNode);
void p4_calculateBigPDecks(p4_node *aNode);
void p4_calculateBigPDecksPart(p4_node *aNode, int pNum);
void p4_calculateBigPDecks_1stD(p4_node *aNode);
void p4_calculateBigPDecks_2ndD(p4_node *aNode);
void p4_calculatePickerDecks(p4_node *aNode);
void p4_setConditionalLikelihoodsOfInteriorNode(p4_node *aNode);
void p4_setConditionalLikelihoodsOfInteriorNodePart(p4_node *aNode, int pNum);
void p4_initializeCL2ToRootComp(p4_node *aNode);
void p4_setCL2Up(p4_node *cl2Node);
void p4_setCL2Down(p4_node *cl2Node, p4_node *aNode);
void p4_dumpCL(double ****theCL);

void p4_calculateExpectedComp(p4_node *aNode);
//void p4_calculateObservedCharFreq(p4_node *aNode, int verbose);


void p4_copyCondLikesFromNodeToNode(p4_node *nA, p4_node *nB);
int p4_verifyCondLikesFromNodeToNode(p4_node *nA, p4_node *nB);
void p4_copyBigPDecksFromNodeToNode(p4_node *nA, p4_node *nB);
int p4_verifyBigPDecksFromNodeToNode(p4_node *nA, p4_node *nB);
