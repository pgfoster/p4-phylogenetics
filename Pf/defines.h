
#define RMATRIX_NONE	0
#define RMATRIX_JC		1
#define RMATRIX_K2P		2  // not used
#define RMATRIX_F81		3
#define RMATRIX_HKY		4  // not used
#define RMATRIX_2P      5  // Used for K2P and HKY
#define RMATRIX_SPECIFIED	20


#define RMATRIX_POISSON	100
#define RMATRIX_CPREV		101
#define RMATRIX_D78		102
#define RMATRIX_JTT		103
#define RMATRIX_MTREV24		104
#define RMATRIX_MTMAM		105
#define RMATRIX_WAG		106
#define RMATRIX_BLOSUM62_A  107 // from MrBayes and phyml, the one I have
#define RMATRIX_BLOSUM62_B  108 // no got
#define RMATRIX_PHAT70    109
#define RMATRIX_RTREV     110
#define RMATRIX_TMJTT94   111
#define RMATRIX_TMLG99   112
#define RMATRIX_LG   113
#define RMATRIX_HIVB 114
#define RMATRIX_MTART 115
#define RMATRIX_MTZOA 116
#define RMATRIX_GCPREV 117
#define RMATRIX_STMTREV 118
#define RMATRIX_PRASREV 121
#define RMATRIX_GNETREV 120
#define RMATRIX_VT 119

#define PIVEC_NONE	0
#define PIVEC_TUPLE	1

#define GAP_CODE       -1
#define QMARK_CODE     -2
#define N_LIKE         -3
#define EQUATES_BASE   -64

#define DATATYPE_NONE  0
#define DATATYPE_STANDARD  1
#define DATATYPE_DNA 2
#define DATATYPE_PROTEIN 3

// #define PINVAR_MIN 0.0
// #define PINVAR_MAX 0.99
// #define KAPPA_MIN 0.000001
// #define KAPPA_MAX 100.0
// #define GAMMA_SHAPE_MIN 0.000001
// #define GAMMA_SHAPE_MAX 300.0
// #define PIVEC_MIN 1.0e-18 // changed from zero jan04, due to a problem with exp(bigQ), see setPrams()
// #define PIVEC_MAX 0.999
// #define RATE_MIN 1.0e-14  // changed from 1e-8 nov 2016
// #define RATE_MAX 0.9999999  // changed from 1.0e8  nov 2016
// #define RELRATE_MIN 1.0e-8
// #define RELRATE_MAX 1.0e8
// #define BRLEN_MIN 1.0e-8  // changed apr02 from 1e-18
// #define BRLEN_MAX 3.0

#define SAME       0
#define DIFFERENT  1
#define NO_ORDER   -10000
