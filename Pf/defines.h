
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

#define PINVAR_MIN 0.0
#define PINVAR_MAX 0.99
#define KAPPA_MIN 0.000001
#define KAPPA_MAX 100.0
#define GAMMA_SHAPE_MIN 0.000001
#define GAMMA_SHAPE_MAX 300.0
#define PIVEC_MIN 1.0e-18 // changed from zero jan04, due to a problem with exp(bigQ), see setPrams()
#define PIVEC_MAX 0.999
#define RATE_MIN 1.0e-8
#define RATE_MAX 1.0e8
#define COVARION_S_MIN 1.0e-15
#define COVARION_S_MAX 1.0e3
#define RELRATE_MIN 1.0e-8
#define RELRATE_MAX 1.0e8
#define BRLEN_MIN 1.0e-8  // changed apr02 from 1e-18
#define BRLEN_MAX 3.0
#define MIXTURE_FREQ_MIN 0.0
#define MIXTURE_FREQ_MAX 1.0
#define MIXTURE_RATE_MIN 1.0e-3
#define MIXTURE_RATE_MAX 1.0e3


#define OPT_PIVEC      0
#define OPT_RMATRIX    1
#define OPT_KAPPA      2
#define OPT_SHAPE      3
#define OPT_COVARION   4
#define OPT_PINVAR     5
#define OPT_RELRATE    6
#define N_OPT          7


#define TUNING_CHAIN_TEMP 0
#define TUNING_PIVEC      1
#define TUNING_RMATRIX    2
#define TUNING_KAPPA      3
#define TUNING_SHAPE      4
#define TUNING_PINVAR     5
#define TUNING_RELRATE    6
#define TUNING_LOCAL      7
#define TUNING_WORM       8
#define TUNING_ROOT_2M    9
#define TUNING_BRLEN     10
#define N_TUNINGS        11

#define CHANGE_PIVEC        0
#define CHANGE_RMATRIX      1
#define CHANGE_KAPPA        2
#define CHANGE_SHAPE        3
#define CHANGE_COVARION     4
#define CHANGE_PINVAR       5
#define CHANGE_RELRATE      6
#define CHANGE_LOCAL        7
#define CHANGE_JPH_LOCAL    8
#define CHANGE_WORM         9
#define CHANGE_MODEL       10
#define CHANGE_ROOT_3      11
#define CHANGE_ROOT_3M     12
#define CHANGE_ROOT_2M     13
#define CHANGE_BRLEN       14
#define N_CHANGES          15

#define TOPOL_LOCAL         1
#define TOPOL_JPH_LOCAL     2
#define TOPOL_WORM          3
#define TOPOL_ROOT_2M       4
#define N_TOPOL             4

#define SAME       0
#define DIFFERENT  1
//#define ALL        0
//#define SOME       1
#define NO_ORDER   -10000

#define PF_MALLOC malloc
#define PF_FREE free

