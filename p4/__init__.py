"""Maximum likelihood and Bayesian phylogenetic analysis.

P4 is a phylogenetic toolkit that does Bayesian and maximum likelihood
phylogenetic analyses of molecular sequences.  It's specialty is that
you can use heterogeneous models, where the model parameters can
differ in different parts of the tree, or over different parts of the
data.

The package comes with-
   - the package itself, which goes where the third-party modules go,
     eg /usr/local/lib/python2.7/site-packages/p4/
   - the p4 script, which goes somewhere in your path.
   - documentation and examples, which go in an appropriate share
     directory, for example in /usr/local/share/doc

You can import the package as usual, by starting python and then
saying

    import p4


but my preferred way to use it is with the p4 script.  Using this
script at the command line starts python and loads in the package.
The advantage of using the p4 script is that you can read data files
from the command line.  To see if it is installed and works, try

    p4 --help

Documentation is at <http://p4.nhm.ac.uk>

"""
from __future__ import print_function
from past.builtins import execfile

import sys
import os
import glob

from .var import var
from . import func
from .func import read  # Make this one top-level, as it is used so often

from .alignment import Alignment
import numpy
from .data import Data
from .model import Model
from .mcmc import Mcmc
from .mcmccheckpointreader import McmcCheckPointReader
from .chain import Chain
from .stmcmc import STMcmc, STMcmcCheckPointReader
from .distancematrix import DistanceMatrix
from .sequencelist import Sequence, SequenceList
from .tree import Tree
from .node import Node, NodeBranch
#from .nexus import Nexus
from .treepartitions import TreePartitions
from .trees import Trees
from .pnumbers import Numbers
from .p4exceptions import P4Error
from .quartet import QuartetSet, Quartet
from .reducedstrictconsensus import Reduced, Intersection
from .triplets import Aho
from .leafsupport import LeafSupport
from .constraints import Constraints
from .quartetjoining import QuartetJoining
from .posteriorsamples import PosteriorSamples
from .var import Var  # To provide access to the property docstrings.
from .geneticcode import GeneticCode


# Read in user-defined p4 stuff on startup.  These might be var
# settings (eg var.verboseRead = 0) or might be user-defined
# functions or methods.  Both RickR and CymonC suggested having an
# P4_STARTUP environment variable.  This is set up, ready to go,
# below, but turned off by default.

# But I've started using a '.p4' directory in the user's home
# directory, that can hold several files.

# Additional p4 directories can be specified by the environment
# variable P4_STARTUP_DIRS, which would be a colon-separated list of
# directory names, eg
# /home/me/lib/p4_startup:/usr/local/lib/p4_startup.  The ~/.p4
# directory remains, and need not be listed in P4_STARTUP_DIRS

if 1:
    verboseStartupFiles = False  # Turn on for debugging...
    if 0:  # If you want it, turn it on.
        try:
            execfile(os.environ['P4_STARTUP'])
            if verboseStartupFiles:
                print('\n\n ***** ...have read p4 config file from $P4_STARTUP *****')
        except KeyError:
            pass
        except IOError:
            pass

    sdd = os.environ.get('P4_STARTUP_DIRS')
    # print("got startup dirs", sdd)
    if sdd:
        sdd = sdd.split(':')
    else:
        sdd = []
    sdd.append(os.path.join(os.path.expanduser('~'), '.p4'))
    fName = None
    pathPat = None
    pyFileNames = None

    for theDir in sdd:
        if os.path.isdir(theDir):
            pathPat = os.path.join(theDir, '*.py')
            pyFileNames = glob.glob(pathPat)
            if pyFileNames:
                for fName in pyFileNames:
                    if verboseStartupFiles:
                        print('...reading %s' % fName)
                    execfile(fName)
    del(fName)
    del(pathPat)
    del(pyFileNames)
    del(theDir)
    del(verboseStartupFiles)
    del(sdd)

if var.excepthookEditor:
    from .interactive import excepthook

