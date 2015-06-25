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

import sys, os, glob
# Check for version 2, then for version 2.3.  Using sys.version_info
# is more convenient than sys.version, but it is only in v2.
#if int(sys.version[0]) < 2:
#    print 'p4 wants Python version 2.3 or better.'
#    sys.exit()
#if sys.version_info[1] < 3:
#    print 'p4 wants Python version 2.3 or better.'
#    sys.exit()

### People have told me that it does not work when they first try it
### after installing, calling p4 from the install directory.  That is
### because the pf.so module is someplace else (ie installed) , but
### python is trying to import the local p4.  No workee, and from the
### users point of view, mysteriously.  So check ...
##try:
##    import pf
##except ImportError:
##    if os.path.isdir('p4') and os.path.isdir('Pf'):
##        #print "cwd is %s" % os.getcwd()
##        # Seems to be in the install directory?
##        if os.path.isdir('build') and not os.path.isfile('p4/pf.so'):
##            try:
##                from Glitch import Glitch
##                gm = ["It looks like you are importing the p4 package from the"]
##                gm.append("install directory, after it has been installed.")
##                gm.append("That does not work.")
##                gm.append("Try going to another directory.")
##                gm.append("If that does not work, maybe there was another installation problem.")
##                raise Glitch, gm
##            except:
##                raise

from Var import var
import func
from func import read # Make this one top-level, as it is used so often
#import p3rlcompleter  # Although pyrepl sounds very interesting ...

from Alignment import Alignment
if var.usePfAndNumpy:
    import numpy
    from Data import Data
    from Model import Model
    from Mcmc import Mcmc
    from McmcCheckPointReader import McmcCheckPointReader
    from Chain import Chain
    from STMcmc import STMcmc,STMcmcCheckPointReader
from DistanceMatrix import DistanceMatrix
from SequenceList import Sequence,SequenceList
from Tree import Tree
from Node import Node,NodeBranch
#from Nexus import Nexus
from TreePartitions import TreePartitions
from Trees import Trees
from Numbers import Numbers
from Glitch import Glitch
from Quartet import QuartetSet, Quartet
from ReducedStrictConsensus import Reduced, Intersection
from Triplets import Aho
from LeafSupport import LeafSupport
from Constraints import Constraints
from QuartetJoining import QuartetJoining
from PosteriorSamples import PosteriorSamples
from Var import Var  # To provide access to the property docstrings.


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
    verboseStartupFiles = False # Turn on for debugging...
    if 0:  # If you want it, turn it on.
        try:
            execfile(os.environ['P4_STARTUP'])
            if verboseStartupFiles:
                print '\n\n ***** ...have read p4 config file from $P4_STARTUP *****'
        except KeyError:
            pass
        except IOError:
            pass

    sdd = os.environ.get('P4_STARTUP_DIRS')
    #print "got startup dirs", sdd
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
                        print '...reading %s' % fName
                    execfile(fName)
    del(fName)
    del(pathPat)
    del(pyFileNames)
    del(theDir)
    del(verboseStartupFiles)
    del(sdd)
