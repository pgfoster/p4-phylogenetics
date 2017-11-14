from __future__ import print_function
import p4.pf as pf
import p4.func
from p4.var import var
import math
import random
import string
import sys
import time
import copy
import os
import pickle
from p4.chain import Chain
from p4.p4exceptions import P4Error
from p4.treepartitions import TreePartitions
from p4.constraints import Constraints
import datetime
import numpy

# for proposal probs
fudgeFactor = {}
fudgeFactor['local'] = 1.0
fudgeFactor['brLen'] = 1.0
fudgeFactor['eTBR'] = 1.0
#fudgeFactor['treeScale'] = 1.0
fudgeFactor['allBrLens'] = 1.0
fudgeFactor['polytomy'] = 1.0
fudgeFactor['rjComp'] = 2.0
fudgeFactor['rjRMatrix'] = 2.0
fudgeFactor['root3'] = 0.02
fudgeFactor['compLocation'] = 0.01
fudgeFactor['rMatrixLocation'] = 0.01
fudgeFactor['gdasrvLocation'] = 0.01
fudgeFactor['allCompsDir'] = 1.0
fudgeFactor['ndch2comp'] = 0.2
fudgeFactor['ndch2alpha'] = 0.04



class McmcTuningsPart(object):

    def __init__(self, partNum):
        object.__setattr__(self, 'num', partNum)
        # It is no longer changed depending on the dim
        object.__setattr__(self, 'comp', 0.3)
        # This would depend on the dim; this is done in Mcmc.__init__()
        object.__setattr__(self, 'compDir', 100.)
        object.__setattr__(self, 'allCompsDir', 500.)
        object.__setattr__(self, 'ndch2_leafCompsDir', 200.)
        object.__setattr__(self, 'ndch2_internalCompsDir', 100.)
        #object.__setattr__(self, 'rjComp', 200.)
        # rMatrix with sliders longer changed depending on the dim (ie size of rMatrix)
        object.__setattr__(self, 'rMatrix', 0.3)
        # rMatrixDir would depend on the dim; this is done in Mcmc.__init__()
        object.__setattr__(self, 'rMatrixDir', 200.)
        object.__setattr__(self, 'twoP', 50.)
        #object.__setattr__(self, 'rjRMatrix', 300.)
        object.__setattr__(self, 'gdasrv', 2.0 * math.log(1.5))  # 0.811
        #object.__setattr__(self, 'gdasrv', 2.0 * math.log(2.0))
        #object.__setattr__(self, 'gdasrvPriorLambda', 0.5)
        object.__setattr__(self, 'pInvar', 0.5)
        object.__setattr__(self, 'compLocation', 0.0)
        object.__setattr__(self, 'rMatrixLocation', 0.0)
        #object.__setattr__(self, 'cmd1_p', 100.0)
        #object.__setattr__(self, 'cmd1_s', 100.0)
        #object.__setattr__(self, 'cmd1_lna', 100.0)
        #object.__setattr__(self, 'cmd1_lnt', 100.0)

    def __setattr__(self, item, val):
        # print "Got request to set %s to %s" % (item, val)
        if item in self.__dict__.keys():
            if not isinstance(val, float):
                gm = ["\nMcmcTuningsPart.__setattr__()  Part %i" % self.num]
                gm.append('Tunings must be floats.')
                raise P4Error(gm)
            # Sanity checking goes here.
            if item in ['comp', 'rMatrix']:
                if val > var.mcmcMaxCompAndRMatrixTuning:
                    gm = ["\nMcmcTuningsPart.__setattr__()  Part %i" %
                          self.num]
                    gm.append("Maximum tuning for '%s' is %s.  Got attempt to set it to %s" % (
                        item, var.mcmcMaxCompAndRMatrixTuning, val))
                    raise P4Error(gm)
            # print "    Part %i, setting tuning of '%s' to %s" % (self.num,
            # item, val)
            object.__setattr__(self, item, val)
        else:
            gm = ["\nMcmcTuningsPart.__setattr__()  Part %i" % self.num]
            gm.append("    Can't set tuning '%s'-- no such tuning." % item)
            raise P4Error(gm)


class McmcTunings(object):

    def __init__(self, nParts):
        object.__setattr__(self, 'nParts', nParts)
        object.__setattr__(self, 'parts', [])
        for pNum in range(nParts):
            self.parts.append(McmcTuningsPart(pNum))
        object.__setattr__(self, 'chainTemp', 0.15)  # was 0.2
        object.__setattr__(self, 'relRate', 0.5)
        # This next tuning is set so that by default the brLens go up or down
        # maximum 10%, ie from 0.909 to 1.1
        object.__setattr__(self, 'local', 2.0 * math.log(1.1))  # 0.1906
        #object.__setattr__(self, 'local', 100.)
        object.__setattr__(self, 'brLen', 2.0 * math.log(2.0))  # 1.386
        # Crux has  2.0 * math.log(1.6) as default
        object.__setattr__(self, 'etbrLambda', 2.0 * math.log(1.6))
        object.__setattr__(self, 'etbrPExt', 0.8)
        object.__setattr__(self, 'brLenPriorLambda', 10.0)
        object.__setattr__(self, 'brLenPriorLambdaForInternals', 1000.0)
        object.__setattr__(self, 'doInternalBrLenPrior', False)
        object.__setattr__(self, 'doPolytomyResolutionClassPrior', False)
        object.__setattr__(self, 'polytomyPriorLogBigC', 0.0)
        object.__setattr__(self, 'brLenPriorType', 'exponential')
        #object.__setattr__(self, 'treeScale', 2.0 * math.log(1.1))
        object.__setattr__(self, 'allBrLens', 2.0 * math.log(1.02))

    def __setattr__(self, item, val):
        # print "Got request to set %s to %s" % (item, val)
        if item in self.__dict__.keys() and item not in ['nParts', 'parts']:
            # Here is where I should do the sanity checking of the new vals.
            if item == 'brLenPriorType':
                assert val in ['exponential', 'uniform']

            # print "    Setting tuning '%s' to %s" % (item, val)
            object.__setattr__(self, item, val)
        else:
            print(self.dump())
            gm = ["\nMcmcTunings.__setattr__()"]
            gm.append("Can't set tuning '%s'-- no such tuning." % item)
            raise P4Error(gm)

    def reprString(self, advice=True):
        lst = ["\nMcmc.tunings:  nParts=%s" % self.nParts]
        spacer = ' ' * 4

        lst.append("%s%15s: %s" % (spacer, 'chainTemp', self.chainTemp))
        lst.append("%s%15s: %7.5f" % (spacer, 'local', self.local))
        lst.append("%s%15s: %5.3f" % (spacer, 'brLen', self.brLen))
        lst.append("%s%15s: %s" %
                   (spacer, 'brLenPriorType', self.brLenPriorType))
        lst.append("%s%15s: %5.3f" % (spacer, 'etbrPExt', self.etbrPExt))
        lst.append("%s%15s: %5.3f" % (spacer, 'etbrLambda', self.etbrLambda))
        lst.append("%s%15s: %.3f" % (spacer, 'relRate', self.relRate))
        #lst.append("%s%15s: %.3f" % (spacer, 'treeScale', self.treeScale))
        lst.append("%s%15s: %.3f" % (spacer, 'allBrLens', self.allBrLens))

        #lst.append("%s%30s: %s" % (spacer, 'chainTemp', self.chainTemp))
        #lst.append("%s%30s: %5.3f" % (spacer, 'local', self.local))
        #lst.append("%s%30s: %5.3f" % (spacer, 'brLen', self.brLen))
        #lst.append("%s%30s: %s" % (spacer, 'relRate', self.relRate))
        lst.append("%s%15s: %5.3f" %
                   (spacer, 'brLenPriorLambda', self.brLenPriorLambda))
        #lst.append("%s%30s: %s" % (spacer, 'doPolytomyResolutionClassPrior', self.doPolytomyResolutionClassPrior))
        #lst.append("%s%30s: %5.3f" % (spacer, 'polytomyPriorLogBigC', self.polytomyPriorLogBigC))

        lst.append("")

        theSig = "%s%22s"
        aLine = theSig % (spacer, 'part-specific tunings')
        for pNum in range(self.nParts):
            aLine += " %10s" % pNum
        lst.append(aLine)

        aLine = theSig % (spacer, '')
        for pNum in range(self.nParts):
            aLine += " %10s" % '------'
        lst.append(aLine)

        aLine = theSig % (spacer, 'comp')
        for pNum in range(self.nParts):
            aLine += " %10.3f" % self.parts[pNum].comp
        lst.append(aLine)

        aLine = theSig % (spacer, 'compDir')
        for pNum in range(self.nParts):
            aLine += " %10.3f" % self.parts[pNum].compDir
        lst.append(aLine)

        aLine = theSig % (spacer, 'allCompsDir')
        for pNum in range(self.nParts):
            aLine += " %10.3f" % self.parts[pNum].allCompsDir
        lst.append(aLine)

        aLine = theSig % (spacer, 'ndch2_leafCompsDir')
        for pNum in range(self.nParts):
            aLine += " %10.3f" % self.parts[pNum].ndch2_leafCompsDir
        lst.append(aLine)

        aLine = theSig % (spacer, 'ndch2_internalCompsDir')
        for pNum in range(self.nParts):
            aLine += " %10.3f" % self.parts[pNum].ndch2_internalCompsDir
        lst.append(aLine)

        #aLine = theSig % (spacer, 'rjComp')
        #for pNum in range(self.nParts):
        #    aLine += " %10.3f" % self.parts[pNum].rjComp
        #lst.append(aLine)

        aLine = theSig % (spacer, 'rMatrix')
        for pNum in range(self.nParts):
            aLine += " %10.3f" % self.parts[pNum].rMatrix
        lst.append(aLine)

        aLine = theSig % (spacer, 'rMatrixDir')
        for pNum in range(self.nParts):
            aLine += " %10.3f" % self.parts[pNum].rMatrixDir
        lst.append(aLine)

        aLine = theSig % (spacer, 'twoP')
        for pNum in range(self.nParts):
            aLine += " %10.3f" % self.parts[pNum].twoP
        lst.append(aLine)

        #aLine = theSig % (spacer, 'rjRMatrix')
        #for pNum in range(self.nParts):
        #    aLine += " %10.3f" % self.parts[pNum].rjRMatrix
        #lst.append(aLine)

        aLine = theSig % (spacer, 'gdasrv')
        for pNum in range(self.nParts):
            #aLine += " %10s" % self.parts[pNum].gdasrv
            aLine += " %10.3f" % self.parts[pNum].gdasrv
        lst.append(aLine)

        aLine = theSig % (spacer, 'pInvar')
        for pNum in range(self.nParts):
            aLine += " %10s" % self.parts[pNum].pInvar
        lst.append(aLine)

        #aLine = theSig % (spacer, 'gdasrvPriorLambda')
        # for pNum in range(self.nParts):
        #    aLine += " %10.2f" % self.parts[pNum].gdasrvPriorLambda
        # lst.append(aLine)

        aLine = theSig % (spacer, 'compLocation')
        for pNum in range(self.nParts):
            aLine += " %10.2f" % self.parts[pNum].compLocation
        lst.append(aLine)

        aLine = theSig % (spacer, 'rMatrixLocation')
        for pNum in range(self.nParts):
            aLine += " %10.2f" % self.parts[pNum].rMatrixLocation
        lst.append(aLine)

        if advice:
            lst.append("\n  To change these settings, do eg")
            lst.append("    yourMcmc.tunings.chainTemp = 0.15")
            lst.append("  or")
            lst.append("    yourMcmc.tunings.parts[0].comp = 0.5")

        return '\n'.join(lst)

    def dump(self, advice=True):
        print(self.reprString(advice))

    def __repr__(self):
        return self.reprString()

# McmcTuningsUsage is used only by autoTune().  It is used to allow
# autoTune() to be able to get proposals easily.  So for example
# mcmc.tuningsUsage.local is set to a Proposal object for local.


class McmcTuningsUsagePart(object):

    def __init__(self, partNum):
        object.__setattr__(self, 'num', partNum)
        object.__setattr__(self, 'comp', [])
        object.__setattr__(self, 'compDir', [])
        object.__setattr__(self, 'allCompsDir', [])
        object.__setattr__(self, 'ndch2_leafCompsDir', [])
        object.__setattr__(self, 'ndch2_internalCompsDir', [])
        #object.__setattr__(self, 'rjComp', [])
        #object.__setattr__(self, 'addComp', [])
        object.__setattr__(self, 'rMatrix', [])
        object.__setattr__(self, 'rMatrixDir', [])
        object.__setattr__(self, 'gdasrv', [])
        object.__setattr__(self, 'pInvar', None)
        object.__setattr__(self, 'compLocation', None)
        object.__setattr__(self, 'rMatrixLocation', None)

    def __setattr__(self, item, val):
        gm = ["\nMcmcTuningsUsagePart.__setattr__()  Part %i" % self.num]
        gm.append("Can't set-- its not allowed.")
        raise P4Error(gm)


class McmcTuningsUsage(object):

    """This class associates tunings with proposals."""

    def __init__(self, nParts):
        object.__setattr__(self, 'nParts', nParts)
        object.__setattr__(self, 'parts', [])
        for pNum in range(nParts):
            self.parts.append(McmcTuningsUsagePart(pNum))
        #object.__setattr__(self, 'chainTemp', None)
        object.__setattr__(self, 'allBrLens', None)
        object.__setattr__(self, 'brLen', None)
        #object.__setattr__(self, 'eTBR', None)
        object.__setattr__(self, 'local', None)
        #object.__setattr__(self, 'polytomy', None)
        object.__setattr__(self, 'relRate', None)
        # root3 has no tuning
        #object.__setattr__(self, 'treeScale', None)

    def __setattr__(self, item, val):
        gm = ["\nMcmcTuningsUsage.__setattr__()"]
        gm.append("Can't set-- its not allowed.")
        raise P4Error(gm)

    def reprString(self):
        nTunings = 0
        lst = ["\nMcmc.tuningsUsage:  nParts=%s" % self.nParts]
        lst.append("Number of proposals used with the various tunings.")
        spacer = ' ' * 4
        #spacer2 = ' ' * 22
        #lst.append("%s%15s: %s" % (spacer, 'chainTemp', self.chainTemp))

        if self.brLen:
            nTunings += 1
            lst.append("%s%15s: 1" % (spacer, 'brLen'))
        else:
            lst.append("%s%15s: None" % (spacer, 'brLen'))

        if self.allBrLens:
            nTunings += 1
            lst.append("%s%15s: 1" % (spacer, 'allBrLens'))
        else:
            lst.append("%s%15s: None" % (spacer, 'allBrLens'))

        if self.local:
            nTunings += 1
            lst.append("%s%15s: 1" % (spacer, 'local'))
        else:
            lst.append("%s%15s: None" % (spacer, 'local'))

        # if self.polytomy:
        #    nTunings += 1
        #    lst.append("%s%15s: 1" % (spacer, 'polytomy'))
        # else:
        #    lst.append("%s%15s: None" % (spacer, 'polytomy'))

        if self.relRate:
            nTunings += 1
            lst.append("%s%15s: 1" % (spacer, 'relRate'))
        else:
            lst.append("%s%15s: None" % (spacer, 'relRate'))

        # if self.treeScale:
        #     nTunings += 1
        #     lst.append("%s%15s: 1" % (spacer, 'treeScale'))
        # else:
        #     lst.append("%s%15s: None" % (spacer, 'treeScale'))
        lst.append("")

        theSig = "%s%15s"
        sig2 = "%s%22s: %i"
        lst.append(theSig % (spacer, 'part-specific tunings'))
        lst.append(theSig % (spacer, '---------------------'))

        for pNum in range(self.nParts):
            lst.append(theSig % (spacer, 'part %i    ' % pNum))
            lst.append(sig2 % (spacer, 'comp', len(self.parts[pNum].comp)))
            lst.append(sig2 %
                       (spacer, 'compDir', len(self.parts[pNum].compDir)))
            lst.append(sig2 %
                       (spacer, 'allCompsDir', len(self.parts[pNum].allCompsDir)))
            lst.append(sig2 %
                       (spacer, 'ndch2_leafCompsDir', len(self.parts[pNum].ndch2_leafCompsDir)))
            lst.append(sig2 %
                       (spacer, 'ndch2_internalCompsDir', len(self.parts[pNum].ndch2_internalCompsDir)))
            lst.append(sig2 %
                       (spacer, 'rMatrix', len(self.parts[pNum].rMatrix)))
            lst.append(sig2 %
                       (spacer, 'rMatrixDir', len(self.parts[pNum].rMatrixDir)))
            lst.append(sig2 % (spacer, 'gdasrv', len(self.parts[pNum].gdasrv)))
            if len(self.parts[pNum].comp):
                nTunings += 1
            if len(self.parts[pNum].rMatrix):
                nTunings += 1
            if len(self.parts[pNum].gdasrv):
                nTunings += 1
            if self.parts[pNum].pInvar:
                lst.append(sig2 % (spacer, 'pInvar', 1))
                nTunings += 1
            else:
                lst.append(sig2 % (spacer, 'pInvar', 0))

            # if self.parts[pNum].compLocation:
            if len(self.parts[pNum].comp) > 1:
                lst.append(sig2 % (spacer, 'compLocation', 1))
                nTunings += 1
            else:
                lst.append(sig2 % (spacer, 'compLocation', 0))

            # if self.parts[pNum].rMatrixLocation:
            if len(self.parts[pNum].rMatrix) > 1:
                lst.append(sig2 % (spacer, 'rMatrixLocation', 1))
                nTunings += 1
            else:
                lst.append(sig2 % (spacer, 'rMatrixLocation', 0))

        lst.append(
            "There are %i tunings for the model (not including the chainTemp tuning, if it exists)." % nTunings)

        return '\n'.join(lst)

    def dump(self):
        print(self.reprString())

    def __repr__(self):
        return self.reprString()


class McmcProposalProbs(dict):

    """User-settable relative proposal probabilities.

    An instance of this class is made as Mcmc.prob, where you can
    do, for example::

        yourMcmc.prob.local = 2.0

    These are relative proposal probs, that do not sum to 1.0, and
    affect the calculation of the final proposal probabilities (ie the
    kind that do sum to 1).  It is a relative setting, and the default
    is 1.0.  Setting it to 0 turns it off.  For small
    probabilities, setting it to 2.0 doubles it.  For bigger
    probabilities, setting it to 2.0 makes it somewhat bigger.

    Check the effect that it has by doing::

        yourMcmc.writeProposalIntendedProbs()
    which prints out the final calculated probabilities. 
    """

    def __init__(self):
        object.__setattr__(self, 'comp', 1.0)
        object.__setattr__(self, 'compDir', 0.0)
        object.__setattr__(self, 'allCompsDir', 0.0)
        object.__setattr__(self, 'ndch2_leafCompsDir', 0.0)
        object.__setattr__(self, 'ndch2_internalCompsDir', 0.0)
        object.__setattr__(self, 'ndch2_leafCompsDirAlpha', 0.0)
        object.__setattr__(self, 'ndch2_internalCompsDirAlpha', 0.0)
        #object.__setattr__(self, 'rjComp', 0.0)
        object.__setattr__(self, 'rMatrix', 1.0)
        object.__setattr__(self, 'rMatrixDir', 0.0)
        #object.__setattr__(self, 'rjRMatrix', 0.0)
        object.__setattr__(self, 'gdasrv', 1.0)
        object.__setattr__(self, 'pInvar', 1.0)
        object.__setattr__(self, 'local', 1.0)
        object.__setattr__(self, 'brLen', 0.0)
        object.__setattr__(self, 'allBrLens', 0.0)
        object.__setattr__(self, 'eTBR', 1.0)
        #object.__setattr__(self, 'treeScale', 0.0)
        object.__setattr__(self, 'polytomy', 0.0)
        object.__setattr__(self, 'root3', 0.0)
        object.__setattr__(self, 'compLocation', 0.0)
        object.__setattr__(self, 'rMatrixLocation', 0.0)
        object.__setattr__(self, 'gdasrvLocation', 0.0)
        object.__setattr__(self, 'relRate', 1.0)
        #object.__setattr__(self, 'cmd1_compDir', 0.0)
        #object.__setattr__(self, 'cmd1_comp0Dir', 0.0)
        #object.__setattr__(self, 'cmd1_allCompDir', 0.0)
        #object.__setattr__(self, 'cmd1_alpha', 0.0)

    def __setattr__(self, item, val):
        # complaintHead = "\nMcmcProposalProbs.__setattr__()"
        gm = ["\nMcmcProposalProbs(). (set %s to %s)" % (item, val)]
        theKeys = self.__dict__.keys()
        if item in theKeys:
            try:
                val = float(val)
                if val < 1e-9:
                    val = 0
                object.__setattr__(self, item, val)
            except:
                gm.append("Should be a float.  Got '%s'" % val)
                raise P4Error(gm)

        else:
            self.dump()
            gm.append("Can't set '%s'-- no such proposal." % item)
            gm.append("See the list of proposal probabilities above for valid proposals.")
            raise P4Error(gm)

    def reprString(self):
        stuff = [
            "\nUser-settable relative proposal probabilities, from yourMcmc.prob"]
        stuff.append("  To change it, do eg ")
        stuff.append("    yourMcmc.prob.comp = 0.0 # turns comp proposals off")
        stuff.append("  Current settings:")
        theKeys = self.__dict__.keys()
        theKeys.sort()
        for k in theKeys:
            stuff.append("        %30s: %s" % (k, getattr(self, k)))
        return '\n'.join(stuff)

    def dump(self):
        print(self.reprString())

    def __repr__(self):
        return self.reprString()


class Proposal(object):

    def __init__(self, theMcmc=None):
        self.name = None
        self.variant = 'gtr'  # only for rMatrix.  2p or gtr
        self.mcmc = theMcmc
        self.nChains = theMcmc.nChains
        self.pNum = -1
        self.mtNum = -1
        self.weight = 1.0
        #self.tuning = None
        self.nProposals = [0] * self.nChains
        self.nAcceptances = [0] * self.nChains
        self.accepted = 0
        self.topologyChanged = 0
        self.nTopologyChangeAttempts = [0] * self.nChains
        self.nTopologyChanges = [0] * self.nChains
        self.doAbort = False
        self.nAborts = [0] * self.nChains

    def dump(self):
        print("proposal name=%-10s pNum=%2i, mtNum=%2i, weight=%5.1f, tuning=%7.2f" % (
            '%s,' % self.name, self.pNum, self.mtNum, self.weight, self.tuning))
        print("    nProposals   by temperature:  %s" % self.nProposals)
        print("    nAcceptances by temperature:  %s" % self.nAcceptances)

    # Some tunings are part-specific, and so are associated with the proposals.
    def _getTuning(self):
        if self.name in ['relRate', 'local', 'brLen', 'treeScale', 'allBrLens']:
            # print "getting tuning for %s, returning %f" % (self.name,
            # getattr(self.mcmc.tunings, self.name))
            return getattr(self.mcmc.tunings, self.name)
        elif self.name in ['eTBR']:
            return getattr(self.mcmc.tunings, 'etbrLambda')
        elif self.name in ['comp', 'compDir', 'allCompsDir', 
                           'ndch2_leafCompsDir', 'ndch2_internalCompsDir',
                           'rjComp', 'rMatrix', 
                           'rMatrixDir', 'rjRMatrix', 'gdasrv', 'pInvar']:
            # print "getting tuning for %s, partNum %i, returning %f" % (
            #    self.name, self.pNum, getattr(self.mcmc.tunings.parts[self.pNum], self.name))
            # the variant attribute is new, and can mess up reading older
            # pickles.
            if self.name in ['rMatrix', 'rMatrixDir'] and self.variant == '2p':
                return getattr(self.mcmc.tunings.parts[self.pNum], 'twoP')
            else:
                return getattr(self.mcmc.tunings.parts[self.pNum], self.name)
        elif self.name in ['compLocation', 'rMatrixLocation']:
            # print "getting tuning for %s, partNum %i, returning %f" % (
            # self.name, self.pNum, getattr(self.mcmc.tunings.parts[self.pNum],
            # self.name))
            if hasattr(self.mcmc.tunings.parts[self.pNum], self.name):
                return getattr(self.mcmc.tunings.parts[self.pNum], self.name)
            else:
                return None
        else:
            return None

    def _setTuning(self, whatever):
        raise P4Error("Can't set tuning this way.")

    def _delTuning(self):
        raise P4Error("Can't del tuning.")

    tuning = property(_getTuning, _setTuning, _delTuning)

class SwapTuner(object):
    """Continuous tuning for swap temperature"""

    def __init__(self, sampleSize):
        self.sampleSize = sampleSize
        self.swaps01 = []
        self.accHi = 0.1    # 10% acceptance
        self.factorHi = 1.333  # Factor if acceptance is > accHi
        self.accLo = 0.02
        self.accB = 0.005
        self.factorB = 1.333   # if acceptance is > accB, but < accLo
        self.factorC = 1.5     # if acceptance is < accB

class Mcmc(object):

    """An MCMC for molecular sequences.

    aTree
                  The tree should have a model and data attached.

    nChains
                  The number of chains in the MCMCMC, default 4

    runNum
                  You may want to do more than one 'run' in the same
                  directory, to facilitate convergence testing
                  (another idea stolen from MrBayes, so thanks to the
                  authors).  The first runNum would be 0, and samples,
                  likelihoods, and checkPoints are written to files
                  with that number.

    sampleInterval
                  Interval at which the cold chain is sampled,
                  including writing a tree, the logLike, and perhaps
                  doing a simulation.

    checkPointInterval
                  Intervals at which the MCMC is checkpointed,
                  meaning that the whole thing is written to a pickle
                  file.  You can re-start from a checkpoint, eg in the
                  event of a crash, or if you just want to make the
                  MCMC run longer.  You can turn off checkPointing by
                  setting it to zero or None.

    simulate
                  a 'binary' number from 1-31 that says what posterior
                  predictive simulation test quantities to collect.
                  This week, we have

                  1.  multinomial log likelihood

                  2.  bigXSquared

                  4.  meanNCharsPerSite

                  8.  c_m, the compStatFromCharFreqs (2 numbers per part--sim and curTree)

                  16. constant sites count (not proportion)

                  So if you say simulate=1, you get multinomial log
                  like.  If you say simulate=2 you get bigXSquared.
                  If you say simulate=6 you get bigXSquared and
                  meanNCharsPerSite.  If you say simulate=31 you get
                  all 5.

    writePrams
                  write changeable model parameters to a file.

    constraints
                  If you want to constrain the topology, say so with
                  a Constraints object here.  A constraints object may
                  enforce multiple and nested constraints.

    To prepare for a run, instantiate an Mcmc object::

        m = Mcmc(t, sampleInterval=1000, checkpointInterval=200000)

    To start it running, do this::

        m.run(1000000) # Tell it the number of generations to do

    You will want to make the last generation land on a
    *checkPointInterval*.  So in this case 1000000 is ok, but 900000 is
    not.

    As it runs, it saves trees and likelihoods at sampleInterval
    intervals (actually whenever the current generation number is
    evenly divisible by the sampleInterval).

    Whenever the current generation number is evenly divisible by the
    checkPointInterval it will write a checkPoint file.  A checkPoint
    file is the whole MCMC without the data, and without C-structs.
    Using a checkPoint, you can re-start an Mcmc from the point you
    left off.  Or, in the event of a crash, you can restart from the
    latest checkPoint.  You can query checkPoints to get information
    about how the chain has been running, and about convergence
    diagnostics.

    In order to restart the MCMC from the end of a previous run, you
    need the data::

        read('yourData.nex')
        d = Data()
        # read the last checkPoint file
        m = func.unPickleMcmc(0, d)  # runNum 0
        m.run(20000)

    Its that easy if your previous run finished properly.  However, if
    your previous run has crashed and you want to restart it from a
    checkPoint, then you will need to repair the sample output files
    to remove samples that were taken after the last checkPoint, but
    before the crash.  Fix the trees, likelihoods, prams, and sims.
    (You probably do not need to beware of confusing gen (eg 9999) and
    gen+1 (eg 10000) issues.)  When you remove trees from the tree
    files to restart, each tree file should have the word ``end;`` at
    the end after the last tree, with a single return after that.
    Other files (prams and so on) would end more simply with a single
    return.

    The checkPoints can help with convergence testing.  To help with
    that, you can use the McmcCheckPointReader class.  It will print
    out a table of average standard deviations of split supports
    between 2 runs, or between 2 checkPoints from the same run.  It
    will print out tables of proposal acceptances to show whether they
    change over the course of the MCMC.

    Another way to monitor convergence of split supports is to use the
    Trees method trackSplitSupports().  That method requires a
    reference tree to tell it the splits that it should track; an
    obvious reference tree to use is the consensus tree, but it need
    not be; and the reference tree need not be fully resolved.  Here
    is how it might be used::

        # read in all 1000 trees in the file
        read('mcmc_trees_0.nex')
        tt = Trees()
        # Make a cons tree only from the last 500 trees
        tt2 = Trees(trees=var.trees[500:])
        tp = TreePartitions(tt2)
        t = tp.consensus()
        # Track the splits in the cons tree
        tt.trackSplitsFromTree(t)

    This makes a verbose output, and also writes the numbers to a file
    in case you want to use them later.

    To make a consensus from the trees from more than one run, you can
    add trees to an existing TreePartitions object, like this::

        tp = TreePartitions('mcmc_trees_0.nex', skip=500)
        tp.read('mcmc_trees_1.nex', skip=500)   # add trees from a second run
        t = tp.consensus()
        t.reRoot(t.node('OutgroupTaxon').parent)
        for n in t.iterInternalsNoRoot():
            n.name = '%.0f' % (100. * n.br.support)
        t.writeNexus('cons.nex')
    """

    def __init__(self, aTree, nChains=4, runNum=0, sampleInterval=100, checkPointInterval=10000, simulate=None, writePrams=True, constraints=None, verbose=True, tuningsFileName=None, swapTuner=150):
        gm = ['Mcmc.__init__()']

        self.verbose = verbose

        if aTree and aTree.model and aTree.data:
            pass
        else:
            gm.append(
                "The tree that you feed to this class should have a model and data attached.")
            raise P4Error(gm)

        if 1:
            if 1:
                if aTree.root.getNChildren() != 3 or not aTree.isFullyBifurcating():
                    gm.append(
                        "Mcmc is not implemented for bifurcating roots, or trees that are not fully bifurcating.")
                    raise P4Error(gm)
            else:
                if aTree.root.getNChildren() < 3:
                    gm.append(
                        "Mcmc is not implemented for roots that have less than 3 children.")
                    raise P4Error(gm)

        if 0:  # Muck with polytomies
            nn = [n for n in aTree.iterInternalsNoRoot()]
            for n in nn:
                r = random.random()
                if r < 0.4:
                    aTree.collapseNode(n)
            aTree.draw()
            # sys.exit()

        self.tree = aTree
        self.constraints = constraints
        if self.constraints:
            self.tree.makeSplitKeys()
            assert isinstance(self.constraints, Constraints)
            mySplitKeys = []
            for n in self.tree.iterInternalsNoRoot():
                mySplitKeys.append(n.br.splitKey)
            # print "input tree splitKeys = %s" % mySplitKeys
            for sk in self.constraints.constraints:
                if sk not in mySplitKeys:
                    # self.tree.draw()
                    gm.append('Constraint %s' %
                              p4.func.getSplitStringFromKey(sk, self.tree.nTax))
                    gm.append('is not in the starting tree.')
                    gm.append(
                        'Maybe you want to make a randomTree with constraints?')
                    raise P4Error(gm)

        try:
            nChains = int(nChains)
        except (ValueError, TypeError):
            gm.append("nChains should be an int, 1 or more.  Got %s" % nChains)
            raise P4Error(gm)
        if nChains < 1:
            gm.append("nChains should be an int, 1 or more.  Got %s" % nChains)
            raise P4Error(gm)
        self.nChains = nChains
        self.chains = []
        self.gen = -1
        self.startMinusOne = -1

        for n in self.tree.iterNodesNoRoot():
            if n.br.len < var.BRLEN_MIN:
                gm.append("node %i brlen (%g)is too short." %
                          (n.nodeNum, n.br.len))
                raise P4Error(gm)
            elif n.br.len > var.BRLEN_MAX:
                gm.append("node %i brlen (%f)is too long." %
                          (n.nodeNum, n.br.len))
                raise P4Error(gm)

        try:
            runNum = int(runNum)
        except (ValueError, TypeError):
            gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
            raise P4Error(gm)
        if runNum < 0:
            gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
            raise P4Error(gm)
        self.runNum = runNum

        # Check that we are not going to over-write good stuff
        ff = os.listdir(os.getcwd())
        hasPickle = False
        for fName in ff:
            if fName.startswith("mcmc_checkPoint_%i." % self.runNum):
                hasPickle = True
                break
        if hasPickle:
            gm.append("runNum is set to %i" % self.runNum)
            gm.append(
                "There is at least one mcmc_checkPoint_%i.xxx file in this directory." % self.runNum)
            gm.append(
                "This is a new Mcmc, and I am refusing to over-write exisiting files.")
            gm.append(
                "Maybe you want to re-start from the latest mcmc_checkPoint_%i file?" % self.runNum)
            gm.append(
                "Otherwise, get rid of the existing mcmc_xxx_%i.xxx files and start again." % self.runNum)
            raise P4Error(gm)

        if var.strictRunNumberChecking:
            # We want to start runs with number 0, so if runNum is more than
            # that, check that there are other runs.
            if self.runNum > 0:
                for runNum2 in range(self.runNum):
                    hasTrees = False
                    for fName in ff:
                        if fName.startswith("mcmc_trees_%i" % runNum2):
                            hasTrees = True
                            break
                    if not hasTrees:
                        gm.append("runNum is set to %i" % self.runNum)
                        gm.append("runNums should go from zero up.")
                        gm.append("There are no mcmc_trees_%i.nex files to show that run %i has been done." % (
                            runNum2, runNum2))
                        gm.append("Set the runNum to that, first.")
                        gm.append(
                            "(To get rid of this requirement, turn off var.strictRunNumberChecking.)")
                        raise P4Error(gm)

        self.sampleInterval = sampleInterval
        self.checkPointInterval = checkPointInterval

        self.proposals = []
        self.proposalsHash = {}
        self.propWeights = []
        self.cumPropWeights = []
        self.totalPropWeights = 0.0

        self.treePartitions = None
        self.likesFileName = "mcmc_likes_%i" % runNum
        self.treeFileName = "mcmc_trees_%i.nex" % runNum
        self.simFileName = "mcmc_sims_%i" % runNum
        self.pramsFileName = "mcmc_prams_%i" % runNum
        self.hypersFileName = "mcmc_hypers_%i" % runNum
        #self.rjKFileName = "mcmc_rjK_%i" % runNum
        self.writePrams = writePrams
        self.writeHypers = True

        self.lastTimeCheck = None

        if simulate:
            try:
                simulate = int(simulate)
            except (ValueError, TypeError):
                gm.append("Arg 'simulate' should be an int, 1-31, inclusive.")
                raise P4Error(gm)
            if simulate <= 0 or simulate > 31:
                gm.append("Arg 'simulate' should be an int, 1-31, inclusive.")
                raise P4Error(gm)
        self.simulate = simulate
        if self.simulate:
            self.simTree = self.tree.dupe()
            self.simTree.data = self.tree.data.dupe()
            self.simTree.calcLogLike(verbose=False)
        else:
            self.simTree = None

        if self.nChains > 1:
            self.swapMatrix = []
            for i in range(self.nChains):
                self.swapMatrix.append([0] * self.nChains)
            if swapTuner:             # a kwarg
                myST = int(swapTuner)
                if myST >= 100:
                    self.swapTuner = SwapTuner(myST)
                else:
                    gm.append("The swapTuner kwarg, the sample size, should be at least 100.  Got %i." % myST)
                    raise P4Error(gm)
            else:
                self.swapTuner = None
        else:
            self.swapMatrix = None
            self.swapTuner = None

        # check the tree
        aTree.calcLogLike(verbose=False)

        if 0:
            # print complaintHead
            print("    logLike of the input tree is %s" % aTree.logLike)

        if not aTree.taxNames:
            gm.append(
                "The tree that you supply should have a 'taxNames' attribute.")
            gm.append("The taxNames should be in the same order as the data.")
            raise P4Error(gm)

        if tuningsFileName:
            if verbose:
                print("Reading tunings from file '%s' ..." % tuningsFileName)
            tf = open(tuningsFileName, 'rb')
            self.tunings = pickle.load(tf)
            tf.close()
            #self.tunings.dump()
        else:
            self.tunings = McmcTunings(self.tree.model.nParts)

        # tuningsUsage is only used by autoTune()
        self.tuningsUsage = McmcTuningsUsage(self.tree.model.nParts)
        self.prob = McmcProposalProbs()

        # Tuning hack --- for old dirichlet proposals

        # Notes on default tunings.  Test with sims.  (XX% acceptance shown)

        # Test with DNA, hetero.
        #   comp tuning of 300 was good (20, 21%, 20, 25, 20, 22).
        #   rMatrix tuning of 10 was too small (36, 46%, 39, 40%, 40, 39).  20 (60, 53).
        #                                       5.0 (17, 18, 28, 26, 36, 45, 15,17)
        #   gdasrv tuning of 2. was good.  (14%, 18%, 20, 7)   3 (16,20, 22, 20)
        #   pInvar tuning of 0.1 was too small.  0.2 too small.  0.5 was ok (40% 28)  0.9 was ok (35% acceptance)
        # Test with protein, mildly hetero comp
        #   comp tuning of 300 too small (1% acceptance)  5000 (46%)  2000. (26%, 27%)  1000. (11%, 14%, 10%, 11%)
        #   gdasrv tuning of 2 was ok (15%, 7%, 9%, 15%, 20%, 23%)   4.0 (26%)  3.0 (20%, 23%)
        #   pInvar tuning of 0.1 was too small (90%), 0.5 (55%, 40%, 50%, 45%)
        # Test with grouped aa's, simulated via protein
        #   comp  300 (8%, 7%)   500 (14, 17, 15,16, 11, 16)  1000 (25, 29%)  700 (24,29, 17,23)
        #   rmatrix  10 (1%)   500 (70%)   50 (33%)   40 (20, 26, 30, 34, 22, 15, 23)
        #   gdasrv   2. (8)     3 (10, 11, 10, 15, 0, 16, 4,6)    4 (9, 15, 25)
        #   pInvar   0.1 (75%, 50%)  0.5 (11., 33, 25, 30, 18,30, 22, 24)

        # The default tunings for comp (300) and rMatrix (300) are
        # appropriate defaults for DNA data, but not good for protein
        # or recoded aa data.  If that is the case, change it.
        #  This is old -- should be checked and then probably deleted.
        if 0:
            for pNum in range(self.tunings.nParts):
                if self.tree.model.parts[pNum].dim == 20:
                    self.tunings.parts[pNum].comp = 2000.
                elif self.tree.model.parts[pNum].dim == 6:
                    self.tunings.parts[pNum].comp = 700.
                    self.tunings.parts[pNum].rMatrix = 1000.
                elif self.tree.model.parts[pNum].dim == 4:  # DNA
                    if self.tree.model.parts[pNum].rMatrices[0].spec == '2p':
                        self.tunings.parts[pNum].rMatrix = 50
                        # print "setting 2p tunings to 50"

        # New tunings, as of May 2010.
        # For real data, with grouped aa's, comp 0.074, rMatrix 0.067
        # Simulated DNA -- 2 comps -- 0.25
        #               -- 2 rmatrices -- 0.167
        # Simulated protein, homog, comp tuning 0.05 -> acceptance 6.7%
        #                           rMatrix tuning 0.005 -> acceptance 31.1%

        # Even newer tunings, as of January 2011.  The sliders for
        # comp and rMatrix now go from approx 0 to approx 1, and so do
        # not depend on the dim (I think).  So leave them at their
        # default.
        if 0:
            for pNum in range(self.tunings.nParts):
                theDim = self.tree.model.parts[pNum].dim
                nRates = ((theDim * theDim) - theDim) / 2
                self.tunings.parts[pNum].comp = 1.0 / theDim
                self.tunings.parts[pNum].rMatrix = 1.0 / nRates

        # Two new tunings --- compDir and rMatrixDir, seem to depend on the dim.
        # And now allCompsDir
        if 1 and not tuningsFileName:
            for pNum in range(self.tunings.nParts):
                theDim = self.tree.model.parts[pNum].dim
                nRates = ((theDim * theDim) - theDim) / 2
                self.tunings.parts[pNum].compDir = 50. * theDim
                self.tunings.parts[pNum].allCompsDir = 100. * theDim
                self.tunings.parts[pNum].ndch2_leafCompsDir = 2000. * theDim
                self.tunings.parts[pNum].ndch2_internalCompsDir = 500. * theDim
                self.tunings.parts[pNum].rMatrixDir = 50. * nRates
            


        # Zap internal node names
        for n in aTree.root.iterInternals():
            if n.name:
                n.name = None

        if self.tree.model.nParts > 1 and self.tree.model.relRatesAreFree:
            self.prob.relRate = 1.0
            if self.verbose:
                print("\nInitiating across-data heterogeneous model...")
                print("\n%23s" % "Additional proposals:")
                print("     relative partition rate = on")
                print("\n  %s" % "[You can turn it off by setting")
                print("  %s" % "yourMcmc.prob.relRate=0.0]")
        else:
            self.prob.relRate = 0.0

        nNodes = len(self.tree.nodes)
        for pNum in range(self.tree.model.nParts):
            mp = self.tree.model.parts[pNum]
            if mp.ndch2:
                if mp.nComps != nNodes:
                    gm.append("Model part %i, ndch2 is on, nNodes is %i, nComps is %i" % (
                        pNum, nNodes, mp.nComps))
                    gm.append("For ndch2 there should be one comp for each node.")
                    raise P4Error(gm)
                if mp.nRMatrices > 1:
                    gm.append("Model part %i, ndch2 is on, nNodes is %i, nRMatrices is %i" % (
                        pNum, nNodes, mp.nRMatrices))
                    gm.append("This week, for ndch2 there should be only one rMatrix.")
                    raise P4Error(gm)
                if mp.nGdasrvs > 1:
                    gm.append("Model part %i, ndch2 is on, nNodes is %i, nGdasrvs is %i" % (
                        pNum, nNodes, mp.nGdasrvs))
                    gm.append("This week, for ndch2 there should be only one gdasrv.")
                    raise P4Error(gm)

                mp.ndch2_globalComp = numpy.array(self.tree.data.parts[pNum].composition())
                while mp.ndch2_globalComp.min()  < var.PIVEC_MIN:
                    for i in range(mp.dim):
                        if mp.ndch2_globalComp[i] < var.PIVEC_MIN:
                            mp.ndch2_globalComp[i] += (1.0 + random.random()) * var.PIVEC_MIN
                    thisSum = mp.ndch2_globalComp.sum()
                    for i in range(mp.dim):
                        mp.ndch2_globalComp[i] /= thisSum

               

        if self.tree.model.isHet:
            props_on = []
            if verbose:
                print("\nInitiating across-tree heterogeneous model...")

            self.tree.setModelThingsNNodes()

            if self.verbose:
                self.tree.summarizeModelThingsNNodes()

            for pNum in range(self.tree.model.nParts):
                mp = self.tree.model.parts[pNum]

                if mp.nComps > 1 and mp.nComps < nNodes:
                    self.prob.compLocation = 1.0
                    if "composition location" not in props_on:
                        props_on.append("composition location")

                if mp.nRMatrices > 1 and mp.nRMatrices < nNodes:
                    self.prob.rMatrixLocation = 1.0
                    if "rate matrix location" not in props_on:
                        props_on.append("rate matrix location")

                if mp.nGdasrvs > 1 and mp.nGdasrvs < nNodes:
                    self.prob.gdasrvLocation = 1.0
                    if "gamma rate location" not in props_on:
                        props_on.append("gamma rate location")

                if mp.ndch2:
                    self.prob.ndch2_leafCompsDir = 1.0
                    thisMString = "ndch2_leafCompsDir"
                    if thisMString not in props_on:
                        props_on.append(thisMString)

                    self.prob.ndch2_internalCompsDir = 1.0
                    thisMString = "ndch2_internalCompsDir"
                    if thisMString not in props_on:
                        props_on.append(thisMString)

                    self.prob.ndch2_leafCompsDirAlpha = 1.0
                    thisMString = "ndch2_leafCompsDirAlpha"
                    if thisMString not in props_on:
                        props_on.append(thisMString)

                    self.prob.ndch2_internalCompsDirAlpha = 1.0
                    thisMString = "ndch2_internalCompsDirAlpha"
                    if thisMString not in props_on:
                        props_on.append(thisMString)

            self.prob.root3 = 1.0

            if verbose:
                props_on.append("root location")
                print("\n%23s" % "Additional proposals:")
                for prop in props_on:
                    print("%30s = on" % prop)
                print("\n  %s" % "[You can of course turn them")
                print("  %s" % "off again before Mcmc.run()]\n")
        else:
            self.prob.root3 = 0.0
            self.prob.compLocation = 0.0
            self.prob.rMatrixLocation = 0.0
            self.prob.gdasrvLocation = 0.0



        # Hidden experimental hacking
        self.doHeatingHack = False
        self.heatingHackTemperature = 5.0
        #self.heatingHackProposalNames = ['local', 'eTBR']
        

        # # Are we using rjComp in any model partitions?
        # # True and False
        # rjCompParts = [mp.rjComp for mp in self.tree.model.parts]
        # rjCompPartNums = [
        #     pNum for pNum in range(self.tree.model.nParts) if rjCompParts[pNum]]
        # # print rjCompParts
        # # print rjCompPartNums
        # if rjCompPartNums:
        #     if verbose:
        #         print("\nInitiating rjComp...")
        #         print("Turning on proposal for rjComp.")

        #     for pNum in rjCompPartNums:
        #         mp = self.tree.model.parts[pNum]
        #         # Do some checking
        #         if mp.nComps <= 1:
        #             gm.append("rjComp is turned on for part %i, but there are %i comps.  Too few." % (
        #                 pNum, mp.nComps))
        #             gm.append("You will want to add more comps.")
        #             raise P4Error(gm)
        #         elif mp.nComps > len(self.tree.nodes):
        #             gm.append("rjComp is turned on for part %i, but there are %i comps.  Too many." % (
        #                 pNum, mp.nComps))
        #             gm.append(
        #                 "That is more than the number of nodes in the tree.")
        #             raise P4Error(gm)

        #         self.prob.rjComp = 1.0

        #         # Calculate the initial pool size, setting rjComp_k
        #         mp.rjComp_k = 0
        #         for comp in mp.comps:
        #             if comp.nNodes:
        #                 mp.rjComp_k += 1
        #                 comp.rj_isInPool = True   # False by default

        #         if verbose:
        #             print("Part %i: nComps %i, pool size (rjComp_k) %i" % (pNum, mp.nComps, mp.rjComp_k))

        #         # This stuff below applies when
        #         # rjCompUniformAllocationPrior is not on, meaning it
        #         # uses the hierarchical allocation prior.

        #         # This next calc of rj_f depends on
        #         # self.tree.setModelThingsNNodes() having been done, above.
        #         mySum = float(len(self.tree.nodes))
        #         for comp in mp.comps:
        #             if comp.nNodes:
        #                 comp.rj_f = comp.nNodes / mySum
        #         # for comp in mp.comps:
        #         #    print comp.num, comp.nNodes, comp.rj_f

        # # Are we using rjRMatrix in any model partitions?
        # # True and False
        # rjRMatrixParts = [mp.rjRMatrix for mp in self.tree.model.parts]
        # rjRMatrixPartNums = [
        #     pNum for pNum in range(self.tree.model.nParts) if rjRMatrixParts[pNum]]
        # # print rjRMatrixParts
        # # print rjRMatrixPartNums
        # if rjRMatrixPartNums:
        #     if verbose:
        #         print("\nInitiating rjRMatrix...")
        #         print("Turning on proposal for rjRMatrix.")

        #     for pNum in rjRMatrixPartNums:
        #         mp = self.tree.model.parts[pNum]
        #         # Do some checking
        #         if mp.nRMatrices <= 1:
        #             gm.append("rjRMatrix is turned on for part %i, but there are %i rMatrices.  Too few." % (
        #                 pNum, mp.nRMatrices))
        #             gm.append("You will want to add more rMatrices.")
        #             raise P4Error(gm)
        #         elif mp.nRMatrices > (len(self.tree.nodes) - 1):
        #             gm.append("rjRMatrix is turned on for part %i, but there are %i rMatrices.  Too many." % (
        #                 pNum, mp.nRMatrices))
        #             gm.append(
        #                 "That is more than the number of branches in the tree.")
        #             raise P4Error(gm)

        #         self.prob.rjRMatrix = 1.0

        #         # Calculate the initial pool size, setting rjRMatrix_k
        #         mp.rjRMatrix_k = 0
        #         for rMatrix in mp.rMatrices:
        #             if rMatrix.nNodes:
        #                 mp.rjRMatrix_k += 1
        #                 rMatrix.rj_isInPool = True   # False by default

        #         if verbose:
        #             print("Part %i: nRMatrices %i, pool size (rjRMatrix_k) %i" % (pNum, mp.nRMatrices, mp.rjRMatrix_k))

        #         # This stuff below applies when
        #         # rjRMatrixUniformAllocationPrior is not on, meaning it
        #         # uses the hierarchical allocation prior.

        #         # This next calc of rj_f depends on
        #         # self.tree.setModelThingsNNodes() having been done, above.
        #         mySum = float(len(self.tree.nodes) - 1.)
        #         for rMatrix in mp.rMatrices:
        #             if rMatrix.nNodes:
        #                 rMatrix.rj_f = rMatrix.nNodes / mySum
        #         # for rMatrix in mp.rMatrices:
        #         #    print rMatrix.num, rMatrix.nNodes, rMatrix.rj_f


        # # Are we doing cmd1 in any model partitions?
        # cmd1Parts = [mp.cmd1 for mp in self.tree.model.parts]  # True and False
        # # empty if there are none that do cmd1
        # cmd1PartNums = [
        #     pNum for pNum in range(self.tree.model.nParts) if cmd1Parts[pNum]]
        # if cmd1PartNums:
        #     if verbose:
        #         print("\nInitiating cmd1 ...")
        #         print("Turning on proposals for cmd1")
        #     self.prob.cmd1_compDir = 1.0
        #     self.prob.cmd1_comp0Dir = 1.0
        #     self.prob.cmd1_allCompDir = 1.0
        #     self.prob.cmd1_alpha = 1.0
        #     for pNum in cmd1PartNums:
        #         mp = self.tree.model.parts[pNum]
        #         mp.cmd1_pi0 = [1.0 / mp.dim] * mp.dim

                

    def _makeProposals(self):
        """Make proposals for the mcmc."""

        gm = ['Mcmc._makeProposals()']

        # The weight determines how often proposals are made.  The
        # weight is the product of 3 things:
        #
        #    1.  User-settable self.prob.something.  Its going to be
        #        1.0 by default, if it is on.
        #
        #    2.  The inherent complexity of the proposal.  Eg if there
        #        is only one parameter, like pInvar, that would have a
        #        complexity of 1, while a DNA composition proposal
        #        would have a complexity of 3, and a protein
        #        composition would have a complexity of 19.
        #
        #    3.  A fudge factor.  Arbitrary witchcraft.
        #
        # So for example, for the 'local' proposal, the inherent
        # complexity is the number of nodes in the tree.  That is the
        # same complexity as 'root3', but I don't think that root3
        # needs to be proposed nearly as often as local, so I upweight
        # local using a fudgeFactor.  So the weight will be
        # p.weight = self.prob.local * len(self.tree.nodes) * fudgeFactor['local']

        # brLen
        if self.prob.brLen:
            p = Proposal(self)
            p.name = 'brLen'
            p.weight = self.prob.brLen * \
                (len(self.tree.nodes) - 1) * fudgeFactor['brLen']
            self.proposals.append(p)
            object.__setattr__(self.tuningsUsage, 'brLen', p)

        # allBrLens
        if self.prob.allBrLens:
            p = Proposal(self)
            p.name = 'allBrLens'
            p.weight = self.prob.allBrLens * \
                (len(self.tree.nodes) - 1) * fudgeFactor['allBrLens']
            self.proposals.append(p)
            object.__setattr__(self.tuningsUsage, 'allBrLens', p)

        # eTBR
        if self.prob.eTBR:
            p = Proposal(self)
            p.name = 'eTBR'
            p.weight = self.prob.eTBR * \
                (len(self.tree.nodes) - 1) * fudgeFactor['eTBR']
            self.proposals.append(p)
            #object.__setattr__(self.tuningsUsage, 'eTBR', p)

        # local
        if self.prob.local:
            p = Proposal(self)
            p.name = 'local'
            p.weight = self.prob.local * \
                (len(self.tree.nodes) - 1) * fudgeFactor['local']
            self.proposals.append(p)
            object.__setattr__(self.tuningsUsage, 'local', p)

        # # treeScale
        # if self.prob.treeScale:
        #     p = Proposal(self)
        #     p.name = 'treeScale'
        #     p.weight = self.prob.treeScale * \
        #         (len(self.tree.nodes) - 1) * fudgeFactor['treeScale']
        #     self.proposals.append(p)
        #     object.__setattr__(self.tuningsUsage, 'treeScale', p)

        # polytomy
        if self.prob.polytomy:
            p = Proposal(self)
            p.name = 'polytomy'
            p.weight = self.prob.polytomy * \
                (len(self.tree.nodes) - 1) * fudgeFactor['polytomy']
            self.proposals.append(p)
            #object.__setattr__(self.tuningsUsage, 'polytomy', p)

        # root3
        if self.prob.root3:
            p = Proposal(self)
            p.name = 'root3'
            if len(self.tree.taxNames) <= 3:
                p.weight = 0.0
            else:
                p.weight = self.prob.root3 * \
                    (len(self.tree.taxNames) - 3) * fudgeFactor['root3']
            self.proposals.append(p)

        # relRate
        if self.prob.relRate:
            if self.tree.model.nParts == 1:
                pass
            if self.tree.model.doRelRates and self.tree.model.relRatesAreFree:
                p = Proposal(self)
                p.name = 'relRate'
                p.weight = self.prob.relRate * self.tree.model.nParts
                self.proposals.append(p)
                object.__setattr__(self.tuningsUsage, 'relRate', p)

        # comp, rMatrix, gdasrv, pInvar, modelThingLocations ...
        for pNum in range(self.tree.model.nParts):
            mp = self.tree.model.parts[pNum]

            # comp
            if self.prob.comp:
                for mtNum in range(mp.nComps):
                    if mp.comps[mtNum].free and not mp.cmd1:
                        p = Proposal(self)
                        p.name = 'comp'
                        p.weight = self.prob.comp * (mp.dim - 1)
                        p.pNum = pNum
                        p.mtNum = mtNum
                        self.proposals.append(p)
                        self.tuningsUsage.parts[pNum].comp.append(p)

            # compDir
            if self.prob.compDir:
                for mtNum in range(mp.nComps):
                    if mp.comps[mtNum].free and not mp.cmd1:
                        p = Proposal(self)
                        p.name = 'compDir'
                        p.weight = self.prob.compDir * (mp.dim - 1)
                        p.pNum = pNum
                        p.mtNum = mtNum
                        self.proposals.append(p)
                        self.tuningsUsage.parts[pNum].compDir.append(p)

            # allCompsDir
            if self.prob.allCompsDir:
                for mtNum in range(mp.nComps):
                    assert mp.comps[mtNum].free
                p = Proposal(self)
                p.name = 'allCompsDir'
                p.weight = self.prob.allCompsDir * (mp.dim - 1) * mp.nComps * fudgeFactor['allCompsDir']
                p.pNum = pNum
                self.proposals.append(p)
                self.tuningsUsage.parts[pNum].allCompsDir.append(p)

            # ndch2_leafCompsDir
            if self.prob.ndch2_leafCompsDir:
                for mtNum in range(mp.nComps):
                    assert mp.comps[mtNum].free
                p = Proposal(self)
                p.name = 'ndch2_leafCompsDir'
                p.weight = self.prob.ndch2_leafCompsDir * (mp.dim - 1) * mp.nComps * fudgeFactor['ndch2comp']
                p.pNum = pNum
                self.proposals.append(p)
                self.tuningsUsage.parts[pNum].ndch2_leafCompsDir.append(p)

            # ndch2_internalCompsDir
            if self.prob.ndch2_internalCompsDir:
                for mtNum in range(mp.nComps):
                    assert mp.comps[mtNum].free
                p = Proposal(self)
                p.name = 'ndch2_internalCompsDir'
                p.weight = self.prob.ndch2_internalCompsDir * (mp.dim - 1) * mp.nComps * fudgeFactor['ndch2comp']
                p.pNum = pNum
                self.proposals.append(p)
                self.tuningsUsage.parts[pNum].ndch2_internalCompsDir.append(p)

            # ndch2_leafCompsDirAlpha
            if self.prob.ndch2_leafCompsDirAlpha:
                p = Proposal(self)
                p.name = 'ndch2_leafCompsDirAlpha'
                p.weight = self.prob.ndch2_leafCompsDirAlpha * (mp.dim - 1) * mp.nComps * fudgeFactor['ndch2alpha']
                p.pNum = pNum
                self.proposals.append(p)
                #self.tuningsUsage.parts[pNum].allCompsDir.append(p)

            # ndch2_internalCompsDirAlpha
            if self.prob.ndch2_internalCompsDirAlpha:
                p = Proposal(self)
                p.name = 'ndch2_internalCompsDirAlpha'
                p.weight = self.prob.ndch2_internalCompsDirAlpha * (mp.dim - 1) * mp.nComps * fudgeFactor['ndch2alpha']
                p.pNum = pNum
                self.proposals.append(p)
                #self.tuningsUsage.parts[pNum].allCompsDir.append(p)


            # # rjComp
            # if self.prob.rjComp:
            #     if mp.rjComp:
            #         for mtNum in range(mp.nComps):
            #             assert mp.comps[mtNum].free
            #         p = Proposal(self)
            #         p.name = 'rjComp'
            #         p.weight = self.prob.rjComp * fudgeFactor['rjComp']
            #         p.pNum = pNum
            #         #p.mtNum = mtNum
            #         self.proposals.append(p)

            # rMatrix
            if self.prob.rMatrix:
                for mtNum in range(mp.nRMatrices):
                    if mp.rMatrices[mtNum].free:
                        p = Proposal(self)
                        p.name = 'rMatrix'
                        if mp.rMatrices[mtNum].spec == '2p':
                            p.weight = self.prob.rMatrix
                            p.variant = '2p'
                        else:
                            p.weight = self.prob.rMatrix * \
                                (((mp.dim * mp.dim) - mp.dim) / 2)
                        p.pNum = pNum
                        p.mtNum = mtNum
                        self.proposals.append(p)
                        self.tuningsUsage.parts[pNum].rMatrix.append(p)

            # rMatrixDir
            if self.prob.rMatrixDir:
                for mtNum in range(mp.nRMatrices):
                    if mp.rMatrices[mtNum].free:
                        p = Proposal(self)
                        p.name = 'rMatrixDir'
                        if mp.rMatrices[mtNum].spec == '2p':
                            p.weight = self.prob.rMatrixDir
                            p.variant = '2p'
                        else:
                            p.weight = self.prob.rMatrixDir * \
                                (((mp.dim * mp.dim) - mp.dim) / 2)
                        p.pNum = pNum
                        p.mtNum = mtNum
                        self.proposals.append(p)
                        self.tuningsUsage.parts[pNum].rMatrixDir.append(p)

            # # rjRMatrix
            # if self.prob.rjRMatrix:
            #     if mp.rjRMatrix:
            #         for mtNum in range(mp.nRMatrices):
            #             assert mp.rMatrices[mtNum].free
            #         p = Proposal(self)
            #         p.name = 'rjRMatrix'
            #         p.weight = self.prob.rjRMatrix * fudgeFactor['rjRMatrix']
            #         p.pNum = pNum
            #         #p.mtNum = mtNum
            #         self.proposals.append(p)

            # gdasrv
            if self.prob.gdasrv:
                for mtNum in range(mp.nGdasrvs):
                    if mp.gdasrvs[mtNum].free:
                        p = Proposal(self)
                        p.name = 'gdasrv'
                        p.weight = self.prob.gdasrv
                        p.pNum = pNum
                        p.mtNum = mtNum
                        self.proposals.append(p)
                        self.tuningsUsage.parts[pNum].gdasrv.append(p)

            # pInvar
            if self.prob.pInvar:
                if mp.pInvar.free:
                    p = Proposal(self)
                    p.name = 'pInvar'
                    p.weight = self.prob.pInvar
                    p.pNum = pNum
                    self.proposals.append(p)
                    object.__setattr__(self.tuningsUsage.parts[pNum], 'pInvar', p)

            # compLocation
            if self.prob.compLocation:
                if mp.nComps > 1 and not mp.cmd1:
                    p = Proposal(self)
                    p.name = 'compLocation'
                    p.weight = self.prob.compLocation * mp.nComps * len(self.tree.nodes) * \
                        fudgeFactor['compLocation']
                    p.pNum = pNum
                    self.proposals.append(p)

            # rMatrixLocation
            if self.prob.rMatrixLocation:
                if mp.nRMatrices > 1:
                    p = Proposal(self)
                    p.name = 'rMatrixLocation'
                    p.weight = self.prob.rMatrixLocation * mp.nRMatrices * (len(self.tree.nodes) - 1) * \
                        fudgeFactor['rMatrixLocation']
                    p.pNum = pNum
                    self.proposals.append(p)

            # gdasrvLocation
            if self.prob.gdasrvLocation:
                if mp.nGdasrvs > 1:
                    p = Proposal(self)
                    p.name = 'gdasrvLocation'
                    p.weight = self.prob.gdasrvLocation * mp.nGdasrvs * (len(self.tree.nodes) - 1) * \
                        fudgeFactor['gdasrvLocation']
                    p.pNum = pNum
                    self.proposals.append(p)

            # # cmd1 stuff
            # if 0 and mp.cmd1:
            #     p = Proposal(self)
            #     p.name = 'cmd1_compDir'
            #     p.weight = self.prob.cmd1_compDir * (mp.dim - 1)
            #     p.pNum = pNum
            #     self.proposals.append(p)

            #     p = Proposal(self)
            #     p.name = 'cmd1_comp0Dir'
            #     p.weight = self.prob.cmd1_comp0Dir * (mp.dim - 1)
            #     p.pNum = pNum
            #     self.proposals.append(p)

            # if 1 and mp.cmd1:
            #     p = Proposal(self)
            #     p.name = 'cmd1_allCompDir'
            #     p.weight = self.prob.cmd1_allCompDir * \
            #         (mp.dim - 1) * len(self.tree.nodes)
            #     p.pNum = pNum
            #     self.proposals.append(p)

            #     p = Proposal(self)
            #     p.name = 'cmd1_alpha'
            #     p.weight = self.prob.cmd1_alpha
            #     p.pNum = pNum
            #     self.proposals.append(p)

        if not self.proposals:
            gm.append("No proposals?")
            raise P4Error(gm)
        self.propWeights = []
        for p in self.proposals:
            # print "%s: %s" % (p.name, p.weight)
            self.propWeights.append(p.weight)
        # print self.propWeights
        self.cumPropWeights = [self.propWeights[0]]
        for i in range(len(self.propWeights))[1:]:
            self.cumPropWeights.append(
                self.cumPropWeights[i - 1] + self.propWeights[i])
        self.totalPropWeights = sum(self.propWeights)
        if self.totalPropWeights < 1e-9:
            gm.append("No proposal weights?")
            raise P4Error(gm)
        for p in self.proposals:
            self.proposalsHash[p.name] = p

    def _refreshProposalProbsAndTunings(self):
        """Adjust proposals after a restart."""

        gm = ['Mcmc._refreshProposalProbsAndTunings()']

        for p in self.proposals:
            # brLen
            if p.name == 'brLen':
                p.weight = self.prob.brLen * \
                    (len(self.tree.nodes) - 1) * fudgeFactor['brLen']

            # allBrLens
            if p.name == 'allBrLens':
                p.weight = self.prob.allBrLens * \
                    (len(self.tree.nodes) - 1) * fudgeFactor['allBrLens']

            # eTBR
            if p.name == 'eTBR':
                p.weight = self.prob.eTBR * \
                    (len(self.tree.nodes) - 1) * fudgeFactor['eTBR']

            # local
            if p.name == 'local':
                p.weight = self.prob.local * \
                    (len(self.tree.nodes) - 1) * fudgeFactor['local']

            # # treeScale
            # if p.name == 'treeScale':
            #     p.weight = self.prob.treeScale * \
            #         (len(self.tree.nodes) - 1) * fudgeFactor['treeScale']

            # polytomy
            if p.name == 'polytomy':
                p.weight = self.prob.polytomy * \
                    (len(self.tree.nodes) - 1) * fudgeFactor['polytomy']

            # root3
            if p.name == 'root3':
                if len(self.tree.taxNames) <= 3:
                    p.weight = 0.0
                else:
                    p.weight = self.prob.root3 * \
                        (len(self.tree.taxNames) - 3) * fudgeFactor['root3']

            # relRate
            if p.name == 'relRate':
                p.weight = self.prob.relRate * self.tree.model.nParts

            # comp
            if p.name == 'comp':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.comp * float(mp.dim - 1)

            # compDir
            if p.name == 'compDir':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.compDir * float(mp.dim - 1)

            # allCompsDir
            if p.name == 'allCompsDir':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.allCompsDir * float(mp.dim - 1) * mp.nComps * fudgeFactor['allCompsDir']

            # ndch2_leafCompsDir
            if p.name == 'ndch2_leafCompsDir':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.ndch2_leafCompsDir * float(mp.dim - 1) * mp.nComps * fudgeFactor['ndch2comp']

            # ndch2_internalCompsDir
            if p.name == 'ndch2_internalCompsDir':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.ndch2_internalCompsDir * float(mp.dim - 1) * mp.nComps * fudgeFactor['ndch2comp']

            # ndch2_leafCompsDirAlpha
            if p.name == 'ndch2_leafCompsDirAlpha':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.ndch2_leafCompsDirAlpha * float(mp.dim - 1) * mp.nComps * fudgeFactor['ndch2alpha']

            # ndch2_internalCompsDirAlpha
            if p.name == 'ndch2_internalCompsDirAlpha':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.ndch2_internalCompsDirAlpha * float(mp.dim - 1) * mp.nComps * fudgeFactor['ndch2alpha']

            # # rjComp
            # if p.name == 'rjComp':
            #     p.weight = self.prob.rjComp * fudgeFactor['rjComp']

            # rMatrix
            if p.name == 'rMatrix':
                mp = self.tree.model.parts[p.pNum]
                if mp.rMatrices[0].spec == '2p':
                    p.weight = self.prob.rMatrix
                else:
                    p.weight = self.prob.rMatrix * \
                        ((((mp.dim * mp.dim) - mp.dim) / 2) - 1)

            # rMatrixDir
            if p.name == 'rMatrixDir':
                mp = self.tree.model.parts[p.pNum]
                if mp.rMatrices[0].spec == '2p':
                    p.weight = self.prob.rMatrixDir
                else:
                    p.weight = self.prob.rMatrixDir * \
                        ((((mp.dim * mp.dim) - mp.dim) / 2) - 1)

            # # rjRMatrix
            # if p.name == 'rjRMatrix':
            #     p.weight = self.prob.rjRMatrix * fudgeFactor['rjRMatrix']

            # gdasrv
            if p.name == 'gdasrv':
                p.weight = self.prob.gdasrv

            # pInvar
            if p.name == 'pInvar':
                p.weight = self.prob.pInvar

            # compLocation
            if p.name == 'compLocation':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.compLocation * mp.nComps * \
                    len(self.tree.nodes) * fudgeFactor['compLocation']

            # rMatrixLocation
            if p.name == 'rMatrixLocation':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.rMatrixLocation * mp.nRMatrices * (len(self.tree.nodes) - 1) * \
                    fudgeFactor['rMatrixLocation']

            # gdasrvLocation
            if p.name == 'gdasrvLocation':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.gdasrvLocation * mp.nGdasrvs * (len(self.tree.nodes) - 1) * \
                    fudgeFactor['gdasrvLocation']

            if p.name.startswith('cmd1_'):
                raise P4Error("fix me!")

        self.propWeights = []
        for p in self.proposals:
            self.propWeights.append(p.weight)
        self.cumPropWeights = [self.propWeights[0]]
        for i in range(len(self.propWeights))[1:]:
            self.cumPropWeights.append(
                self.cumPropWeights[i - 1] + self.propWeights[i])
        self.totalPropWeights = sum(self.propWeights)
        if self.totalPropWeights < 1e-9:
            gm.append("No proposal weights?")
            raise P4Error(gm)

    def writeProposalAcceptances(self):
        """Pretty-print the proposal acceptances."""

        if (self.gen - self.startMinusOne) <= 0:
            print("\nwriteProposalAcceptances()  There is no info in memory. ")
            print(" Maybe it was just emptied after writing to a checkpoint?  ")
            print("If so, read the checkPoint and get the proposalAcceptances from there.")
        else:

            spacer = ' ' * 8
            print("\nProposal acceptances, run %i, for %i gens, from gens %i to %i, inclusive." % (
                self.runNum, (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen))
            print("%s %30s %5s %5s %10s %13s%8s" % (spacer, 'proposal', 'part', 'num', 'nProposals', 'acceptance(%)', 'tuning'))
            for p in self.proposals:
                print("%s" % spacer, end=' ')
                print("%30s" % p.name, end=' ')
                if p.pNum != -1:
                    print(" %3i " % p.pNum, end=' ')
                else:
                    print("   - ", end=' ')
                if p.mtNum != -1:
                    print(" %3i " % p.mtNum, end=' ')
                else:
                    print("   - ", end=' ')
                print("%10i" % p.nProposals[0], end=' ')

                if p.nProposals[0]:  # Don't divide by zero
                    print("       %5.1f " % (100.0 * float(p.nAcceptances[0]) / float(p.nProposals[0])), end=' ')
                else:
                    print("           - ", end=' ')

                if p.tuning == None:
                    print("      -", end=' ')
                elif p.tuning < 2.0:
                    print("  %5.3f" % p.tuning, end=' ')
                else:
                    print("%7.1f" % p.tuning, end=' ')
                print()

            # Tabulate topology changes for 'local', if any were attempted.
            doTopol = 0
            p = None
            try:
                p = self.proposalsHash['local']
            except KeyError:
                pass
            if p:
                for tNum in range(self.nChains):
                    if p.nTopologyChangeAttempts[tNum]:
                        doTopol = 1
                        break
                if doTopol:
                    p = self.proposalsHash['local']
                    print("'Local' proposal-- attempted topology changes")
                    print("%s tempNum   nProps nAccepts percent nTopolChangeAttempts nTopolChanges percent" % spacer)
                    for tNum in range(self.nChains):
                        print("%s" % spacer, end=' ')
                        print("%4i " % tNum, end=' ')
                        print("%9i" % p.nProposals[tNum], end=' ')
                        print("%8i" % p.nAcceptances[tNum], end=' ')
                        print("  %5.1f" % (100.0 * float(p.nAcceptances[tNum]) / float(p.nProposals[tNum])), end=' ')
                        print("%20i" % p.nTopologyChangeAttempts[tNum], end=' ')
                        print("%13i" % p.nTopologyChanges[tNum], end=' ')
                        print("  %5.1f" % (100.0 * float(p.nTopologyChanges[tNum]) / float(p.nTopologyChangeAttempts[tNum])))
                else:
                    print("%sFor the 'local' proposals, there were no attempted" % spacer)
                    print("%stopology changes in any of the chains." % spacer)

            # do the same for eTBR
            doTopol = 0
            p = None
            try:
                p = self.proposalsHash['eTBR']
            except KeyError:
                pass
            if p:
                for tNum in range(self.nChains):
                    if p.nTopologyChangeAttempts[tNum]:
                        doTopol = 1
                        break
                if doTopol:
                    p = self.proposalsHash['eTBR']
                    print("'eTBR' proposal-- attempted topology changes")
                    print("%s tempNum   nProps nAccepts percent nTopolChangeAttempts nTopolChanges percent" % spacer)
                    for tNum in range(self.nChains):
                        print("%s" % spacer, end=' ')
                        print("%4i " % tNum, end=' ')
                        print("%9i" % p.nProposals[tNum], end=' ')
                        print("%8i" % p.nAcceptances[tNum], end=' ')
                        print("  %5.1f" % (100.0 * float(p.nAcceptances[tNum]) / float(p.nProposals[tNum])), end=' ')
                        print("%20i" % p.nTopologyChangeAttempts[tNum], end=' ')
                        print("%13i" % p.nTopologyChanges[tNum], end=' ')
                        print("  %5.1f" % (100.0 * float(p.nTopologyChanges[tNum]) / float(p.nTopologyChangeAttempts[tNum])))
                else:
                    print("%sFor the 'eTBR' proposals, there were no attempted" % spacer)
                    print("%stopology changes in any of the chains." % spacer)

            # Check for aborts.
            p = None
            try:
                p = self.proposalsHash['local']
            except KeyError:
                pass
            if p:
                if hasattr(p, 'nAborts'):
                    if p.nAborts[0]:
                        print("The 'local' proposal had %i aborts. (Not counted in nProps above.)" % p.nAborts[0])
                        print("(Aborts might be due to brLen proposals too big or too small)")
                        if self.constraints:
                            print("(Or, more likely, due to violated constraints.)")
                    else:
                        print("The 'local' proposal had no aborts (either due to brLen proposals")
                        print("too big or too small, or due to violated constraints).")
            p = None
            try:
                p = self.proposalsHash['eTBR']
            except KeyError:
                pass
            if p:
                if hasattr(p, 'nAborts'):
                    if p.nAborts[0]:
                        print("The 'eTBR' proposal had %i aborts.  (Not counted in nProps above.)" % p.nAborts[0])
                        assert self.constraints
                        print("(Aborts due to violated constraints)")
                    else:
                        if self.constraints:
                            print("The 'eTBR' proposal had no aborts (due to violated constraints).")

            for pN in ['polytomy', 'compLocation', 'rMatrixLocation', 'gdasrvLocation']:
                p = None
                try:
                    p = self.proposalsHash[pN]
                except KeyError:
                    pass
                if p:
                    if hasattr(p, 'nAborts'):
                        print("The %15s proposal had %5i aborts." % (p.name, p.nAborts[0]))

    def writeSwapMatrix(self):
        print("\nChain swapping, for %i gens, from gens %i to %i, inclusive." % (
            (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen))
        print("    Swaps are presented as a square matrix, nChains * nChains.")
        print("    Upper triangle is the number of swaps proposed between two chains.")
        print("    Lower triangle is the percent swaps accepted.")
        print("    The current tunings.chainTemp is %5.3f\n" % self.tunings.chainTemp)
        print(" " * 10, end=' ')
        for i in range(self.nChains):
            print("%7i" % i, end=' ')
        print()
        print(" " * 10, end=' ')
        for i in range(self.nChains):
            print("   ----", end=' ')
        print()
        for i in range(self.nChains):
            print(" " * 7, "%2i" % i, end=' ')
            for j in range(self.nChains):
                if i < j:  # upper triangle
                    print("%7i" % self.swapMatrix[i][j], end=' ')
                elif i == j:
                    print("      -", end=' ')
                else:
                    if self.swapMatrix[j][i] == 0:  # no proposals
                        print("      -", end=' ')
                    else:
                        print("  %5.1f" % (100.0 * float(self.swapMatrix[i][j]) / float(self.swapMatrix[j][i])), end=' ')
            print()

    def _makeChainsAndProposals(self):
        """Make chains and proposals."""

        gm = ['Mcmc._makeChainsAndProposals()']

        # random.seed(0)

        # Make chains, if needed
        if not self.chains:
            self.chains = []
            for chNum in range(self.nChains):
                aChain = Chain(self)
                # Temperature.  Set this way to start, but it changes.
                aChain.tempNum = chNum
                self.chains.append(aChain)
        if not self.proposals:
            self._makeProposals()

            # If we are going to be doing the resolution class prior
            # in the polytomy move, we want to pre-compute the logs of
            # T_{n,m}.  Its a vector with indices (ie m) from zero to
            # nTax-2 inclusive.
            if self.proposalsHash.has_key('polytomy') and self.tunings.doPolytomyResolutionClassPrior:
                p = self.proposalsHash['polytomy']
                bigT = p4.func.nUnrootedTreesWithMultifurcations(self.tree.nTax)
                p.logBigT = [0.0] * (self.tree.nTax - 1)
                for i in range(1, self.tree.nTax - 1):
                    p.logBigT[i] = math.log(bigT[i])
                # print p.logBigT

    def _setOutputTreeFile(self):
        """Setup the (output) tree file for the mcmc."""

        gm = ['Mcmc._setOutputTreeFile()']

        # Write the preamble for the trees outfile.
        self.treeFile = open(self.treeFileName, 'w')
        self.treeFile.write('#nexus\n\n')
        self.treeFile.write('begin taxa;\n')
        self.treeFile.write('  dimensions ntax=%s;\n' % self.tree.nTax)
        self.treeFile.write('  taxlabels')
        for tN in self.tree.taxNames:
            self.treeFile.write(' %s' % p4.func.nexusFixNameIfQuotesAreNeeded(tN))
        self.treeFile.write(';\nend;\n\n')

        self.treeFile.write('begin trees;\n')
        self.translationHash = {}
        i = 1
        for tName in self.tree.taxNames:
            self.translationHash[tName] = i
            i += 1

        self.treeFile.write('  translate\n')
        for i in range(self.tree.nTax - 1):
            self.treeFile.write('    %3i %s,\n' % (
                i + 1, p4.func.nexusFixNameIfQuotesAreNeeded(self.tree.taxNames[i])))
        self.treeFile.write('    %3i %s\n' % (
            self.tree.nTax, p4.func.nexusFixNameIfQuotesAreNeeded(self.tree.taxNames[-1])))
        self.treeFile.write('  ;\n')

        # write the models comment
        if self.tree.model.isHet:
            if (not self.tree.model.parts[0].ndch2) or (self.tree.model.parts[0].ndch2 and self.tree.model.parts[0].ndch2_writeComps):
                self.treeFile.write('  [&&p4 models p%i' % self.tree.model.nParts)
                for pNum in range(self.tree.model.nParts):
                    self.treeFile.write(
                        ' c%i.%i' % (pNum, self.tree.model.parts[pNum].nComps))
                    self.treeFile.write(
                        ' r%i.%i' % (pNum, self.tree.model.parts[pNum].nRMatrices))
                    self.treeFile.write(
                        ' g%i.%i' % (pNum, self.tree.model.parts[pNum].nGdasrvs))
                self.treeFile.write(']\n')
        self.treeFile.write('  [Tree numbers are gen+1]\n')
        self.treeFile.close()

        if 0:
            self.prob.dump()
            self.tunings.dump()
            self.writeProposalProbs()
        if 0:
            for p in self.proposals:
                p.dump()
            # return

    def run(self, nGensToDo, verbose=True):
        """Start the Mcmc running."""

        gm = ['Mcmc.run()']

        # Keep track of the first gen of this call to run(), maybe restart
        firstGen = self.gen + 1

        # Hidden experimental hack
        if self.doHeatingHack:
            print("Heating hack is turned on.")
            assert self.nChains == 1, "MCMCMC does not work with the heating hack"
            print("Heating hack temperature is %.2f" % self.heatingHackTemperature)
            #print("Heating hack affects proposals %s" % self.heatingHackProposalNames)

        if self.checkPointInterval:
            # We want a couple of things:
            #  1.  The last gen should be on checkPointInterval.  For
            #      example, if the checkPointInterval is 200, then doing
            #      100 or 300 generations will not be allowed cuz the
            #      chain would continue past the checkPoint-- bad.  Or if
            #      you re-start after 500 gens and change to a
            #      checkPointInterval of 200, then you won't be allowed to
            #      do 500 gens.
            # if ((self.gen + 1) + nGensToDo) % self.checkPointInterval == 0:
            if nGensToDo % self.checkPointInterval == 0:
                pass
            else:
                gm.append(
                    "With the current settings, the last generation won't be on a checkPointInterval.")
                gm.append("self.gen+1=%i, nGensToDo=%i, checkPointInterval=%i" % ((self.gen + 1),
                                                                                  nGensToDo, self.checkPointInterval))
                raise P4Error(gm)
            #  2.  We also want the checkPointInterval to be evenly
            #      divisible by the sampleInterval.
            if self.checkPointInterval % self.sampleInterval == 0:
                pass
            else:
                gm.append(
                    "The checkPointInterval (%i) should be evenly divisible" % self.checkPointInterval)
                gm.append("by the sampleInterval (%i)." % self.sampleInterval)
                raise P4Error(gm)

        if self.proposals:
            # Its either a re-start, or it has been thru autoTune().
            # I can tell the difference by self.gen, which is -1 after
            # autoTune()
            if self.gen == -1:
                self._makeChainsAndProposals()
                self._setOutputTreeFile()
                if self.simulate:
                    self._writeSimFileHeader(self.tree)
            # The probs and tunings may have been changed by the user.
            self._refreshProposalProbsAndTunings()

            # This stuff below should be the same as is done after pickling,
            # see below.
            self.startMinusOne = self.gen

            # Start the tree partitions over.
            self.treePartitions = None
            # Zero the proposal counts
            for p in self.proposals:
                p.nProposals = [0] * self.nChains
                p.nAcceptances = [0] * self.nChains
                p.nTopologyChangeAttempts = [0] * self.nChains
                p.nTopologyChanges = [0] * self.nChains
            # Zero the swap matrix
            if self.nChains > 1:
                self.swapMatrix = []
                for i in range(self.nChains):
                    self.swapMatrix.append([0] * self.nChains)

        else:
            self._makeChainsAndProposals()
            self._setOutputTreeFile()
            if self.simulate:
                self._writeSimFileHeader(self.tree)
        if verbose:
            self.writeProposalIntendedProbs()
            sys.stdout.flush()

        coldChainNum = 0

        # If polytomy is turned on, then it is possible to get a star
        # tree, in which case local will not work.  So if we have both
        # polytomy and local proposals, we should also have brLen.
        if self.proposalsHash.has_key("polytomy") and self.proposalsHash.has_key("local"):
            if not self.proposalsHash.has_key('brLen'):
                gm.append(
                    "If you have polytomy and local proposals, you should have a brLen proposal as well.")
                gm.append(
                    "It can have a low proposal probability, but it needs to be there.")
                gm.append("Turn it on by eg yourMcmc.prob.brLen = 0.001")
                raise P4Error(gm)

        # # Are we using rjComp in any model partitions?
        # rjCompParts = [mp.rjComp for mp in self.chains[
        #     coldChainNum].curTree.model.parts]  # True and False
        # rjCompPartNums = [pNum for pNum in range(
        #     self.chains[coldChainNum].curTree.model.nParts) if rjCompParts[pNum]]
        # # print rjCompParts
        # # print rjCompPartNums

        # # Are we using rjRMatrix in any model partitions?
        # rjRMatrixParts = [mp.rjRMatrix for mp in self.chains[
        #     coldChainNum].curTree.model.parts]  # True and False
        # rjRMatrixPartNums = [pNum for pNum in range(
        #     self.chains[coldChainNum].curTree.model.nParts) if rjRMatrixParts[pNum]]
        # # print rjRMatrixParts
        # # print rjRMatrixPartNums
        # # print self.chains[0].curTree.model.parts[1].rjRMatrix_k

        if self.gen > -1:
            # it is a re-start, so we need to back over the "end;" in the tree
            # files.
            f2 = open(self.treeFileName, 'a+')
            pos = -1
            while 1:
                f2.seek(pos, 2)
                c = f2.read(1)
                if c == ';':
                    break
                pos -= 1
            # print "pos now %i" % pos
            pos -= 3  # end;
            f2.seek(pos, 2)
            c = f2.read(4)
            # print "got c = '%s'" % c
            if c != "end;":
                gm.append(
                    "Mcmc.run().  Failed to find and remove the 'end;' at the end of the tree file.")
                raise P4Error(gm)
            else:
                f2.seek(pos, 2)
                f2.truncate()
            f2.close()

            if verbose:
                print()
                print("Re-starting the MCMC run %i from gen=%i" % (self.runNum, self.gen))
                print("Set to do %i more generations." % nGensToDo)
                if self.writePrams:
                    if self.chains[0].curTree.model.nFreePrams == 0:
                        print("There are no free prams in the model, so I am turning writePrams off.")
                        self.writePrams = False
                sys.stdout.flush()

            self.startMinusOne = self.gen
        else:
            if verbose:
                if self.nChains > 1:
                    print("Using Metropolis-coupled MCMC, with %i chains.  Temperature %.2f" % (self.nChains, self.tunings.chainTemp))
                else:
                    print("Not using Metropolis-coupled MCMC.")
                print("Starting the MCMC %s run %i" % ((self.constraints and "(with constraints)" or ""), self.runNum))
                print("Set to do %i generations." % nGensToDo)

            if self.writePrams:
                if self.chains[0].curTree.model.nFreePrams == 0:
                    if verbose:
                        print("There are no free prams in the model, so I am turning writePrams off.")
                    self.writePrams = False
                else:
                    pramsFile = open(self.pramsFileName, 'a')
                    self.chains[0].curTree.model.writePramsProfile(pramsFile, self.runNum)
                    pramsFile.write("genPlus1")
                    self.chains[0].curTree.model.writePramsHeaderLine(pramsFile)
                    pramsFile.close()
            if self.writeHypers:
                if not self.tree.model.parts[0].ndch2:     # and therefore all model parts, this week
                    self.writeHypers = False
                else:
                    hypersFile = open(self.hypersFileName, 'a')
                    hypersFile.write('genPlus1')
                    self.chains[0].curTree.model.writeHypersHeaderLine(hypersFile)
                    hypersFile.close()

            # if 0 and rjCompPartNums:
            #     rjKFile = open(self.rjKFileName, 'w')
            #     rjKFile.write(
            #         "# k_comp_max, a constant, is the number of comp vectors in a part in total\n")
            #     rjKFile.write(
            #         "#             (both in the pool and not in the pool).\n")
            #     rjKFile.write(
            #         "# ck, aka ModelPart.rjComp_k, is the number of comp vectors in the 'pool' (for each part)\n")
            #     rjKFile.write(
            #         "# k_0 for comps (ck0 below) is the number of comp vectors on the tree (for each part)\n")
            #     rjKFile.write("#\n")
            #     for pNum in rjCompPartNums:
            #         rjKFile.write("# part%i: " % pNum)
            #         rjKFile.write(" k_comp_max = %i\n" % self.chains[
            #                       coldChainNum].curTree.model.parts[pNum].nComps)
            #     rjKFile.write("#\n")

            # if 0 and rjRMatrixPartNums:
            #     if not rjCompPartNums:
            #         rjKFile = open(self.rjKFileName, 'w')
            #     rjKFile.write(
            #         "# k_rMatrix_max, a constant, is the number of rMatrices in a part in total\n")
            #     rjKFile.write(
            #         "#             (both in the pool and not in the pool).\n")
            #     rjKFile.write(
            #         "# rk, aka ModelPart.rjRMatrix_k, is the number of rMatrices in the 'pool' (for each part)\n")
            #     rjKFile.write(
            #         "# k_0 for rMatrices (rk0 below) is the number of rMatrices on the tree (for each part)\n")
            #     rjKFile.write("#\n")
            #     for pNum in rjRMatrixPartNums:
            #         rjKFile.write("# part%i: " % pNum)
            #         rjKFile.write(" k_rMatrix_max = %i\n" % self.chains[
            #                       coldChainNum].curTree.model.parts[pNum].nRMatrices)
            #     rjKFile.write("#\n")

            # if 0 and rjCompPartNums or rjRMatrixPartNums:
            #     rjKFile.write("# %10s " % 'genPlus1')
            #     for pNum in rjCompPartNums:
            #         rjKFile.write("%7s " % 'p%i_ck' % pNum)
            #         rjKFile.write("%8s " % 'p%i_ck0' % pNum)
            #     for pNum in rjRMatrixPartNums:
            #         rjKFile.write("%7s " % 'p%i_rk' % pNum)
            #         rjKFile.write("%8s " % 'p%i_rk0' % pNum)
            #     rjKFile.write("\n")
            #     rjKFile.close()

        if verbose:
            print("Sampling every %i." % self.sampleInterval)
            if self.checkPointInterval:
                print("CheckPoints written every %i." % self.checkPointInterval)
            if nGensToDo <= 20000:
                print("One dot is 100 generations.")
            else:
                print("One dot is 1000 generations.")
            sys.stdout.flush()

        self.treePartitions = None
        realTimeStart = time.time()
        self.lastTimeCheck = time.time()
        for ch in self.chains:
            # 1 means do all
            pf.p4_copyCondLikes(ch.curTree.cTree, ch.propTree.cTree, 1)
            # 1 means do all
            pf.p4_copyBigPDecks(ch.curTree.cTree, ch.propTree.cTree, 1)
            ch.verifyIdentityOfTwoTreesInChain()

        abortableProposals = ['local', 'polytomy', 'compLocation',
                              'rMatrixLocation', 'gdasrvLocation', 'rjComp', 'rjRMatrix']

        for gNum in range(nGensToDo):
            self.gen += 1
            # Do an initial time estimate based on 100 gens
            if nGensToDo > 100 and self.gen - firstGen == 100:
                diff_secs = time.time() - realTimeStart
                total_secs = (float(nGensToDo) / float(100)) * float(diff_secs)
                deltaTime = datetime.timedelta(seconds=int(round(total_secs)))
                print("Estimated completion time: %s days, %s" % (
                    deltaTime.days, time.strftime("%H:%M:%S", time.gmtime(deltaTime.seconds))))

            # Above is a list of proposals where it is possible to abort.
            # When a gen(aProposal) is made, below, aProposal.doAbort
            # might be set, in which case we want to skip it for this
            # gen.  But we want to start each 'for chNum' loop with
            # doAborts all turned off.

            for chNum in range(self.nChains):
                failure = True
                nAttempts = 0
                while failure:
                    # Get the next proposal
                    gotIt = False
                    safety = 0
                    while not gotIt:
                        theRan = random.uniform(0.0, self.totalPropWeights)
                        for i in range(len(self.cumPropWeights)):
                            if theRan < self.cumPropWeights[i]:
                                break
                        aProposal = self.proposals[i]
                        gotIt = True
                        if aProposal.name == 'local':
                            # Can't do local on a star tree.
                            if self.chains[chNum].curTree.nInternalNodes == 1:
                                aProposal = self.proposalsHash['brLen']
                        elif aProposal.name == 'root3':
                            # Can't do root3 on a star tree.
                            if self.chains[chNum].curTree.nInternalNodes == 1:
                                gotIt = False
                        # elif aProposal.name in ['comp', 'compDir']:
                        #     # If we do RJ on this part, make sure the comp is
                        #     # actually in the RJ pool.
                        #     mp = self.chains[chNum].curTree.model.parts[
                        #         aProposal.pNum]
                        #     if mp.rjComp:
                        #         mt = mp.comps[aProposal.mtNum]
                        #         if not mt.rj_isInPool:
                        #             gotIt = False
                        elif aProposal.name in ['rMatrix']:
                            # If we do RJ on this part, make sure the rMatrix
                            # is actually in the RJ pool.
                            mp = self.chains[chNum].curTree.model.parts[
                                aProposal.pNum]
                            if mp.rjRMatrix:
                                mt = mp.rMatrices[aProposal.mtNum]
                                if not mt.rj_isInPool:
                                    gotIt = False
                        if aProposal.doAbort:
                            gotIt = False
                        # if aProposal.name == 'gdasrv':
                        #    gotIt = False
                        if aProposal.name in ['comp', 'compDir']:
                            if self.chains[chNum].curTree.model.parts[aProposal.pNum].nComps > 1:
                                nNodes = self.chains[chNum].curTree.model.parts[
                                    aProposal.pNum].comps[aProposal.mtNum].nNodes
                                if nNodes == 0:
                                    gotIt = False
                        safety += 1
                        if safety > 1000:
                            gm.append(
                                "Could not find a proposal after %i attempts." % safety)
                            gm.append("Possibly a programming error.")
                            gm.append(
                                "Or possibly it is just a pathologically frustrating Mcmc.")
                            raise P4Error(gm)

                    # if gNum % 2:
                    #    aProposal = self.proposalsHash['brLen']
                    # else:
                    #    aProposal = self.proposalsHash['comp']

                    if 0:
                        print("==== gNum=%i, chNum=%i, aProposal=%s (part %i)" % (
                            gNum, chNum, aProposal.name, aProposal.pNum), end=' ')
                        sys.stdout.flush()
                        # print gNum,

                    # success returns None
                    failure = self.chains[chNum].gen(aProposal)

                    if 0:
                        if failure:
                            print("    failure")
                        else:
                            print()

                    nAttempts += 1
                    if nAttempts > 1000:
                        gm.append(
                            "Was not able to do a successful generation after %i attempts." % nAttempts)
                        raise P4Error(gm)
                # print "   Mcmc.run(). finished a gen on chain %i" % (chNum)
                for pr in abortableProposals:
                    if self.proposalsHash.has_key(pr):
                        self.proposalsHash[pr].doAbort = False

            # Do swap, if there is more than 1 chain.
            if self.nChains == 1:
                coldChain = 0
            else:
                # Chain swapping stuff was lifted from MrBayes.  Thanks again.
                chain1, chain2 = random.sample(self.chains, 2)

                thisCh1Temp = None
                thisCh2Temp = None
                # Use the upper triangle of swapMatrix for nProposed's
                if chain1.tempNum < chain2.tempNum:
                    self.swapMatrix[chain1.tempNum][chain2.tempNum] += 1
                    thisCh1Temp = chain1.tempNum
                    thisCh2Temp = chain2.tempNum
                else:
                    self.swapMatrix[chain2.tempNum][chain1.tempNum] += 1
                    thisCh1Temp = chain2.tempNum
                    thisCh2Temp = chain1.tempNum
                    

                lnR = (1.0 / (1.0 + (self.tunings.chainTemp * chain1.tempNum))
                        ) * chain2.curTree.logLike
                lnR += (1.0 / (1.0 + (self.tunings.chainTemp * chain2.tempNum))
                        ) * chain1.curTree.logLike
                lnR -= (1.0 / (1.0 + (self.tunings.chainTemp * chain1.tempNum))
                        ) * chain1.curTree.logLike
                lnR -= (1.0 / (1.0 + (self.tunings.chainTemp * chain2.tempNum))
                        ) * chain2.curTree.logLike

                if lnR < -100.0:
                    r = 0.0
                elif lnR >= 0.0:
                    r = 1.0
                else:
                    r = math.exp(lnR)

                acceptSwap = 0
                if random.random() < r:
                    acceptSwap = 1

                # for continuous temperature tuning with self.swapTuner
                if self.swapTuner and thisCh1Temp == 0 and thisCh2Temp == 1:
                    if acceptSwap:
                        self.swapTuner.swaps01.append(1.0)
                    else:
                        self.swapTuner.swaps01.append(0.0)
                    if len(self.swapTuner.swaps01) > self.swapTuner.sampleSize:
                        self.swapTuner.swaps01.pop(0)   # fifo
                    if len(self.swapTuner.swaps01) == self.swapTuner.sampleSize:
                        thisSum = sum(self.swapTuner.swaps01)
                        thisAccepted = thisSum / self.swapTuner.sampleSize
                        
                        
                        if thisAccepted > self.swapTuner.accHi:   # acceptance too high; temperature too low
                            oldTemp = self.tunings.chainTemp
                            self.tunings.chainTemp *= self.swapTuner.factorHi
                            #print("swap accepted %.2f, increase temp from %.3f to %.3f" % (thisAccepted, oldTemp, self.tunings.chainTemp))
                            self.swapTuner.swaps01 = []
                        elif thisAccepted < self.swapTuner.accLo:  # acceptance too low; temperature too high
                            oldTemp = self.tunings.chainTemp
                            if thisAccepted > self.swapTuner.accB:
                                self.tunings.chainTemp /= self.swapTuner.factorB
                            else:
                                self.tunings.chainTemp /= self.swapTuner.factorC
                            #print("swap accepted %.3f, decrease temp from %.3f to %.3f" % (thisAccepted, oldTemp, self.tunings.chainTemp))
                            self.swapTuner.swaps01 = []
                        
                        

                if acceptSwap:
                    # Use the lower triangle of swapMatrix to keep track of
                    # nAccepted's
                    if chain1.tempNum < chain2.tempNum:
                        self.swapMatrix[chain2.tempNum][chain1.tempNum] += 1
                    else:
                        self.swapMatrix[chain1.tempNum][chain2.tempNum] += 1

                    # Do the swap
                    chain1.tempNum, chain2.tempNum = chain2.tempNum, chain1.tempNum

                # Find the cold chain, the one where tempNum is 0
                coldChainNum = -1
                for i in range(len(self.chains)):
                    if self.chains[i].tempNum == 0:
                        coldChainNum = i
                        break
                if coldChainNum == -1:
                    gm.append(
                        "Unable to find which chain is the cold chain.  Bad.")
                    raise P4Error(gm)

            # If it is a writeInterval, write stuff
            if (self.gen + 1) % self.sampleInterval == 0:
                if 1:
                    likesFile = open(self.likesFileName, 'a')
                    likesFile.write(
                        '%11i %f\n' % (self.gen + 1, self.chains[coldChainNum].curTree.logLike))
                    likesFile.close()

                    # Check the likelihood every write interval
                    if 0:
                        oldLike = self.chains[coldChainNum].curTree.logLike
                        print("gen+1 %11i  %f  " % (
                            self.gen+1, 
                            self.chains[coldChainNum].curTree.logLike), end=' ')
                        self.chains[coldChainNum].curTree.calcLogLike(verbose=False)
                        newLike = self.chains[coldChainNum].curTree.logLike
                        print("%f" % self.chains[coldChainNum].curTree.logLike, end=' ')
                        likeDiff = math.fabs(oldLike - newLike)
                        if likeDiff > 1e-14:
                            print("%f" % likeDiff)
                        else:
                            print()
                                                                  

                    treeFile = open(self.treeFileName, 'a')
                    treeFile.write("  tree t_%i = [&U] " % (self.gen + 1))
                    if self.tree.model.parts[0].ndch2:     # and therefore all model parts
                        if self.tree.model.parts[0].ndch2_writeComps:
                            self.chains[coldChainNum].curTree.writeNewick(treeFile,
                                                                          withTranslation=1,
                                                                          translationHash=self.translationHash,
                                                                          doMcmcCommandComments=True)
                        else:
                            self.chains[coldChainNum].curTree.writeNewick(treeFile,
                                                                          withTranslation=1,
                                                                          translationHash=self.translationHash,
                                                                          doMcmcCommandComments=False)

                    else:
                        self.chains[coldChainNum].curTree.writeNewick(treeFile,
                                                                      withTranslation=1,
                                                                      translationHash=self.translationHash,
                                                                      doMcmcCommandComments=self.tree.model.isHet)
                    treeFile.close()

                if 0 and rjCompPartNums:  # we made rjCompPartNums above
                    rjKFile = open(self.rjKFileName, 'a')
                    rjKFile.write("%12i " % (self.gen + 1))
                    for pNum in rjCompPartNums:
                        mp = self.chains[
                            coldChainNum].curTree.model.parts[pNum]
                        assert mp.rjComp
                        # k is the number of comp vectors in the pool
                        # k_0 is the number of comp vectors on the tree
                        # k_max is the number of comp vectors in total
                        # mp.rjComp_k is the number of comp vectors in the pool
                        k = 0
                        k_0 = 0
                        for cNum in range(mp.nComps):
                            theComp = mp.comps[cNum]
                            if theComp.nNodes:
                                k_0 += 1
                            if theComp.rj_isInPool:
                                k += 1
                        assert k == mp.rjComp_k
                        rjKFile.write("%6i " % k)
                        rjKFile.write("%7i " % k_0)
                    if not rjRMatrixPartNums:
                        rjKFile.write("\n")
                        rjKFile.close()

                if 0 and rjRMatrixPartNums:  # we made rjRMatrixPartNums above
                    if not rjCompPartNums:
                        rjKFile = open(self.rjKFileName, 'a')
                        rjKFile.write("%12i " % (self.gen + 1))
                    for pNum in rjRMatrixPartNums:
                        # print "doing pNum %i" % pNum
                        mp = self.chains[
                            coldChainNum].curTree.model.parts[pNum]
                        assert mp.rjRMatrix
                        # k is the number of rMatrices in the pool
                        # k_0 is the number of rMatrices on the tree
                        # k_max is the number of rMatrices in total
                        # mp.rjRMatrix_k is the number of comp vectors in the
                        # pool
                        k = 0
                        k_0 = 0
                        for rNum in range(mp.nRMatrices):
                            theRMatrix = mp.rMatrices[rNum]
                            if theRMatrix.nNodes:
                                k_0 += 1
                            if theRMatrix.rj_isInPool:
                                k += 1
                        if k != mp.rjRMatrix_k:
                            gm.append("k=%i, mp.rjRMatrix_k=%i" %
                                      (k, mp.rjRMatrix_k))
                            raise P4Error(gm)
                        rjKFile.write("%6i " % k)
                        rjKFile.write("%7i " % k_0)
                    rjKFile.write("\n")
                    rjKFile.close()

                if self.writePrams:
                    pramsFile = open(self.pramsFileName, 'a')
                    #pramsFile.write("%12i " % (self.gen + 1))
                    pramsFile.write("%12i" % (self.gen + 1))
                    self.chains[coldChainNum].curTree.model.writePramsLine(pramsFile)
                    pramsFile.close()

                if self.writeHypers:
                    hypersFile = open(self.hypersFileName, 'a')
                    hypersFile.write("%12i" % (self.gen + 1))
                    self.chains[coldChainNum].curTree.model.writeHypersLine(hypersFile)
                    hypersFile.close()

                # Do a simulation
                if self.simulate:
                    # print "about to simulate..."
                    self._doSimulate(self.chains[coldChainNum].curTree)
                    # print "...finished simulate."

                # Do other stuff.
                if hasattr(self, 'hook'):
                    self.hook(self.chains[coldChainNum].curTree)

                if 0 and self.constraints:
                    print("Mcmc x1c")
                    print(self.chains[0].verifyIdentityOfTwoTreesInChain())
                    print("b checking curTree ..")
                    self.chains[0].curTree.checkSplitKeys()
                    print("b checking propTree ...")
                    self.chains[0].propTree.checkSplitKeys()
                    print("Mcmc xxx")

                # Add curTree to treePartitions
                if self.treePartitions:
                    self.treePartitions._getSplitsFromTree(
                        self.chains[coldChainNum].curTree)
                else:
                    self.treePartitions = TreePartitions(
                        self.chains[coldChainNum].curTree)
                # After _getSplitsFromTree, need to follow, at some point,
                # with _finishSplits().  Do that when it is pickled, or at the
                # end of the run.

                # Checking and debugging constraints
                if 0 and self.constraints:
                    print("Mcmc x1d")
                    print(self.chains[coldChainNum].verifyIdentityOfTwoTreesInChain())
                    print("c checking curTree ...")
                    self.chains[coldChainNum].curTree.checkSplitKeys()
                    print("c checking propTree ...")
                    self.chains[coldChainNum].propTree.checkSplitKeys()
                    # print "c checking that all constraints are present"
                    #theSplits = [n.br.splitKey for n in self.chains[0].curTree.iterNodesNoRoot()]
                    # for sk in self.constraints.constraints:
                    #    if sk not in theSplits:
                    #        gm.append("split %i is not present in the curTree." % sk)
                    #        raise P4Error(gm)
                    print("Mcmc zzz")

                # Check that the curTree has all the constraints
                if self.constraints:
                    splitsInCurTree = [
                        n.br.splitKey for n in self.chains[coldChainNum].curTree.iterInternalsNoRoot()]
                    for sk in self.constraints.constraints:
                        if sk not in splitsInCurTree:
                            gm.append("Programming error.")
                            gm.append(
                                "The current tree (the last tree sampled) does not contain constraint")
                            gm.append(
                                "%s" % p4.func.getSplitStringFromKey(sk, self.tree.nTax))
                            raise P4Error(gm)

                # If it is a checkPointInterval, pickle
                if self.checkPointInterval and (self.gen + 1) % self.checkPointInterval == 0:
                    self._checkPoint()

                    # The stuff below needs to be done in a re-start as well.
                    # See above "if self.proposals:"
                    self.startMinusOne = self.gen

                    # Start the tree partitions over.
                    self.treePartitions = None
                    # Zero the proposal counts
                    for p in self.proposals:
                        p.nProposals = [0] * self.nChains
                        p.nAcceptances = [0] * self.nChains
                        p.nTopologyChangeAttempts = [0] * self.nChains
                        p.nTopologyChanges = [0] * self.nChains
                        p.nAborts = [0] * self.nChains
                    # Zero the swap matrix
                    if self.nChains > 1:
                        self.swapMatrix = []
                        for i in range(self.nChains):
                            self.swapMatrix.append([0] * self.nChains)

                    

            # Reassuring pips ...
            # We want to skip the first gen of every call to run()
            if firstGen != self.gen:
                if nGensToDo <= 20000:
                    if (self.gen - firstGen) % 1000 == 0:
                        if verbose:
                            deltaTime = self._doTimeCheck(
                                nGensToDo, firstGen, 1000)
                            if deltaTime.days:
                                timeString = "%s days, %s" % (
                                    deltaTime.days, time.strftime("%H:%M:%S", time.gmtime(deltaTime.seconds)))
                            else:
                                timeString = time.strftime(
                                    "%H:%M:%S", time.gmtime(deltaTime.seconds))
                            print("%10i - %s" % (self.gen, timeString))

                        else:
                            sys.stdout.write(".")
                            sys.stdout.flush()
                    elif (self.gen - firstGen) % 100 == 0:
                        sys.stdout.write(".")
                        sys.stdout.flush()
                else:
                    if (self.gen - firstGen) % 50000 == 0:
                        if verbose:
                            deltaTime = self._doTimeCheck(
                                nGensToDo, firstGen, 50000)
                            if deltaTime.days:
                                timeString = "%s days, %s" % (
                                    deltaTime.days, time.strftime("%H:%M:%S", time.gmtime(deltaTime.seconds)))
                            else:
                                timeString = time.strftime(
                                    "%H:%M:%S", time.gmtime(deltaTime.seconds))
                            print("%10i - %s" % (self.gen, timeString))
                        else:
                            sys.stdout.write(".")
                            sys.stdout.flush()
                    elif (self.gen - firstGen) % 1000 == 0:
                        sys.stdout.write(".")
                        sys.stdout.flush()

        # Gens finished.  Clean up.
        print()
        if verbose:
            print("Finished %s generations." % nGensToDo)

        treeFile = open(self.treeFileName, 'a')
        treeFile.write('end;\n\n')
        treeFile.close()

    def _doTimeCheck(self, nGensToDo, firstGen, genInterval):
        """Time check 

        firstGen is the first generation of this call to Mcmc.run() else
        timing fails on restart"""
        nowTime = time.time()
        diff_secs = nowTime - self.lastTimeCheck
        total_secs = (float(nGensToDo - (self.gen - firstGen)) /
                      float(genInterval)) * float(diff_secs)
        deltaTime = datetime.timedelta(seconds=int(round(total_secs)))
        self.lastTimeCheck = nowTime
        return deltaTime

    def _writeSimFileHeader(self, curTree):
        simFile = open(self.simFileName, 'a')

        simFile.write(" genPlus1")
        # If self.simulate contains a 1, do unconstrained log like
        if 1 & self.simulate:
            for pNum in range(self.simTree.data.nParts):
                simFile.write(' uncLike%i' % pNum)
        if 2 & self.simulate:  # If self.simulate contains a 2, do bigX^2
            for pNum in range(self.simTree.model.nParts):
                simFile.write(' bigXSq%i' % pNum)
        # If self.simulate contains a 4, do meanNCharPerSite
        if 4 & self.simulate:
            for pNum in range(self.simTree.model.nParts):
                simFile.write(' meanNCharsPerSite%i' % pNum)
        # If self.simulate contains an 8, do c_m, the compStatFromCharFreqs
        if 8 & self.simulate:
            for pNum in range(self.simTree.model.nParts):
                simFile.write(' c_mSim%i  c_mOrig%i' % (pNum, pNum))
        # If self.simulate contains a 16, do constant sites count
        if 16 & self.simulate:
            for pNum in range(self.simTree.model.nParts):
                simFile.write(' nConstSites%i' % pNum)
        simFile.write('\n')
        simFile.close()

    def _doSimulate(self, curTree):
        curTree.copyToTree(self.simTree)
        curTree.model.copyValsTo(self.simTree.model)
        self.simTree.simulate()
        simFile = open(self.simFileName, 'a')
        simFile.write(" %11i" % (self.gen + 1))
        # If self.simulate contains a 1, do unconstrained log like
        if 1 & self.simulate:
            for p in self.simTree.data.parts:
                simFile.write(' %f' % pf.getUnconstrainedLogLike(p.cPart))
        if 2 & self.simulate:  # If self.simulate contains a 2, do bigX^2
            #ret2 = self.simTree.data.compoChiSquaredTest(verbose=0, skipColumnZeros=True)
            # for pNum in range(self.simTree.model.nParts):
            #    simFile.write(' %f' % ret[pNum][0])
            ret = self.simTree.data.simpleBigXSquared()
            for pNum in range(self.simTree.model.nParts):
                simFile.write(' %f' % ret[pNum])
            # for i in range(len(ret)):
            #    if math.fabs(ret[i][0] - ret2[i]) > 0.000001:
            # print "The two methods of bigXSquared calculations differed.  %f
            # and %f" % (ret[i], ret2[i])
        # If self.simulate contains a 4, do meanNCharPerSite
        if 4 & self.simulate:
            ret = self.simTree.data.meanNCharsPerSite()
            # ret is a list, one number per part
            for pNum in range(self.simTree.model.nParts):
                simFile.write(' %f' % ret[pNum])
        # If self.simulate contains an 8, do c_m, the compStatFromCharFreqs
        if 8 & self.simulate:
            ret = self.simTree.compStatFromCharFreqs()
            ret2 = curTree.compStatFromCharFreqs()
            # ret is a list, one number per part
            for pNum in range(self.simTree.model.nParts):
                simFile.write(' %f  %f' % (ret[pNum], ret2[pNum]))
                # print ' compStatFromCharFreqs: %f  %f' % (ret[pNum],
                # ret2[pNum])
        # If self.simulate contains a 16, do constant sites count
        if 16 & self.simulate:
            ret = self.simTree.data.simpleConstantSitesCount()
            # ret is a list, one number per part
            for pNum in range(self.simTree.model.nParts):
                simFile.write(' %i' % ret[pNum])
        simFile.write('\n')
        simFile.close()

    def _checkPoint(self):

        if 0:
            for chNum in range(self.nChains):
                ch = self.chains[chNum]
                print("chain %i ==================" % chNum)
                ch.curTree.summarizeModelThingsNNodes()

        #print("Before removing data", end=' ')
        #self.chains[0].curTree.calcLogLike(verbose=True, resetEmpiricalComps=False)

        # Make a copy of self, but with no cStuff.
        # But we don't want to copy data.  So detach it.
        savedData = self.tree.data
        self.tree.data = None
        if self.simulate:
            savedSimData = self.simTree.data
            self.simTree.data = None
        for chNum in range(self.nChains):
            ch = self.chains[chNum]
            ch.curTree.data = None
            ch.propTree.data = None

        theCopy = copy.deepcopy(self)

        # Re-attach data to self.
        self.tree.data = savedData
        self.tree.calcLogLike(verbose=False, resetEmpiricalComps=False)
        if self.simulate:
            self.simTree.data = savedSimData
            self.simTree.calcLogLike(verbose=False, resetEmpiricalComps=False)
        for chNum in range(self.nChains):
            ch = self.chains[chNum]
            ch.curTree.data = savedData
            #print("After restoring data", end=' ')
            ch.curTree.calcLogLike(verbose=False, resetEmpiricalComps=False)
            ch.propTree.data = savedData
            ch.propTree.calcLogLike(verbose=False, resetEmpiricalComps=False)

        # Get rid of c-pointers in the copy
        theCopy.tree.deleteCStuff()
        theCopy.tree.data = None
        if self.simulate:
            theCopy.simTree.deleteCStuff()
            theCopy.simTree.data = None
        theCopy.treePartitions._finishSplits()
        theCopy.likesFile = None
        theCopy.treeFile = None
        for chNum in range(theCopy.nChains):
            ch = theCopy.chains[chNum]
            ch.curTree.deleteCStuff()
            #ch.curTree.data = None
            ch.propTree.deleteCStuff()
            #ch.propTree.data = None

        # Pickle it.
        fName = "mcmc_checkPoint_%i.%i" % (self.runNum, self.gen + 1)
        f = open(fName, 'wb')
        pickle.dump(theCopy, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    def autoTune(self, gensPerProposal=500, verbose=True, giveUpAfter=10, writeTunings=True, carryOn=False):
        """Attempt to tune the Mcmc automatically.  A bit of a hack.

        Here we let the Mcmc run for a while, and then examine the
        proposal acceptances to see if they are ok, and if they are
        not, make adjustments and do it again.

        We let the chain run in cycles for the number of proposals
        times gensPerProposal (default 500) gens.  The proposals are
        made randomly but with expected equal frequency, to counter
        the effects that proposals being common and rare might
        otherwise have.  At the end of a cycle, the proposals that are
        affected by tunings are examined to see if the acceptances are
        within range.  This week, I am aiming for acceptances from
        10-70% (advice from the authors of MrBayes-- thanks again),
        but to be safe I make adjustments if the acceptances are
        outside 15-60%.  (Although I accept as low as 5% acceptance
        for local).  If adjustments are made, then another cycle is
        done.  When all the acceptances are ok, then the operation
        stops.  If arg 'giveUpAfter' (by default 10) cycles complete
        without getting it right, it gives up.

        The carryOn arg is set to False by default, meaning that if it gives up
        after so many cycles then it dies with a P4Error.  However, if you set
        carryOn to True then it does not die, even though it is not tuned.  This
        may be useful for difficult tunings, as a partially tuned chain may be
        better than completely untuned.

        It is complicated a little because for tree-hetero models,
        some proposals use the same tuning.  For example, if there is
        more than one rMatrix in a given partition, all the rMatrix
        proposals in that partition will use the same tuning.
        Sometimes you might see one rMatrix with too low an acceptance
        rate and its neighbor with too high an acceptance rate-- in
        that case the tuning is left alone.

        The chainTemp is also tested.  I test the acceptance between the cold
        chain and the first heated chain.  If acceptance is less than 1% then
        the temperature is deemed too high and so is lowered, and if the
        acceptance is more than 10% then the temperature is deemed too low and
        raised.

        It is a bit of a hack, so you might see a tuning adjusted on
        one cycle, and then that adjustment is reversed on another
        cycle.

        If you follow this method directly with Mcmc.run(), it uses
        the chains in the state they are left in by this method, which
        might save some burn-in time.

        News: you can pickle the tunings from this method by turning on the arg
        *writeTunings*, which writes a pickle file.  The name of the pickle file
        incorporates the ``runNum``, eg ``mcmc_tunings_0.pickle`` for runNum 0.
        You can then read it by::

            tf = open(tuningsFileName, 'rb')
            tunings = pickle.load(tf)
            tf.close()
            tunings.dump()
        
        and you can then apply the autoTune tuning values to another Mcmc, like
        this::

            m = Mcmc(t, tuningsFileName='myTuningsFile.pickle')


        """

        gm = ['Mcmc.autoTune()']

        if writeTunings:
            tuningsFileName = "mcmc_tunings_%i.pickle" % self.runNum
            if os.path.isfile(tuningsFileName):
                gm.append("Arg 'writeTunings' is on")
                gm.append("File '%s' already exists." % tuningsFileName)
                gm.append("I'm refusing to over-write.  Delete it or move it.")
                raise P4Error(gm)

        if self.proposals:  # Its a re-start
            self.gen = -1
            self.startMinusOne = -1
            #self.proposals = []
            #self.propWeights = []
            #self.cumPropWeights = []
            #self.totalPropWeights = 0.0
            self.treePartitions = None
            # Zero the proposal counts
            for p in self.proposals:
                p.nProposals = [0] * self.nChains
                p.nAcceptances = [0] * self.nChains
                p.nTopologyChangeAttempts = [0] * self.nChains
                p.nTopologyChanges = [0] * self.nChains
            # Zero the swap matrix
            if self.nChains > 1:
                self.swapMatrix = []
                for i in range(self.nChains):
                    self.swapMatrix.append([0] * self.nChains)

        if not self.proposals:
            self._makeChainsAndProposals()

        #coldChainNum = 0
        nGensToDo = gensPerProposal * len(self.proposals)
        if verbose:
            print("Starting the MCMC autoTune()")
            print("There are %i proposals." % len(self.proposals))
            print("Set to do %i samples." % nGensToDo)
            print("One dot is 100 generations.")

        if 0:
            for pr in self.proposals:
                print("pNum = %2i  mtNum = %2i   %s" % (pr.pNum, pr.mtNum, pr.name))  
            print()
            print(self.tuningsUsage)

        for ch in self.chains:
            ch.verifyIdentityOfTwoTreesInChain()

        print("Before autoTune() ...", end=' ')
        self.tunings.dump(advice=False)
        # return
        needsToBeTuned = True  # To start.
        roundCounter = 0

        while needsToBeTuned:
            if verbose:
                print("================ autoTune() round %i ================" % roundCounter)

            # self.chains[0].curTree.model.dump()
            needsToBeTuned = False
            for gNum in range(nGensToDo):
                self.gen += 1
                for chNum in range(self.nChains):
                    # Get the next proposal
                    gotIt = False
                    safety = 0
                    while not gotIt:
                        aProposal = random.choice(self.proposals)
                        gotIt = True
                        if aProposal.name == 'local':
                            # Can't do local on a star tree.
                            if self.chains[chNum].curTree.nInternalNodes == 1:
                                #aProposal = self.proposalsHash['brLen']
                                gotIt = False
                        elif aProposal.name == 'root3':
                            # Can't do root3 on a star tree.
                            if self.chains[chNum].curTree.nInternalNodes == 1:
                                gotIt = False
                        safety += 1
                        if safety > 100:
                            gm.append(
                                "I've been unable to find a suitable proposal after 100 tries.")
                            gm.append(
                                "Its probably a star tree, and none of the proposals can use star trees.")
                            raise P4Error(gm)

                    # print aProposal.name
                    self.chains[chNum].gen(aProposal)

                # Do swap, if there is more than 1 chain.
                if self.nChains == 1:
                    #coldChain = 0
                    pass
                else:
                    # Chain swapping stuff was lifted from MrBayes.  Thanks
                    # again.
                    chain1, chain2 = random.sample(self.chains, 2)

                    # Use the upper triangle of swapMatrix for nProposed's
                    if chain1.tempNum < chain2.tempNum:
                        self.swapMatrix[chain1.tempNum][chain2.tempNum] += 1
                    else:
                        self.swapMatrix[chain2.tempNum][chain1.tempNum] += 1

                    lnR = (
                        1.0 / (1.0 + (self.tunings.chainTemp * chain1.tempNum))) * chain2.curTree.logLike
                    lnR += (1.0 / (1.0 + (self.tunings.chainTemp *
                                          chain2.tempNum))) * chain1.curTree.logLike
                    lnR -= (1.0 / (1.0 + (self.tunings.chainTemp *
                                          chain1.tempNum))) * chain1.curTree.logLike
                    lnR -= (1.0 / (1.0 + (self.tunings.chainTemp *
                                          chain2.tempNum))) * chain2.curTree.logLike

                    if lnR < -100.0:
                        r = 0.0
                    elif lnR >= 0.0:
                        r = 1.0
                    else:
                        r = math.exp(lnR)

                    acceptSwap = 0
                    if random.random() < r:
                        acceptSwap = 1

                    if acceptSwap:
                        # Use the lower triangle of swapMatrix to keep track of
                        # nAccepted's
                        if chain1.tempNum < chain2.tempNum:
                            self.swapMatrix[chain2.tempNum][
                                chain1.tempNum] += 1
                        else:
                            self.swapMatrix[chain1.tempNum][
                                chain2.tempNum] += 1

                        # Do the swap
                        temporary = chain1.tempNum
                        chain1.tempNum = chain2.tempNum
                        chain2.tempNum = temporary

                # Checking and debugging constraints
                if 0 and self.constraints:
                    print("Mcmc x1d")
                    print(self.chains[0].verifyIdentityOfTwoTreesInChain())
                    print("c checking curTree ...")
                    self.chains[0].curTree.checkSplitKeys()
                    print("c checking propTree ...")
                    self.chains[0].propTree.checkSplitKeys()
                    print("c checking that all constraints are present")
                    theSplits = [
                        n.br.splitKey for n in self.chains[0].curTree.iterNodesNoRoot()]
                    for sk in self.constraints.constraints:
                        if sk not in theSplits:
                            gm.append(
                                "split %i is not present in the curTree." % sk)
                            raise P4Error(gm)
                    print("Mcmc zzz")

                # Reassuring pips ...
                # if self.gen and self.gen % 1000 == 0:
                #    print "%10i" % self.gen
                # elif self.gen and self.gen % 100 == 0:
                if self.gen and self.gen % 100 == 0:
                    sys.stdout.write(".")
                    sys.stdout.flush()

            if verbose:
                print()

            atLeast = 100
            for i in range(len(self.proposals)):
                p = self.proposals[i]
                if p.nProposals[0] < atLeast:
                    self.writeProposalProbs()
                    gm.append(
                        "nProposals for proposal %i (%s) is only %i." % (i, p.name, p.nProposals[0]))
                    gm.append(
                        "We want at least %i samples per proposal." % atLeast)
                    gm.append("The sample size is not big enough.")
                    raise P4Error(gm)
            # self.writeProposalAcceptances()
            if 0:
                for p in self.proposals:
                    accepted = float(
                        p.nAcceptances[0]) / float(p.nProposals[0])
                    print("%25s  %5.3f" % (p.name, accepted))

            # Here is where we go over each tuning and ask whether the
            # proposal acceptance is within limits.  There might be
            # more than one proposal for a tuning, so we take that
            # into account.  The limits are 0.1 to 0.7, as suggested
            # by MrBayes.  However, we don't want it too close to the
            # border, so we use 'safe' limits.
            safeLower = 0.15
            safeUpper = 0.60
            safeMultiUpper = 0.40  # For allBrLens, compDir, allCompsDir, 
            safeMultiLower = 0.05
            # It appears that branch length lower limits should be
            # very low.  Say 5%.
            brLenLower = 0.05

            if verbose:
                print("Acceptances for the tunings:")

            theSig = "%25s  %5.3f"
            sig2 = " %-10s"
            sig3 = "%40s %s"

            if self.tuningsUsage.brLen:
                p = self.tuningsUsage.brLen
                accepted = float(p.nAcceptances[0]) / float(p.nProposals[0])
                if verbose:
                    print(theSig % ("brLen", accepted), end=' ')
                if accepted < brLenLower:
                    if verbose:
                        print(sig2 % "too small", end=' ')
                    oldTuning = self.tunings.brLen
                    self.tunings.brLen /= 2.0
                    if verbose:
                        print("tuning currently %5.3f; halve it to %5.3f" % (oldTuning, self.tunings.brLen))
                    needsToBeTuned = True
                elif accepted > safeUpper:
                    if verbose:
                        print(sig2 % "too big", end=' ')
                    oldTuning = self.tunings.brLen
                    self.tunings.brLen *= 2.0
                    if verbose:
                        print("tuning currently %5.3f; double it to %5.3f" % (oldTuning, self.tunings.brLen))
                    needsToBeTuned = True
                else:
                    if verbose:
                        print(sig2 % "ok")

            if self.tuningsUsage.allBrLens:
                p = self.tuningsUsage.allBrLens
                accepted = float(p.nAcceptances[0]) / float(p.nProposals[0])
                if verbose:
                    print(theSig % ("allBrLens", accepted), end=' ')
                if accepted < brLenLower:                     # good?
                    if verbose:
                        print(sig2 % "too small", end=' ')
                    oldTuning = self.tunings.allBrLens
                    self.tunings.allBrLens /= 2.0
                    if verbose:
                        print("tuning currently %5.3f; halve it to %5.3f" % (oldTuning, self.tunings.allBrLens))
                    needsToBeTuned = True
                elif accepted > safeMultiUpper:
                    if verbose:
                        print(sig2 % "too big", end=' ')
                    oldTuning = self.tunings.allBrLens
                    self.tunings.allBrLens *= 2.0
                    if verbose:
                        print("tuning currently %5.3f; double it to %5.3f" % (oldTuning, self.tunings.allBrLens))
                    needsToBeTuned = True
                else:
                    if verbose:
                        print(sig2 % "ok")

            if self.tuningsUsage.local:
                p = self.tuningsUsage.local
                accepted = float(p.nAcceptances[0]) / float(p.nProposals[0])
                if verbose:
                    print(theSig % ("local", accepted), end=' ')
                if accepted < brLenLower:
                    if verbose:
                        print(sig2 % "too small", end=' ')
                    oldTuning = self.tunings.local
                    self.tunings.local /= 2.0
                    #self.tuningsUsage.local.tuning = self.tunings.local
                    if verbose:
                        print("tuning currently %5.3f; halve it to %5.3f" % (oldTuning, self.tunings.local))
                    needsToBeTuned = True
                elif accepted > safeUpper:
                    if verbose:
                        print(sig2 % "too big", end=' ')
                    oldTuning = self.tunings.local
                    self.tunings.local *= 2.0
                    #self.tuningsUsage.local.tuning = self.tunings.local
                    if verbose:
                        print("tuning currently %5.3f; double it to %5.3f" % (oldTuning, self.tunings.local))
                    needsToBeTuned = True
                else:
                    if verbose:
                        print(sig2 % "ok")

            if self.tuningsUsage.relRate:
                p = self.tuningsUsage.relRate
                accepted = float(p.nAcceptances[0]) / float(p.nProposals[0])
                if verbose:
                    print(theSig % ("relRate", accepted), end=' ')
                if accepted < safeLower:
                    if verbose:
                        print(sig2 % "too small", end=' ')
                    oldTuning = self.tunings.relRate
                    self.tunings.relRate /= 2.0
                    #self.tuningsUsage.relRate.tuning = self.tunings.relRate
                    if verbose:
                        print("tuning currently %5.3f; halve it to %5.3f" % (oldTuning, self.tunings.relRate))
                    needsToBeTuned = True
                elif accepted > safeUpper:
                    if verbose:
                        print(sig2 % "too big", end=' ')
                    oldTuning = self.tunings.relRate
                    self.tunings.relRate *= 2.0
                    #self.tuningsUsage.relRate.tuning = self.tunings.relRate
                    if verbose:
                        print("tuning currently %5.3f; double it to %5.3f" % (oldTuning, self.tunings.relRate))
                    needsToBeTuned = True
                else:
                    if verbose:
                        print(sig2 % "ok")

            # if self.tuningsUsage.treeScale:
            #     p = self.tuningsUsage.treeScale
            #     accepted = float(p.nAcceptances[0]) / float(p.nProposals[0])
            #     if verbose:
            #         print theSig % ("treeScale", accepted),
            #     if accepted < safeLower: 
            #         if verbose:
            #             print sig2 % "too small",
            #         oldTuning = self.tunings.treeScale
            #         self.tunings.treeScale /= 2.0
            #         if verbose:
            #             print "tuning currently %5.3f; halve it to %5.3f" % (oldTuning, self.tunings.treeScale)
            #         needsToBeTuned = True
            #     elif accepted > safeUpper:
            #         if verbose:
            #             print sig2 % "too big",
            #         oldTuning = self.tunings.treeScale
            #         self.tunings.treeScale *= 2.0
            #         if verbose:
            #             print "tuning currently %5.3f; double it to %5.3f" % (oldTuning, self.tunings.treeScale)
            #         needsToBeTuned = True
            #     else:
            #         if verbose:
            #             print sig2 % "ok"



            for pNum in range(self.tuningsUsage.nParts):
                if verbose:
                    print("%15s %i" % ("part", pNum))

                # comp
                if self.tuningsUsage.parts[pNum].comp:
                    isTooBig = 0
                    isTooSmall = 0
                    #isVerySmall = 0
                    for p in self.tuningsUsage.parts[pNum].comp:
                        accepted = float(
                            p.nAcceptances[0]) / float(p.nProposals[0])
                        if verbose:
                            print(theSig % ("comp", accepted), end=' ')
                        # if accepted < 0.01:
                        #    if verbose:
                        #        print sig2 % "very small"
                        #    isVerySmall += 1
                        #    isTooSmall += 1
                        if accepted < safeLower:
                            if verbose:
                                print(sig2 % "too small")
                            isTooSmall += 1
                        elif accepted > safeUpper:
                            if verbose:
                                print(sig2 % "too big")
                            isTooBig += 1
                        else:
                            if verbose:
                                print(sig2 % "ok")
                    if isTooBig and isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Combination of too big and too small."))
                            print(sig3 % (' ', "-> leaving it alone."))
                    elif isTooBig:
                        alreadyAtMax = False
                        if math.fabs(self.tunings.parts[pNum].comp - var.mcmcMaxCompAndRMatrixTuning) < 0.0000001:
                            alreadyAtMax = True
                        if verbose:
                            print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].comp))
                        if not alreadyAtMax:
                            myNewVal = self.tunings.parts[pNum].comp * 1.5
                            newlyAtMax = False
                            if myNewVal >= var.mcmcMaxCompAndRMatrixTuning:
                                self.tunings.parts[
                                    pNum].comp = var.mcmcMaxCompAndRMatrixTuning
                                newlyAtMax = True
                            else:
                                self.tunings.parts[pNum].comp = myNewVal
                            # Transfer the new tuning to the proposals
                            # for mt in self.tuningsUsage.parts[pNum].comp:
                            #    mt.tuning = self.tunings.parts[pNum].comp
                            if verbose:
                                print(sig3 % (' ', "-> increase it to %.3f" % self.tunings.parts[pNum].comp))
                                if newlyAtMax:
                                    print(sig3 % (' ', "(now at the maximum)"))
                            needsToBeTuned = True
                        else:
                            if verbose:
                                print(sig3 % (' ', "already at max %.3f" % self.tunings.parts[pNum].comp))
                    elif isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].comp))
                        # if accepted is very small, decrease the tuning by a lot
                        # if isVerySmall:
                        #    self.tunings.parts[pNum].comp /= 3.
                        # else:
                        self.tunings.parts[pNum].comp /= 1.5
                        # Transfer the new tuning to the proposals
                        # for mt in self.tuningsUsage.parts[pNum].comp:
                        #    mt.tuning = self.tunings.parts[pNum].comp
                        if verbose:
                            print(sig3 % (' ', "-> decrease it to %.3f" % self.tunings.parts[pNum].comp))
                        needsToBeTuned = True

                # compDir
                if self.tuningsUsage.parts[pNum].compDir:
                    isTooBig = 0
                    isTooSmall = 0
                    #isVerySmall = 0
                    for p in self.tuningsUsage.parts[pNum].compDir:
                        accepted = float(p.nAcceptances[0]) / float(p.nProposals[0])
                        if verbose:
                            print(theSig % ("compDir", accepted), end=' ')
                        # if accepted < 0.01:
                        #    if verbose:
                        #        print sig2 % "very small"
                        #    isVerySmall += 1
                        #    isTooSmall += 1
                        if accepted < safeMultiLower:
                            if verbose:
                                print(sig2 % "too small")
                            isTooSmall += 1
                        elif accepted > safeMultiUpper:
                            if verbose:
                                print(sig2 % "too big")
                            isTooBig += 1
                        else:
                            if verbose:
                                print(sig2 % "ok")
                    if isTooBig and isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Combination of too big and too small."))
                            print(sig3 % (' ', "-> leaving it alone."))
                    elif isTooBig:
                        if verbose:
                            print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].compDir))
                        self.tunings.parts[pNum].compDir /= 1.5
                        if verbose:
                            print(sig3 % (' ', "-> decrease it to %.3f" % self.tunings.parts[pNum].compDir))
                        needsToBeTuned = True
                    elif isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].compDir))
                        self.tunings.parts[pNum].compDir *= 1.5
                        if verbose:
                            print(sig3 % (' ', "-> increase it to %.3f" % self.tunings.parts[pNum].compDir))
                        needsToBeTuned = True

                # allCompsDir
                if self.tuningsUsage.parts[pNum].allCompsDir:
                    isTooBig = 0
                    isTooSmall = 0
                    #isVerySmall = 0
                    for p in self.tuningsUsage.parts[pNum].allCompsDir:
                        accepted = float(p.nAcceptances[0]) / float(p.nProposals[0])
                        if verbose:
                            print(theSig % ("allCompsDir", accepted), end=' ')
                        if accepted < safeMultiLower:
                            if verbose:
                                print(sig2 % "too small")
                            isTooSmall += 1
                        elif accepted > safeMultiUpper:
                            if verbose:
                                print(sig2 % "too big")
                            isTooBig += 1
                        else:
                            if verbose:
                                print(sig2 % "ok")
                    if isTooBig and isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Combination of too big and too small."))
                            print(sig3 % (' ', "-> leaving it alone."))
                    elif isTooBig:
                        if verbose:
                            print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].allCompsDir))
                        self.tunings.parts[pNum].allCompsDir /= 1.5
                        if verbose:
                            print(sig3 % (' ', "-> decrease it to %.3f" % self.tunings.parts[pNum].allCompsDir))
                        needsToBeTuned = True
                    elif isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].allCompsDir))
                        self.tunings.parts[pNum].allCompsDir *= 1.5
                        if verbose:
                            print(sig3 % (' ', "-> increase it to %.3f" % self.tunings.parts[pNum].allCompsDir))
                        needsToBeTuned = True

                # ndch2_leafCompsDir
                if self.tuningsUsage.parts[pNum].ndch2_leafCompsDir:
                    isTooBig = 0
                    isTooSmall = 0
                    #isVerySmall = 0
                    for p in self.tuningsUsage.parts[pNum].ndch2_leafCompsDir:
                        accepted = float(p.nAcceptances[0]) / float(p.nProposals[0])
                        if verbose:
                            print(theSig % ("ndch2_leafCompsDir", accepted), end=' ')
                        if accepted < safeMultiLower:
                            if verbose:
                                print(sig2 % "too small")
                            isTooSmall += 1
                        elif accepted > safeMultiUpper:
                            if verbose:
                                print(sig2 % "too big")
                            isTooBig += 1
                        else:
                            if verbose:
                                print(sig2 % "ok")
                    if isTooBig and isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Combination of too big and too small."))
                            print(sig3 % (' ', "-> leaving it alone."))
                    elif isTooBig:
                        if verbose:
                            print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].ndch2_leafCompsDir))
                        self.tunings.parts[pNum].ndch2_leafCompsDir /= 2.0
                        if verbose:
                            print(sig3 % (' ', "-> decrease it to %.3f" % self.tunings.parts[pNum].ndch2_leafCompsDir))
                        needsToBeTuned = True
                    elif isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].ndch2_leafCompsDir))
                        self.tunings.parts[pNum].ndch2_leafCompsDir *= 1.5
                        if verbose:
                            print(sig3 % (' ', "-> increase it to %.3f" % self.tunings.parts[pNum].ndch2_leafCompsDir))
                        needsToBeTuned = True

                # ndch2_internalCompsDir
                if self.tuningsUsage.parts[pNum].ndch2_internalCompsDir:
                    isTooBig = 0
                    isTooSmall = 0
                    #isVerySmall = 0
                    for p in self.tuningsUsage.parts[pNum].ndch2_internalCompsDir:
                        accepted = float(p.nAcceptances[0]) / float(p.nProposals[0])
                        if verbose:
                            print(theSig % ("ndch2_internalCompsDir", accepted), end=' ')
                        if accepted < safeMultiLower:
                            if verbose:
                                print(sig2 % "too small")
                            isTooSmall += 1
                        elif accepted > safeMultiUpper:
                            if verbose:
                                print(sig2 % "too big")
                            isTooBig += 1
                        else:
                            if verbose:
                                print(sig2 % "ok")
                    if isTooBig and isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Combination of too big and too small."))
                            print(sig3 % (' ', "-> leaving it alone."))
                    elif isTooBig:
                        if verbose:
                            print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].ndch2_internalCompsDir))
                        self.tunings.parts[pNum].ndch2_internalCompsDir /= 2.0
                        if verbose:
                            print(sig3 % (' ', "-> decrease it to %.3f" % self.tunings.parts[pNum].ndch2_internalCompsDir))
                        needsToBeTuned = True
                    elif isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].ndch2_internalCompsDir))
                        self.tunings.parts[pNum].ndch2_internalCompsDir *= 1.5
                        if verbose:
                            print(sig3 % (' ', "-> increase it to %.3f" % self.tunings.parts[pNum].ndch2_internalCompsDir))
                        needsToBeTuned = True

                # rMatrix
                if self.tuningsUsage.parts[pNum].rMatrix:
                    isTooBig = 0
                    isTooSmall = 0
                    #isVerySmall = 0
                    for p in self.tuningsUsage.parts[pNum].rMatrix:
                        accepted = float(
                            p.nAcceptances[0]) / float(p.nProposals[0])
                        if verbose:
                            print(theSig % ("rMatrix", accepted), end=' ')
                        # if accepted < 0.01:
                        #    if verbose:
                        #        print sig2 % "very small"
                        #    isTooSmall += 1
                        #    isVerySmall += 1
                        if accepted < safeLower:
                            if verbose:
                                print(sig2 % "too small")
                            isTooSmall += 1
                        elif accepted > safeUpper:
                            if verbose:
                                print(sig2 % "too big")
                            isTooBig += 1
                        else:
                            if verbose:
                                print(sig2 % "ok")
                    if isTooBig and isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Combination of too big and too small."))
                            print(sig3 % (' ', "-> leaving it alone."))
                    else:
                        if self.tree.model.parts[pNum].rMatrices[0].spec == '2p':
                            if isTooBig:
                                if verbose:
                                    print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].twoP))
                                self.tunings.parts[pNum].twoP /= 2.0
                                if verbose:
                                    print(sig3 % (' ', "-> halve it to %.3f" % self.tunings.parts[pNum].twoP))
                                needsToBeTuned = True

                            elif isTooSmall:
                                if verbose:
                                    print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].twoP))
                                self.tunings.parts[pNum].twoP *= 2.0
                                if verbose:
                                    print(sig3 % (' ', "-> double it to %.3f" % self.tunings.parts[pNum].twoP))
                                needsToBeTuned = True

                        else:
                            if isTooBig:
                                alreadyAtMax = False
                                if math.fabs(self.tunings.parts[pNum].rMatrix - var.mcmcMaxCompAndRMatrixTuning) < 0.0000001:
                                    alreadyAtMax = True
                                if verbose:
                                    print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].rMatrix))
                                if not alreadyAtMax:
                                    myNewVal = self.tunings.parts[
                                        pNum].rMatrix * 1.5
                                    newlyAtMax = False
                                    if myNewVal >= var.mcmcMaxCompAndRMatrixTuning:
                                        self.tunings.parts[
                                            pNum].rMatrix = var.mcmcMaxCompAndRMatrixTuning
                                        newlyAtMax = True
                                    else:
                                        self.tunings.parts[
                                            pNum].rMatrix = myNewVal
                                    # Transfer the new tuning to the proposals
                                    # for mt in self.tuningsUsage.parts[pNum].rMatrix:
                                    #    mt.tuning = self.tunings.parts[pNum].rMatrix
                                    if verbose:
                                        print(sig3 % (' ', "-> increase it to %.3f" % self.tunings.parts[pNum].rMatrix))
                                        if newlyAtMax:
                                            print(sig3 % (' ', "(now at the maximum)"))
                                    needsToBeTuned = True
                                else:
                                    if verbose:
                                        print(sig3 % (' ', "already at max %.3f" % self.tunings.parts[pNum].rMatrix))

                            elif isTooSmall:
                                if verbose:
                                    print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].rMatrix))
                                # if accepted is very small, decrease the tuning by a lot
                                # if isVerySmall:
                                #    self.tunings.parts[pNum].rMatrix /= 3.0
                                # else:
                                self.tunings.parts[pNum].rMatrix /= 1.5
                                # Transfer the new tuning to the proposals
                                # for mt in self.tuningsUsage.parts[pNum].rMatrix:
                                #    mt.tuning = self.tunings.parts[pNum].rMatrix
                                if verbose:
                                    print(sig3 % (' ', "-> decrease it to %.3f" % self.tunings.parts[pNum].rMatrix))
                                needsToBeTuned = True

                # rMatrixDir
                if self.tuningsUsage.parts[pNum].rMatrixDir:
                    isTooBig = 0
                    isTooSmall = 0
                    #isVerySmall = 0
                    for p in self.tuningsUsage.parts[pNum].rMatrixDir:
                        accepted = float(
                            p.nAcceptances[0]) / float(p.nProposals[0])
                        if verbose:
                            print(theSig % ("rMatrixDir", accepted), end=' ')
                        # if accepted < 0.01:
                        #    if verbose:
                        #        print sig2 % "very small"
                        #    isTooSmall += 1
                        #    isVerySmall += 1
                        if accepted < safeMultiLower:
                            if verbose:
                                print(sig2 % "too small")
                            isTooSmall += 1
                        elif accepted > safeMultiUpper:
                            if verbose:
                                print(sig2 % "too big")
                            isTooBig += 1
                        else:
                            if verbose:
                                print(sig2 % "ok")
                    if isTooBig and isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Combination of too big and too small."))
                            print(sig3 % (' ', "-> leaving it alone."))
                    else:
                        if self.tree.model.parts[pNum].rMatrices[0].spec == '2p':
                            if isTooBig:
                                if verbose:
                                    print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].twoP))
                                self.tunings.parts[pNum].twoP /= 2.0
                                if verbose:
                                    print(sig3 % (' ', "-> halve it to %.3f" % self.tunings.parts[pNum].twoP))
                                needsToBeTuned = True

                            elif isTooSmall:
                                if verbose:
                                    print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].twoP))
                                self.tunings.parts[pNum].twoP *= 2.0
                                if verbose:
                                    print(sig3 % (' ', "-> double it to %.3f" % self.tunings.parts[pNum].twoP))
                                needsToBeTuned = True

                        else:
                            if isTooBig:
                                if verbose:
                                    print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].rMatrixDir))
                                self.tunings.parts[pNum].rMatrixDir = self.tunings.parts[pNum].rMatrixDir / 1.5
                                if verbose:
                                    print(sig3 % (' ', "-> decrease it to %.3f" % self.tunings.parts[pNum].rMatrixDir))
                                needsToBeTuned = True

                            elif isTooSmall:
                                if verbose:
                                    print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].rMatrixDir))
                                self.tunings.parts[pNum].rMatrixDir *= 1.5
                                if verbose:
                                    print(sig3 % (' ', "-> increase it to %.3f" % self.tunings.parts[pNum].rMatrixDir))
                                needsToBeTuned = True

                # gdasrv
                if self.tuningsUsage.parts[pNum].gdasrv:
                    isTooBig = 0
                    isTooSmall = 0
                    for p in self.tuningsUsage.parts[pNum].gdasrv:
                        accepted = float(
                            p.nAcceptances[0]) / float(p.nProposals[0])
                        if verbose:
                            print(theSig % ("gdasrv", accepted), end=' ')
                        if accepted < safeLower:
                            if verbose:
                                print(sig2 % "too small")
                            isTooSmall += 1
                        elif accepted > safeUpper:
                            if verbose:
                                print(sig2 % "too big")
                            isTooBig += 1
                        else:
                            if verbose:
                                print(sig2 % "ok")
                    if isTooBig and isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Combination of too big and too small."))
                            print(sig3 % (' ', "-> leaving it alone."))
                    elif isTooBig:
                        if verbose:
                            print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].gdasrv))
                        self.tunings.parts[pNum].gdasrv *= 2.0
                        # Transfer the new tuning to the proposals
                        # for mt in self.tuningsUsage.parts[pNum].gdasrv:
                        #    mt.tuning = self.tunings.parts[pNum].gdasrv
                        if verbose:
                            print(sig3 % (' ', "-> double it to %.3f" % self.tunings.parts[pNum].gdasrv))
                        needsToBeTuned = True

                    elif isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].gdasrv))
                        self.tunings.parts[pNum].gdasrv /= 2.0
                        # Transfer the new tuning to the proposals
                        # for mt in self.tuningsUsage.parts[pNum].gdasrv:
                        #    mt.tuning = self.tunings.parts[pNum].gdasrv
                        if verbose:
                            print(sig3 % (' ', "-> halve it to %.3f" % self.tunings.parts[pNum].gdasrv))
                        needsToBeTuned = True

                # pInvar
                if self.tuningsUsage.parts[pNum].pInvar:
                    isTooBig = 0
                    isTooSmall = 0
                    p = self.tuningsUsage.parts[pNum].pInvar
                    accepted = float(
                        p.nAcceptances[0]) / float(p.nProposals[0])
                    if verbose:
                        print(theSig % ("pInvar", accepted), end=' ')
                    if accepted < safeLower:
                        if verbose:
                            print(sig2 % "too small")
                        isTooSmall += 1
                    elif accepted > safeUpper:
                        if verbose:
                            print(sig2 % "too big")
                        isTooBig += 1
                    else:
                        if verbose:
                            print(sig2 % "ok")
                    if isTooBig and isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Combination of too big and too small."))
                            print(sig3 % (' ', "-> leaving it alone."))
                    elif isTooBig:
                        if verbose:
                            print(sig3 % (' ', "Current tuning is %.3f" % self.tunings.parts[pNum].pInvar))
                        self.tunings.parts[pNum].pInvar *= 2.0
                        # Transfer the new tuning to the proposals
                        #self.tuningsUsage.parts[pNum].pInvar.tuning = self.tunings.parts[pNum].pInvar
                        if verbose:
                            print(sig3 % (' ', "-> double it to %.3f" % self.tunings.parts[pNum].pInvar))
                        needsToBeTuned = True

                    elif isTooSmall:
                        if verbose:
                            print(sig3 % (' ', "Current tuning is %.1f" % self.tunings.parts[pNum].pInvar))
                        self.tunings.parts[pNum].pInvar /= 2.0
                        # Transfer the new tuning to the proposals
                        #self.tuningsUsage.parts[pNum].pInvar.tuning = self.tunings.parts[pNum].pInvar
                        if verbose:
                            print(sig3 % (' ', "-> halve it to %.3f" % self.tunings.parts[pNum].pInvar))
                        needsToBeTuned = True

            if 1 and self.nChains > 1:
                # Try to autoTune the swaps, by adjusting the chainTemp.  New
                # arbitrary rules, to encourage a high temperature.  I think
                # that the most important number is between the cold chain and
                # the next chain, and it should be low.  I think 1-10% should be
                # OK.

                # First check that there were enough proposals to make a valid
                # calculation.
                tooFews = 0
                for i in range(self.nChains)[:-1]:
                    for j in range(i + 1, self.nChains):
                        if self.swapMatrix[i][j] < atLeast:
                            tooFews += 1
                if tooFews:
                    self.writeSwapMatrix()
                    gm.append(
                        "checking the swap matrix, but there were too few samples taken.")
                    gm.append("I want at least %i samples, but %i swaps had fewer than that." % (
                        atLeast, tooFews))
                    raise P4Error(gm)

                
                # diagonal = []
                # for i in range(self.nChains)[:-1]:
                #     j = i + 1
                #     accepted = float(
                #         self.swapMatrix[j][i]) / float(self.swapMatrix[i][j])
                #     diagonal.append(accepted)
                # print(diagonal)
                # meanDiagonal = p4.func.mean(diagonal)
                # isOK = True
                # if meanDiagonal > 0.7:  # chainTemp is too low
                #     isOK = False
                #     if verbose:
                #         self.writeSwapMatrix()
                #         print("    Current chainTemp is %5.3f; too low" % self.tunings.chainTemp)
                #     self.tunings.chainTemp *= 1.3333
                #     if verbose:
                #         print("    Mean of acceptances on the diagonal is %5.3f" % meanDiagonal)
                #         print("      -> raising the chainTemp by one third, to %5.3f" % self.tunings.chainTemp)
                #     needsToBeTuned = True
                # elif meanDiagonal < 0.5:  # maybe too high
                #     accepted = float(
                #         self.swapMatrix[self.nChains - 1][0]) / float(self.swapMatrix[0][self.nChains - 1])
                #     # print "coldest minus hottest acceptance is %5.3f" %
                #     # accepted
                #     if accepted < 0.01:
                #         isOK = False
                #         if verbose:
                #             self.writeSwapMatrix()
                #             print("    Current chainTemp is %5.3f; too high" % self.tunings.chainTemp)
                #             print("    Mean of acceptances on the diagonal is %5.3f" % meanDiagonal)
                #         if meanDiagonal > 0.25:
                #             self.tunings.chainTemp /= 1.3333
                #             if verbose:
                #                 print("      -> lowering the chainTemp by one quarter, to %5.3f" % self.tunings.chainTemp)
                #         elif meanDiagonal > 0.1:
                #             self.tunings.chainTemp /= 2.0
                #             if verbose:
                #                 print("      -> lowering the chainTemp by half, to %5.3f" % self.tunings.chainTemp)
                #         else:
                #             self.tunings.chainTemp /= 3.0
                #             if verbose:
                #                 print("      -> dividing the chainTemp by 3, to %5.3f" % self.tunings.chainTemp)

                #         needsToBeTuned = True
                myAccepted = float(self.swapMatrix[1][0]) / float(self.swapMatrix[0][1])
                isOK = True
                if myAccepted > 0.1:   # temperature too low
                    isOK = False
                    if verbose:
                        self.writeSwapMatrix()
                        print("    Current chainTemp is %5.3f; too low" % self.tunings.chainTemp)
                    self.tunings.chainTemp *= 1.3333
                    if verbose:
                        print("    Acceptance of chain 0 with chain 1 is is %5.3f" % myAccepted)
                        print("      -> raising the chainTemp by one third, to %5.3f" % self.tunings.chainTemp)
                    needsToBeTuned = True
                elif myAccepted < 0.01:  # temperature too high
                    isOK = False
                    if verbose:
                        self.writeSwapMatrix()
                        print("    Current chainTemp is %5.3f; too high" % self.tunings.chainTemp)
                        print("    Acceptance of chain 0 with chain 1 is %5.3f" % myAccepted)
                    if myAccepted > 0.005:
                        self.tunings.chainTemp /= 1.3333
                        if verbose:
                            print("      -> lowering the chainTemp by one quarter, to %5.3f" % self.tunings.chainTemp)
                    elif myAccepted > 0.001:
                        self.tunings.chainTemp /= 2.0
                        if verbose:
                            print("      -> lowering the chainTemp by half, to %5.3f" % self.tunings.chainTemp)
                    else:
                        self.tunings.chainTemp /= 3.0
                        if verbose:
                            print("      -> dividing the chainTemp by 3, to %5.3f" % self.tunings.chainTemp)

                    needsToBeTuned = True
                        
                if isOK:
                    if verbose:
                        print("    Chain temp appears to be ok.")

            roundCounter += 1
            # if roundCounter >= 1:
            #    needsToBeTuned = False
            if needsToBeTuned and roundCounter >= giveUpAfter:
                self.tunings.dump(advice=False)
                myMessage = "autoTune() has gone thru %i rounds, and it still needs tuning." % giveUpAfter
                if not carryOn:
                    gm.append(myMessage)
                    gm.append("Giving up.  Do it by hand?  Or set carryOn to not give up?")
                    if writeTunings:
                        print("Writing tunings to pickle file '%s'" % tuningsFileName)
                        self.pickleTunings(tuningsFileName)
                    raise P4Error(gm)
                else:
                    # carry on
                    if verbose:
                        print(myMessage)
                        print('carryOn is set, so continuing anyway ...')
                    break
                    
                

            # self.writeProposalProbs()

            # This stuff below should be the same as is done after pickling
            self.startMinusOne = self.gen

            # Start the tree partitions over.
            self.treePartitions = None
            # Zero the proposal counts
            for p in self.proposals:
                p.nProposals = [0] * self.nChains
                p.nAcceptances = [0] * self.nChains
                p.nTopologyChangeAttempts = [0] * self.nChains
                p.nTopologyChanges = [0] * self.nChains
            # Zero the swap matrix
            if self.nChains > 1:
                self.swapMatrix = []
                for i in range(self.nChains):
                    self.swapMatrix.append([0] * self.nChains)

            # print "End of round %i, needsToBeTuned = %s" % (roundCounter - 1, needsToBeTuned)
            # break

        print("\nAfter autoTune() ...", end=' ')
        self.tunings.dump(advice=False)
        self.gen = -1
        self.startMinusOne = -1
        #self.proposals = []
        #self.propWeights = []
        #self.cumPropWeights = []
        #self.totalPropWeights = 0.0
        self.treePartitions = None
        # Zero the proposal counts
        for p in self.proposals:
            p.nProposals = [0] * self.nChains
            p.nAcceptances = [0] * self.nChains
            p.nTopologyChangeAttempts = [0] * self.nChains
            p.nTopologyChanges = [0] * self.nChains
        if self.nChains > 1:
            self.swapMatrix = []
            for i in range(self.nChains):
                self.swapMatrix.append([0] * self.nChains)
        else:
            self.swapMatrix = None
        #self.chains = []

        if writeTunings:
            print("Writing tunings to pickle file '%s'" % tuningsFileName)
            self.pickleTunings(tuningsFileName)

    def pickleTunings(self, tuningsFileName):
        f = open(tuningsFileName, 'wb')
        pickle.dump(self.tunings, f, pickle.HIGHEST_PROTOCOL)
        f.close()
        
    def writeProposalProbs(self, makeDict=False):
        """(Another) Pretty-print the proposal probabilities.

        See also Mcmc.writeProposalAcceptances().
        """

        nProposals = len(self.proposals)
        if not nProposals:
            print("Mcmc.writeProposalProbs().  No proposals (yet?).")
            return
        #intended = self.propWeights[:]
        # for i in range(len(intended)):
        #    intended[i] /= self.totalPropWeights
        # if math.fabs(sum(intended) - 1.0 > 1e-15):
        #    raise P4Error("bad sum of intended proposal probs. %s" % sum(intended))

        nAttained = [0] * nProposals
        nAccepted = [0] * nProposals
        for i in range(nProposals):
            nAttained[i] = self.proposals[i].nProposals[0]
            nAccepted[i] = self.proposals[i].nAcceptances[0]
        sumAttained = float(sum(nAttained))  # should be zero or nGen
        if not sumAttained:
            print("Mcmc.writeProposalProbs().  No proposals have been made.")
            print("Possibly, due to it being a checkPoint interval, nProposals have all been set to zero.")
            return
        # assert int(sumAttained) == self.gen + 1, "sumAttained is %i, should be gen+1, %i." % (
        #    int(sumAttained), self.gen + 1)
        probAttained = []
        for i in range(len(nAttained)):
            probAttained.append(100.0 * float(nAttained[i]) / sumAttained)
        if math.fabs(sum(probAttained) - 100.0 > 1e-13):
            raise P4Error(
                "bad sum of attained proposal probs. %s" % sum(probAttained))

        if not makeDict:
            spacer = ' ' * 4
            print("\nProposal probabilities (%)")
            # print "There are %i proposals" % len(self.proposals)
            print("For %i gens, from gens %i to %i, inclusive." % (
                (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen))
            print("%2s %11s %11s  %11s %10s %18s %5s %5s" % ('', 'nProposals', 'proposed(%)',
                                                             'accepted(%)', 'tuning', 'proposal', 'part', 'num'))
            for i in range(len(self.proposals)):
                print("%2i" % i, end=' ')
                p = self.proposals[i]
                print("   %7i " % self.proposals[i].nProposals[0], end=' ')
                print("   %5.1f    " % probAttained[i], end=' ')
                if nAttained[i]:
                    print("   %5.1f   " % (100.0 * float(nAccepted[i]) / float(nAttained[i])), end=' ')
                else:
                    print("       -   ", end=' ')
                if p.tuning == None:
                    print("       -    ", end=' ')
                elif p.tuning < 2.0:
                    print("   %7.3f  " % p.tuning, end=' ')
                else:
                    print("   %7.1f  " % p.tuning, end=' ')
                print(" %15s" % p.name, end=' ')
                if p.pNum != -1:
                    print(" %3i " % p.pNum, end=' ')
                else:
                    print("   - ", end=' ')
                if p.mtNum != -1:
                    print(" %3i " % p.mtNum, end=' ')
                else:
                    print("   - ", end=' ')
                print()
        else:
            rd = {}
            for i in range(len(self.proposals)):
                p = self.proposals[i]
                if p.pNum != -1:
                    if p.mtNum != -1:
                            pname = p.name + "_%i_%i" % (p.pNum, p.mtNum)
                    else:
                        pname = p.name + "_%i" % p.pNum
                else:
                    if p.mtNum != -1:
                        pname = p.name + "_%i" % p.mtNum
                    else:
                        pname = p.name
                rd[pname] = []
                rd[pname].append(self.proposals[i].nProposals[0])
                rd[pname].append(probAttained[i])
                if nAttained[i]:
                    rd[pname].append((100.0 * float(nAccepted[i]) / float(nAttained[i])))
                else:
                    rd[pname].append(None)
                if p.tuning == None:
                    rd[pname].append(None)
                elif p.tuning < 2.0:
                    rd[pname].append(p.tuning)
                else:
                    rd[pname].append(p.tuning)
            return rd

    def writeProposalIntendedProbs(self):
        """Tabulate the intended proposal probabilities.
        """

        nProposals = len(self.proposals)
        if not nProposals:
            print("Mcmc.writeProposalIntendedProbs().  No proposals (yet?).")
            return
        intended = self.propWeights[:]
        for i in range(len(intended)):
            intended[i] /= self.totalPropWeights
        if math.fabs(sum(intended) - 1.0 > 1e-14):
            raise P4Error("bad sum of intended proposal probs. %s" % sum(intended))

        spacer = ' ' * 4
        print("\nIntended proposal probabilities (%)")
        # print "There are %i proposals" % len(self.proposals)
        print("%2s %11s %30s %5s %5s %12s" % ('', 'intended(%)', 'proposal', 'part', 'num', 'tuning'))
        for i in range(len(self.proposals)):
            print("%2i" % i, end=' ')
            p = self.proposals[i]
            print("   %6.2f    " % (100. * intended[i]), end=' ')

            print(" %27s" % p.name, end=' ')

            if p.pNum != -1:
                print(" %3i " % p.pNum, end=' ')
            else:
                print("   - ", end=' ')
            if p.mtNum != -1:
                print(" %3i " % p.mtNum, end=' ')
            else:
                print("   - ", end=' ')

            if p.tuning == None:
                print(" %12s "% '    -   ', end=' ')
            else:
                if p.tuning < 0.1:
                    print(" %12.4g" % p.tuning, end=' ')
                elif p.tuning < 1.0:
                    print(" %12.4f" % p.tuning, end=' ')
                else:
                    print(" %12s " % p.tuning, end=' ')
            print()
