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
from p4.node import Node
from p4.pnumbers import Numbers
from p4.simtemp import SimTemp,SimTempTemp
import datetime
import numpy
import logging
import statistics
from collections import deque

# for proposal probs
fudgeFactor = {}
fudgeFactor['local'] = 1.0
fudgeFactor['brLen'] = 1.0
fudgeFactor['eTBR'] = 1.0
fudgeFactor['allBrLens'] = 1.0
fudgeFactor['polytomy'] = 1.0
fudgeFactor['root3'] = 0.05
fudgeFactor['root3n'] = 0.1
fudgeFactor['root2'] = 2.0
fudgeFactor['compLocation'] = 0.01
fudgeFactor['rMatrixLocation'] = 0.01
fudgeFactor['gdasrvLocation'] = 0.01
fudgeFactor['allCompsDir'] = 1.0
fudgeFactor['allRMatricesDir'] = 1.0
fudgeFactor['ndch2comp'] = 0.2
fudgeFactor['ndch2priorRefComp'] = 1.0
fudgeFactor['ndch2alpha'] = 0.01
fudgeFactor['ndrh2rates'] = 0.2
fudgeFactor['ndrh2alpha'] = 0.01
fudgeFactor['ndrh2priorRefRMatrix'] = 1.0


class McmcTuningsPart(object):

    def __init__(self, partNum):

        self.num = partNum
        self.default = {}
        # It is no longer changed depending on the dim
        self.default['comp'] = 0.3
        # This would depend on the dim; this is done in Mcmc.__init__()
        self.default['compDir'] = 100.
        self.default['allCompsDir'] = 500.
        self.default['ndch2_leafCompsDir'] = 200.
        self.default['ndch2_internalCompsDir'] = 100.
        self.default['ndch2_leafCompsDirAlpha'] = 2.0 * math.log(1.2)
        self.default['ndch2_internalCompsDirAlpha'] = 2.0 * math.log(3.0)
        self.default['ndch2_priorRefComp'] = 100.
        # rMatrix with sliders no longer changed depending on the dim (ie size of rMatrix)
        self.default['rMatrix'] = 0.3
        # rMatrixDir would depend on the dim; this is done in Mcmc.__init__()
        self.default['rMatrixDir'] = 200.
        self.default['allRMatricesDir'] = 500.
        self.default['twoP'] = 50.

        self.default['ndrh2_leafRatesDir'] = 200.
        self.default['ndrh2_internalRatesDir'] = 100.
        self.default['ndrh2_leafRatesDirAlpha'] = 2.0 * math.log(1.2)
        self.default['ndrh2_internalRatesDirAlpha'] = 2.0 * math.log(3.0)
        self.default['ndrh2_priorRefRMatrix'] = 100.


        self.default['gdasrv'] = 2.0 * math.log(1.5)  # 0.811
        self.default['pInvar'] = 0.5
        #self.default['rMatrixLocation'] = 0.0


class McmcTunings(object):

    def __init__(self, nParts):
        self.nParts = nParts
        self.parts = []
        for pNum in range(nParts):
            self.parts.append(McmcTuningsPart(pNum))
        self.default = {}
        self.default['relRate'] = 0.5
        # This next tuning is set so that by default the brLens go up or down
        # maximum 10%, ie from 0.909 to 1.1
        self.default['local'] = 2.0 * math.log(1.1)  # 0.1906
        self.default['brLen'] = 2.0 * math.log(2.0)  # 1.386
        # Crux has  2.0 * math.log(1.6) as self.default
        self.default['etbrLambda'] = 2.0 * math.log(1.6)
        self.default['etbrPExt'] = 0.8
        self.default['brLenPriorLambda'] = 10.0
        self.default['brLenPriorLambdaForInternals'] = 1000.0
        self.default['doInternalBrLenPrior'] = False
        self.default['brLenPriorType'] = 'exponential'
        self.default['allBrLens'] = 2.0 * math.log(1.02)   # 0.0396
        self.default['root2'] = 0.1
        #self.default['root3'] = 10.0
        #self.default['root3n'] = 10.0



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
        #object.__setattr__(self, 'comp', 1.0)
        object.__setattr__(self, 'compDir', 0.0)
        object.__setattr__(self, 'allCompsDir', 1.0)
        object.__setattr__(self, 'ndch2_leafCompsDir', 0.0)
        object.__setattr__(self, 'ndch2_internalCompsDir', 0.0)
        object.__setattr__(self, 'ndch2_leafCompsDirAlpha', 0.0)
        object.__setattr__(self, 'ndch2_internalCompsDirAlpha', 0.0)
        object.__setattr__(self, 'ndch2_priorRefCompDir', 0.0)
        #object.__setattr__(self, 'rMatrix', 1.0)
        object.__setattr__(self, 'rMatrixDir', 0.0)
        object.__setattr__(self, 'allRMatricesDir', 1.0)

        object.__setattr__(self, 'ndrh2_leafRatesDir', 0.0)
        object.__setattr__(self, 'ndrh2_internalRatesDir', 0.0)
        object.__setattr__(self, 'ndrh2_leafRatesDirAlpha', 0.0)
        object.__setattr__(self, 'ndrh2_internalRatesDirAlpha', 0.0)
        object.__setattr__(self, 'ndrh2_priorRefRMatrixDir', 0.0)

        object.__setattr__(self, 'gdasrv', 1.0)
        object.__setattr__(self, 'pInvar', 1.0)
        object.__setattr__(self, 'local', 1.0)
        object.__setattr__(self, 'brLen', 0.0)
        object.__setattr__(self, 'allBrLens', 1.0)
        object.__setattr__(self, 'eTBR', 1.0)
        object.__setattr__(self, 'polytomy', 0.0)
        object.__setattr__(self, 'root3', 0.0)
        object.__setattr__(self, 'root3n', 0.0)
        object.__setattr__(self, 'root2', 0.0)
        object.__setattr__(self, 'compLocation', 0.0)
        object.__setattr__(self, 'rMatrixLocation', 0.0)
        object.__setattr__(self, 'relRate', 1.0)

    def __setattr__(self, item, val):
        # complaintHead = "\nMcmcProposalProbs.__setattr__()"
        gm = ["\nMcmcProposalProbs(). (set %s to %s)" % (item, val)]
        theKeys = self.__dict__.keys()
        if item in theKeys:
            try:
                val = float(val)
                if val < 1e-15:
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
        theKeys = list(self.__dict__.keys())
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
        #self.mtNum = -1
        self.weight = 1.0
        self.tuning = None
        self.tuningLimitHi = None
        self.tuningLimitLo = None
        self.tunings = {}
        self.nProposals = [0] * theMcmc.nChains
        self.nAcceptances = [0] * theMcmc.nChains
        self.accepted = 0
        self.topologyChanged = 0
        self.nTopologyChangeAttempts = [0] * theMcmc.nChains
        self.nTopologyChanges = [0] * theMcmc.nChains
        self.doAbort = False
        self.nAborts = [0] * theMcmc.nChains

        self.tnSampleSize = 250
        self.tnNSamples = [0] * theMcmc.nChains
        self.tnNAccepts = [0] * theMcmc.nChains
        self.tnAccVeryHi = None
        self.tnAccHi = None
        self.tnAccLo = None
        self.tnAccVeryLo = None
        self.tnFactorVeryHi = None
        self.tnFactorHi = None
        self.tnFactorLo = None
        self.tnFactorVeryLo = None

        self.prior = None

    def dump(self):
        print("proposal name=%-10s pNum=%2s, weight=%s, tuning=%s" % (
            '%s,' % self.name, self.pNum, self.weight, self.tuning))
        #print("proposal name=%-10s pNum=%2i, weight=%5.1f, tuning=%7.2f" % (
        #    '%s,' % self.name, self.pNum, self.weight, self.tuning))
        #print("    nProposals   by temperature:  %s" % self.nProposals)
        #print("    nAcceptances by temperature:  %s" % self.nAcceptances)


    def tune(self, tempNum):
        assert self.tnSampleSize >= 100.
        assert self.tnNSamples[tempNum] >= self.tnSampleSize
        acc = self.tnNAccepts[tempNum] / self.tnNSamples[tempNum]
        doMessage = False
        if acc > self.tnAccHi:
            oldTn = self.tuning[tempNum]
            if acc > self.tnAccVeryHi:
                self.tuning[tempNum] *= self.tnFactorVeryHi
            else:
                self.tuning[tempNum] *= self.tnFactorHi
            doMessage = True
        elif acc < self.tnAccLo:
            oldTn = self.tuning[tempNum]
            if acc < self.tnAccVeryLo:
                self.tuning[tempNum] *= self.tnFactorVeryLo
            else:
                self.tuning[tempNum] *= self.tnFactorLo
            doMessage = True
        extraMessage = None
        if doMessage and self.tuningLimitHi and (self.tuning[tempNum] > self.tuningLimitHi):
            self.tuning[tempNum] = self.tuningLimitHi
            extraMessage = " (tuningLimitHi)"
        self.tnNSamples[tempNum] = 0
        self.tnNAccepts[tempNum] = 0
        if var.mcmc_logTunings and doMessage:
            message = "%s tune  gen=%i tempNum=%i acceptance=%.3f " % (self.name, self.mcmc.gen, tempNum, acc)
            message += "(target %.3f -- %.3f) " % (self.tnAccLo, self.tnAccHi)
            message += "Adjusting tuning from %g to %g" % (oldTn, self.tuning[tempNum])
            if extraMessage:
                message += extraMessage
            #print(message)
            self.mcmc.logger.info(message)



class Proposals(object):
    def __init__(self):
        self.proposals = []
        self.topologyProposalsDict = {}
        self.propWeights = []
        self.cumPropWeights = []
        self.totalPropWeights = 0.0
        self.intended = None
        self.pDict = {}

    def summary(self):
        print("There are %i proposals" % len(self.proposals))
        for p in self.proposals:
            print("proposal name=%-10s pNum=%2s, weight=%s, tuning=%s" % (
                '%s,' % p.name, p.pNum, p.weight, p.tuning))
            
    def calculateWeights(self):
        gm = ["Proposals.calculateWeights()"]
        self.propWeights = []
        for p in self.proposals:
            #print("%s: %s" % (p.name, p.weight))
            self.propWeights.append(p.weight)
        #print(self.propWeights)
        self.cumPropWeights = [self.propWeights[0]]
        for i in range(len(self.propWeights))[1:]:
            self.cumPropWeights.append(
                self.cumPropWeights[i - 1] + self.propWeights[i])
        self.totalPropWeights = sum(self.propWeights)
        if self.totalPropWeights < 1e-9:
            gm.append("No proposal weights?")
            raise P4Error(gm)
        self.intended = self.propWeights[:]
        for i in range(len(self.intended)):
            self.intended[i] /= self.totalPropWeights
        if math.fabs(sum(self.intended) - 1.0 > 1e-14):
            raise P4Error("bad sum of intended proposal probs. %s" % sum(self.intended))
        #print(self.intended)

    def chooseProposal(self, equiProbableProposals):
        if equiProbableProposals:
            return random.choice(self.proposals)
        else:
            theRan = random.uniform(0.0, self.totalPropWeights)
            for i in range(len(self.cumPropWeights)):
                if theRan < self.cumPropWeights[i]:
                    break
            return self.proposals[i]
        

    def writeProposalIntendedProbs(self):
        """Tabulate the intended proposal probabilities"""

        spacer = ' ' * 4
        print("\nIntended proposal probabilities (%)")
        print("There are %i proposals" % len(self.proposals))
        print("%2s %11s %35s %5s %12s" % ('', 'intended(%)', 'proposal', 'part', 'tuning'))
        for i in range(len(self.proposals)):
            print("%2i" % i, end=' ')
            p = self.proposals[i]
            print("   %6.2f    " % (100. * self.intended[i]), end=' ')

            print(" %32s" % p.name, end=' ')

            if p.pNum != -1:
                print(" %3i " % p.pNum, end=' ')
            else:
                print("   - ", end=' ')

            if p.tuning == None:
                print(" %12s "% '    -   ', end=' ')
            else:
                if p.tuning[0] < 0.1:
                    print(" %12.4g" % p.tuning[0], end=' ')
                elif p.tuning[0] < 1.0:
                    print(" %12.4f" % p.tuning[0], end=' ')
                elif p.tuning[0] < 10.0:
                    print(" %12.3f" % p.tuning[0], end=' ')
                elif p.tuning[0] < 1000.0:
                    print(" %12.1f" % p.tuning[0], end=' ')
                else:
                    print(" %12.2g " % p.tuning[0], end=' ')
            print()

    def writeTunings(self):
        print("Proposal tunings:")
        print("%20s %12s" % ("proposal name", "tuning"))
        for p in self.proposals:
            print("%20s" % p.name, end=' ')
            if p.tuning:
                # if p.tuning < 10.0:
                #     print("%12.3f" % p.tuning, end=' ')
                # else:
                #     print("%12.1f" % p.tuning, end=' ')
                print(p.tuning)
            else:
                print("    %4s    " % '-', end=' ')
            print()



class SwapTuner(object):
    """Continuous tuning for swap temperature, matrix version"""

    def __init__(self, sampleSize):
        assert sampleSize >= 100
        self.sampleSize = sampleSize
        self.swaps01_nAttempts = 0
        self.swaps01_nSwaps = 0

        self.tnAccVeryHi = 0.18
        self.tnAccHi = 0.12
        self.tnAccLo = 0.04
        self.tnAccVeryLo = 0.01
        self.tnFactorVeryHi = 1.4
        self.tnFactorHi = 1.2
        self.tnFactorLo = 0.9
        self.tnFactorVeryLo = 0.6
        self.tnFactorZero = 0.4


    def tune(self, theMcmc):
        """A bad idea, at least as implemented.  It is unstable."""
        assert self.swaps01_nAttempts >= self.sampleSize
        acc = float(self.swaps01_nSwaps) / self.swaps01_nAttempts    # float() for Py2
        #print("SwapTuner.tune() nSwaps %i, nAttemps %i, acc %s" % (
        #    self.swaps01_nSwaps, self.swaps01_nAttempts, acc))
        doMessage = False
        direction = None
        if acc > self.tnAccHi:
            oldTn = theMcmc.chainTemp
            if acc > self.tnAccVeryHi:
                theMcmc.chainTemp *= self.tnFactorVeryHi
            else:
                theMcmc.chainTemp *= self.tnFactorHi
            doMessage = True
            direction = 'Increase'
        elif acc < self.tnAccLo:
            oldTn = theMcmc.chainTemp
            if acc == 0.0:   # no swaps at all
                theMcmc.chainTemp *= self.tnFactorZero
            elif acc < self.tnAccVeryLo:
                theMcmc.chainTemp *= self.tnFactorVeryLo
            else:
                theMcmc.chainTemp *= self.tnFactorLo
            doMessage = True
            direction = 'Decrease'
        self.swaps01_nAttempts = 0
        self.swaps01_nSwaps = 0
        if doMessage:
            message = "%s tune  gen=%i acceptance=%.3f " % ('chainTemp', theMcmc.gen, acc)
            message += "(target %.3f -- %.3f) " % (self.tnAccLo, self.tnAccHi)
            message += "%s chainTemp from %g to %g" % (direction, oldTn, theMcmc.chainTemp)
            #print(message)
            theMcmc.logger.info(message)

class SwapTunerV(object):
    """Continuous tuning for swap temperature, vector version"""

    def __init__(self, theMcmc):
        assert theMcmc.swapTunerSampleSize >= 100
        self.mcmc = theMcmc
        self.nChains = theMcmc.nChains

        # These are for adjacent pairs. Eg for attempts between chains 0 and 1,
        # we increment self.nAttempts[0], ie it is indexed with the lower number
        # in the pair.
        self.nAttempts = [0] * self.nChains
        self.nSwaps = [0] * self.nChains

        self.tnAccVeryHi = 0.30
        self.tnAccHi = 0.25
        self.tnAccLo = 0.10
        self.tnAccVeryLo = 0.05
        self.tnFactorVeryHi = 1.4
        self.tnFactorHi = 1.2
        self.tnFactorLo = 0.9
        self.tnFactorVeryLo = 0.6
        self.tnFactorZero = 0.4


    def tune(self, theTempNum):
        assert self.nAttempts[theTempNum] >= self.mcmc.swapTunerSampleSize
        acc = float(self.nSwaps[theTempNum]) / self.nAttempts[theTempNum]    # float() for Py2
        # print("SwapTunerV.tune() theTempNum %i, nSwaps %i, nAttemps %i, acc %s" % (
        #     theTempNum, self.nSwaps[theTempNum], self.nAttempts[theTempNum], acc))
        # print("tempDiffs %s" % self.mcmc.chainTempDiffs)
        # print("temps     %s" % self.mcmc.chainTemps)

        tnLimitHi = 1.0e5
        doMessage = False
        direction = None
        oldTn = self.mcmc.chainTempDiffs[theTempNum]
        if acc > self.tnAccHi:
            if self.mcmc.chainTempDiffs[theTempNum] >= tnLimitHi:
                message = f"chainTemp tune  gen={self.mcmc.gen} tempNum={theTempNum}  diff %g >= tnLimitHi %g --- no change" % (
                    self.mcmc.chainTempDiffs[theTempNum], tnLimitHi)
                self.mcmc.logger.info(message)
                direction = "no change"
            else:
                if acc > self.tnAccVeryHi:
                    self.mcmc.chainTempDiffs[theTempNum] *= self.tnFactorVeryHi
                else:
                    self.mcmc.chainTempDiffs[theTempNum] *= self.tnFactorHi
                doMessage = True
                direction = 'Increase'
        elif acc < self.tnAccLo:
            #if self.mcmc.chainTempDiffs[theTempNum] <= self.tnLimitLo:
            #    direction = "no change"
            #else:
            if acc == 0.0:   # no swaps at all
                self.mcmc.chainTempDiffs[theTempNum] *= self.tnFactorZero
            elif acc < self.tnAccVeryLo:
                self.mcmc.chainTempDiffs[theTempNum] *= self.tnFactorVeryLo
            else:
                self.mcmc.chainTempDiffs[theTempNum] *= self.tnFactorLo
            doMessage = True
            direction = 'Decrease'
        self.nAttempts[theTempNum] = 0
        self.nSwaps[theTempNum] = 0
        if direction != "no change":
            if doMessage:
                message = "%s tune  gen=%i tempNum=%i acceptance=%.3f " % ('chainTemp', self.mcmc.gen, theTempNum, acc)
                message += "(target %.3f -- %.3f) " % (self.tnAccLo, self.tnAccHi)
                message += "%s chainTempDiff from %g to %g" % (direction, oldTn, self.mcmc.chainTempDiffs[theTempNum])
                #print(message)
                self.mcmc.logger.info(message)
        # Make chainTemps from chainTempDiffs
        if self.mcmc.doHeatingHack:
            self.mcmc.chainTemps = [self.mcmc.heatingHackTemperature]
        else:
            self.mcmc.chainTemps = [0.0]
        for dNum in range(self.mcmc.nChains - 1):
            self.mcmc.chainTemps.append(self.mcmc.chainTempDiffs[dNum] + self.mcmc.chainTemps[-1])
        if doMessage:
            message = "new chainTemps gen=%i " % (self.mcmc.gen)
            for cT in self.mcmc.chainTemps:
                message += "%10.2f" % cT
            self.mcmc.logger.info(message)




class Mcmc(object):

    """An MCMC for molecular sequences.

    aTree
                  The tree should have a model and data attached.

    nChains
                  The number of chains in the MCMCMC, default 4

    runNum
                  You may want to do more than one 'run' in the same directory,
                  to facilitate convergence testing.  The first runNum would be
                  0, and samples, likelihoods, and checkPoints are written to
                  files with that number.

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

        m = Mcmc(t, sampleInterval=1000, checkPointInterval=200000)

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

    def __init__(self, aTree, nChains=4, runNum=0, sampleInterval=100, checkPointInterval=10000, simulate=None, writePrams=True, constraints=None, verbose=True, simTempNTemps=None, simTempMax=10.0):
        gm = ['Mcmc.__init__()']
        self.verbose = verbose

        if aTree and aTree.model and aTree.data:
            pass
        else:
            gm.append("The tree that you feed to this class should have a model and data attached.")
            raise P4Error(gm)


        self.isBiRoot = False
        rootNChildren = aTree.root.getNChildren()
        if var.mcmc_allowUnresolvedStartingTree:
            if rootNChildren == 2:
                self.isBiRoot = True
                # Add a rooter node to aTree.nodes, not "in" the tree
                # I did that, but both local and eTBR caused grief. 
        else:
            if rootNChildren == 2:
                ret = aTree.isFullyBifurcating(verbose=True, biRoot=True)
                if not ret:
                    gm.append("The tree has a bifurcating root, but otherwise is not fully bifurcating.")
                    raise P4Error(gm)
                self.isBiRoot = True
                # Add a rooter node to aTree.nodes, not "in" the tree
                n = Node()
                n.isLeaf = 1
                n.name = 'tempRooter'
                n.nodeNum = var.NO_ORDER
                aTree.nodes.append(n)
                aTree.setPreAndPostOrder()

            elif rootNChildren == 3:
                ret = aTree.isFullyBifurcating(verbose=True, biRoot=False)
                if not ret:
                    gm.append("The tree has a trifurcating root, but otherwise is not fully bifurcating.")
                    raise P4Error(gm)
                self.isBiRoot = False
            elif rootNChildren == 1:
                ret = aTree.isFullyBifurcating(verbose=True, biRoot=False)
                if not ret:
                    gm.append("The tree has a monofurcating root, but otherwise is not fully bifurcating.")
                    raise P4Error(gm)
                self.isBiRoot = False
            else:
                gm.append("Mcmc only allows trifurcating or bifurcating roots.  This tree root has %i children" % rootNChildren)
                raise P4Error(gm)
            

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
                    gm.append('Constraint %i %s' % (sk, p4.func.getSplitStringFromKey(sk, self.tree.nTax)))
                    gm.append('is not in the starting tree.')
                    gm.append('Maybe you want to make a randomTree with constraints?')
                    raise P4Error(gm)
            if self.constraints.rooting:
                ret = self.constraints.areConsistentWithTreeRoot(self.tree)
                if not ret:
                    gm.append("The starting tree is not consistent with the constraints.cTree")
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
        self.chainTemp = 1.0
        self.chainTemps = []
        self.simTemp = None
        self.simTemp_doTuneTempsAndLogPi = True
        self.simTemp_doLogTemps = True
        self.simTemp_tempsLogFName = "mcmc_simTempTempsLog"
        self.simTemp_tempsFlob = None     # Open file

        # If we are doing simulated tempering ...
        if simTempNTemps:
            
            try:
                simTempNTemps = int(simTempNTemps)
            except (ValueError, TypeError):
                gm.append("If set, simTempNTemps should be an int, 2 or more.  Got %s" % simTempNTemps)
                raise P4Error(gm)
            if simTempNTemps < 2:
                gm.append("If set, simTempNTemps should be an int, 2 or more.  Got %s" % simTempNTemps)
                raise P4Error(gm)
            if self.nChains != 1:
                gm.append("If simTempNTemps is set, nChains should only be 1.  Got %i" % self.nChains)
            self.simTemp = SimTemp(self)
            self.simTemp.nTemps = simTempNTemps

            try:
                thisSimTempMax = float(simTempMax)
            except (ValueError, TypeError):
                gm.append("If doing simTemp, simTempMax should be set to a float more than zero.  Got  %s" % simTempMax)
                raise P4Error(gm)
            if thisSimTempMax <= 0.0:
                gm.append("If doing simTemp, simTempMax should be set to a float more than zero.  Got  %s" % simTempMax)
                raise P4Error(gm)

            #self.simTemp_maxTemp = thisSimTempMax

            if 0:
                # This makes the simTemp temperatures evenly spaced.
                stepSize = self.simTemp.maxTemp / (self.simTemp.nTemps - 1)
                self.simTemp.temps = [SimTempTemp(0.0)]
                for tNum in range(1,self.simTemp.nTemps):
                    tmp = SimTempTemp(tNum * stepSize)
                    self.simTemp.temps.append(tmp)
            if 1:
                # This makes more small temperatures and fewer big temperatures
                # The bigger the logBase, the more curvey the curve.  Smaller is more linear
                logBase = var.mcmc_simTemp_tempCurveLogBase
                factor = thisSimTempMax / (math.pow(logBase, self.simTemp.nTemps - 1) - 1.0)
                tmp = SimTempTemp()
                tmp.temp = 0.0
                tmp.tempNum = 0
                self.simTemp.temps = [tmp]
                for tNum in range(1,self.simTemp.nTemps):
                    tmp = SimTempTemp()
                    tmp.tempNum = tNum
                    tmp.temp = (math.pow(logBase, tNum) - 1.) * factor
                    self.simTemp.temps.append(tmp)

            self.simTemp.tNumSampleSize = var.mcmc_simTemp_tNumSampleSizePerTemp * self.simTemp.nTemps
            self.simTemp.setTempDiffsFromTemps()

            if 1:
                print("Initial settings of temperatures in Mcmc.__init__()")
                for tNum,tmp in enumerate(self.simTemp.temps):
                    print("%2i  %10.3f" % (tNum, tmp.temp))

            self.simTemp.tNumSample = deque([-1] * self.simTemp.tNumSampleSize)
        
            
        # Check that branch lengths are neither too short nor too long
        for n in self.tree.iterNodesNoRoot():
            if n.br.len < var.BRLEN_MIN:
                gm.append("Mcmc.__init__()  node %i brlen (%g)is too short." %
                          (n.nodeNum, n.br.len))
                raise P4Error(gm)
            elif n.br.len > var.BRLEN_MAX:
                gm.append("Mcmc.__init__()  node %i brlen (%f)is too long." %
                          (n.nodeNum, n.br.len))
                raise P4Error(gm)

        # Get the run number
        try:
            runNum = int(runNum)
        except (ValueError, TypeError):
            gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
            raise P4Error(gm)
        if runNum < 0:
            gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
            raise P4Error(gm)
        self.runNum = runNum

        self.loggerPrinter = None
        self.logger = None
        self._setLogger()

        # Check that we are not going to over-write good stuff
        ff = os.listdir(os.getcwd())
        hasPickle = False
        for fName in ff:
            if fName.startswith("mcmc_checkPoint_%i." % self.runNum):
                hasPickle = True
                break
        if hasPickle:
            gm.append("runNum is set to %i" % self.runNum)
            gm.append("There is at least one mcmc_checkPoint_%i.xxx file in this directory." % self.runNum)
            gm.append("This is a new Mcmc, and I am refusing to over-write exisiting files.")
            gm.append("Maybe you want to re-start from the latest mcmc_checkPoint_%i file?" % self.runNum)
            gm.append("Otherwise, get rid of the existing mcmc_xxx_%i.xxx files and start again." % self.runNum)
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
                        gm.append("(To get rid of this requirement, turn off var.strictRunNumberChecking.)")
                        raise P4Error(gm)

        self.sampleInterval = sampleInterval
        self.checkPointInterval = checkPointInterval

        self.props = Proposals()
        self.tunableProps = """allBrLens allCompsDir brLen compDir 
                    gdasrv local ndch2_internalCompsDir 
                    ndch2_internalCompsDirAlpha ndch2_leafCompsDir 
                    ndch2_leafCompsDirAlpha ndch2_priorRefCompDir
                    pInvar rMatrixDir allRMatricesDir ndrh2_priorRefRMatrixDir
                    ndrh2_internalRatesDir ndrh2_internalRatesDirAlpha ndrh2_leafRatesDir 
                    ndrh2_leafRatesDirAlpha relRate """.split()

        # Maybetunablebutnotnow  compLocation eTBR polytomy root3 root3n rMatrixLocation root2

        self.treePartitions = None
        self.likesFileName = "mcmc_likes_%i" % runNum
        self.treeFileName = "mcmc_trees_%i.nex" % runNum
        self.simFileName = "mcmc_sims_%i" % runNum
        self.pramsFileName = "mcmc_prams_%i" % runNum
        self.hypersFileName = "mcmc_hypers_%i" % runNum
        self.simTempFileName = "mcmc_simTemp_%i" % runNum
        self.ssLikesFileName = "mcmc_ssLikes_%i" % runNum
        # self.siteLikesFileName = "mcmc_siteLikes_%i" % runNum

        self.likesFile = None
        self.treeFile = None
        self.simFile = None
        self.pramsFile = None
        self.hypersFile = None
        self.ssLikesFile = None
        # self.siteLikesFile = None

        self.writePrams = writePrams
        # self.writeSiteLikes = False
        self.doCpo = False
        self.cpo_startGen = None
        self.cpo_sumsOfInverseSiteLikes = None
        self.cpo_nSamples = None
        self.writeHypers = True
        self.coldChainNum = -1

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
        self.simTree = None

        if self.nChains > 1:
            self.swapMatrix = []
            for i in range(self.nChains):
                self.swapMatrix.append([0] * self.nChains)
        else:
            self.swapMatrix = None

        self.swapTuner = None
        self.swapTunerDoTuning = True
        self.swapTunerSampleSize = 250
        self.stickyRootComp = False

        # check the tree, and tree+model+data
        if not aTree.taxNames:
            gm.append("The tree that you supply should have a 'taxNames' attribute.")
            gm.append("The taxNames should be in the same order as the data.")
            raise P4Error(gm)
        aTree.calcLogLike(verbose=False)

        if 0:
            # print complaintHead
            print("    logLike of the input tree is %s" % aTree.logLike)

        # Default tunings
        self._tunings = McmcTunings(self.tree.model.nParts)
        self.polytomyUseResolutionClassPrior = False
        self.polytomyPriorLogBigC = 0.0

        self.prob = McmcProposalProbs()

        # New tunings --- compDir and rMatrixDir, seem to depend on the dim.
        # And now allCompsDir
        for pNum in range(self._tunings.nParts):
            theDim = self.tree.model.parts[pNum].dim
            nRates = ((theDim * theDim) - theDim) / 2
            self._tunings.parts[pNum].default['compDir'] = 50. * theDim
            self._tunings.parts[pNum].default['allCompsDir'] = 100. * theDim
            self._tunings.parts[pNum].default['ndch2_leafCompsDir'] = 2000. * theDim
            self._tunings.parts[pNum].default['ndch2_internalCompsDir'] = 500. * theDim
            self._tunings.parts[pNum].default['ndch2_priorRefCompDir'] = 500. * theDim
            self._tunings.parts[pNum].default['rMatrixDir'] = 50. * nRates
            self._tunings.parts[pNum].default['allRMatricesDir'] = 100. * nRates
            self._tunings.parts[pNum].default['ndrh2_leafRatesDir'] = 50. * nRates
            self._tunings.parts[pNum].default['ndrh2_internalRatesDir'] = 50. * nRates
            self._tunings.parts[pNum].default['ndrh2_priorRefRMatrixDir'] = 50. * nRates
            

        # Zap internal node names
        for n in aTree.root.iterInternals():
            if n.name:
                n.name = None

        # If we need relRate, turn it on, and say so.
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

        nNodes = len(list(self.tree.iterNodes()))
        for pNum in range(self.tree.model.nParts):
            mp = self.tree.model.parts[pNum]
            dp = self.tree.data.parts[pNum]
            if mp.ndch2:
                if mp.nComps != nNodes:
                    gm.append("Model part %i, ndch2 is on, nNodes is %i, nComps is %i" % (
                        pNum, nNodes, mp.nComps))
                    gm.append("For ndch2 there should be one comp for each node.")
                    raise P4Error(gm)
                if mp.nGdasrvs > 1:
                    gm.append("Model part %i, ndch2 is on, nNodes is %i, nGdasrvs is %i" % (
                        pNum, nNodes, mp.nGdasrvs))
                    gm.append("This week, for ndch2 there should be only one gdasrv.")
                    raise P4Error(gm)

                for n in self.tree.iterNodes():
                    assert n.nodeNum == n.parts[pNum].compNum

                if not mp.ndch2_priorRefComp:
                    mp.ndch2_priorRefComp = numpy.array(dp.composition())

                mp.ndch2_priorRefComp = numpy.array(mp.ndch2_priorRefComp, dtype=float) # if not already a numpy.array
                assert len(mp.ndch2_priorRefComp) == mp.dim
                theMin = numpy.min(mp.ndch2_priorRefComp)
                if theMin < var.PIVEC_MIN:
                    gm.append("ndch2_priorRefComp min value is too low (%s)" % theMin)
                    gm.append("Set it 'by hand'")
                    raise P4Error(gm)
                mp.ndch2_priorRefComp /= theMin
                while mp.ndch2_priorRefComp.min()  < var.PIVEC_MIN:
                    for i in range(mp.dim):
                        if mp.ndch2_priorRefComp[i] < var.PIVEC_MIN:
                            mp.ndch2_priorRefComp[i] += (1.0 + (1.1 * random.random())) * var.PIVEC_MIN
                    mp.ndch2_priorRefComp /= mp.ndch2_priorRefComp.sum()

                # print(f"Mcmc init() mp.ndch2_priorRefComp is {mp.ndch2_priorRefComp}")

                # ususal comp proposals should not be on if we are doing ndch2
                self.prob.compDir = 0.0
                self.prob.allCompsDir = 0.0
                self.prob.compLocation = 0.0

            if mp.ndrh2:
                if mp.nRMatrices != nNodes - 1:
                    gm.append("Model part %i, ndrh2 is on, nNodes is %i, nRMatices is %i" % (
                        pNum, nNodes, mp.nRMatrices))
                    gm.append("For ndrh2 there should be one comp for each node except the root.")
                    raise P4Error(gm)

                for n in self.tree.iterNodesNoRoot():
                    assert n.nodeNum == n.br.parts[pNum].rMatrixNum + 1

                thisRMatrixLen = int(((mp.dim * mp.dim) - mp.dim) / 2)
                if not mp.ndrh2_priorRefRMatrix:
                    mp.ndrh2_priorRefRMatrix = numpy.ones(thisRMatrixLen)

                mp.ndrh2_priorRefRMatrix = numpy.array(mp.ndrh2_priorRefRMatrix, dtype=float) # if not already a numpy.array
                assert len(mp.ndrh2_priorRefRMatrix) == thisRMatrixLen
                mp.ndrh2_priorRefRMatrix /= numpy.min(mp.ndrh2_priorRefRMatrix)

                while mp.ndrh2_priorRefRMatrix.min()  < var.RATE_MIN:
                    for i in range(thisRMatrixLen):
                        if mp.ndrh2_priorRefRMatrix[i] < var.RATE_MIN:
                            mp.ndrh2_priorRefRMatrix[i] += (1.0 + (1.1 * random.random())) * var.RATE_MIN
                    mp.ndrh2_priorRefRMatrix /= mp.ndrh2_priorRefRMatrix.sum()
                # print(f"Mcmc init() mp.ndrh2_priorRefRMatrix is {mp.ndrh2_priorRefRMatrix}")
                

                # Turn off non-ndrh2 rmatrix moves
                self.prob.rMatrixDir = 0.0
                self.prob.allRMatricesDir = 0.0
                self.prob.rMatrixLocation = 0.0

        if self.tree.model.isHet:
            props_on = []
            if verbose:
                print("\nInitiating across-tree heterogeneous model...")

            self.tree.setModelComponentsNNodes()

            if self.verbose:
                self.tree.summarizeModelComponentsNNodes()

            for pNum in range(self.tree.model.nParts):
                mp = self.tree.model.parts[pNum]

                if mp.nComps > 1 and mp.nComps < nNodes:
                    self.prob.compLocation = 1.0
                    if "composition location" not in props_on:
                        props_on.append("composition location")

                if mp.nRMatrices > 1 and mp.nRMatrices < (nNodes - 1):
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

                    self.prob.ndch2_priorRefCompDir = 1.0
                    thisMString = "ndch2_priorRefCompDir"
                    if thisMString not in props_on:
                        props_on.append(thisMString)

                if mp.ndrh2:
                    self.prob.ndrh2_leafRatesDir = 1.0
                    thisMString = "ndrh2_leafRatesDir"
                    if thisMString not in props_on:
                        props_on.append(thisMString)

                    self.prob.ndrh2_internalRatesDir = 1.0
                    thisMString = "ndrh2_internalRatesDir"
                    if thisMString not in props_on:
                        props_on.append(thisMString)

                    self.prob.ndrh2_leafRatesDirAlpha = 1.0
                    thisMString = "ndrh2_leafRatesDirAlpha"
                    if thisMString not in props_on:
                        props_on.append(thisMString)

                    self.prob.ndrh2_internalRatesDirAlpha = 1.0
                    thisMString = "ndrh2_internalRatesDirAlpha"
                    if thisMString not in props_on:
                        props_on.append(thisMString)

                    self.prob.ndrh2_priorRefRMatrixDir = 1.0
                    thisMString = "ndrh2_priorRefRMatrixDir"
                    if thisMString not in props_on:
                        props_on.append(thisMString)

            if self.constraints and self.constraints.rooting: # a root constraint
                self.prob.root2 = 0.0
                self.prob.root3 = 0.0
                self.prob.root3n = 0.0
                print("There is a constraint on the root position, so turning root proposals off.")
            else:
                print("root proposals are now turned off by default")
                self.prob.root2 = 0.0
                self.prob.root3 = 0.0
                self.prob.root3n = 0.0

                # if self.isBiRoot:
                #     self.prob.root2 = 1.0
                #     thisMString = "root2 (root location)"
                #     if thisMString not in props_on:
                #         props_on.append(thisMString)
                # else:
                #     self.prob.root3 = 1.0
                #     thisMString = "root3 (root location)"
                #     if thisMString not in props_on:
                #         props_on.append(thisMString)
                #     self.prob.root3n = 1.0
                #     thisMString = "root3n (root location, neighbours)"
                #     if thisMString not in props_on:
                #         props_on.append(thisMString)
            
            if verbose:
                print("\n%23s" % "Additional proposals:")
                for prop in props_on:
                    print("%30s = on" % prop)
                print("\n  %s" % "[All set to a relative proposal probability of 1.0. You can adjust those,")
                print("  %s" % " or turn them off, eg by setting m.prob.theProposal=0.0, before Mcmc.run()]\n")
        else:
            self.prob.root3 = 0.0
            self.prob.root3n = 0.0
            self.prob.root2 = 0.0
            self.prob.compLocation = 0.0
            self.prob.rMatrixLocation = 0.0
            #self.prob.gdasrvLocation = 0.0

        splash = p4.func.splash2(verbose=False)
        for aLine in splash:
            print(aLine)
            self.logger.info(aLine)

        # Check the data for blank sequences, partition by partition
        hasBlanks = False
        blankSeqNums = []
        for partNum in range(self.tree.data.nParts):
            p = self.tree.data.parts[partNum]
            partBlankSeqNums = []
            for taxNum in range(self.tree.data.nTax):
                nSites = pf.partSequenceSitesCount(p.cPart, taxNum)  # no gaps, no missings
                # print(f"Mcmc.__init__()  partNum {partNum} taxNum {taxNum} nSites {nSites}")
                if not nSites:
                    partBlankSeqNums.append(taxNum)
            if partBlankSeqNums:
                hasBlanks = True
            blankSeqNums.append(partBlankSeqNums)
        if hasBlanks:
            self.logger.info("Blank sequences were found.  For each partition, the sequence numbers are ---")
            self.logger.info(f"{blankSeqNums}")
            self.logger.info("These will be skipped in bigXSq simulations.")
            self.blankSeqNums = blankSeqNums
        else:
            self.blankSeqNums = None

        self.swapVector = var.mcmc_swapVector   # Usually True
        if self.swapVector:
            if self.nChains > 1:
                print("%16s: %s" % ('swapVector', "on"))
                self.swapTuner = SwapTunerV(self)
                print("%16s: %s" % ('swapTuner', "on"))
        else:
            if self.nChains > 1:
                print("%16s: %s" % ('swapVector', "off"))


        # Hidden experimental hacking
        self.doHeatingHack = False
        self.heatingHackTemperature = 5.0
        self.originalHeatingHackTemperature = 5.0
        self.heatingHackRamper = None

        # Whether logging from the Pf module is turned on.
        # When it is turned on, a callback is set up to self._logFromPfModule()
        self.setupPfLogging = False

        # Whether to do root3 and root3n tuning or not
        # self.doRoot3Tuning = False
        # self.doRoot3nTuning = False

        if not var.gsl_rng:
            var.gsl_rng = pf.gsl_rng_get()
            pf.gsl_rng_set(var.gsl_rng, int(time.time()))
        the_gsl_rng_size = pf.gsl_rng_size(var.gsl_rng) # size of the state
        # A place to store the state.  Empty to start.  It is stored during a checkpoint.
        self.gsl_rng_state_ndarray = numpy.array(['0'] * the_gsl_rng_size, numpy.dtype('B'))  # B is unsigned byte

        self.randomState = None    # for python random library

        # For re-starts, if these values have changed or have been forgotten, allow checking and warning in Mcmc.run()
        self.init_PIVEC_MIN = var.PIVEC_MIN
        self.init_RATE_MIN = var.RATE_MIN
        self.init_BRLEN_MIN = var.BRLEN_MIN
        self.init_GAMMA_SHAPE_MIN = var.GAMMA_SHAPE_MIN

        # For stepping stone power posteriors
        self.doSteppingStone = False
        self.ssBeta = None

    def _setLogger(self):
        """Make a logger."""

        # Need to reset, or else basicConfig will be ignored
        log = logging.getLogger()  # root logger
        for hdlr in log.handlers[:]:  # remove all old handlers
            log.removeHandler(hdlr)

        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s %(message)s',
                            datefmt='[%Y-%m-%d %H:%M]',
                            filename="mcmc_log_%i" % self.runNum,
                            filemode='a')

        # This logger only logs to the file, not to stderr.
        self.logger = logging.getLogger("logFileOnly")
        
    def _logFromPfModule(self, treeAddress, message):
        
        tinfo = None
        for chNum,ch in enumerate(self.chains):
            if ch.curTree.cTree == treeAddress:
                tinfo = "gen %i, chain %i, curTree; " % (self.gen, chNum)
            elif ch.propTree.cTree == treeAddress:
                tinfo = "gen %i, chain %i, propTree; " % (self.gen, chNum)
        if not tinfo:
            tinfo = "unknown tree; "
        message = tinfo + message
        #print(treeAddress, message)
        self.logger.info(message)
        

    def _makeProposals(self):
        """Make proposals for the mcmc."""

        gm = ['Mcmc._makeProposals()']

        # The weight determines how often proposals are made.  The
        # weight is the product of 4 things:
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
        #    3.  In the NDCH or NDRH model, the number of comps or rMatrices in
        #        the partition.  So two comps doubles the weight compared to one
        #        comp.
        #
        #    4.  A fudge factor.  Arbitrary witchcraft.
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
            p.tuning = [self._tunings.default[p.name]] * self.nChains
            p.brLenPriorType = self._tunings.default['brLenPriorType']
            p.brLenPriorLambda = self._tunings.default['brLenPriorLambda']
            p.weight = self.prob.brLen * \
                (len(self.tree.nodes) - 1) * fudgeFactor['brLen']

            p.tnAccVeryHi = 0.7
            p.tnAccHi = 0.6
            p.tnAccLo = 0.1
            p.tnAccVeryLo = 0.05

            p.tnFactorVeryHi = 1.6
            p.tnFactorHi = 1.2
            p.tnFactorLo = 0.8
            p.tnFactorVeryLo = 0.7

            self.props.proposals.append(p)

        # allBrLens
        if self.prob.allBrLens:
            p = Proposal(self)
            p.name = 'allBrLens'
            p.tuning = [self._tunings.default[p.name]] * self.nChains
            p.tuningLimitHi = 0.2
            p.brLenPriorType = self._tunings.default['brLenPriorType']
            p.brLenPriorLambda = self._tunings.default['brLenPriorLambda']
            p.weight = self.prob.allBrLens * \
                (len(self.tree.nodes) - 1) * fudgeFactor['allBrLens']

            p.tnAccVeryHi = 0.4
            p.tnAccHi = 0.15
            p.tnAccLo = 0.05
            p.tnAccVeryLo = 0.03

            p.tnFactorVeryHi = 1.6
            p.tnFactorHi = 1.2
            p.tnFactorLo = 0.8
            p.tnFactorVeryLo = 0.7

            self.props.proposals.append(p)



        # eTBR
        if self.prob.eTBR:
            p = Proposal(self)
            p.name = 'eTBR'
            p.etbrLambda = self._tunings.default['etbrLambda']
            p.etbrPExt = self._tunings.default['etbrPExt']
            p.brLenPriorType = self._tunings.default['brLenPriorType']
            p.brLenPriorLambda = self._tunings.default['brLenPriorLambda']
            p.weight = self.prob.eTBR * \
                (len(self.tree.nodes) - 1) * fudgeFactor['eTBR']
            self.props.proposals.append(p)

        # local
        if self.prob.local:
            p = Proposal(self)
            p.name = 'local'
            p.tuning = [self._tunings.default[p.name]] * self.nChains
            p.brLenPriorType = self._tunings.default['brLenPriorType']
            p.brLenPriorLambda = self._tunings.default['brLenPriorLambda']
            p.weight = self.prob.local * \
                (len(self.tree.nodes) - 1) * fudgeFactor['local']

            p.tnAccVeryHi = 0.7
            p.tnAccHi = 0.6
            p.tnAccLo = 0.05
            p.tnAccVeryLo = 0.03

            p.tnFactorVeryHi = 1.6
            p.tnFactorHi = 1.2
            p.tnFactorLo = 0.8
            p.tnFactorVeryLo = 0.7

            self.props.proposals.append(p)

        # polytomy
        if self.prob.polytomy:
            p = Proposal(self)
            p.name = 'polytomy'
            p.brLenPriorType = self._tunings.default['brLenPriorType']
            p.brLenPriorLambda = self._tunings.default['brLenPriorLambda']
            p.polytomyUseResolutionClassPrior = self.polytomyUseResolutionClassPrior
            p.polytomyPriorLogBigC = self.polytomyPriorLogBigC
            p.weight = self.prob.polytomy * \
                (len(self.tree.nodes) - 1) * fudgeFactor['polytomy']
            self.props.proposals.append(p)

        # root3
        if self.prob.root3:
            if len(self.tree.taxNames) <= 3:
                pass
            else:
                p = Proposal(self)
                p.name = 'root3'
                #p.tuning = [self._tunings.default[p.name]] * self.nChains
                p.weight = self.prob.root3 * \
                           self.tree.nInternalNodes * fudgeFactor['root3']

                p.tnAccVeryHi = 0.4
                p.tnAccHi = 0.3
                p.tnAccLo = 0.05
                p.tnAccVeryLo = 0.03

                p.tnFactorVeryHi = 0.7
                p.tnFactorHi = 0.8
                p.tnFactorLo = 1.2
                p.tnFactorVeryLo = 1.6

                self.props.proposals.append(p)

        # root3n
        if self.prob.root3n:
            if len(self.tree.taxNames) <= 3:
                pass
            else:
                p = Proposal(self)
                p.name = 'root3n'
                #p.tuning = [self._tunings.default[p.name]] * self.nChains
                p.weight = self.prob.root3n * \
                           self.tree.nInternalNodes * fudgeFactor['root3n']

                p.tnAccVeryHi = 0.5
                p.tnAccHi = 0.4
                p.tnAccLo = 0.1
                p.tnAccVeryLo = 0.06

                p.tnFactorVeryHi = 0.7
                p.tnFactorHi = 0.8
                p.tnFactorLo = 1.2
                p.tnFactorVeryLo = 1.6

                self.props.proposals.append(p)

        # root2
        if self.prob.root2:
            if len(self.tree.taxNames) <= 3:
                pass
            else:
                p = Proposal(self)
                p.name = 'root2'
                p.tuning = [self._tunings.default[p.name]] * self.nChains
                p.weight = self.prob.root2 * self.tree.nInternalNodes * fudgeFactor['root3n']

                p.tnAccVeryHi = 0.5
                p.tnAccHi = 0.4
                p.tnAccLo = 0.1
                p.tnAccVeryLo = 0.06

                p.tnFactorVeryHi = 1.6
                p.tnFactorHi = 1.2
                p.tnFactorLo = 0.8
                p.tnFactorVeryLo = 0.7

                self.props.proposals.append(p)

        # relRate
        if self.prob.relRate:
            if self.tree.model.nParts == 1:
                pass
            if self.tree.model.doRelRates and self.tree.model.relRatesAreFree:
                p = Proposal(self)
                p.name = 'relRate'
                p.tuning = [self._tunings.default[p.name]] * self.nChains
                p.weight = self.prob.relRate * self.tree.model.nParts

                p.tnAccVeryHi = 0.7
                p.tnAccHi = 0.6
                p.tnAccLo = 0.1
                p.tnAccVeryLo = 0.06

                p.tnFactorVeryHi = 1.6
                p.tnFactorHi = 1.2
                p.tnFactorLo = 0.8
                p.tnFactorVeryLo = 0.7

                self.props.proposals.append(p)

        # comp, rMatrix, gdasrv, pInvar, modelComponentLocations ...
        for pNum in range(self.tree.model.nParts):
            mp = self.tree.model.parts[pNum]

            # # comp
            # if self.prob.comp:
            #     if mp.comps and mp.comps[0].free:
            #         for mtNum in range(mp.nComps):
            #             assert mp.comps[mtNum].free
            #         p = Proposal(self)
            #         p.name = 'comp'
            #         p.tuning = self._tunings.parts[pNum].default[p.name]
            #         p.weight = self.prob.comp * (mp.dim - 1) * mp.nComps
            #         p.pNum = pNum

            #         p.tnAccVeryHi = 0.7
            #         p.tnAccHi = 0.6
            #         p.tnAccLo = 0.1
            #         p.tnAccVeryLo = 0.06

            #         p.tnFactorVeryHi = 1.6
            #         p.tnFactorHi = 1.2
            #         p.tnFactorLo = 0.8
            #         p.tnFactorVeryLo = 0.7

            #         self.props.proposals.append(p)

            # compDir
            if self.prob.compDir:
                if mp.comps and mp.comps[0].free:
                    for mtNum in range(mp.nComps):
                        assert mp.comps[mtNum].free
                    p = Proposal(self)
                    p.name = 'compDir'
                    p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                    p.weight = self.prob.compDir * (mp.dim - 1) * mp.nComps
                    p.pNum = pNum

                    p.tnAccVeryHi = 0.4
                    p.tnAccHi = 0.15
                    p.tnAccLo = 0.05
                    p.tnAccVeryLo = 0.03

                    p.tnFactorVeryHi = 0.7
                    p.tnFactorHi = 0.8
                    p.tnFactorLo = 1.2
                    p.tnFactorVeryLo = 1.6

                    self.props.proposals.append(p)

            # allCompsDir
            if self.prob.allCompsDir:
                if mp.comps and mp.comps[0].free:
                    for mtNum in range(mp.nComps):
                        assert mp.comps[mtNum].free
                    p = Proposal(self)
                    p.name = 'allCompsDir'
                    p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                    p.weight = self.prob.allCompsDir * (mp.dim - 1) * mp.nComps * fudgeFactor['allCompsDir']
                    p.pNum = pNum

                    p.tnAccVeryHi = 0.4
                    p.tnAccHi = 0.15
                    p.tnAccLo = 0.05
                    p.tnAccVeryLo = 0.03

                    p.tnFactorVeryHi = 0.7
                    p.tnFactorHi = 0.8
                    p.tnFactorLo = 1.2
                    p.tnFactorVeryLo = 1.4

                    self.props.proposals.append(p)

            # ndch2_leafCompsDir
            if self.prob.ndch2_leafCompsDir:
                for mtNum in range(mp.nComps):
                    assert mp.comps[mtNum].free
                p = Proposal(self)
                p.name = 'ndch2_leafCompsDir'
                p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                p.weight = self.prob.ndch2_leafCompsDir * (mp.dim - 1) * mp.nComps * fudgeFactor['ndch2comp']
                p.pNum = pNum

                p.tnAccVeryHi = 0.4
                p.tnAccHi = 0.15
                p.tnAccLo = 0.05
                p.tnAccVeryLo = 0.03

                p.tnFactorVeryHi = 0.7
                p.tnFactorHi = 0.8
                p.tnFactorLo = 1.2
                p.tnFactorVeryLo = 1.4

                self.props.proposals.append(p)

            # ndch2_internalCompsDir
            if self.prob.ndch2_internalCompsDir:
                for mtNum in range(mp.nComps):
                    assert mp.comps[mtNum].free
                p = Proposal(self)
                p.name = 'ndch2_internalCompsDir'
                p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                p.weight = self.prob.ndch2_internalCompsDir * (mp.dim - 1) * mp.nComps * fudgeFactor['ndch2comp']
                p.pNum = pNum

                p.tnAccVeryHi = 0.4
                p.tnAccHi = 0.15
                p.tnAccLo = 0.05
                p.tnAccVeryLo = 0.03

                p.tnFactorVeryHi = 0.7
                p.tnFactorHi = 0.8
                p.tnFactorLo = 1.2
                p.tnFactorVeryLo = 1.4

                self.props.proposals.append(p)


            # ndch2_leafCompsDirAlpha
            if self.prob.ndch2_leafCompsDirAlpha:
                p = Proposal(self)
                p.name = 'ndch2_leafCompsDirAlpha'
                p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                p.weight = self.prob.ndch2_leafCompsDirAlpha * (mp.dim - 1) * mp.nComps * fudgeFactor['ndch2alpha']
                p.pNum = pNum

                p.tnAccVeryHi = 0.7
                p.tnAccHi = 0.6
                p.tnAccLo = 0.1
                p.tnAccVeryLo = 0.06

                p.tnFactorVeryHi = 1.6
                p.tnFactorHi = 1.2
                p.tnFactorLo = 0.8
                p.tnFactorVeryLo = 0.7

                self.props.proposals.append(p)

            # ndch2_internalCompsDirAlpha
            if self.prob.ndch2_internalCompsDirAlpha:
                p = Proposal(self)
                p.name = 'ndch2_internalCompsDirAlpha'
                p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                p.weight = self.prob.ndch2_internalCompsDirAlpha * (mp.dim - 1) * mp.nComps * fudgeFactor['ndch2alpha']
                p.pNum = pNum

                p.tnAccVeryHi = 0.7
                p.tnAccHi = 0.6
                p.tnAccLo = 0.1
                p.tnAccVeryLo = 0.06

                p.tnFactorVeryHi = 1.6
                p.tnFactorHi = 1.2
                p.tnFactorLo = 0.8
                p.tnFactorVeryLo = 0.7

                self.props.proposals.append(p)

            # ndch2_priorRefCompDir
            if self.prob.ndch2_priorRefCompDir:
                p = Proposal(self)
                p.name = 'ndch2_priorRefCompDir'
                p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                p.weight = self.prob.ndch2_priorRefCompDir * (mp.dim - 1)  * fudgeFactor['ndch2priorRefComp']
                p.pNum = pNum

                p.tnAccVeryHi = 0.4
                p.tnAccHi = 0.15
                p.tnAccLo = 0.05
                p.tnAccVeryLo = 0.03

                p.tnFactorVeryHi = 0.7
                p.tnFactorHi = 0.8
                p.tnFactorLo = 1.2
                p.tnFactorVeryLo = 1.4

                self.props.proposals.append(p)




            # rMatrixDir
            if self.prob.rMatrixDir:
                if mp.rMatrices and mp.rMatrices[0].free:
                    for mtNum in range(mp.nRMatrices):
                        assert mp.rMatrices[mtNum].free
                    p = Proposal(self)
                    p.name = 'rMatrixDir'
                    if mp.rMatrices[mtNum].spec == '2p':
                        p.weight = self.prob.rMatrixDir         # Is this enough?
                        p.variant = '2p'
                        p.tuning = [self._tunings.parts[pNum].default['twoP']] * self.nChains
                    else:
                        p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                        p.weight = self.prob.rMatrixDir * mp.nRMatrices * \
                            ((((mp.dim * mp.dim) - mp.dim) / 2) - 1)
                    p.pNum = pNum

                    p.tnAccVeryHi = 0.4
                    p.tnAccHi = 0.15
                    p.tnAccLo = 0.05
                    p.tnAccVeryLo = 0.03

                    p.tnFactorVeryHi = 0.7
                    p.tnFactorHi = 0.8
                    p.tnFactorLo = 1.2
                    p.tnFactorVeryLo = 1.4


                    self.props.proposals.append(p)

            # allRMatricesDir
            if self.prob.allRMatricesDir:
                if mp.rMatrices and mp.rMatrices[0].free:
                    for mtNum in range(mp.nRMatrices):
                        assert mp.rMatrices[mtNum].free
                    p = Proposal(self)
                    p.name = 'allRMatricesDir'
                    if mp.rMatrices[mtNum].spec == '2p':
                        raise P4Error("Proposal allRMatricesDir is not yet working for 2p.  Fixme.")
                    else:
                        p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                        p.weight = self.prob.allRMatricesDir * mp.nRMatrices * \
                            ((((mp.dim * mp.dim) - mp.dim) / 2) - 1) * fudgeFactor['allRMatricesDir']
                    p.pNum = pNum

                    p.tnAccVeryHi = 0.4
                    p.tnAccHi = 0.15
                    p.tnAccLo = 0.05
                    p.tnAccVeryLo = 0.03

                    p.tnFactorVeryHi = 0.7
                    p.tnFactorHi = 0.8
                    p.tnFactorLo = 1.2
                    p.tnFactorVeryLo = 1.4

                    self.props.proposals.append(p)


            # ndrh2_leafRatesDir
            if self.prob.ndrh2_leafRatesDir:
                for mtNum in range(mp.nRMatrices):
                    assert mp.rMatrices[mtNum].free
                p = Proposal(self)
                p.name = 'ndrh2_leafRatesDir'
                p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                p.weight = self.prob.ndrh2_leafRatesDir * \
                    ((((mp.dim * mp.dim) - mp.dim) / 2.) - 1.) * \
                    mp.nRMatrices * fudgeFactor['ndrh2rates']
                p.pNum = pNum

                p.tnAccVeryHi = 0.4
                p.tnAccHi = 0.15
                p.tnAccLo = 0.05
                p.tnAccVeryLo = 0.03

                p.tnFactorVeryHi = 0.7
                p.tnFactorHi = 0.8
                p.tnFactorLo = 1.2
                p.tnFactorVeryLo = 1.4

                self.props.proposals.append(p)

            # ndrh2_internalRatesDir
            if self.prob.ndrh2_internalRatesDir:
                for mtNum in range(mp.nRMatrices):
                    assert mp.rMatrices[mtNum].free
                p = Proposal(self)
                p.name = 'ndrh2_internalRatesDir'
                p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                p.weight = self.prob.ndrh2_internalRatesDir * \
                    ((((mp.dim * mp.dim) - mp.dim) / 2.) - 1.) * \
                    mp.nRMatrices * fudgeFactor['ndrh2rates']
                p.pNum = pNum

                p.tnAccVeryHi = 0.4
                p.tnAccHi = 0.15
                p.tnAccLo = 0.05
                p.tnAccVeryLo = 0.03

                p.tnFactorVeryHi = 0.7
                p.tnFactorHi = 0.8
                p.tnFactorLo = 1.2
                p.tnFactorVeryLo = 1.4

                self.props.proposals.append(p)

            # ndrh2_leafRatesDirAlpha
            if self.prob.ndrh2_leafRatesDirAlpha:
                p = Proposal(self)
                p.name = 'ndrh2_leafRatesDirAlpha'
                p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                p.weight = self.prob.ndrh2_leafRatesDirAlpha * \
                    ((((mp.dim * mp.dim) - mp.dim) / 2.) - 1.) * \
                    mp.nRMatrices * fudgeFactor['ndrh2alpha']
                p.pNum = pNum

                p.tnAccVeryHi = 0.7
                p.tnAccHi = 0.6
                p.tnAccLo = 0.1
                p.tnAccVeryLo = 0.06

                p.tnFactorVeryHi = 1.6
                p.tnFactorHi = 1.2
                p.tnFactorLo = 0.8
                p.tnFactorVeryLo = 0.7

                self.props.proposals.append(p)

            # ndrh2_internalRatesDirAlpha
            if self.prob.ndrh2_internalRatesDirAlpha:
                p = Proposal(self)
                p.name = 'ndrh2_internalRatesDirAlpha'
                p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                p.weight = self.prob.ndrh2_internalRatesDirAlpha * \
                    ((((mp.dim * mp.dim) - mp.dim) / 2.) - 1.) * \
                    mp.nRMatrices * fudgeFactor['ndrh2alpha']
                p.pNum = pNum

                p.tnAccVeryHi = 0.7
                p.tnAccHi = 0.6
                p.tnAccLo = 0.1
                p.tnAccVeryLo = 0.06

                p.tnFactorVeryHi = 1.6
                p.tnFactorHi = 1.2
                p.tnFactorLo = 0.8
                p.tnFactorVeryLo = 0.7

                self.props.proposals.append(p)

            # ndrh2_priorRefRMatrixDir
            if self.prob.ndrh2_priorRefRMatrixDir:
                for mtNum in range(mp.nRMatrices):
                    assert mp.rMatrices[mtNum].free
                p = Proposal(self)
                p.name = 'ndrh2_priorRefRMatrixDir'
                p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                p.weight = self.prob.ndrh2_priorRefRMatrixDir * \
                    ((((mp.dim * mp.dim) - mp.dim) / 2.) - 1.) * \
                    fudgeFactor['ndrh2priorRefRMatrix']
                p.pNum = pNum

                p.tnAccVeryHi = 0.4
                p.tnAccHi = 0.15
                p.tnAccLo = 0.05
                p.tnAccVeryLo = 0.03

                p.tnFactorVeryHi = 0.7
                p.tnFactorHi = 0.8
                p.tnFactorLo = 1.2
                p.tnFactorVeryLo = 1.4

                self.props.proposals.append(p)



            # gdasrv
            if self.prob.gdasrv:
                if mp.nGdasrvs and mp.gdasrvs[0].free:
                    assert mp.nGdasrvs == 1
                    p = Proposal(self)
                    p.name = 'gdasrv'
                    p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                    p.weight = self.prob.gdasrv
                    p.pNum = pNum

                    p.tnAccVeryHi = 0.7
                    p.tnAccHi = 0.6
                    p.tnAccLo = 0.15
                    p.tnAccVeryLo = 0.10

                    p.tnFactorVeryHi = 2.0
                    p.tnFactorHi = 1.5
                    p.tnFactorLo = 0.75
                    p.tnFactorVeryLo = 0.5

                    self.props.proposals.append(p)

            # pInvar
            if self.prob.pInvar:
                if mp.pInvar.free:
                    p = Proposal(self)
                    p.name = 'pInvar'
                    p.tuning = [self._tunings.parts[pNum].default[p.name]] * self.nChains
                    p.weight = self.prob.pInvar
                    p.pNum = pNum

                    p.tnAccVeryHi = 0.7
                    p.tnAccHi = 0.6
                    p.tnAccLo = 0.1
                    p.tnAccVeryLo = 0.06

                    p.tnFactorVeryHi = 1.6
                    p.tnFactorHi = 1.2
                    p.tnFactorLo = 0.8
                    p.tnFactorVeryLo = 0.7

                    self.props.proposals.append(p)

            # compLocation
            if self.prob.compLocation:
                if mp.nComps > 1:
                    p = Proposal(self)
                    p.name = 'compLocation'
                    p.tuning = None #self._tunings.parts[pNum].default[p.name]
                    p.weight = self.prob.compLocation * mp.nComps * len(self.tree.nodes) * \
                        fudgeFactor['compLocation']
                    p.pNum = pNum
                    self.props.proposals.append(p)

            # rMatrixLocation
            if self.prob.rMatrixLocation:
                if mp.nRMatrices > 1:
                    p = Proposal(self)
                    p.name = 'rMatrixLocation'
                    p.tuning = None #self._tunings.parts[pNum].default[p.name]
                    p.weight = self.prob.rMatrixLocation * mp.nRMatrices * (len(self.tree.nodes) - 1) * \
                        fudgeFactor['rMatrixLocation']
                    p.pNum = pNum
                    self.props.proposals.append(p)

            # # gdasrvLocation
            # if self.prob.gdasrvLocation:
            #     if mp.nGdasrvs > 1:
            #         p = Proposal(self)
            #         p.name = 'gdasrvLocation'
            #         p.weight = self.prob.gdasrvLocation * mp.nGdasrvs * (len(self.tree.nodes) - 1) * \
            #             fudgeFactor['gdasrvLocation']
            #         p.pNum = pNum
            #         self.props.proposals.append(p)


        if not self.props.proposals:
            gm.append("No proposals?")
            raise P4Error(gm)
        for p in self.props.proposals:
            if p.name in ['local', 'eTBR', 'polytomy']:
                self.props.topologyProposalsDict[p.name] = p
        self.props.calculateWeights()

    def _refreshProposalProbs(self):
        """Adjust proposals after a restart."""

        gm = ['Mcmc._refreshProposalProbs()']

        for p in self.props.proposals:
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
                p.weight = self.prob.root3 * \
                    self.tree.nInternalNodes * fudgeFactor['root3']

            # root3n
            if p.name == 'root3n':
                p.weight = self.prob.root3n * \
                    self.tree.nInternalNodes * fudgeFactor['root3n']

            # relRate
            if p.name == 'relRate':
                p.weight = self.prob.relRate * self.tree.model.nParts

            # # comp
            # if p.name == 'comp':
            #     mp = self.tree.model.parts[p.pNum]
            #     p.weight = self.prob.comp * float(mp.dim - 1) * mp.nComps

            # compDir
            if p.name == 'compDir':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.compDir * (mp.dim - 1) * mp.nComps

            # allCompsDir
            if p.name == 'allCompsDir':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.allCompsDir * (mp.dim - 1) * mp.nComps * fudgeFactor['allCompsDir']

            # ndch2_leafCompsDir
            if p.name == 'ndch2_leafCompsDir':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.ndch2_leafCompsDir * (mp.dim - 1) * mp.nComps * fudgeFactor['ndch2comp']

            # ndch2_internalCompsDir
            if p.name == 'ndch2_internalCompsDir':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.ndch2_internalCompsDir * (mp.dim - 1) * mp.nComps * fudgeFactor['ndch2comp']

            # ndch2_leafCompsDirAlpha
            if p.name == 'ndch2_leafCompsDirAlpha':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.ndch2_leafCompsDirAlpha * (mp.dim - 1) * mp.nComps * fudgeFactor['ndch2alpha']

            # ndch2_internalCompsDirAlpha
            if p.name == 'ndch2_internalCompsDirAlpha':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.ndch2_internalCompsDirAlpha * (mp.dim - 1) * mp.nComps * fudgeFactor['ndch2alpha']

            # ndch2_priorRefCompDir
            if p.name == 'ndch2_priorRefCompDir':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.ndch2_priorRefCompDir * (mp.dim - 1)  * fudgeFactor['ndch2priorRefComp']

            # # rMatrix
            # if p.name == 'rMatrix':
            #     mp = self.tree.model.parts[p.pNum]
            #     if mp.rMatrices[0].spec == '2p':
            #         p.weight = self.prob.rMatrix
            #     else:
            #         p.weight = self.prob.rMatrix * mp.nRMatrices * \
            #             ((((mp.dim * mp.dim) - mp.dim) / 2) - 1)

            # rMatrixDir
            if p.name == 'rMatrixDir':
                mp = self.tree.model.parts[p.pNum]
                if mp.rMatrices[0].spec == '2p':
                    p.weight = self.prob.rMatrixDir
                else:
                    p.weight = self.prob.rMatrixDir * mp.nRMatrices * \
                        ((((mp.dim * mp.dim) - mp.dim) / 2.) - 1)

            # allRMatricesDir
            if p.name == 'allRMatricesDir':
                mp = self.tree.model.parts[p.pNum]
                #if mp.rMatrices[0].spec == '2p':
                #    p.weight = self.prob.rMatrixDir
                #else:
                p.weight = self.prob.allRMatricesDir * mp.nRMatrices * \
                        ((((mp.dim * mp.dim) - mp.dim) / 2.) - 1.) * fudgeFactor['allRMatricesDir']

            # ndrh2_leafRatesDir
            if p.name == 'ndrh2_leafRatesDir':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.ndrh2_leafRatesDir * \
                    ((((mp.dim * mp.dim) - mp.dim) / 2.) - 1.) * \
                    mp.nRMatrices * fudgeFactor['ndrh2rates']

            # ndrh2_internalRatesDir
            if p.name == 'ndrh2_internalRatesDir':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.ndrh2_internalRatesDir * \
                    ((((mp.dim * mp.dim) - mp.dim) / 2.) - 1.) * \
                    mp.nRMatrices * fudgeFactor['ndrh2rates']

            # ndrh2_leafRatesDirAlpha
            if p.name == 'ndrh2_leafRatesDirAlpha':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.ndrh2_leafRatesDirAlpha * \
                    ((((mp.dim * mp.dim) - mp.dim) / 2.) - 1.) * \
                    mp.nRMatrices * fudgeFactor['ndrh2alpha']

            # ndrh2_internalRatesDirAlpha
            if p.name == 'ndrh2_internalRatesDirAlpha':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.ndrh2_internalRatesDirAlpha * \
                    ((((mp.dim * mp.dim) - mp.dim) / 2.) - 1.) * \
                    mp.nRMatrices * fudgeFactor['ndrh2alpha']

            # ndrh2_priorRefRMatrix
            if p.name == 'ndrh2_priorRefRMatrixDir':
                mp = self.tree.model.parts[p.pNum]
                p.weight = self.prob.ndrh2_priorRefRMatrixDir * \
                    ((((mp.dim * mp.dim) - mp.dim) / 2.) - 1.) * \
                    fudgeFactor['ndrh2priorRefRMatrix']




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

            # # gdasrvLocation
            # if p.name == 'gdasrvLocation':
            #     mp = self.tree.model.parts[p.pNum]
            #     p.weight = self.prob.gdasrvLocation * mp.nGdasrvs * (len(self.tree.nodes) - 1) * \
            #         fudgeFactor['gdasrvLocation']


        self.props.calculateWeights()


    def writeProposalAcceptances(self):
        """Pretty-print the proposal acceptances."""

        if (self.gen - self.startMinusOne) <= 0:
            print("\nwriteProposalAcceptances()  There is no info in memory. ")
            print(" Maybe it was just emptied after writing to a checkpoint?  ")
            print("If so, read the checkPoint and get the proposalAcceptances from there.")
            return

        spacer = ' ' * 8
        print("\nProposal acceptances, run %i, for %i gens, from gens %i to %i, inclusive." % (
            self.runNum, (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen))
        print("%s %30s %5s %10s %13s%10s" % (spacer, 'proposal', 'part', 'nProposals', 'acceptance(%)', 'tuning'))
        for p in self.props.proposals:
            print("%s" % spacer, end=' ')
            print("%30s" % p.name, end=' ')
            if p.pNum != -1:
                print(" %3i " % p.pNum, end=' ')
            else:
                print("   - ", end=' ')
            print("%10i" % p.nProposals[0], end=' ')

            if p.nProposals[0]:  # Don't divide by zero
                print("       %5.1f " % (100.0 * float(p.nAcceptances[0]) / float(p.nProposals[0])), end=' ')
            else:
                print("           - ", end=' ')

            if p.tuning == None:
                print("      -", end=' ')
            elif p.tuning[0] < 2.0:
                print("  %8.4f" % p.tuning[0], end=' ')
            elif p.tuning[0] < 20.0:
                print("  %8.3f" % p.tuning[0], end=' ')
            elif p.tuning[0] < 200.0:
                print("  %8.1f" % p.tuning[0], end=' ')
            else:
                print("  %8.3g" % p.tuning[0], end=' ')
            print()

        # Tabulate topology changes for 'local', if any were attempted.
        doTopol = 0
        p = self.props.topologyProposalsDict.get('local')
        if p:
            for tNum in range(self.nChains):
                if p.nTopologyChangeAttempts[tNum]:
                    doTopol = 1
                    break
            if doTopol:
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
        p = self.props.topologyProposalsDict.get('eTBR')
        if p:
            for tNum in range(self.nChains):
                if p.nTopologyChangeAttempts[tNum]:
                    doTopol = 1
                    break
            if doTopol:
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
        p = self.props.topologyProposalsDict.get('local')
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

        p = self.props.topologyProposalsDict.get('eTBR')
        if p:
            if hasattr(p, 'nAborts'):
                if p.nAborts[0]:
                    print("The 'eTBR' proposal had %i aborts.  (Not counted in nProps above.)" % p.nAborts[0])
                    assert self.constraints
                    print("(Aborts due to violated constraints)")
                else:
                    if self.constraints:
                        print("The 'eTBR' proposal had no aborts (due to violated constraints).")

        for pN in ['polytomy']:
            for p in self.props.proposals:
                if p.name == pN:
                    if hasattr(p, 'nAborts'):
                        print("The %15s proposal had %5i aborts." % (p.name, p.nAborts[0]))
        for pN in ['compLocation', 'rMatrixLocation', 'gdasrvLocation']:
            for p in self.props.proposals:
                if p.name == pN:
                    if hasattr(p, 'nAborts'):
                        print("The %15s proposal (part %i) had %5i aborts." % (p.name, p.pNum, p.nAborts[0]))

        if self.nChains > 1:
            print("\n\nAcceptances and tunings by temperature")
            print("%s %30s %5s %5s %10s %13s%10s" % (
                spacer, 'proposal', 'part', 'tempNum', 'nProposals', 'acceptance(%)', 'tuning'))
            for p in self.props.proposals:
                for tempNum in range(self.nChains):
                    print("%s" % spacer, end=' ')
                    print("%30s" % p.name, end=' ')
                    if p.pNum != -1:
                        print(" %3i " % p.pNum, end=' ')
                    else:
                        print("   - ", end=' ')
                    print(" %3i " % tempNum, end=' ')
                    print("%10i" % p.nProposals[tempNum], end=' ')

                    if p.nProposals[tempNum]:  # Don't divide by zero
                        print("       %5.1f " % (
                            100.0 * float(p.nAcceptances[tempNum]) / float(p.nProposals[tempNum])), end=' ')
                    else:
                        print("           - ", end=' ')

                    if p.tuning == None:
                        print("      -", end=' ')
                    elif p.tuning[tempNum] < 2.0:
                        print("  %8.4f" % p.tuning[tempNum], end=' ')
                    elif p.tuning[tempNum] < 20.0:
                        print("  %8.3f" % p.tuning[tempNum], end=' ')
                    elif p.tuning[tempNum] < 200.0:
                        print("  %8.1f" % p.tuning[tempNum], end=' ')
                    else:
                        print("  %8.3g" % p.tuning[tempNum], end=' ')
                    print()

    def writeSwapMatrix(self):
        gm = ["Mcmc.writeSwapMatrix()"]
        if self.swapVector:
            gm.append("yourMcmc.swapVector is turned on, so you should use Mcmc.writeSwapVector() instead.")
            gm.append("or yourMcmcCheckPointReader.writeSwapVectors()")
            raise P4Error(gm)
        # swapVector is default now, so this would not be used.
        print("\nChain swapping, for %i gens, from gens %i to %i, inclusive." % (
            (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen))
        #print("    Swaps are presented as a square matrix, nChains * nChains.")
        print("    Upper triangle is the number of swaps proposed between two chains.")
        print("    Lower triangle is the percent swaps accepted.")
        #print("    The chainTemp is %s\n" % self.chainTemp)
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

    def writeSwapVector(self):
        assert self.swapVector
        assert self.nChains > 1
        # swapVector is default now.  However, self.swapMatrix is used and has the info.
        print("\nMcmc.writeSwapVector(). Chain swapping, for %i gens, from gens %i to %i, inclusive." % (
            (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen))
        print("    swapVector is on, so swaps occur only between adjacent chains.")
        chTmpString = '  '.join([f"{it:6.3f}" for it in self.chainTemps])
        print(f"    last chainTemps      {chTmpString}")
        chTmpDiffsString = '  '.join([f"{it:6.3f}" for it in self.chainTempDiffs[:-1]])
        print(f"    last chainTempDiffs  {chTmpDiffsString}")
        # This may be set differently in the run, but it is not saved (still true?),
        # so this could be invalid.  Fix this!
        #print(f"    self.swapTunerSampleSize is {self.swapTunerSampleSize}")
        print()

        print(f"{'chains':10}", end=' ')
        for i in range(1, self.nChains):
            s = f"{i-1}--{i}"
            print(f"{s:>9s}", end=' ')
        print()

        print(f"{' ':10}", end=' ')
        for i in range(1, self.nChains):
            print(f"{'------':>9s}", end=' ')
        print()
        
        print(f"{'proposed':10}", end=' ')
        for i in range(1, self.nChains):
            j = i - 1
            print(f"{self.swapMatrix[j][i]:9}", end=' ')
        print()

        print(f"{'accepted':10}", end=' ')
        for i in range(1, self.nChains):
            j = i - 1
            print(f"{self.swapMatrix[i][j]:9}", end=' ')
        print()

        print(f"{'accepted':10}", end=' ')
        for i in range(1, self.nChains):
            j = i - 1
            if self.swapMatrix[j][i] == 0:  # no proposals
                myProportion = "-"
            else:
                myProportion = f"{self.swapMatrix[i][j] / self.swapMatrix[j][i]:.6f}"
            print(f"{myProportion:>9}", end=' ')
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
        if not self.props.proposals:
            self._makeProposals()

            # If we are going to be doing the resolution class prior
            # in the polytomy move, we want to pre-compute the logs of
            # T_{n,m}.  Its a vector with indices (ie m) from zero to
            # nTax-2 inclusive.
            p = self.props.topologyProposalsDict.get('polytomy')
            if p and self.polytomyUseResolutionClassPrior:
                bigT = p4.func.nUnrootedTreesWithMultifurcations(self.tree.nTax)
                p.logBigT = [0.0] * (self.tree.nTax - 1)
                for i in range(1, self.tree.nTax - 1):
                    p.logBigT[i] = math.log(bigT[i])
                # print p.logBigT

            # Make a dictionary.  If there is more than one partition,
            # some proposals may have the same name, so it is a double
            # index pDict[proposal.pNum][proposal.name].
            # Some props (eg topology moves) have a pNum of -1, so that is extra.
            self.props.pDict = {}
            self.props.pDict[-1] = {}
            for pNum in range(self.tree.data.nParts):
                self.props.pDict[pNum] = {}
            for pr in self.props.proposals:
                self.props.pDict[pr.pNum][pr.name] = pr



    def _setOutputTreeFile(self):
        """Setup the (output) tree file for the mcmc."""

        gm = ['Mcmc._setOutputTreeFile()']

        # Check if it exists.
        if os.path.isfile(self.treeFileName):
            gm.append(f"The file '{self.treeFileName}' already exists.")
            raise P4Error(gm)

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
        # if self.tree.model.isHet:
        #     writeModelComment = True
        #     if self.tree.model.parts[0].ndch2:
        #         writeModelComment = False
        #         if self.tree.model.parts[0].ndch2_writeComps:
        #             writeModelComment = True
        #     if self.tree.model.parts[0].ndrh2:
        #         writeModelComment = False
        #         if self.tree.model.parts[0].ndrh2_writeRMatrices:
        #             writeModelComment = True
        writeModelComment = False
        for mp in self.tree.model.parts:
            if mp.nComps > 1 or mp.nRMatrices > 1 or mp.nGdasrvs > 1:
                writeModelComment = True
                break
                    
        if writeModelComment:
            self.treeFile.write('  [&&p4 models p%i' % self.tree.model.nParts)
            for pNum in range(self.tree.model.nParts):
                self.treeFile.write(' c%i.%i' % (pNum, self.tree.model.parts[pNum].nComps))
                self.treeFile.write(' r%i.%i' % (pNum, self.tree.model.parts[pNum].nRMatrices))
                self.treeFile.write(' g%i.%i' % (pNum, self.tree.model.parts[pNum].nGdasrvs))
            self.treeFile.write(']\n')
        self.treeFile.write('  [Tree numbers are gen+1]\nend;\n')

        if 0:
            self.prob.dump()
            self.tunings.dump()
            self.writeProposalProbs()
        if 0:
            for p in self.proposals:
                p.dump()
            # return

        

        

    def run(self, nGensToDo, verbose=True, equiProbableProposals=False, writeSamples=True):
        """Start the Mcmc running."""

        gm = ['Mcmc.run()']

        # Keep track of the first gen of this call to run(), maybe restart
        firstGen = self.gen + 1

        # Don't forget to set the PIVEC_MIN etc in a restart.
        # The hasattr() part is for older checkpoints that do not have self.init_PIVEC_MIN
        if var.mcmc_doCheck_PIVEC_MIN_etc and hasattr(self, "init_PIVEC_MIN"):
            doWarning = False
            if var.PIVEC_MIN != self.init_PIVEC_MIN:
                gm.append(f"Initial var.PIVEC_MIN for this run was {self.init_PIVEC_MIN}, but now is {var.PIVEC_MIN}")
                doWarning = True
            if var.RATE_MIN != self.init_RATE_MIN:
                gm.append(f"Initial var.RATE_MIN for this run was {self.init_RATE_MIN}, but now is {var.RATE_MIN}")
                doWarning = True
            if var.BRLEN_MIN != self.init_BRLEN_MIN:
                gm.append(f"Initial var.BRLEN_MIN for this run was {self.init_BRLEN_MIN}, but now is {var.BRLEN_MIN}")
                doWarning = True
            if var.GAMMA_SHAPE_MIN != self.init_GAMMA_SHAPE_MIN:
                gm.append(f"Initial var.GAMMA_SHAPE_MIN for this run was {self.init_GAMMA_SHAPE_MIN}, but now is {var.GAMMA_SHAPE_MIN}")
                doWarning = True
            if doWarning:
                gm.append("Were these values forgotten in a restart?  If so, add them to the restart script.")
                gm.append("Are these values intentionally different?  If so, you can turn off checking by setting, in your script ---")
                gm.append("var.mcmc_doCheck_PIVEC_MIN_etc = False")
                raise P4Error(gm)

        # Hidden experimental hack
        if self.doHeatingHack:
            print("Heating hack is turned on.")
            print("Heating hack temperature is %.2f" % self.heatingHackTemperature)
            #print("Heating hack affects proposals %s" % self.heatingHackProposalNames)
            self.originalHeatingHackTemperature = self.heatingHackTemperature

        if self.simTemp:
            if self.doHeatingHack:
                gm.append("If we are doing simTemp, then doHeatingHack should be off.")
                raise P4Error(gm)
        
            if not self.simTemp.temps:
                gm.append("If we are doing simTemp, we should have SimTempTemp objects in m.simTemp.temps.")
                raise P4Error(gm)
            if len(self.simTemp.temps) != self.simTemp.nTemps:
                gm.append("If we are doing simTemp, we should have %i SimTempTemp objects in m.simTemp.temps." % \
                          self.simTemp)
                raise P4Error(gm)
            for it in self.simTemp.temps:
                if not isinstance(it, SimTempTemp):
                    gm.append("Each item in m.simTemp_temps should be a SimTempTemp object.")
                    raise P4Error(gm)

            self.simTemp.nTempChangeProposals = 0
            self.simTemp.nTempChangeAccepts = 0
            self.simTemp.tNums = []
            for tmp in self.simTemp.temps:
                tmp.occupancy = 0
                tmp.nProposalsUp = 0
                tmp.nAcceptedUp = 0
                tmp.nProposalsDn = 0
                tmp.nAcceptedDn = 0

            if self.simTemp_doLogTemps:
                self.simTemp.writeTempsOneLineHeader()

        if self.prob.polytomy:
            for pNum in range(self.tree.model.nParts):
                mp = self.tree.model.parts[pNum]
                if mp.ndch2:
                    gm.append("Part %i uses ndch2" % pNum)
                    gm.append("Ndch2 does not work with the polytomy move.")
                    raise P4Error(gm)
                if mp.ndrh2:
                    gm.append("Part %i uses ndrh2" % pNum)
                    gm.append("Ndrh2 does not work with the polytomy move.")
                    raise P4Error(gm)
                
        if self.doSteppingStone:
            # self.logger.info("doSteppingStone is on, so turn writing checkpoints off")
            self.checkPointInterval = False

        if writeSamples and self.checkPointInterval:
            # Check that checkPointInterval makes sense.
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
                gm.append("With the current settings, the last generation won't be on a checkPointInterval.")
                gm.append("self.gen+1=%i, nGensToDo=%i, checkPointInterval=%i" % ((self.gen + 1),
                                                                                  nGensToDo, self.checkPointInterval))
                raise P4Error(gm)
            #  2.  We also want the checkPointInterval to be evenly
            #      divisible by the sampleInterval.
            if not self.simTemp:
                if self.checkPointInterval % self.sampleInterval == 0:
                    pass
                else:
                    gm.append("The checkPointInterval (%i) should be evenly divisible" % self.checkPointInterval)
                    gm.append("by the sampleInterval (%i)." % self.sampleInterval)
                    raise P4Error(gm)


        if self.props.proposals:  # It is a re-start

            # The probs and tunings may have been changed by the user.
            self._refreshProposalProbs()

            # This stuff below should be the same as is done after pickling,
            # see below.
            self.startMinusOne = self.gen

            # Start the tree partitions over.
            self.treePartitions = None
            # Zero the proposal counts
            for p in self.props.proposals:
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
            # The swap vector is just the diagonal of the swap matrix
            if self.swapVector and self.nChains > 1:
                if self.swapTunerDoTuning:
                    print("======= Resetting chainTempDiffs and chainTemps")
                    # These are differences in temperatures between adjacent chains.  The last one is not used.
                    self.chainTempDiffs = [self.chainTemp] * self.nChains
                    # These are cumulative, summed over the diffs.  This needs to be done whenever the diffs change
                    self.chainTemps = [0.0]
                    for dNum in range(self.nChains - 1):
                        self.chainTemps.append(self.chainTempDiffs[dNum] + self.chainTemps[-1])
                else:
                    # doTuning has been turned off
                    assert self.chainTemps, "swapTunerDoTuning is off.  Please specify chainTemps."
                    assert len(self.chainTemps) == self.nChains
                    assert self.chainTemps[0] == 0.0

                    # set chainTempDiffs based on chainTemps, ie the reverse of above
                    self.chainTempDiffs = []
                    for dNum in range(1, self.nChains):
                        self.chainTempDiffs.append(self.chainTemps[dNum] - self.chainTemps[dNum - 1])
                    self.chainTempDiffs.append(0.0) # with reason? -- maybe not.

            if self.simulate:
                self.simTree = self.tree.dupe()
                self.simTree.data = self.tree.data.dupe()
                self.simTree.calcLogLike(verbose=False)

            self._makeChainsAndProposals()
            self._setOutputTreeFile()
            if self.simulate:
                self._writeSimFileHeader(self.tree)

        if self.setupPfLogging:
            for ch in self.chains:
                pf.setMcmcTreeCallback(ch.curTree.cTree, self._logFromPfModule)
                pf.setMcmcTreeCallback(ch.propTree.cTree, self._logFromPfModule)
                # Should I do self.tree and self.simTree as well?

        if self.doCpo:
            if self.cpo_sumsOfInverseSiteLikes:
                pass
            else:
                nSites = sum([p.nChar for p in self.tree.data.parts])
                self.cpo_sumsOfInverseSiteLikes = [0.0] * nSites
                self.cpo_nSamples = 0

        if verbose:
            self.props.writeProposalIntendedProbs()
            sys.stdout.flush()

            if self.stickyRootComp:
                print("stickyRootComp is turned on.")
                # compLocation should be turned off.  Is it?
                if self.prob.compLocation > 0.0:
                    gm.append("If stickyRootComp is turned on, then compLocation should be off.")
                    gm.append("Turn it off, by eg (for Mcmc object m)")
                    gm.append("m.prob.compLocation = 0.0")
                    raise P4Error(gm)
                    
        # Find the cold chain, the one where tempNum is 0
        self.coldChainNum = -1
        for i in range(len(self.chains)):
            if self.chains[i].tempNum == 0:
                self.coldChainNum = i
                break
        if self.coldChainNum == -1:
            gm.append("Unable to find which chain is the cold chain.  That is Bad.")
            raise P4Error(gm)


        # # If polytomy is turned on, then it is possible to get a star
        # # tree, in which case local will not work.  So if we have both
        # # polytomy and local proposals, we should also have brLen.
        # if 'polytomy' in self.proposalsHash and 'local' in self.proposalsHash:
        #     if 'brLen' not in self.proposalsHash:
        #         gm.append("If you have polytomy and local proposals, you should have a brLen proposal as well.")
        #         gm.append("It can have a low proposal probability, but it needs to be there.")
        #         gm.append("Turn it on by eg yourMcmc.prob.brLen = 0.001")
        #         raise P4Error(gm)

        if writeSamples and not self.doSteppingStone:
            # it is a re-start, so we need to back over the "end;" in the tree
            # files.
            if self.treeFile:
                if not self.treeFile.closed:
                    self.treeFile.close()

            f2 = open(self.treeFileName, 'r+b')
            # print(f2, type(f2), f2.tell())   
            # <_io.BufferedRandom name='mcmc_trees_0.nex'> <class '_io.BufferedRandom'> 0
            pos = -1
            while 1:
                f2.seek(pos, 2)
                c = f2.read(1)
                if c == b';':
                    break
                pos -= 1
            # print "pos now %i" % pos
            pos -= 3  # end;
            f2.seek(pos, 2)
            c = f2.read(4)
            # print "got c = '%s'" % c
            if c != b"end;":
                gm.append("Mcmc.run().  Failed to find and remove the 'end;' at the end of the tree file.")
                raise P4Error(gm)
            else:
                f2.seek(pos, 2)
                f2.truncate()
            f2.close()

        if self.gen > -1:

            self.logger.info("Re-starting the MCMC run %i from gen=%i" % (self.runNum, self.gen))

            if verbose:
                print()
                print("Re-starting the MCMC run %i from gen=%i" % (self.runNum, self.gen))
                if not writeSamples:
                    print("Arg 'writeSamples' is off" )
                print("Set to do %i more generations." % nGensToDo)
                if writeSamples and self.writePrams:
                    if self.chains[0].curTree.model.nFreePrams == 0:
                        print("There are no free prams in the model, so I am turning writePrams off.")
                        self.writePrams = False
                sys.stdout.flush()

            self.startMinusOne = self.gen
        else:
            self.logger.info("Starting the MCMC %s run %i" % (
                (self.constraints and "(with constraints)" or ""), self.runNum))
            if self.simTemp:
                self.logger.info("Using simulated tempering MCMC, with %i temperatures" % (
                    self.simTemp.nTemps))
            if verbose:
                if self.nChains > 1:
                    print("Using Metropolis-coupled MCMC, with %i chains." % self.nChains)
                    if self.swapTunerDoTuning:
                        print("Temperatures are tuned (self.swapTunerDoTuning is on)")
                    else:
                        print("self.swapTunerDoTuning is turned off, so temperatures will not be tuned.")
                else:
                    print("Single chain --- not using Metropolis-coupled MCMC.")
                if self.simTemp:
                    print("Using simulated tempering MCMC, with %i temperatures" % (
                        self.simTemp.nTemps))
                print("Starting the MCMC %s run %i" % (
                    (self.constraints and "(with constraints)" or ""), self.runNum))
                print("Set to do %i generations." % nGensToDo)

            if self.writePrams:
                if self.chains[0].curTree.model.nFreePrams == 0:
                    if verbose:
                        print("There are no free prams in the model, so I am turning writePrams off.")
                    self.writePrams = False
                else:
                    self.pramsFile = open(self.pramsFileName, 'w')
                    self.chains[0].curTree.model.writePramsProfile(self.pramsFile, self.runNum)
                    self.pramsFile.write("genPlus1")
                    self.chains[0].curTree.model.writePramsHeaderLine(self.pramsFile)
            if self.writeHypers:
                if not self.tree.model.parts[0].ndch2 and not self.tree.model.parts[0].ndrh2:     
                    # and therefore all model parts, this week  -- lame!
                    self.writeHypers = False
                else:
                    self.hypersFile = open(self.hypersFileName, 'w')
                    self.hypersFile.write('genPlus1')
                    self.chains[0].curTree.model.writeHypersHeaderLine(self.hypersFile)

        if not writeSamples:
            self.logger.info("Mcmc.run() arg 'writeSamples' is off, so samples are not being written")
        if equiProbableProposals:
            self.logger.info("Mcmc.run() arg 'equiProbableProposals' is turned on")
        if self.doSteppingStone:
            self.logger.info("Mcmc.run() doSteppingStone is turned on")
            self.ssLikesFile = open(self.ssLikesFileName, 'a')



        if verbose:
            if writeSamples:
                print("Sampling every %i." % self.sampleInterval)
            else:
                print("Arg 'writeSamples' is off, so samples are not being written")
            if equiProbableProposals:
                print("Arg 'equiProbableProposals' is turned on")
            if self.checkPointInterval:
                print("CheckPoints written every %i." % self.checkPointInterval)
            else:
                print("CheckPoints are not scheduled to be written")
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

        abortableProposals = ['local', 'polytomy', 'compLocation', 'eTBR'
                              'rMatrixLocation', 'gdasrvLocation']

        ##################################################
        ############### Main loop ########################
        ##################################################

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
                    # Choose a proposal
                    gotIt = False
                    safety = 0
                    while not gotIt:
                        # equiProbableProposals is True or False.  Usually False.
                        aProposal = self.props.chooseProposal(equiProbableProposals)
                        if aProposal:
                            gotIt = True

                        if aProposal.name == 'local':
                            # Can't do local on a star tree.
                            if self.chains[chNum].curTree.nInternalNodes == 1:
                                #aProposal = self.proposalsHash['brLen']
                                gotIt = False

                        elif aProposal.name in ['root3', 'root3n']:
                            # Can't do root3 or root3n on a star tree.
                            if self.chains[chNum].curTree.nInternalNodes == 1:
                                gotIt = False

                        if aProposal.doAbort:
                            gotIt = False

                        safety += 1
                        if safety > 1000:
                            gm.append("Could not find a proposal after %i attempts." % safety)
                            gm.append("Possibly a programming error.")
                            gm.append("Or possibly it is just a pathologically frustrating Mcmc.")
                            raise P4Error(gm)

                    if 0:
                        print("==== Mcmc.run() gNum=%i, chNum=%i, aProposal=%s (part %i)" % (
                            gNum, chNum, aProposal.name, aProposal.pNum), end='\n')
                        #sys.stdout.flush()

                    # success returns None
                    failure = self.chains[chNum].gen(aProposal)
                    #failure = None

                    if failure:
                        myWarn = "Mcmc.run() main loop.  Proposal %s generated a 'failure'. Should be debugged." % aProposal.name
                        self.logger.warning(myWarn)

                    nAttempts += 1
                    if nAttempts > 1000:
                        gm.append("Was not able to do a successful generation after %i attempts." % nAttempts)
                        raise P4Error(gm)

                    # Continuous tuning.  We have a tuning, and propose/accept
                    # tallies for each temperature, kept separately.  Note that
                    # since chNum does not equal tempNum, the most recently
                    # incremented values (and the ones we want to tune now) will
                    # be aProposal.tnNSamples[tempNum] and
                    # aProposal.tnNAccepts[tempNum], where tempNum will most
                    # likely not be chNum.  So we get the tempNum from this
                    # chNum, and tune it.
                    
                    # From init
                    # self.tunableProps = """allBrLens allCompsDir brLen compDir 
                    # gdasrv local ndch2_internalCompsDir 
                    # ndch2_internalCompsDirAlpha ndch2_leafCompsDir 
                    # ndch2_leafCompsDirAlpha ndch2_priorRefCompDir
                    # pInvar rMatrixDir allRMatricesDir ndrh2_priorRefRMatrix
                    # ndrh2_internalRatesDir ndrh2_internalRatesDirAlpha ndrh2_leafRatesDir 
                    # ndrh2_leafRatesDirAlpha relRate """.split()

                    # maybeTunablesButNotNow  compLocation eTBR polytomy root3 rMatrixLocation root2

                    if aProposal.name in self.tunableProps:
                        tempNum = self.chains[chNum].tempNum
                        if aProposal.tnNSamples[tempNum] >= aProposal.tnSampleSize:
                            aProposal.tune(tempNum)
                    # elif aProposal.name == 'root3' and self.doRoot3Tuning:
                    #     tempNum = self.chains[chNum].tempNum
                    #     if aProposal.tnNSamples[tempNum] >= aProposal.tnSampleSize:
                    #         aProposal.tune(tempNum)
                    # elif aProposal.name == 'root3n' and self.doRoot3nTuning:
                    #     tempNum = self.chains[chNum].tempNum
                    #     if aProposal.tnNSamples[tempNum] >= aProposal.tnSampleSize:
                    #         aProposal.tune(tempNum)

                # print "   Mcmc.run(). finished a gen on chain %i" % (chNum)
                for p in self.props.proposals:
                    if p.name in abortableProposals:
                        p.doAbort = False


            if self.nChains == 1:
                if self.simTemp:
                    self.simTemp.temps[self.simTemp.curTemp].occupancy += 1
                    # if (self.gen + 1) % self.simTemp.tempChangeProposeFreq == 0:
                    self.simTemp.tNums.append(self.simTemp.curTemp)

                    # this is a deque.  Append on the left and pop on the right
                    self.simTemp.tNumSample.appendleft(self.simTemp.curTemp)
                    self.simTemp.tNumSample.pop() # from the right, of course

                    self.simTemp.proposeTempChange()

                    if self.simTemp_doTuneTempsAndLogPi:
                        if self.gen >= self.simTemp.tNumSampleSize:
                            if self.gen % self.simTemp.tNumSampleSize == 0:                    
                                self.simTemp.tuneTempsAndLogPi()




            else:
                # Propose swap, if there is more than 1 chain.
                self._proposeSwapChainsInMcmcmc()


            # =====================================
            # If it is a writeInterval, write stuff
            # =====================================

            doWrite = False
            if self.simTemp:
                if self.simTemp.curTemp == 0:
                    # Thinning the simTemp chain is awkward.  Here is a hack.
                    # var.mcmc_simTemp_thinning is a list of digit strings, by default ['0']
                    if var.mcmc_simTemp_thinning:
                        genString = f"{self.gen + 1}"
                        for ending in var.mcmc_simTemp_thinning:
                            if genString.endswith(ending):
                                doWrite = True
                    else:
                        doWrite = True
            else:
                if ((self.gen + 1) % self.sampleInterval) == 0:
                    doWrite = True
            
            if doWrite and self.doSteppingStone:
                if writeSamples:
                    self.ssLikesFile.write('%f\n' % self.chains[self.coldChainNum].curTree.logLike)


            if doWrite and not self.doSteppingStone:
                if writeSamples:
                    self._writeSample()


                # Do a simulation
                if self.simulate:
                    # print "about to simulate..."
                    self._doSimulate(self.chains[self.coldChainNum].curTree)
                    # print "...finished simulate."

                # Do other stuff.
                if hasattr(self, 'hook'):
                    self.hook(self.chains[self.coldChainNum].curTree)

                if 0 and self.constraints:
                    print("Mcmc x1c")
                    print(self.chains[0].verifyIdentityOfTwoTreesInChain())
                    print("b checking curTree ..")
                    self.chains[0].curTree.checkSplitKeys()
                    print("b checking propTree ...")
                    self.chains[0].propTree.checkSplitKeys()
                    print("Mcmc xxx")

                if 0:
                    print("Mcmc.run() xxx curTree propTree")
                    a = self.chains[0].curTree
                    b = self.chains[0].propTree
                    for nNum in range(len(a.nodes)):
                        if a.nodes[nNum] != a.root:
                            print(f"{nNum:3} {a.nodes[nNum].br.rawSplitKey} {b.nodes[nNum].br.rawSplitKey}")
                        else:
                            print(f"{nNum:3} is root")

                # Add curTree to treePartitions
                if self.treePartitions:
                    self.treePartitions.getSplitsFromTree(
                        self.chains[self.coldChainNum].curTree)
                else:
                    self.treePartitions = TreePartitions(
                        self.chains[self.coldChainNum].curTree)

                if 0:
                    print("Mcmc.run() zzz curTree propTree")
                    a = self.chains[0].curTree
                    b = self.chains[0].propTree
                    for nNum in range(len(a.nodes)):
                        if a.nodes[nNum] != a.root:
                            print(f"{nNum:3} {a.nodes[nNum].br.rawSplitKey} {b.nodes[nNum].br.rawSplitKey}")
                        else:
                            print(f"{nNum:3} is root")

                # After _getSplitsFromTree, need to follow, at some point,
                # with _finishSplits().  Do that when it is pickled, or at the
                # end of the run.

                # Checking and debugging constraints
                if 0 and self.constraints:
                    print("Mcmc x1d")
                    print(self.chains[self.coldChainNum].verifyIdentityOfTwoTreesInChain())
                    print("c checking curTree ...")
                    self.chains[self.coldChainNum].curTree.checkSplitKeys()
                    print("c checking propTree ...")
                    self.chains[self.coldChainNum].propTree.checkSplitKeys()
                    # print "c checking that all constraints are present"
                    #theSplits = [n.br.splitKey for n in self.chains[0].curTree.iterNodesNoRoot()]
                    # for sk in self.constraints.constraints:
                    #    if sk not in theSplits:
                    #        gm.append("split %i is not present in the curTree." % sk)
                    #        raise P4Error(gm)
                    print("Mcmc zzz")

                # Double check that the curTree has all the constraints
                # The proposals will have done this already, so this is a double check to find errors.
                if self.constraints:
                    #print("Mcmc.run() double checking constraints")
                    for sk in self.constraints.constraints:
                        ret = self.chains[self.coldChainNum].curTree.nodeForSplitKeyDict.get(sk)
                        if not ret:
                            gm.append("Programming error.")
                            gm.append(f"The current tree (the last tree sampled) does not contain constraint sk")
                            gm.append("%s" % p4.func.getSplitStringFromKey(sk, self.tree.nTax))
                            raise P4Error(gm)
                    if self.constraints.rooting:
                        ret = self.constraints.areConsistentWithTreeRoot(self.chains[self.coldChainNum].curTree)
                        if not ret:
                            gm.append("Programming error.")
                            gm.append(f"The current tree (the last tree sampled) is not consistent with root constraint")
                            raise P4Error(gm)



            
            # Do checkpoints
            doCheckPoint = False
            if not writeSamples:
                doCheckPoint = False
            elif self.doSteppingStone:
                doCheckPoint = False
            else:
                if self.checkPointInterval and (gNum + 1) % self.checkPointInterval == 0:
                    doCheckPoint = True                
            
            if doCheckPoint:

                if self.simTemp:
                    fout = open(self.simTempFileName, 'a')
                    print("-" * 50, file=fout)
                    print("gen+1 %11i" % (self.gen + 1), file=fout)
                    self.simTemp.writeTempsTable(flob=fout, zeroCountsAfter=True)
                    fout.close()

                # print(f"writing checkpoint at self.gen {self.gen}, gNum {gNum}")
                self.checkPoint()

                # The stuff below needs to be done in a re-start as well.
                # See above "if self.proposals:"
                self.startMinusOne = self.gen

                # Start the tree partitions over.
                self.treePartitions = None
                # Zero the proposal counts
                for p in self.props.proposals:
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

        # End of the Main loop.  Gens finished.  Clean up.
        print()
        if verbose:
            print("Finished %s generations." % nGensToDo)

        if writeSamples and not self.doSteppingStone:
            # Write "end;" regardless of writeSamples, or else subsequent attempts to remove it will fail.
            # No, that does not work for simTemp trialAndError.  Why not?
            if not self.treeFile or self.treeFile.closed:
                self.treeFile = open(self.treeFileName, 'a')
            self.treeFile.write('end;\n\n')

        if self.likesFile:
            self.likesFile.close()
        if self.treeFile:
            self.treeFile.close()
        if self.simFile:
            self.simFile.close()
        if self.pramsFile:
            self.pramsFile.close()
        if self.hypersFile:
            self.hypersFile.close()
        if self.ssLikesFile:
            self.ssLikesFile.close()
        #if self.siteLikesFile:
        #    self.siteLikesFile.close()




    def _writeSample(self):
        try:
            self.likesFile.write('%11i %f\n' % (self.gen + 1, self.chains[self.coldChainNum].curTree.logLike))
        except (AttributeError, ValueError):
            self.likesFile = open(self.likesFileName, 'a')
            self.likesFile.write('%11i %f\n' % (self.gen + 1, self.chains[self.coldChainNum].curTree.logLike))

        # Check the likelihood every write interval
        if 0:
            oldLike = self.chains[self.coldChainNum].curTree.logLike
            print("gen+1 %11i  %f  " % (
                self.gen+1, 
                self.chains[self.coldChainNum].curTree.logLike), end=' ')
            self.chains[self.coldChainNum].curTree.calcLogLike(verbose=False)
            newLike = self.chains[self.coldChainNum].curTree.logLike
            print("%f" % self.chains[self.coldChainNum].curTree.logLike, end=' ')
            likeDiff = math.fabs(oldLike - newLike)
            if likeDiff > 1e-14:
                print("%f" % likeDiff)
            else:
                print()

        try:
            self.treeFile.write("  tree t_%i = [&U] " % (self.gen + 1))
        except (AttributeError, ValueError):
            self.treeFile = open(self.treeFileName, 'a')
            self.treeFile.write("  tree t_%i = [&U] " % (self.gen + 1))

        # if self.tree.model.parts[0].ndch2 or self.tree.model.parts[0].ndrh2: # lazy programming -- and therefore all model parts
        #     if self.tree.model.parts[0].ndch2_writeComps  or self.tree.model.parts[0].ndrh2_writeRMatrices:
        #         self.chains[self.coldChainNum].curTree.writeNewick(self.treeFile,
        #                                                       withTranslation=1,
        #                                                       translationHash=self.translationHash,
        #                                                       doMcmcCommandComments=True)
        #     else:
        #         self.chains[self.coldChainNum].curTree.writeNewick(self.treeFile,
        #                                                       withTranslation=1,
        #                                                       translationHash=self.translationHash,
        #                                                       doMcmcCommandComments=False)

        # else:

        self.chains[self.coldChainNum].curTree.writeNewick(self.treeFile,
                                                          withTranslation=1,
                                                          translationHash=self.translationHash,
                                                          doMcmcCommandComments=self.tree.model.isHet)

        if self.writePrams:
            try:
                self.pramsFile.write("%12i" % (self.gen + 1))
            except (AttributeError, ValueError):
                self.pramsFile = open(self.pramsFileName, 'a')
                self.pramsFile.write("%12i" % (self.gen + 1))
            self.chains[self.coldChainNum].curTree.model.writePramsLine(self.pramsFile)

        if self.writeHypers:
            try:
                self.hypersFile.write("%12i" % (self.gen + 1))
            except (AttributeError, ValueError):
                self.hypersFile = open(self.hypersFileName, 'a')
                self.hypersFile.write("%12i" % (self.gen + 1))
            self.chains[self.coldChainNum].curTree.model.writeHypersLine(self.hypersFile)

        # if self.writeSiteLikes:
        #     try:
        #         pf.p4_treeLogLike(self.chains[self.coldChainNum].curTree.cTree, 1)
        #         self.siteLikesFile.write(f"{self.gen + 1:11} ")
        #     except (AttributeError, ValueError):
        #         self.siteLikesFile = open(self.siteLikesFileName, 'a')
        #         self.siteLikesFile.write(f"{self.gen + 1:11} ")
        #     siteLikes = []
        #     for p in self.tree.data.parts:
        #         siteLikes += pf.getSiteLikes(p.cPart)
        #     for sL in siteLikes:
        #         self.siteLikesFile.write("%.17g " % sL)
        #     self.siteLikesFile.write("\n")

        if self.doCpo:
            if self.gen + 1 > self.cpo_startGen:
                siteNum = 0
                pf.p4_treeLogLike(self.chains[self.coldChainNum].curTree.cTree, 1)
                for dp in self.chains[self.coldChainNum].curTree.data.parts:
                    siteLikes = pf.getSiteLikes(dp.cPart)
                    for sL in siteLikes:
                        self.cpo_sumsOfInverseSiteLikes[siteNum] += 1./sL
                        siteNum += 1
                self.cpo_nSamples += 1

    def _proposeSwapChainsInMcmcmc(self):
        if self.swapVector:
            rTempNum1 = random.randrange(self.nChains - 1)
            rTempNum2 = rTempNum1 + 1
            chain1 = None
            chain2 = None
            for ch in self.chains:
                if ch.tempNum == rTempNum1:
                    chain1 = ch
                elif ch.tempNum == rTempNum2:
                    chain2 = ch
            assert chain1 and chain2

            # Use the upper triangle of swapMatrix for nAttempts
            self.swapMatrix[chain1.tempNum][chain2.tempNum] += 1

            lnR = (1.0 / (1.0 + (self.chainTemps[chain1.tempNum]))
                    ) * chain2.curTree.logLike
            lnR += (1.0 / (1.0 + (self.chainTemps[chain2.tempNum]))
                    ) * chain1.curTree.logLike
            lnR -= (1.0 / (1.0 + (self.chainTemps[chain1.tempNum]))
                    ) * chain1.curTree.logLike
            lnR -= (1.0 / (1.0 + (self.chainTemps[chain2.tempNum]))
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
            if self.swapTuner:
                # Index the nAttempts and nSwaps with the lower of the two tempNum's, which would be chain1.tempNum
                self.swapTuner.nAttempts[chain1.tempNum] += 1
                if acceptSwap:
                    self.swapTuner.nSwaps[chain1.tempNum] += 1
                if self.swapTunerDoTuning:
                    if self.swapTuner.nAttempts[chain1.tempNum] >= self.swapTunerSampleSize:
                        self.swapTuner.tune(chain1.tempNum)
                        # tune() zeros nAttempts and nSwaps counters

            if acceptSwap:
                # Use the lower triangle of swapMatrix to keep track of
                # nAccepted's
                assert chain1.tempNum < chain2.tempNum
                self.swapMatrix[chain2.tempNum][chain1.tempNum] += 1

                # Do the swap
                chain1.tempNum, chain2.tempNum = chain2.tempNum, chain1.tempNum

        else:    # swap matrix
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


            lnR = (1.0 / (1.0 + (self.chainTemp * chain1.tempNum))
                    ) * chain2.curTree.logLike
            lnR += (1.0 / (1.0 + (self.chainTemp * chain2.tempNum))
                    ) * chain1.curTree.logLike
            lnR -= (1.0 / (1.0 + (self.chainTemp * chain1.tempNum))
                    ) * chain1.curTree.logLike
            lnR -= (1.0 / (1.0 + (self.chainTemp * chain2.tempNum))
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

            # # for continuous temperature tuning with self.swapTuner
            # if self.swapTuner and thisCh1Temp == 0 and thisCh2Temp == 1:
            #     self.swapTuner.swaps01_nAttempts += 1
            #     if acceptSwap:
            #         self.swapTuner.swaps01_nSwaps += 1
            #     if self.swapTuner.swaps01_nAttempts >= self.swapTuner.sampleSize:
            #         self.swapTuner.tune(self)
            #         # tune() zeros nAttempts and nSwaps counters

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
        self.coldChainNum = -1
        for i in range(len(self.chains)):
            if self.chains[i].tempNum == 0:
                self.coldChainNum = i
                break
        if self.coldChainNum == -1:
            gm.append("Unable to find which chain is the cold chain.  Bad.")
            raise P4Error(gm)


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
        self.simFile = open(self.simFileName, 'a')

        self.simFile.write(" genPlus1")
        # If self.simulate contains a 1, do unconstrained log like
        if 1 & self.simulate:
            for pNum in range(self.simTree.data.nParts):
                self.simFile.write(' uncLike%i' % pNum)
        if 2 & self.simulate:  # If self.simulate contains a 2, do bigX^2
            for pNum in range(self.simTree.model.nParts):
                self.simFile.write(' bigXSq%i' % pNum)
        # If self.simulate contains a 4, do meanNCharPerSite
        if 4 & self.simulate:
            for pNum in range(self.simTree.model.nParts):
                self.simFile.write(' meanNCharsPerSite%i' % pNum)
        # If self.simulate contains an 8, do c_m, the compStatFromCharFreqs
        if 8 & self.simulate:
            for pNum in range(self.simTree.model.nParts):
                self.simFile.write(' c_mSim%i  c_mOrig%i' % (pNum, pNum))
        # If self.simulate contains a 16, do constant sites count
        if 16 & self.simulate:
            for pNum in range(self.simTree.model.nParts):
                self.simFile.write(' nConstSites%i' % pNum)
        self.simFile.write('\n')

    def _doSimulate(self, curTree):
        curTree.copyToTree(self.simTree)
        curTree.model.copyValsTo(self.simTree.model)
        self.simTree.simulate()
        try:
            self.simFile.write(" %11i" % (self.gen + 1))
        except (AttributeError, ValueError):
            self.simFile = open(self.simFileName, 'a')
            self.simFile.write(" %11i" % (self.gen + 1))

        # If self.simulate contains a 1, do unconstrained log like
        if 1 & self.simulate:
            for p in self.simTree.data.parts:
                self.simFile.write(' %f' % pf.getUnconstrainedLogLike(p.cPart))

        if 2 & self.simulate:  # If self.simulate contains a 2, do bigX^2
            if self.blankSeqNums == None:
                ret = self.simTree.data.simpleBigXSquared()
                for pNum in range(self.simTree.model.nParts):
                    self.simFile.write(' %f' % ret[pNum])
            else:
                # We do not want sims corresponding to blank seqs to contribute to the bigXSq
                ret2 = self.simTree.data.compoChiSquaredTest(verbose=0, skipTaxNums=self.blankSeqNums, skipColumnZeros=True)
                for pNum in range(self.simTree.model.nParts):
                    self.simFile.write(' %f' % ret2[pNum][0])

            # check ...
            # for i in range(len(ret)):
            #    if math.fabs(ret[i][0] - ret2[i]) > 0.000001:
            # print "The two methods of bigXSquared calculations differed.  %f
            # and %f" % (ret[i], ret2[i])

        # If self.simulate contains a 4, do meanNCharPerSite
        if 4 & self.simulate:
            ret = self.simTree.data.meanNCharsPerSite()
            # ret is a list, one number per part
            for pNum in range(self.simTree.model.nParts):
                self.simFile.write(' %f' % ret[pNum])
        # If self.simulate contains an 8, do c_m, the compStatFromCharFreqs
        if 8 & self.simulate:
            ret = self.simTree.compStatFromCharFreqs()
            ret2 = curTree.compStatFromCharFreqs()
            # ret is a list, one number per part
            for pNum in range(self.simTree.model.nParts):
                self.simFile.write(' %f  %f' % (ret[pNum], ret2[pNum]))
                # print ' compStatFromCharFreqs: %f  %f' % (ret[pNum],
                # ret2[pNum])
        # If self.simulate contains a 16, do constant sites count
        if 16 & self.simulate:
            ret = self.simTree.data.simpleConstantSitesCount()
            # ret is a list, one number per part
            for pNum in range(self.simTree.model.nParts):
                self.simFile.write(' %i' % ret[pNum])
        self.simFile.write('\n')

    def checkPoint(self):

        if 0:
            for chNum,ch in enumerate(self.chains):
                print("chain %i ==================" % chNum)
                ch.curTree.summarizeModelComponentsNNodes()

        #print("Before removing data", end=' ')
        #self.chains[0].curTree.calcLogLike(verbose=True, resetEmpiricalComps=False)

        # Make a copy of self, but with no cStuff.
        # But we don't want to copy data or logger.  So detach them.
        savedData = self.tree.data
        self.tree.data = None

        # The logger does not pickle
        savedLogger = self.logger
        self.logger = None
        savedLoggerPrinter = self.loggerPrinter
        self.loggerPrinter = None

        #print(f"checkPoint() self.simTemp_tempsFlob is {self.simTemp_tempsFlob}", end=' ')
        #print(f"type {type(self.simTemp_tempsFlob)}")  # type <class '_io.TextIOWrapper'>
        # This does not pickle
        self.simTemp_tempsFlob = None


        # gsl_rng state
        the_gsl_rng_size = pf.gsl_rng_size(var.gsl_rng) # size of the state
        assert self.gsl_rng_state_ndarray.shape[0] == the_gsl_rng_size
        pf.gsl_rng_getstate(var.gsl_rng, self.gsl_rng_state_ndarray)

        # random (python module) state
        self.randomState = random.getstate()

        if self.setupPfLogging:
            # Get rid of the mcmcTreeCallbacks
            for ch in self.chains:
                pf.unsetMcmcTreeCallback(ch.curTree.cTree)
                pf.unsetMcmcTreeCallback(ch.propTree.cTree)

        if self.simulate:
            savedSimData = self.simTree.data
            self.simTree.data = None
        for ch in self.chains:
            ch.curTree.data = None
            ch.curTree.savedLogLike = ch.curTree.logLike
            ch.propTree.data = None

        # Open files don't copy.  Furthermore, merely closing a file
        # is not enough --- nor is it enough to set the variable to
        # None (although that should be enough!  --- it works for
        # small examples elsewhere!  Why not here?!?).  It needs to be
        # deleted, and then set to None.  Extra: closing a closed file
        # is not an error.
        for myf in [
                self.likesFile,
                self.treeFile,
                self.simFile,
                self.pramsFile,
                self.hypersFile,
                # self.siteLikesFile,
        ]:
            if myf:
                if not myf.closed:
                    myf.close()
            del(myf)

        self.likesFile = None
        self.treeFile = None
        self.simFile = None
        self.pramsFile = None
        self.hypersFile = None
        # self.siteLikesFile = None

        theCopy = copy.deepcopy(self)

        # Re-attach data and logger to self.
        self.tree.data = savedData
        self.logger = savedLogger
        self.loggerPrinter = savedLoggerPrinter
        self.tree.calcLogLike(verbose=False, resetEmpiricalComps=False)
        if self.simulate:
            self.simTree.data = savedSimData
            self.simTree.calcLogLike(verbose=False, resetEmpiricalComps=False)
        for chNum,ch in enumerate(self.chains):
            ch.curTree.data = savedData
            #print("After restoring data", end=' ')
            ch.curTree.calcLogLike(verbose=False, resetEmpiricalComps=False)
            theDiff = math.fabs(ch.curTree.savedLogLike - ch.curTree.logLike)
            if theDiff > 0.01:
                theMessage = "Mcmc.checkPoint(), chainNum %i. " % chNum
                theMessage += "Bad likelihood calculation just before writing the checkpoint. "
                theMessage += "Old curTree.logLike %g, new curTree.logLike %g, diff %g" % (ch.curTree.savedLogLike, ch.curTree.logLike, theDiff)
                self.logger.info(theMessage)
                print()
                print(theMessage)
                sys.stdout.flush()
            ch.propTree.data = savedData
            #ch.propTree.dump(node=True)
            #ch.propTree.draw()
            ch.propTree.calcLogLike(verbose=False, resetEmpiricalComps=False)
        if self.setupPfLogging:
            for ch in self.chains:
                pf.setMcmcTreeCallback(ch.curTree.cTree, self._logFromPfModule)
                pf.setMcmcTreeCallback(ch.propTree.cTree, self._logFromPfModule)

        # Get rid of c-pointers in the copy
        theCopy.tree.deleteCStuff()
        theCopy.tree.data = None
        if self.simulate:
            theCopy.simTree.deleteCStuff()
            theCopy.simTree.data = None
        #theCopy.treePartitions._finishSplits()
        theCopy.likesFile = None
        theCopy.treeFile = None
        for ch in theCopy.chains:
            ch.curTree.deleteCStuff()
            #ch.curTree.data = None
            ch.propTree.deleteCStuff()
            #ch.propTree.data = None

        # Pickle it.
        fName = "mcmc_checkPoint_%i.%i" % (self.runNum, self.gen + 1)
        f = open(fName, 'wb')
        pickle.dump(theCopy, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    def writeProposalProbs(self, makeDict=False):
        """(Another) Pretty-print the proposal probabilities.

        See also Mcmc.writeProposalAcceptances().
        """

        nProposals = len(self.props.proposals)
        if not nProposals:
            print("Mcmc.writeProposalProbs().  No proposals (yet?).")
            return

        nAttained = [0] * nProposals
        nAccepted = [0] * nProposals
        for i in range(nProposals):
            nAttained[i] = self.props.proposals[i].nProposals[0]
            nAccepted[i] = self.props.proposals[i].nAcceptances[0]
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
            raise P4Error("bad sum of attained proposal probs. %s" % sum(probAttained))

        if not makeDict:
            spacer = ' ' * 4
            print("\nProposal probabilities (%)")
            # print "There are %i proposals" % len(self.proposals)
            print("For %i gens, from gens %i to %i, inclusive." % (
                (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen))
            print("%2s %11s %11s %11s  %11s %10s %30s %5s" % ('', 'nProposals', 
                                                         'intended(%)', 'proposed(%)',
                                                         'accepted(%)', 'tuning', 
                                                         'proposal', 'part'))
            for i in range(len(self.props.proposals)):
                print("%2i" % i, end=' ')
                p = self.props.proposals[i]
                print("   %7i " % self.props.proposals[i].nProposals[0], end=' ')
                print("   %5.1f   " % (100.0 * self.props.intended[i]), end=' ')
                print("   %5.1f    " % probAttained[i], end=' ')
                if nAttained[i]:
                    print("   %5.1f   " % (100.0 * float(nAccepted[i]) / float(nAttained[i])), end=' ')
                else:
                    print("       -   ", end=' ')
                if p.tuning == None:
                    print("       -    ", end=' ')
                elif p.tuning[0] < 2.0:
                    print("   %7.3f  " % p.tuning[0], end=' ')
                else:
                    print(" %9.3g  " % p.tuning[0], end=' ')
                print(" %27s" % p.name, end=' ')
                if p.pNum != -1:
                    print(" %3i " % p.pNum, end=' ')
                else:
                    print("   - ", end=' ')
                print()
        else:
            rd = {}
            for i in range(len(self.props.proposals)):
                p = self.props.proposals[i]
                if p.pNum != -1:
                    pname = p.name + "_%i" % p.pNum
                else:
                    pname = p.name
                rd[pname] = []
                rd[pname].append(p.nProposals[0])
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


    def lpml(self, verbose=True):
        """Calculate log pseudo marginal likelihood from CPO"""
        gm = ["Mcmc.lpml()"]
        if not self.doCpo:
            gm.append("CPO was not turned on, via m.doCpo = True")
            raise P4Error(gm)
        if not self.cpo_nSamples:
            gm.append("No cpo samples?")
            raise P4Error(gm)
        a = numpy.array(self.cpo_sumsOfInverseSiteLikes)
        a /= self.cpo_nSamples
        a = -numpy.log(a)
        theSum = a.sum()
        if verbose:
            print(f"cpo_nSamples: {self.cpo_nSamples}, lpml: {theSum}")
        return theSum
