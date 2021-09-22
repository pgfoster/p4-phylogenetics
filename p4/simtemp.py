import math
import statistics
import random
from collections import deque
import sys

class SimTempTemp(object):
    def __init__(self):
        self.tempNum = 0
        self.temp = 0.0
        self.tempDiff = None
        self.logPiDiff = None
        self.occupancy = 0
        self.nProposalsUp = 0
        self.nAcceptedUp = 0
        self.nProposalsDn = 0
        self.nAcceptedDn = 0

        self.short_nProposalsUp = 0
        self.short_nAcceptedUp = 0
        #self.short_nProposalsDn = 0
        #self.short_nAcceptedDn = 0

class SimTemp(object):
    def __init__(self, mcmc):
        self.mcmc = mcmc
        self.nTemps = 2
        self.temps = []
        self.curTemp = 0
        self.nTempChangeProposals = 0
        self.nTempChangeAccepts = 0
        self.tNums = []
        # self.tempChangeProposeFreq = 1
        self.tNumSampleSize = None
        self.tNumSample = None

    def setTempDiffsFromTemps(self):
        for tNum in range(self.nTemps - 1):
            tmpA = self.temps[tNum]
            tmpB = self.temps[tNum + 1]
            tmpA.tempDiff = tmpB.temp - tmpA.temp
            assert tmpA.tempDiff > 0.0

    def setTempsFromTempDiffs(self):
        assert self.temps[0].temp == 0.0
        for tNum in range(self.nTemps - 1):
            tmpA = self.temps[tNum]
            tmpB = self.temps[tNum + 1]
            assert tmpA.tempDiff > 0.0
            tmpB.temp = tmpA.temp + tmpA.tempDiff

    def dumpTemps(self, flob=sys.stdout, zeroCountsAfter=True):
        """Make a nice table of SimTempTemp objects

        Arg flob should be an open file-like object.
        """

        #print(f"Temperatures at gen {self.mcmc.gen}")
        print("%4s %10s %10s %12s %10s %10s %10s %10s %10s" % (
            "indx", "temp", "tempDiff", "logPiDiff", "occupancy", "nPropsUp", "accptUp", "nPropsDn", "accptDn"), 
              file=flob)
 
        for tNum, tmp in enumerate(self.temps):

            # If it is the hottest temp, don't print tempDiff or logPiDiff
            if tNum == self.nTemps - 1:
                print("%4s %10.3f %10s %12s %10i" % (
                    tNum,
                    tmp.temp, 
                    "-",       # tempDiff
                    "-",       # logPiDiff 
                    tmp.occupancy), file=flob, end='')
                # No proposals up
                print(" %10s %10s" % ("-", "-"), file=flob, end='')
            else:
                # print(f"tmp.tempNum, {tmp.tempNum}")
                # print(f"tmp.temp, {tmp.temp}")
                # print(f"tmp.tempDiff, {tmp.tempDiff}")
                # print(f"tmp.logPiDiff, {tmp.logPiDiff}")
                

                print("%4s %10.3f %10.3f %12.4f %10i" % (
                    tNum,
                    tmp.temp,
                    tmp.tempDiff,
                    tmp.logPiDiff,
                    tmp.occupancy), file=flob, end='')
                if tmp.nProposalsUp:
                    print(" %10i %10.4f" % (tmp.nProposalsUp, (tmp.nAcceptedUp/tmp.nProposalsUp)), file=flob, end='')
                else:  # don't divide by zero
                    print(" %10s %10s" % ("-", "-"), file=flob, end='')

            if tNum == 0:  # no props down
                print(" %10s %10s" % ("-", "-"), file=flob, end='')
            else:
                if tmp.nProposalsDn:
                    print(" %10i %10.4f" % (tmp.nProposalsDn, (tmp.nAcceptedDn/tmp.nProposalsDn)), file=flob, end='')
                else:  # don't divide by zero
                    print(" %10s %10s" % ("-", "-"), file=flob, end='')
            print(file=flob)

        # Zero
        if zeroCountsAfter:
            for tmp in self.temps:
                tmp.occupancy = 0
                tmp.nProposalsUp = 0
                tmp.nAcceptedUp = 0
                tmp.nProposalsDn = 0
                tmp.nAcceptedDn = 0

    def dumpTempsOneLine(self, flob=sys.stdout, zeroCountsAfter=True):
        """Single-line output of SimTempTemp objects

        Arg flob should be an open file-like object.
        """

        #print(f"Temperatures at gen {self.mcmc.gen}")
        print("%4s %10s %10s %12s %10s %10s %10s %10s %10s" % (
            "indx", "temp", "tempDiff", "logPiDiff", "occupancy", "nPropsUp", "accptUp", "nPropsDn", "accptDn"), 
              file=flob)
 
        for tNum, tmp in enumerate(self.temps):

            # If it is the hottest temp, don't print tempDiff or logPiDiff
            if tNum == self.nTemps - 1:
                print("%4s %10.3f %10s %12s %10i" % (
                    tNum,
                    tmp.temp, 
                    "-",       # tempDiff
                    "-",       # logPiDiff 
                    tmp.occupancy), file=flob, end='')
                # No proposals up
                print(" %10s %10s" % ("-", "-"), file=flob, end='')
            else:
                # print(f"tmp.tempNum, {tmp.tempNum}")
                # print(f"tmp.temp, {tmp.temp}")
                # print(f"tmp.tempDiff, {tmp.tempDiff}")
                # print(f"tmp.logPiDiff, {tmp.logPiDiff}")
                

                print("%4s %10.3f %10.3f %12.4f %10i" % (
                    tNum,
                    tmp.temp,
                    tmp.tempDiff,
                    tmp.logPiDiff,
                    tmp.occupancy), file=flob, end='')
                if tmp.nProposalsUp:
                    print(" %10i %10.4f" % (tmp.nProposalsUp, (tmp.nAcceptedUp/tmp.nProposalsUp)), file=flob, end='')
                else:  # don't divide by zero
                    print(" %10s %10s" % ("-", "-"), file=flob, end='')

            if tNum == 0:  # no props down
                print(" %10s %10s" % ("-", "-"), file=flob, end='')
            else:
                if tmp.nProposalsDn:
                    print(" %10i %10.4f" % (tmp.nProposalsDn, (tmp.nAcceptedDn/tmp.nProposalsDn)), file=flob, end='')
                else:  # don't divide by zero
                    print(" %10s %10s" % ("-", "-"), file=flob, end='')
            print(file=flob)


    def tuneTempsAndLogPi(self):

        # We get the occupancy sample from self.tNumSample
        # We get the acceptance up from tmp.short_nProposalsUp and tmp.short_nAcceptedUp

        target = self.tNumSampleSize / self.nTemps
        occs = [self.tNumSample.count(tNum)/target for tNum in range(self.nTemps)]
        # print("occs ", occs)


        # print(f"gen:{self.mcmc.gen} props:{self.temps[0].short_nProposalsUp}", end=" ")
        # print("accepts:{self.temps[0].short_nAcceptedUp}", end=" ")

        acceptsUp = []
        for tmp in self.temps[:-1]:
            if tmp.short_nProposalsUp > (target * 0.5):
                aU = tmp.short_nAcceptedUp / tmp.short_nProposalsUp
                acceptsUp.append(aU)
            else:
                acceptsUp.append(None)
            tmp.short_nAcceptedUp = 0
            tmp.short_nProposalsUp = 0
        #print("acceptsUp: ", acceptsUp)


        if self.mcmc.simTemp_doLogTemps:
            print(f"{self.mcmc.gen:9d}", file=self.mcmc.simTemp_tempsFlob, end='')
            for tmp in self.temps[:-1]:
                print(f"{tmp.tempDiff:13.2f}", file=self.mcmc.simTemp_tempsFlob, end='')
            print(f"{self.temps[-1].temp:13.2f}", file=self.mcmc.simTemp_tempsFlob, end='')
            for tmp in self.temps[:-1]:
                print(f"{tmp.logPiDiff:13.3f}", file=self.mcmc.simTemp_tempsFlob, end='')
            for occ in occs:
                print(f"{occ:10.3f}", file=self.mcmc.simTemp_tempsFlob, end='')
            
            # for tmp in self.temps[:-1]:
            #     if tmp.nProposalsUp:
            #         acc = tmp.nAcceptedUp/tmp.nProposalsUp
            #         print(f"{acc:4.2f}", file=self.mcmc.simTemp_tempsFlob, end='')
            #     else:
            #         print(f"{"-":4s}", file=self.mcmc.simTemp_tempsFlob, end='')                    

            print("", file=self.mcmc.simTemp_tempsFlob, end='\n')


        # Rules:
        # - if the hot occupancy is too high, raise the logPiDiff
        # - if the hot occupancy is too low, lower the logPiDiff

        # Does this work? -- Not done yet
        # If the acceptance up (of a particular tmp) is too high, raise the temp diff of that tmp
        # If the acceptance (of a particular tmp) is too low, lower the temp diff of that tmp

        # - If the number of temps is more than 2 ---
        #     - If the cold occupancy is too high, raise all middle temps  (including the highest, so it is not fixed)
        #     - If the cold occupancy is too low, lower the middle temps  (including the highest)
        #   - If middle occupancies are off, adjust the temperature individually
        #     - If occupancy is too high, lower the temperature
        #     - if occupancy is too low, raise the temperature


        # Now in detail ---
        # - if the hot occupancy is too high, raise the logPiDiff
        # - if the hot occupancy is  too low, lower the logPiDiff
        fudge = 0.1
        factor = occs[-1]
        if factor > 2.0:
            factor = 2.0
        elif factor < 0.5:
            factor = 0.5
        factor -= 1.0     # center on zero
        factor *= fudge   # scale it
        factor += 1.0     # center it on 1.
        for tNum in range(self.nTemps - 1):
            tmp = self.temps[tNum]
            old = tmp.logPiDiff
            tmp.logPiDiff *= factor
        if 1:
            if occs[-1] < 1.0:
                print(f"&-  hotOccs too low {occs[-1]:4.2f}; multiply all logPiDiffs by factor {factor:4.2f}")
            else:
                print(f"&+  hotOccs too high{occs[-1]:4.2f}; multiply all logPiDiffs by factor {factor:4.2f}")

        # TODO
        # If the acceptance is too high, raise the temp
        # If the acceptance is too low, lower the temp

        # - If the number of temps is more than 2 ---
        if self.nTemps > 2:
            #   - If the cold occupancy is too high, raise all middle temps (including the highest, so it is not fixed)
            #   - If the cold occupancy is too low, lower all the middle temps (including the highest)
            if 1:
                # This adjusts the first tempDiff, so it affects all the temps (except the coldest of course, which is zero) 
                originalMaxTemp = self.temps[-1].temp
                tmp = self.temps[0]
                thisFudge = 0.1
                factor = occs[0]
                if factor > 2.0:
                    factor = 2.0
                elif factor < 0.5:
                    factor = 0.5
                factor -= 1.0
                factor *= thisFudge
                factor += 1.0
                old = tmp.tempDiff
                tmp.tempDiff *= factor
                self.setTempsFromTempDiffs()

                # Fixed max temp?
                if 0:
                    if originalMaxTemp > self.temps[-2].temp:
                        #print(f"aaA tempDiff[-2] {self.temps[-2].tempDiff:8.2f}  {self.temps[-1].temp:8.2f} ")
                        self.temps[-1].temp = originalMaxTemp
                        self.setTempDiffsFromTemps()
                        #print(f"aaB tempDiff[-2] {self.temps[-2].tempDiff:8.2f}  {self.temps[-1].temp:8.2f} ")
                if 1:
                    if occs[0] < 1.0:
                        print("%- cold occupancy is too low. ", end=' ')
                    else:
                        print("%+ cold occupancy is too high. ", end=' ')
                
                    print(f"coldOcc:{occs[0]:5.3f}; ", end=' ') 
                    print(f"old tempDiff[0]:{old:7.3f} factor:{factor:4.2f} new tempDiff[0]:{tmp.tempDiff:7.3f}")

            # - If middle occupancies are off, adjust the temperature individually
            #   - If occupancy is too high, lower the temperature ( by lowering the tNum-1 tempDiff)
            #   - if occupancy is too low, raise the temperature
            if 1:
                originalMaxTemp = self.temps[-1].temp
                thisFudge = 0.1
                for tNum in range(1, self.nTemps - 1):
                    factor = occs[tNum]
                    if factor > 2.0:
                        factor = 2.0
                    elif factor < 0.5:
                        factor = 0.5
                    factor -= 1.0
                    factor *= thisFudge
                    factor += 1.0
                    rfactor = 1.0/factor
                    old = self.temps[tNum].temp
                    td = self.temps[tNum - 1].tempDiff
                    
                    td *= rfactor
                    assert td > 0.0
                    #print(f"xxx tNum:{tNum} {self.temps[-1].tempNum} {self.temps[-2].tempDiff:8.2f} {self.temps[-1].temp:8.2f} ", end='')
                    self.temps[tNum - 1].tempDiff = td
                    self.setTempsFromTempDiffs()
                    #print(f"yyy {self.temps[-1].temp} ")

                    if 1:
                        if occs[tNum] < 1.0:
                            print(f"=+ Temp {tNum} occupation too low (rfactor:{rfactor:5.2f}), so raise ", end='')
                        else:
                            print(f"=- Temp {tNum} occupation too high (rfactor:{rfactor:5.2f}), so lower ", end='')
                        print(f"(old tempDiff:{self.temps[tNum - 1].tempDiff:5.2f} new tempDiff {td:5.2f}", end='')
                        print(f"temp from {old:5.2f} to {self.temps[tNum].temp:5.2f}")

                # Fixed max temp?
                if 0:
                    if originalMaxTemp > self.temps[-2].temp:
                        #print(f"bbA tempDiff[-2] {self.temps[-2].tempDiff:8.2f}  {self.temps[-1].temp:8.2f} ")
                        self.temps[-1].temp = originalMaxTemp
                        self.setTempDiffsFromTemps()
                        #print(f"bbB tempDiff[-2] {self.temps[-2].tempDiff:8.2f}  {self.temps[-1].temp:8.2f} ")


    def proposeTempChange(self):

        # Method and symbols are from:

        # Geyer, C. J., and Thompson, E. A. (1995). Annealing Markov
        # chain Monte Carlo with applications to ancestral
        # inference. Journal of the American Statistical Association,
        # 90, 909–920.

        i = self.curTemp   # index
        nTemps = self.nTemps
        ch = self.mcmc.chains[0]

        if i == 0:
            j = 1
            qij = 1.
            if nTemps == 2:
                qji = 1.0
            else:
                qji = 0.5
        elif i == nTemps - 1:
            j = nTemps - 2
            qij = 1.
            if nTemps == 2:
                qji = 1.0
            else:
                qji = 0.5
        else:
            if random.random() < 0.5:
                j = i - 1
            else:
                j = i + 1
            qij = 0.5
            if j == 0:
                qji = 1.0
            elif j == nTemps - 1:
                qji = 1.0
            else:
                qji = 0.5

        # tmp's are SimTempTemp objects
        tmpI = self.temps[i]
        tmpJ = self.temps[j]

        beta_i = 1.0 / (1.0 + tmpI.temp)
        beta_j = 1.0 / (1.0 + tmpJ.temp)

        # h_i_Atx is notation from Geyer and Thompson.  I want the log form
        # print("curr log like %f, beta_i %f, beta_j %f" % (ch.curTree.logLike, beta_i, beta_j))
        log_h_i_Atx = ch.curTree.logLike * beta_i 
        log_h_j_Atx = ch.curTree.logLike * beta_j
        logTempRatio = log_h_j_Atx - log_h_i_Atx
        #logTempRatio = ch.log_density_cur * (beta_j - beta_i)

        #logPseudoPriorRatio =  math.log(tmpJ.pi / tmpI.pi)
        if i < j:
            logPseudoPriorRatio = -tmpI.logPiDiff
        else:
            logPseudoPriorRatio = tmpJ.logPiDiff

        #logPseudoPriorRatio =  tmpJ.logPi - tmpI.logPi
        # print("to calc logPseudoPriorRatio:", tmpJ.logPi, "-", tmpI.logPi, "=", logPseudoPriorRatio)
        logHastingsRatio = math.log(qji / qij)

        if 0:
            print("hastingsRatio gen %i: i=%i, j=%i, hastingsRatio %f, logHastingsRatio %f" % (
                self.gen, i, j, (qji/qij), logHastingsRatio))


        # Geyer and Thompson use "r".  I want the log form
        log_r = logTempRatio + logPseudoPriorRatio + logHastingsRatio


        acceptMove = False
        if log_r < -100.0:
            acceptMove = False
        elif log_r >= 0.0:
            acceptMove = True
        else:
            r = math.exp(log_r)
            if random.random() < r:
                acceptMove = True
        #print("log_r %f, acceptMove %s" % (log_r, acceptMove))

        if acceptMove:
            self.curTemp = j
            self.nTempChangeAccepts += 1
        self.nTempChangeProposals += 1
        if j > i:
            tmpI.nProposalsUp += 1
            tmpI.short_nProposalsUp += 1
            if acceptMove:
                tmpI.nAcceptedUp += 1
                tmpI.short_nAcceptedUp += 1
        else:
            tmpI.nProposalsDn += 1
            if acceptMove:
                tmpI.nAcceptedDn += 1

