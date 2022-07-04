import p4.pf as pf
import sys
import random
import math
import p4.func
from p4.var import var
from p4.p4exceptions import P4Error
import numpy


class BigQAndEig(object):  # not used

    def __init__(self, dim, comp, rMatrix):
        self.dim = dim
        self.comp = comp
        self.rMatrix = rMatrix
        self.bigR = numpy.zeros((dim, dim), numpy.float)
        self.bigQ = numpy.zeros((dim, dim), numpy.float)
        self.eval = None
        self.evec = None
        self.inv_evec = None

        if not comp.val:
            raise P4Error("comp.val should be set at this point.")

        for i in range(self.dim):
            if self.comp.val[i] < var.PIVEC_MIN:
                print("bad comp, %f" % self.comp.val[i])
        self.setBigR()
        self.setBigQ()
        self.eig()

    def setBigR(self):
        i = 0
        k = 0
        while i < self.dim - 2:
            j = i + 1
            while j < self.dim:
                # print i,j,k
                self.bigR[i, j] = self.rMatrix.val[k]
                self.bigR[j, i] = self.rMatrix.val[k]
                k += 1
                j += 1
            i += 1
        self.bigR[self.dim - 2, self.dim - 1] = 1.0
        self.bigR[self.dim - 1, self.dim - 2] = 1.0
        # print self.bigR

    def setBigQ(self):
        for i in range(self.dim):
            rowSum = 0.0
            for j in range(self.dim):
                if i != j:
                    self.bigQ[i, j] = self.bigR[i, j] * self.comp.val[j]
                    rowSum += self.bigQ[i, j]
            self.bigQ[i, i] = -rowSum
        self.bigQ /= -sum(diagonal(self.bigQ) * self.comp.val)
        # print self.bigQ

    def eig(self):
        self.eval, self.evec = numpy.linalg.eig(self.bigQ)
        self.evec = numpy.transpose(self.evec)
        self.inv_evec = numpy.linalg.inv(self.evec)

    def calcBigP(self, mat, brLen, nCat, rates):

        # Making a new array of zeros is faster than passing a
        # pre-allocated result (which needs zeroing), cuz zeroing
        # it takes too long.
        #result = zeros((self.dim, self.dim), numpy.float)
        for catNum in range(nCat):
            for i in range(self.dim):
                for j in range(self.dim):
                    mat[catNum, i, j] = 0.0

        #xx = exp(self.eval * brLen)
        # for i in range(self.dim):
        #    for j in range(self.dim):
        #        for k in range(self.dim):
        #            mat[i,j] += self.evec[i,k] * self.inv_evec[k,j] * xx[k]

        for catNum in range(nCat):
            if rates:
                xx = numpy.exp(self.eval * brLen * rates[catNum])
            else:
                xx = numpy.exp(self.eval * brLen)
            for i in range(self.dim):
                for j in range(self.dim):
                    for k in range(self.dim):
                        mat[catNum, i, j] += self.evec[i, k] * \
                            self.inv_evec[k, j] * xx[k]
        # print mat
        #import sys; sys.exit()

    def calcBigPCat(self, mat, brLenTimesRate):
        for i in range(self.dim):
            for j in range(self.dim):
                mat[i, j] = 0.0

        xx = numpy.exp(self.eval * brLenTimesRate)
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    mat[i, j] += self.evec[i, k] * self.inv_evec[k, j] * xx[k]
        # print mat
        #import sys; sys.exit()


class ModelPart(object):

    """Model description for one data partition."""

    def __init__(self):
        self.num = -1
        self.dim = 0
        self.comps = []
        self.rMatrices = []
        self.gdasrvs = []
        self.nGammaCat = 1
        self.pInvar = None  # only one, not a list
        self.relRate = 1.0
        self.isHet = 0
        self.symbols = None
        #self.bigQAndEigArray = None
        self.bQETneedsReset = None

        self.nCat = 0      # This is set by Tree.modelSanityCheck()

        self.ndch2 = False
        self.ndch2_leafAlpha = 2.0      # Leaf
        self.ndch2_internalAlpha = 2.0  # Internal
        self.ndch2_priorRefComp = None
        self.ndch2_writeComps = True

        self.ndrh2 = False
        self.ndrh2_leafAlpha = 2.0      # Leaf
        self.ndrh2_internalAlpha = 2.0  # Internal
        self.ndrh2_priorRefRMatrix = None
        self.ndrh2_writeRMatrices = True
        
    @property
    def nComps(self):
        return len(self.comps)

    @property
    def nRMatrices(self):
        return len(self.rMatrices)

    @property
    def nGdasrvs(self):
        return len(self.gdasrvs)

    def setCStuff(self, theModel):
        gm = ['ModelPart.setCStuff()']
        #print(gm[0])
        for mNum in range(self.nComps):
            #print(mNum)
            mt = self.comps[mNum]
            #print("comp %i = %s" % (mNum, mt.val))
            if 1:
                theSum = numpy.sum(mt.val)
                theDiff = math.fabs(1.0 - theSum)
                if theDiff > 1e-15:
                    print("Model.setCStuff().  1.0 - sum(comp.val) = %g" % theDiff)
            # for i in range(self.dim):
            #     if mt.val[i] < var.PIVEC_MIN:
            #         gm = ['ModelPart.setCStuff()']
            #         gm.append("Part %i, comp %i, val %i is too small. %g" % (
            #             self.num, mNum, i, mt.val[i]))
            #         gm.append("mt.val = %s" % mt.val)
            #         gm.append("Programming error!  This should not happen.")
            #         raise P4Error(gm)
            #     # no longer needed, as comp.val is a numpy array
            #     #pf.p4_setCompVal(theModel.cModel, self.num, mNum, i, mt.val[i])
            assert numpy.min(mt.val) >= var.PIVEC_MIN

        # Do the rMatrices
        for mNum in range(self.nRMatrices):
            mt = self.rMatrices[mNum]
            # print "Model.setCStuff()   rMatrix.spec is %s" % mt.spec
            if mt.spec == '2p':
                pf.p4_setKappa(theModel.cModel, self.num, mNum, mt.val[0])
            elif mt.free or mt.spec == 'specified':
                k = 0
                if var.rMatrixNormalizeTo1:
                    for i in range(self.dim - 1):
                        for j in range(i + 1, self.dim):
                            # print i, j, k
                            pf.p4_setRMatrixBigR(theModel.cModel, self.num, mNum, i, j, mt.val[k])
                            k += 1
                else:
                    for i in range(self.dim - 2):
                        for j in range(i + 1, self.dim):
                            # print i, j, k
                            pf.p4_setRMatrixBigR(theModel.cModel, self.num, mNum, i, j, mt.val[k])
                            k += 1

        # Do the gdasrvs
        # if self.nGammaCat > 1:
        #    for mNum in range(self.nGdasrvs):
        #        mt = self.gdasrvs[mNum]
        # pf.p4_setGdasrvVal(theModel.cModel, self.num, mNum, mt.val)   # no
        # longer needed.

        # Do the pInvar and relRate
        pf.p4_setPInvarVal(theModel.cModel, self.num, self.pInvar.val)
        pf.p4_setRelRateVal(theModel.cModel, self.num, self.relRate)

    def copyValsTo(self, otherModelPart):
        sp = self
        op = otherModelPart

        # comps
        for mtNum in range(sp.nComps):
            smt = sp.comps[mtNum]
            if smt.free:
                omt = op.comps[mtNum]
                for i in range(sp.dim):
                    omt.val[i] = smt.val[i]
            if sp.isHet:
                omt = op.comps[mtNum]
                omt.nNodes = smt.nNodes

        # rMatrices
        for mtNum in range(sp.nRMatrices):
            smt = sp.rMatrices[mtNum]
            if smt.free:
                omt = op.rMatrices[mtNum]
                if smt.spec == '2p':
                    omt.val[0] = smt.val[0]
                else:
                    if var.rMatrixNormalizeTo1:
                        for i in range(((sp.dim * sp.dim) - sp.dim) // 2):
                            omt.val[i] = smt.val[i]
                    else:
                        for i in range((((sp.dim * sp.dim) - sp.dim) // 2) - 1):
                            omt.val[i] = smt.val[i]
                if sp.isHet:
                    omt.nNodes = smt.nNodes

        # gdasrvs
        for mtNum in range(sp.nGdasrvs):
            smt = sp.gdasrvs[mtNum]
            if smt.free:
                omt = op.gdasrvs[mtNum]
                omt.val[0] = smt.val[0]

        # pInvar
        if sp.pInvar.free:
            op.pInvar.val = sp.pInvar.val

        # relRate
        op.relRate = sp.relRate

        # ndch2
        if sp.ndch2:
            op.ndch2_leafAlpha = sp.ndch2_leafAlpha
            op.ndch2_internalAlpha = sp.ndch2_internalAlpha
            for i in range(sp.dim):
                op.ndch2_priorRefComp[i] = sp.ndch2_priorRefComp[i]

        # ndrh2
        if sp.ndrh2:
            op.ndrh2_leafAlpha = sp.ndrh2_leafAlpha
            op.ndrh2_internalAlpha = sp.ndrh2_internalAlpha
            rLen = int(((sp.dim * sp.dim) - sp.dim) / 2)
            for i in range(rLen):
                op.ndrh2_priorRefRMatrix[i] = sp.ndrh2_priorRefRMatrix[i]

    def copyBQETneedsResetTo(self, otherModelPart):
        sp = self
        op = otherModelPart
        # if hasattr(sp.bQETneedsReset, 'size'):  # Can't simply ask 'if
        # sp.bQETneedsReset'
        if 1:
            for i in range(sp.nComps):
                for j in range(sp.nRMatrices):
                    op.bQETneedsReset[i][j] = sp.bQETneedsReset[i][j]

    def copyNNodesTo(self, otherModelPart):
        sp = self
        op = otherModelPart

        # comps
        for mtNum in range(sp.nComps):
            op.comps[mtNum].nNodes = sp.comps[mtNum].nNodes

        # rMatrices
        for mtNum in range(sp.nRMatrices):
            op.rMatrices[mtNum].nNodes = sp.rMatrices[mtNum].nNodes

        # gdasrvs
        for mtNum in range(sp.nGdasrvs):
            op.gdasrvs[mtNum].nNodes = sp.gdasrvs[mtNum].nNodes

    def writeEmpiricalProteinModelInPAMLFormat(self, compNum, rMatrixNum, outFileName):
        gm = ["ModelPart.writeEmpiricalProteinModelInPAMLFormat()"]

        if self.dim != 20: 
            gm.append(f"dim should be 20, got {self.dim}")
            raise P4Error(gm)
 
        # Get comp
        cVal = self.comps[compNum].val
        assert isinstance(cVal, numpy.ndarray)
        assert cVal.shape == (20,)

        # Get rMatrix
        r = self.rMatrices[rMatrixNum]
        if r.spec in var.rMatrixProteinSpecs:
            rVal = p4.func.getProteinEmpiricalModelRMatrix(r.spec, upperTriangle=False) # full r matrix
        else:
            # print(r.val)
            assert isinstance(r.val, numpy.ndarray)
            if r.val.shape not in [(190,)]:
                gm.append(f"r.val.shape is {r.val.shape}, expecting (190,)")
                raise P4Error(gm)
            rVal = numpy.zeros((20,20))
            counter = 0
            for row in range(0,20):
                for col in range(row+1,20):
                    rVal[row][col] = r.val[counter]
                    rVal[col][row] = r.val[counter]
                    counter += 1
            assert counter == 190
        assert rVal.shape == (20,20)

        # write it
        fout = open(outFileName, "w")
        for row in range(1,20):
            for col in range(0,row):   # lower triangle
                print(rVal[row][col], end=" ", file=fout)
            print(file=fout)

        print(file=fout)
        for it in range(20):
            print(cVal[it], end=" ", file=fout)
        print(file=fout)
        fout.close()




class Model(object):

    def __init__(self, nParts):
        self.parts = []
        self.cModel = None
        self.doRelRates = 0

        self.relRatesAreFree = 0
        """Boolean that says whether relRates for each data partition are free parameters."""

        self.nFreePrams = 0
        self.isHet = 0
        for i in range(nParts):
            mp = ModelPart()
            mp.num = i
            self.parts.append(mp)

    @property
    def nParts(self):
        return len(self.parts)

    def __del__(self, freeModel=pf.p4_freeModel):
        if self.cModel:
            freeModel(self.cModel)
            self.cModel = None

    def dump(self):
        """Print a summary of self.  Part by part."""

        # def _writeCharFreqToOpenFile(theCharFreq, dim, symbols, openFile, offset=23):
        # def _writeRMatrixTupleToOpenFile(theTuple, dim, openFile, offset=23):

        print("\nModel.dump().  nParts=%s" % self.nParts)

        if self.nParts == 0:
            print("nParts is zero.")
            return
        if self.doRelRates:
            if self.relRatesAreFree:
                print("Relative rates are set and free.")
            else:
                print("Relative rates are set but not free.")

        for pNum in range(self.nParts):
            mp = self.parts[pNum]
            print("\n==== Part %s, dim=%s" % (pNum, mp.dim))

            if self.doRelRates:
                print("\n  ---- relRate = %s" % mp.relRate)

            print("\n  ---- Comp Info")
            if mp.nComps == 0:
                print("    No comps.")
            else:
                print()
                print("%4s %6s %9s %6s %4s" % (
                    '', 'num', 'spec', 'free', 'symb'))
                for i in range(mp.nComps):
                    c = mp.comps[i]
                    print("%4s" % '', end=' ')
                    print("%6s" % c.num, end=' ')
                    print("%9s" % c.spec, end=' ')
                    print("%6s" % c.free, end=' ')
                    print("%4s" % c.symbol, end=' ')
                    print('')
                print('')
                for i in range(mp.nComps):
                    c = mp.comps[i]
                    print('%6s part %i, num %i' % ('', pNum, i))
                    if c.val is not None:
                        print(" " * 14, end='')
                        if mp.symbols:
                            theseSymbols = mp.symbols
                        else:
                            theseSymbols = '?' * mp.dim
                        p4.func._writeCharFreqToOpenFile(
                            c.val, mp.dim, theseSymbols, sys.stdout, offset=15)
                        print(" %s]" % c.spec)
                    else:
                        print(' ' * 14, 'No values defined.')

            #import sys; sys.exit()
            print("\n  ---- rMatrix Info")
            if mp.nRMatrices == 0:
                print("    No rMatrices.")
            else:
                print()
                print("%4s %6s %9s %6s %6s" % (
                    '', 'num', 'spec', 'free', 'symb'))

                for i in range(mp.nRMatrices):
                    c = mp.rMatrices[i]
                    print("%4s" % '', end=' ')
                    print("%6s" % c.num, end=' ')

                    if c.spec in var.rMatrixSpecs:
                        theSpec = c.spec
                    else:
                        theSpec = 'unknown'

                    print("%9s" % theSpec, end=' ')
                    print("%6s" % c.free, end=' ')
                    print("%6s" % c.symbol)
                print('')
                for i in range(mp.nRMatrices):
                    c = mp.rMatrices[i]
                    if theSpec in var.rMatrixProteinSpecs:
                        pass
                    elif theSpec == '2p':
                        print("%6s part %i, num %i" % ('', pNum, i))
                        print(" " * 15, end='')
                        print("kappa = %f" % c.val[0])
                    else:
                        print("%6s part %i, num %i" % ('', pNum, i))
                        print(" " * 15, end='')
                        if mp.dim > 2:
                            p4.func._writeRMatrixTupleToOpenFile(
                                c.val, mp.dim, sys.stdout, offset=15)
                        elif mp.dim == 2:
                            print("dim = 2.  No free values")
                        else:
                            raise P4Error(
                                "dim=%s  How did *that* happen?" % mp.dim)

            print("\n  ---- gdasrv Info.  nGammaCat=%s" % mp.nGammaCat)
            if mp.nGdasrvs == 0:
                print("    No gdasrvs.")
            else:
                print()
                print("%4s %6s %6s %6s %8s" % (
                    '', 'num', 'free', 'symb', 'val'))
                for i in range(mp.nGdasrvs):
                    c = mp.gdasrvs[i]
                    print("%4s" % '', end=' ')
                    print("%6s" % c.num, end=' ')
                    print("%6s" % c.free, end=' ')
                    print("%6s" % c.symbol, end=' ')
                    if c.val:
                        print("    %6.3f" % c.val)
                    else:
                        print("%6s" % 'None')  # This should never happen

            print("\n  ---- pInvar Info")
            c = mp.pInvar
            if c:
                print(" %6s %6s %8s" % ('', 'free', 'val'))
                print("%6s" % '', end=' ')  # indent
                print("%6s" % c.free, end=' ')
                if c.val or c.val == 0.0:
                    print("    %6.3f" % c.val)
                else:
                    print("%6s" % 'None')  # This should never happen
            else:
                print("    No pInvar.")


# def allocBigQAndEig(self):
# for pNum in range(self.nParts):
##            mp = self.parts[pNum]
##            qe = []
# for cNum in range(mp.nComps):
##                qe.append([None] * mp.nRMatrices)
# for cNum in range(mp.nComps):
# for rNum in range(mp.nRMatrices):
##                    qe[cNum][rNum] = BigQAndEig(mp.dim, mp.comps[cNum], mp.rMatrices[rNum])
##            mp.bigQAndEigArray = qe

    def writePramsProfile(self, flob, runNum):
        """Write commented lines as a key to the model prams."""

        flob.write("# Model.writePramsProfile() runNum %i\n" % runNum)
        flob.write("# \n")
        flob.write("# There are %i free parameters in this model.\n" %
                   self.nFreePrams)
        flob.write("# The parameters sampled at a given gen+1 are written on one line.\n")
        flob.write("# The first number is gen+1, followed by any changeable model parameters.\n")
        flob.write("# The number of changeable model parameters may be more than the \n")
        flob.write("#  number of free model parameters, eg for DNA composition, there are\n")
        flob.write("#  3 free parameters, but 4 changeable parameters.\n")
        flob.write("# The column number or column range for each parameter or\n")
        flob.write("#  parameter set is given twice-- the first is zero-based \n")
        flob.write("#  numbering, and the second is one-based numbering.\n")
        flob.write("# Numbers in square brackets are the number of parameters listed.\n")
        flob.write("# \n")
        nPrams = 0
        spacer1 = ' ' * 23
        spacer2 = ' ' * 12
        flob.write("#   0-based 1-based\n")
        flob.write("#   ------- -------\n")
        flob.write("#%s%i data partitions\n" % (spacer1, self.nParts))
        if self.doRelRates:
            if self.relRatesAreFree:
                flob.write("#%sRelative rates are set and free.\n" % spacer1)
            else:
                flob.write(
                    "#%sRelative rates are set but not free.\n" % spacer1)

        pramsList = []
        zeroBasedColNum = 1
        oneBasedColNum = 2
        for pNum in range(self.nParts):
            pramsList.append([])
            flob.write("#%sData Partition %i\n" % (spacer1, pNum))
            mp = self.parts[pNum]
            if self.doRelRates and self.relRatesAreFree:
                flob.write("#   %7i %7i " % (zeroBasedColNum, oneBasedColNum))
                flob.write("%srelRate[1]\n" % spacer2)
                pramsList[pNum].append(['relRate', 1])
                nPrams += 1
                zeroBasedColNum += 1
                oneBasedColNum += 1
            if mp.nComps:
                if mp.ndch2 and not mp.ndch2_writeComps:
                    pass
                else:
                    for i in range(mp.nComps):
                        if mp.comps[i].free:
                            flob.write("#   ")
                            begin = zeroBasedColNum
                            end = zeroBasedColNum + (mp.dim - 1)
                            theRangeString = "%s-%s" % (begin, end)
                            flob.write("%7s " % theRangeString)
                            begin = oneBasedColNum
                            end = oneBasedColNum + (mp.dim - 1)
                            theRangeString = "%s-%s" % (begin, end)
                            flob.write("%7s " % theRangeString)
                            flob.write("%scomp[%i]\n" % (spacer2, mp.dim))
                            pramsList[pNum].append(['comp', mp.dim])
                            nPrams += mp.dim
                            zeroBasedColNum += mp.dim
                            oneBasedColNum += mp.dim

            if mp.nRMatrices:
                if mp.ndrh2 and not mp.ndrh2_writeRMatrices:
                    pass
                else:
                    for i in range(mp.nRMatrices):
                        mt = mp.rMatrices[i]
                        if mt.free:
                            if mt.spec == '2p':
                                flob.write("#   %7i %7i " %
                                           (zeroBasedColNum, oneBasedColNum))
                                flob.write("%srMatrix[1]\n" % spacer2)
                                pramsList[pNum].append(['rMatrix', 1])
                                nPrams += 1
                                zeroBasedColNum += 1
                                oneBasedColNum += 1
                            else:
                                lenMtVal = len(mt.val)
                                flob.write("#   ")
                                begin = zeroBasedColNum
                                end = zeroBasedColNum + (lenMtVal - 1)
                                theRangeString = "%s-%s" % (begin, end)
                                flob.write("%7s " % theRangeString)
                                begin = oneBasedColNum
                                end = oneBasedColNum + (lenMtVal - 1)
                                theRangeString = "%s-%s" % (begin, end)
                                flob.write("%7s " % theRangeString)
                                flob.write("%srMatrix[%i]\n" % (spacer2, lenMtVal))
                                pramsList[pNum].append(['rMatrix', lenMtVal])
                                nPrams += lenMtVal
                                zeroBasedColNum += lenMtVal
                                oneBasedColNum += lenMtVal
            if mp.nGdasrvs:
                for i in range(mp.nGdasrvs):
                    mt = mp.gdasrvs[i]
                    if mt.free:
                        flob.write("#   %7i %7i " %
                                   (zeroBasedColNum, oneBasedColNum))
                        flob.write("%sgdasrv[1]\n" % spacer2)
                        pramsList[pNum].append(['gdasrv', 1])
                        nPrams += 1
                        zeroBasedColNum += 1
                        oneBasedColNum += 1
            if mp.pInvar and mp.pInvar.free:
                flob.write("#   %7i %7i " % (zeroBasedColNum, oneBasedColNum))
                flob.write("%spInvar[1]\n" % spacer2)
                pramsList[pNum].append(['pInvar', 1])
                nPrams += 1
                zeroBasedColNum += 1
                oneBasedColNum += 1

        flob.write("# \n")
        flob.write("# %i changeable prams in all.\n" % nPrams)
        flob.write("# \n")

        f = open("mcmc_pramsProfile_%i.py" % runNum, 'w')
        f.write(
            "# This file, 'mcmc_pramsProfile_%i.py', is used by func.summarizeMcmcPrams() and PosteriorSamples\n" % runNum)
        f.write("pramsProfile = %s\n" % pramsList)
        f.write("nPrams = %i\n" % nPrams)
        f.close()

    def writePramsLine(self, flob):
        """Write a line of model parameters for mcmc output."""

        profile1 = "\t%12.6f"
        profile2 = "\t%10.8f"
        for pNum in range(self.nParts):
            mp = self.parts[pNum]
            if self.doRelRates and self.relRatesAreFree:
                flob.write(profile1 % mp.relRate)
            if mp.nComps:
                if mp.ndch2 and not mp.ndch2_writeComps:
                    pass
                else:
                    for i in range(mp.nComps):
                        mt = mp.comps[i]
                        if mt.free:
                            for j in mt.val:
                                flob.write(profile2 % j)
            if mp.nRMatrices:
                if mp.ndrh2 and not mp.ndrh2_writeRMatrices:
                    pass
                else:
                    for i in range(mp.nRMatrices):
                        mt = mp.rMatrices[i]
                        if mt.free:
                            for j in mt.val:
                                flob.write(profile2 % j)
            if mp.nGdasrvs:
                for i in range(mp.nGdasrvs):
                    mt = mp.gdasrvs[i]
                    if mt.free:
                        flob.write(profile1 % mt.val[0])
            if mp.pInvar and mp.pInvar.free:
                flob.write(profile2 % mp.pInvar.val)
        flob.write("\n")

    def writeHypersLine(self, flob):
        """Write a line of model hyperparameters for mcmc output."""

        profile1 = "\t%12.6f"
        profile2 = "\t%10.8f"
        for pNum in range(self.nParts):
            mp = self.parts[pNum]
            if mp.ndch2:
                flob.write(profile1 % mp.ndch2_leafAlpha)
                flob.write(profile1 % mp.ndch2_internalAlpha)
            if mp.ndrh2:
                flob.write(profile1 % mp.ndrh2_leafAlpha)
                flob.write(profile1 % mp.ndrh2_internalAlpha)
        flob.write("\n")

    def writePramsHeaderLine(self, flob):
        """Write a header line of model parameters for mcmc output."""

        for pNum in range(self.nParts):
            mp = self.parts[pNum]
            if self.doRelRates and self.relRatesAreFree:
                flob.write('\trelRate.%i' % pNum)
            if mp.nComps:
                if mp.ndch2 and not mp.ndch2_writeComps:
                    pass
                else:
                    for i in range(mp.nComps):
                        mt = mp.comps[i]
                        if mt.free:
                            for j in range(len(mt.val)):
                                flob.write('\tcomp.%i.%i.%i' % (pNum, i, j))
            if mp.nRMatrices:
                if mp.ndrh2 and not mp.ndrh2_writeRMatrices:
                    pass
                else:
                    for i in range(mp.nRMatrices):
                        mt = mp.rMatrices[i]
                        if mt.free:
                            for j in range(len(mt.val)):
                                flob.write('\trMatrix.%i.%i.%i' % (pNum, i, j))
            if mp.nGdasrvs:
                for i in range(mp.nGdasrvs):
                    mt = mp.gdasrvs[i]
                    if mt.free:
                        flob.write('\tgdasrv.%i.%i' % (pNum, i))
            if mp.pInvar and mp.pInvar.free:
                flob.write('\tpInvar.%i' % pNum)
        flob.write("\n")

    def writeHypersHeaderLine(self, flob):
        """Write a header line of model hyperparameters for mcmc output."""

        for pNum in range(self.nParts):
            mp = self.parts[pNum]
            if mp.ndch2:
                flob.write('\tndch2_leafAlpha.%i' % pNum)
                flob.write('\tndch2_internalAlpha.%i' % pNum)
            if mp.ndrh2:
                flob.write('\tndrh2_leafAlpha.%i' % pNum)
                flob.write('\tndrh2_internalAlpha.%i' % pNum)
        flob.write("\n")

    def allocCStuff(self):
        complaintHead = '\nModel.allocCStuff()'
        if 0:
            print("nParts = %s %s" % (self.nParts, type(self.nParts)))
            print("doRelRates = %s %s" % (self.doRelRates, type(self.doRelRates)))
            print("relRatesAreFree = %s %s" % (self.relRatesAreFree, type(self.relRatesAreFree)))
            print("nFreePrams = %s  %s" % (self.nFreePrams, type(self.nFreePrams)))                 # a float ?!
            print("isHet = %s %s" % (self.isHet, type(self.isHet)))
            print("_rMatrixNormaizeTo1 = %s %s" % (var._rMatrixNormalizeTo1, type(var._rMatrixNormalizeTo1)))
        if self.cModel:
            gm = [complaintHead]
            gm.append("About to alloc cModel, but cModel already exists.(%s)" % self.cModel)
            raise P4Error(gm)
        self.cModel = pf.p4_newModel(self.nParts, 
                                     self.doRelRates, 
                                     self.relRatesAreFree, 
                                     int(self.nFreePrams), 
                                     self.isHet,
                                     var._rMatrixNormalizeTo1,
                                     var._PINVAR_MIN,
                                     var._PINVAR_MAX,
                                     var._KAPPA_MIN,
                                     var._KAPPA_MAX,
                                     var._GAMMA_SHAPE_MIN,
                                     var._GAMMA_SHAPE_MAX,
                                     var._PIVEC_MIN,
                                     var._PIVEC_MAX,
                                     var._RATE_MIN,
                                     var._RATE_MAX,
                                     var._RELRATE_MIN,
                                     var._RELRATE_MAX,
                                     var._BRLEN_MIN,
                                     var._BRLEN_MAX)
        #print("                                          ==== new cModel %s"  % self.cModel)
        for pNum in range(self.nParts):
            mp = self.parts[pNum]
            if 0:
                print(self.cModel, end=' ')
                print(pNum, end=' ')
                print(mp.dim, end=' ')
                print(mp.nComps, end=' ')
                print(mp.nRMatrices, end=' ')
                print(mp.nGdasrvs, end=' ')
                print(mp.nGammaCat, end=' ')
                print(mp.pInvar.free, end=' ')

            nGammaCat = mp.nGammaCat

            mp.bQETneedsReset = numpy.ones(
                (mp.nComps, mp.nRMatrices), numpy.int32)

            pf.p4_newModelPart(self.cModel,
                               pNum, mp.dim, mp.nComps, mp.nRMatrices, mp.nGdasrvs, nGammaCat,
                               mp.pInvar.free, mp.bQETneedsReset)
            for mNum in range(mp.nComps):
                mt = mp.comps[mNum]
                assert mt.val is not None  # mt.val is a numpy array
                pf.p4_newComp(self.cModel, pNum, mNum, mt.free, mt.val)
            for mNum in range(mp.nRMatrices):
                mt = mp.rMatrices[mNum]
                theSpec = None
                if mt.spec == 'ones':
                    theSpec = 100
                elif mt.spec == 'specified' or mt.spec == 'optimized':
                    theSpec = 20
                elif mt.spec == '2p':
                    theSpec = 5
<<<<<<< HEAD
                elif mt.spec == 'cpREV':
                    theSpec = 101
                elif mt.spec == 'd78':
                    theSpec = 102
                elif mt.spec == 'jtt':
                    theSpec = 103
                elif mt.spec == 'mtREV24':
                    theSpec = 104
                elif mt.spec == 'mtmam':
                    theSpec = 105
                elif mt.spec == 'wag':
                    theSpec = 106
                # the one I have, from MrBayes.  Same as phyml.
                elif mt.spec == 'blosum62':
                    theSpec = 107
                # elif mt.spec == 'blosum62b':
                #    theSpec = 108
                # elif mt.spec == 'phat70':
                #    theSpec = 109
                elif mt.spec == 'rtRev':
                    theSpec = 110
                elif mt.spec == 'tmjtt94':
                    theSpec = 111
                elif mt.spec == 'tmlg99':
                    theSpec = 112
                elif mt.spec == 'lg':
                    theSpec = 113
                elif mt.spec == 'hivb':
                    theSpec = 114
                elif mt.spec == 'mtart':
                    theSpec = 115
                elif mt.spec == 'mtzoa':
                    theSpec = 116
                elif mt.spec == 'gcpREV':
                    theSpec = 117
                elif mt.spec == 'stmtREV':
                    theSpec = 118
                elif mt.spec == 'vt':
                    theSpec = 119
                elif mt.spec == 'pmb':
                    theSpec = 120
                elif mt.spec == 'prasREV':
                    theSpec = 121
                elif mt.spec == 'gnetREV':
                    theSpec = 122
                elif mt.spec in var.rMatrixProteinSpecs:
                    theSpec = var.rMatrixProteinSpecNumberForNameDict[mt.spec]
                else:
                    gm = [complaintHead]
                    gm.append("Programming error.")
                    gm.append("Bad rMatrix spec '%s'.  Part %i" %
                              (mt.spec, pNum))
                    gm.append("Should be one of %s" % var.rMatrixSpecs)
                    raise P4Error(gm)
                # This next function sets the values of the rMatrices.
                # For protein matrices, the values are set only in c.
                # For ones and specified matrices, the rmatrix is set
                # to all ones-- any specifed values are poked in
                # later, in ModelPart.setCStuff()
                pf.p4_newRMatrix(self.cModel, pNum, mNum, mt.free, theSpec)
            for mNum in range(mp.nGdasrvs):
                mt = mp.gdasrvs[mNum]
                # print 'about to pf.p4_newGdasrv(), mt.val=%f, mt.val[0]=%f' %
                # (mt.val, mt.val[0])
                mt.c = pf.p4_newGdasrv(
                    self.cModel, pNum, mNum, mt.nGammaCat, mt.free, mt.val, mt.freqs, mt.rates)
                mt.calcRates()
        # print "finished Model.allocCStuff()"

    def setCStuff(self, partNum=None):
        complaintHead = '\nModel.setCStuff()'

        if partNum == None:
            for pNum in range(self.nParts):
                self.parts[pNum].setCStuff(self)
        else:
            self.parts[partNum].setCStuff(self)

    def restoreFreePrams(self, prams):

        complaintHead = '\nModel.restoreFreePrams()'
        if 0:
            print(complaintHead)
            print("nFreePrams = %s" % self.nFreePrams)
            print("prams=%s" % prams)
            sys.exit()
        pos = 0
        for pNum in range(self.nParts):
            mp = self.parts[pNum]

            # Do the comps
            for mNum in range(mp.nComps):
                mt = mp.comps[mNum]
                if mt.free:
                    mt.spec = 'optimized'
                    # Restore all but the last val
                    for i in range(mp.dim - 1):
                        #mt.val[i] = prams[pos]
                        pos += 1
                    # # Calculate the last val
                    # mt.val[mp.dim - 1] = 1.0 - sum(mt.val[:-1])
                    # # Make sure the vals are not too small
                    # needsNormalizing = 0
                    # theSum = 0.0
                    # for i in range(mp.dim):
                    #     if mt.val[i] < var.PIVEC_MIN:
                    #         mt.val[i] = var.PIVEC_MIN + (var.PIVEC_MIN * 0.2) + (var.PIVEC_MIN * random.random())
                    #         needsNormalizing = 1
                    #     theSum += mt.val[i]
                    # if needsNormalizing or theSum != 1.0:
                    #     for i in range(mp.dim):
                    #         mt.val[i] /= theSum
                    assert math.fabs(1.0 - numpy.sum(mt.val)) < 1.e-15
                    if numpy.min(mt.val) < var.PIVEC_MIN:
                        gm = [complaintHead]
                        gm.append("var.PIVEC_MIN is currently %g" % var.PIVEC_MIN)
                        gm.append("Got comp %i val %g ; too small" % (mNum, numpy.min(mt.val)))
                        raise P4Error(gm)
                    if 0:
                        theSum = sum(mt.val)
                        print("restoreFreePrams(). pNum %i, comp %i, sum=%g, 1.0 - sum = %g" % (
                            pNum, mNum, theSum, (1.0 - theSum)))

            # Do the rMatrices.
            for mNum in range(mp.nRMatrices):
                mt = mp.rMatrices[mNum]
                if mt.free:
                    if mt.spec == '2p':
                        mt.val[0] = prams[pos]
                        pos += 1
                    else:
                        mt.spec = 'optimized'
                        k = 0
                        for i in range(mp.dim - 2):
                            for j in range(i + 1, mp.dim):
                                # print i, j, k, pos
                                mt.val[k] = prams[pos]
                                k += 1
                                pos += 1
                        if var.rMatrixNormalizeTo1:
                            # Calculate the last value
                            mt.val[-1] = 1.0 - sum(mt.val[:-1])
                            needsNormalizing = 0
                            theSum = 0.0
                            # print mt.val, sum(mt.val)
                            for i in range(len(mt.val)):
                                if mt.val[i] < var.RATE_MIN:
                                    mt.val[i] = var.RATE_MIN + \
                                        (var.RATE_MIN * random.random())
                                    needsNormalizing = 1
                                theSum += mt.val[i]
                            if needsNormalizing or theSum != 1.0:
                                mt.val /= theSum

            # Do the gdasrvs
            if mp.nGammaCat > 1:
                for mNum in range(mp.nGdasrvs):
                    mt = mp.gdasrvs[mNum]
                    if mt.free:
                        # print 'restoreFreePrams().  mt.val=%f, mt.val[0]=%f, prams[pos]=%f' % (
                        # mt.val, mt.val[0], prams[pos])
                        if (mt.val - prams[pos]) > 1.e-10:
                            raise P4Error(
                                "restoreFreePrams. bad gdasrv non-restore.")
                        #mt.val = prams[pos]
                        pos += 1

            # Do the pInvar
            if mp.pInvar.free:
                mp.pInvar.val = prams[pos]
                pos += 1

        # Do the relRates, after the loop
        if self.relRatesAreFree:
            for pNum in range(self.nParts - 1):
                mp = self.parts[pNum]
                mp.relRate = prams[pos]
                pos += 1
            # Get the last relRate
            mp = self.parts[self.nParts - 1]
            mp.relRate = pf.p4_getRelRate(self.cModel, self.nParts - 1)

    def copyValsTo(self, otherModel):
        #complaintHead = '\nModel.copyValsTo()'

        for pNum in range(self.nParts):
            sp = self.parts[pNum]
            op = otherModel.parts[pNum]
            sp.copyValsTo(op)

    def copyBQETneedsResetTo(self, otherModel):
        for pNum in range(self.nParts):
            self.parts[pNum].copyBQETneedsResetTo(otherModel.parts[pNum])

    def copyNNodesTo(self, otherModel):
        #complaintHead = '\nModel.copyNNodesTo()'

        if self.isHet:
            for pNum in range(self.nParts):
                sp = self.parts[pNum]
                if sp.isHet:
                    op = otherModel.parts[pNum]
                    sp.copyNNodesTo(op)

    def verifyValsWith(self, otherModel):
        complaintHead = '\nModel.verifyValsWith()'

        isBad = 0
        epsilon1 = 1.e-15

        for pNum in range(self.nParts):
            sp = self.parts[pNum]
            op = otherModel.parts[pNum]

            # comps
            for mtNum in range(sp.nComps):
                smt = sp.comps[mtNum]
                if smt.free:
                    omt = op.comps[mtNum]
                    for i in range(sp.dim):
                        if math.fabs(smt.val[i] - omt.val[i]) > epsilon1:
                            print("Model.verifyValsWith()  comp vals differ.")
                            isBad = 1
                            break
            if self.isHet:
                for mtNum in range(sp.nComps):
                    if op.comps[mtNum].nNodes != sp.comps[mtNum].nNodes:
                        print("Model.verifyValsWith()  nNodes differ.")
                        isBad = 1
                        break

            # rMatrices
            for mtNum in range(sp.nRMatrices):
                smt = sp.rMatrices[mtNum]
                if smt.free:
                    omt = op.rMatrices[mtNum]
                    if smt.spec == '2p':
                        if math.fabs(smt.val[0] - omt.val[0]) > epsilon1:
                            isBad = 1
                            break
                    else:
                        # for i in range((((sp.dim * sp.dim) - sp.dim) // 2) - 1):
                        for i in range(len(smt.val)):
                            if math.fabs(smt.val[i] - omt.val[i]) > epsilon1:
                                isBad = 1
                                break
            if self.isHet:
                for mtNum in range(sp.nRMatrices):
                    if op.rMatrices[mtNum].nNodes != sp.rMatrices[mtNum].nNodes:
                        isBad = 1
                        break

            # gdasrvs
            for mtNum in range(sp.nGdasrvs):
                smt = sp.gdasrvs[mtNum]
                if smt.free:
                    omt = op.gdasrvs[mtNum]
                    if math.fabs(smt.val[0] - omt.val[0]) > epsilon1:
                        isBad = 1
                        break

            # pInvar
            if sp.pInvar.free:
                if math.fabs(sp.pInvar.val - op.pInvar.val) > epsilon1:
                    isBad = 1
            if isBad:
                break

            # relRate
            if self.doRelRates and self.relRatesAreFree:
                if math.fabs(sp.relRate - op.relRate) > epsilon1:
                    isBad = 1
            if isBad:
                break

            if sp.ndch2:
                if math.fabs(sp.ndch2_leafAlpha - op.ndch2_leafAlpha) > epsilon1:
                    isBad = 1
                    break
                if math.fabs(sp.ndch2_internalAlpha - op.ndch2_internalAlpha) > epsilon1:
                    isBad = 1
                    break
                for i in range(sp.dim):
                    if math.fabs(sp.ndch2_priorRefComp[i] - op.ndch2_priorRefComp[i]) > epsilon1:
                        isBad = 1
                        break
            if sp.ndrh2:
                if math.fabs(sp.ndrh2_leafAlpha - op.ndrh2_leafAlpha) > epsilon1:
                    isBad = 1
                    break
                if math.fabs(sp.ndrh2_internalAlpha - op.ndrh2_internalAlpha) > epsilon1:
                    isBad = 1
                    break
                rLen = int(((sp.dim * sp.dim) - sp.dim) / 2)
                for i in range(rLen):
                    if math.fabs(sp.ndrh2_priorRefRMatrix[i] - op.ndrh2_priorRefRMatrix[i]) > epsilon1:
                        isBad = 1
                        break
                

            # bQETneedsReset
            # if hasattr(sp.bQETneedsReset, 'size'):  # Can't simply ask 'if
            # sp.bQETneedsReset'
            if 0:
                for i in range(sp.nComps):
                    for j in range(sp.nRMatrices):
                        # integers
                        if op.bQETneedsReset[i][j] != sp.bQETneedsReset[i][j]:
                            print(complaintHead)
                            print("Mis-matched bQETneedsReset.")
                            print("self.bQETneedsReset is")
                            print(sp.bQETneedsReset)
                            print("other.bQETneedsReset is")
                            print(op.bQETneedsReset)
                            isBad = 1
                            break
            if 1:
                # bQETneedsReset is a numpy array of ints
                assert isinstance(op.bQETneedsReset, numpy.ndarray)
                # ret is an array of booleans.
                ret = op.bQETneedsReset != sp.bQETneedsReset
                if numpy.any(ret):
                    print(complaintHead)
                    print("Mis-matched bQETneedsReset.  Part %i" % pNum)
                    print("self.bQETneedsReset is")
                    print(sp.bQETneedsReset)
                    print("other.bQETneedsReset is")
                    print(op.bQETneedsReset)
                    print("positions differing:")
                    print(ret)
                    isBad = 1
                    break

        if isBad:
            print(complaintHead)
            print("    Mis-matched model prams.")
            return var.DIFFERENT
        return var.SAME

    def getBigQ(self, pNum=0, compNum=0, rMatrixNum=0):
        """Returns a dim * dim numpy array with the bigQ """

        mp = self.parts[pNum]
        c = mp.comps[compNum]
        r = mp.rMatrices[rMatrixNum]
        a = numpy.zeros((mp.dim, mp.dim), numpy.float)
        pf.getBigQ(self.cModel, mp.dim, pNum, compNum, rMatrixNum, a)
        return a
