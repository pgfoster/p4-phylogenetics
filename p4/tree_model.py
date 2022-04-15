import sys
import string
import math
import copy
import os
import p4.func
import time
import glob
from p4.var import var
from p4.p4exceptions import P4Error
from p4.node import Node, NodePart, NodeBranch, NodeBranchPart
from p4.distancematrix import DistanceMatrix

import numpy
import p4.pf as pf
from p4.model import Model
from p4.data import Data
from p4.alignment import Part
import random

if True:
    @property
    def data(self):
        """(property) The data object"""
        return self._data

    @data.setter
    def data(self, theData):
        """Sets self.data, and self.nParts"""

        # Two cases.  Either
        #   1.  Self already has a data.
        #                - deleteCStuff(), and replace the old data with the new data.
        #   2.  Self has no data.
        #     a.  Has a model.  So just check for compatibility.
        #     b.  Has no model, so has never seen a data before.

        complaintHead = 'Tree.data (property setter)'
        if self.name:
            complaintHead += ", tree '%s'" % self.name
        if self.fName:
            complaintHead += ", file %s" % self.fName
        gm = [complaintHead]

        if isinstance(theData, Data) or theData == None:
            pass
        else:
            gm.append("Set data only to a Data object, or None, ok?")
            raise P4Error(gm)

        if self.data or self.model:
            # Normally there won't be anything to delete, but you never know...
            self.deleteCStuff()

        if not theData:
            self._data = None
            return

        if self.model:

            # We have seen a data object before (otherwise we would not
            # have been able to set the model).  Check for compatibility.

            # print "_setData() here.   self.model exits."

            if not self.taxNames:
                gm.append(
                    "Self has model, but no taxNames.  Programming error.")
                raise P4Error(gm)

            # Check for same number of taxa
            treeNTax = len(self.taxNames)
            dataNTax = len(theData.taxNames)
            if treeNTax != dataNTax:
                gm.append("The number of taxa in the tree (%s)" % treeNTax)
                gm.append("is not the same as in the data (%s)" % dataNTax)
                raise P4Error(gm)

            # Check for mis-matched taxNames
            isBad = 0
            for tn in self.taxNames:
                if tn not in theData.taxNames:
                    isBad = 1
                    break
            for tn in theData.taxNames:
                if tn not in self.taxNames:
                    isBad = 1
                    break
            if isBad:
                gm.append("TaxName mismatch between the tree and the data.")
                gm.append("Here they are, sorted to show mis-matches.")
                gm.append("    %25s    %25s" % ('data', 'tree'))
                self.taxNames.sort()
                theData.taxNames.sort()
                for i in range(len(self.taxNames)):
                    if theData.taxNames[i] == self.taxNames[i]:
                        gm.append("    %25s    %25s" %
                                  (theData.taxNames[i], self.taxNames[i]))
                    else:
                        gm.append("    %25s    %25s  ***" %
                                  (theData.taxNames[i], self.taxNames[i]))
                raise P4Error(gm)

            # Same number of parts
            if len(theData.parts) != self.model.nParts:
                gm.append("nParts mis-match.  len(theData.parts)=%i, model.nParts=%i" % (
                    len(theData.parts), self.model.nParts))
                raise P4Error(gm)

            # Check dims and symbols in the parts
            for pNum in range(self.model.nParts):
                if theData.parts[pNum].dim != self.model.parts[pNum].dim:
                    gm.append("Parts dim mis-match.")
                    raise P4Error(gm)
                if theData.parts[pNum].symbols != self.model.parts[pNum].symbols:
                    gm.append("Parts symbols mis-match.")
                    raise P4Error(gm)

            # Set seqNum
            for n in self.iterLeavesNoRoot():
                if n.seqNum != self.taxNames.index(n.name):
                    gm.append("seqNums do not match up with taxNames.")
                    raise P4Error(gm)

            self._data = theData

        else:
            # When you do this method, _setData(), self gets a suitable
            # model.  We have here no model, so we may have never seen a
            # data before.  Or we might have just lost the model.

            # Check for same number of taxa
            treeNTax = 0
            treeTaxNames = []
            for n in self.iterNodes():
                if n.isLeaf:
                    treeNTax += 1
                    treeTaxNames.append(n.name)
            dataNTax = len(theData.taxNames)
            if treeNTax != dataNTax:
                gm.append("The number of taxa in the tree (%s)" % treeNTax)
                gm.append("is not the same as in the data (%s)" % dataNTax)
                raise P4Error(gm)

            # Check for mis-matched taxNames
            isBad = 0
            for tn in treeTaxNames:
                if tn not in theData.taxNames:
                    isBad = 1
                    break
            for tn in theData.taxNames:
                if tn not in treeTaxNames:
                    isBad = 1
                    break
            if isBad:
                gm.append("TaxName mismatch between the tree and the data.")
                gm.append("Here they are, sorted to show mis-matches.")
                gm.append("    %25s  %25s" % ('data', 'tree'))
                treeTaxNames.sort()
                theData.taxNames.sort()
                for i in range(len(treeTaxNames)):
                    if theData.taxNames[i] == treeTaxNames[i]:
                        gm.append("    %25s  %25s" %
                                  (theData.taxNames[i], treeTaxNames[i]))
                    else:
                        gm.append("*** %25s  %25s" %
                                  (theData.taxNames[i], treeTaxNames[i]))
                raise P4Error(gm)

            # attach
            self.taxNames = theData.taxNames
            self._data = theData
            #self.nParts = len(theData.parts)

            # Now that nParts is known ...
            # print "_setData.  len(theData.parts) = %s" % len(theData.parts)
            # calls self.deleteCStuff()
            self.model = Model(len(theData.parts))
            for n in self.iterNodes():
                if n.parts:
                    n.parts = []
                for i in range(self.model.nParts):
                    n.parts.append(NodePart())
            for n in self.iterNodes():
                if n != self.root:
                    n.br.parts = []
                    for i in range(self.model.nParts):
                        n.br.parts.append(NodeBranchPart())

            # Set modelPart dims and symbols
            for pNum in range(self.model.nParts):
                self.model.parts[pNum].dim = theData.parts[pNum].dim
                self.model.parts[pNum].symbols = theData.parts[pNum].symbols

            # There is intentionally no default pInvar, forcing the user to be
            # explicit.
            for p in self.model.parts:
                p.pInvar = None

            # Set seqNum
            for n in self.iterNodes():
                if n.isLeaf:
                    n.seqNum = self.taxNames.index(n.name)

    @property
    def model(self):
        """(property) The model object"""
        return self._model

    @model.setter
    def model(self, theModel):
        gm = ['Tree.model (property setter)']
        # print gm[0]
        # print "    Got '%s'" % theModel
        if isinstance(theModel, Model) or theModel == None:
            pass
        else:
            gm.append("Attempt to set Tree.model to '%s'.  " % theModel)
            gm.append(
                "Don't set the model to anything other than 'None' or a Model, ok?  ")
            gm.append("(And generally the user only sets it to None.)  ")
            raise P4Error(gm)
        if self.model or self.data:
            self.deleteCStuff()
            # print 'Tree._setModel()  finished deleteCStuff()'
        self._model = theModel

    @model.deleter
    def model(self):
        gm = ['Tree.model']
        gm.append("Caught an attempt to delete self.model, but")
        gm.append("self.model is a property, so you shouldn't delete it.")
        gm.append("But you can set it to None if you like.")
        raise P4Error(gm)


    def _checkModelThing(self, partNum, symbol, complaintHead):
        gm = [complaintHead]
        if not self.data:
            gm.append("No data.  Set the data first.")
            raise P4Error(gm)

        if not self.model:
            # When you _setData(), a model object of suitable dimensions
            # is made and attached to self.  If we have got here, it is
            # because the model has subsequently been lost.  So just
            # re-instate it.
            #self._setData(self.data)
            gm.append("No model.  Set the model first.  Or FIXME?")
            raise P4Error(gm)

        if partNum < 0 or partNum >= self.model.nParts:
            gm.append("Zero-based partNum (%s) is out of range (of %s parts)" %
                      (partNum, self.model.nParts))
            raise P4Error(gm)

        if symbol:
            if not isinstance(symbol, str) or len(symbol) != 1:
                gm.append("Symbols must be 1-length strings.")
                raise P4Error(gm)
            if symbol == '?':
                gm.append("Got assigned text drawing symbol '?'.")
                gm.append("Don't use it-- it is reserved for missing model components")
                raise P4Error(gm)

    def newComp(self, partNum=0, free=0, spec='empirical', val=None, symbol=None):
        """Make, attach, and return a new Comp object.

        The arg *spec* should be a string, one of::

          'equal'          no val
          'empirical'      no val
          'specified'      val=[aList]
          'wag', etc       no val
             (ie one of the empirical protein models, including
             cpREV, d78, jtt, mtREV24, mtmam, wag, etc)

        If spec='specified', then you specify dim or dim-1 values in a
        list as the 'val' arg.

        This method returns a Comp object, which you can ignore if it is a
        tree-homogeneous model.  However, if it is a tree-hetero model
        then you may want to get that Comp object so that you can place
        it on the tree explicitly with setModelComponentOnNode(), like this::

            c0 = newComp(partNum=0, free=1, spec='empirical')
            c1 = newComp(partNum=0, free=1, spec='empirical')
            myTree.setModelComponentOnNode(c0, node=myTree.root, clade=1)
            myTree.setModelComponentOnNode(c1, node=5, clade=1)
            myTree.setModelComponentOnNode(c1, node=18, clade=0)

        Alternatively, you can simply let p4 place them randomly::

            newComp(partNum=0, free=1, spec='empirical')
            newComp(partNum=0, free=1, spec='empirical')
            myTree.setModelComponentsOnNodesRandomly()

        Calculation of probability matrices for likelihood calcs etc are
        wrong when there are any comp values that are zero, so that is not
        allowed.  Any zeros are converted to var.PIVEC_MIN, which is 1e-13
        this week.  Hopefully close enough to zero for you.
        """

        gm = ['Tree.newComp()']

        self._checkModelThing(partNum, symbol, gm[0])
        if self.model.cModel:
            self.deleteCStuff()
        mt = Comp()
        mt.partNum = partNum
        mt.free = free

        # spec
        if spec not in var.compSpecs:
            gm.append("The spec should be one of %s" % var.compSpecs)
            raise P4Error(gm)
        mt.spec = spec

        mt.num = len(self.model.parts[partNum].comps)
        if symbol:
            mt.symbol = symbol
        else:
            try:
                mt.symbol = var.modelSymbols[mt.num]
            except IndexError:
                gm.append("You have asked for modelSymbol number %i" % mt.num)
                gm.append("There are only %i in var.modelSymbols." % len(var.modelSymbols))
                raise P4Error(gm)

        self.model.parts[partNum].comps.append(mt)

        # assign val
        dim = self.model.parts[partNum].dim
        if spec == 'equal':
            mt.val = numpy.ones(dim, float) / dim
        elif spec == 'empirical':
            assert mt.val is None
        elif spec == 'specified':
            # print(f"Got val {val}, type {type(val)}")
            # if val == None or val == []:
            #     gm.append("Specified comp, but no val.")
            #     raise P4Error(gm)
            try:
                val = list(val)
            except TypeError:
                gm.append("Comp is 'specified', but bad 'val' arg.")
                gm.append("The 'val' arg should be a list or tuple.")
                raise P4Error(gm)
            if len(val) == dim or len(val) == dim - 1:
                pass
            else:
                gm.append("Bad length for val arg (%i).  Should be dim or dim-1 long." % len(val))
                gm.append("(Dim for this part is %i)" % dim)
                raise P4Error(gm)

            # I allow val's of dim or dim-1.
            if len(val) == dim - 1:
                lastVal = 1.0 - sum(val)
                if lastVal > 0.0:
                    val = val + [1.0 - sum(val)]
                else:
                    gm.append("Bad comp vals %s" % val)
                    gm.append("sum to 1.0 or more.")
                    raise P4Error(gm)
            else:  # len = dim
                theSum = sum(val)
                theDiff = math.fabs(theSum - 1.0)

                # How big to make the delta?  With reasonably good, normalized
                # protein comps (where all the values had just been divided by
                # the total, so it should have summed to 1.0 at that point) I
                # kept getting 1.1e-16.  So make it 5.e-16.  No, too small.
                # With protein I am getting diffs of 8.88178e-16, so make the
                # cutoff 9e-16

                if theDiff > 9.e-16:  # 1e-17 was too small for protein
                    gm.append("Bad comp vals %s" % val)
                    gm.append("does not sum to 1.0")
                    gm.append("The sum = %f" % theSum)
                    gm.append("abs(1.0 - theSum) = %g" % theDiff)
                    raise P4Error(gm)

            # Are any specified values less than PIVEC_MIN?
            needsNormalizing = 0
            for i in range(len(val)):
                thisVal = val[i]
                if thisVal < var.PIVEC_MIN:
                    print(gm[0])
                    print("    Specifying a comp of zero for a character is not allowed.")
                    print("    The minimum is %g" % var.PIVEC_MIN)
                    myVal = (1.5 + random.random()) * var.PIVEC_MIN
                    print("    Re-setting to %g" % myVal)
                    val[i] = myVal
                    needsNormalizing = 1

            if needsNormalizing:
                theSum = sum(val)
                for i in range(len(val)):
                    val[i] /= theSum
                if math.fabs(sum(val) - 1.0) > 5.e-16:
                    gm.append("Bad comp vals %s" % val)
                    gm.append("does not sum to 1.0")
                    raise P4Error(gm)
            # print "sum(val) - 1.0 = %f (%g)" % (sum(val) - 1.0, sum(val) -
            # 1.0)
            mt.val = val

        # Empirical protein comps are from the dat files in PAML.  Thanks,
        # Ziheng!
        elif spec == 'cpREV':
            mt.val = [0.0755, 0.0621, 0.0410, 0.0371, 0.0091,
                      0.0382, 0.0495, 0.0838, 0.0246, 0.0806,
                      0.1011, 0.0504, 0.0220, 0.0506, 0.0431,
                      0.0622, 0.0543, 0.0181, 0.0307, 0.0660]
            #theSum = sum(mt.val)
            # for i in range(len(mt.val)):
            #    mt.val[i] /= theSum
        elif spec == 'd78':
            # These first values have a couple more decimal places.  I
            # think I got these from some obscure code in a back corner of
            # the NCBI ftp site.  I believe it is obtainable by raising
            # the d78 matrix to a high power.  It is a more precise comp,
            # but is not used here because it is not standard.
            # mt.val = [0.08713, 0.04090, 0.04043, 0.04687, 0.03347,
            #          0.03826, 0.04953, 0.08861, 0.03362, 0.03689,
            #          0.08536, 0.08048, 0.01475, 0.03977, 0.05068,
            #          0.06958, 0.05854, 0.01049, 0.02992, 0.06472]

            # This next set of values is from her paper, and is the set
            # that everybody uses.
            # mt.val = [0.087, 0.041, 0.040, 0.047, 0.033,
            #          0.038, 0.05, 0.089, 0.034, 0.037,
            #          0.085, 0.08, 0.015, 0.04, 0.051,
            #          0.07, 0.058, 0.01, 0.03, 0.065]

            # These values are from Goldman's recommendations (Kosiol &
            # Goldman)
            mt.val = [0.087127, 0.040904, 0.040432, 0.046872, 0.033474,
                      0.038255, 0.049530, 0.088612, 0.033619, 0.036886,
                      0.085357, 0.080481, 0.014753, 0.039772, 0.050680,
                      0.069577, 0.058542, 0.010494, 0.029916, 0.064718]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'jtt':
            # mt.val = [0.077,0.051, 0.043, 0.052, 0.02,
            #          0.041, 0.062, 0.074, 0.023, 0.052,
            #          0.091, 0.059, 0.024, 0.04, 0.051,
            #          0.069, 0.059, 0.014, 0.032, 0.066]

            # again, a Goldman recommendation
            mt.val = [0.076862, 0.051057, 0.042546, 0.051269, 0.020279,
                      0.041061, 0.061820, 0.074714, 0.022983, 0.052569,
                      0.091111, 0.059498, 0.023414, 0.040530, 0.050532,
                      0.068225, 0.058518, 0.014336, 0.032303, 0.066374]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'mtREV24':
            mt.val = [0.072, 0.019, 0.039, 0.019, 0.006,
                      0.025, 0.024, 0.056, 0.028, 0.088,
                      0.168, 0.023, 0.054, 0.061, 0.054,
                      0.072, 0.086, 0.029, 0.033, 0.043]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'mtmam':
            mt.val = [0.0692, 0.0184, 0.0400, 0.0186, 0.0065,
                      0.0238, 0.0236, 0.0557, 0.0277, 0.0905,
                      0.1675, 0.0221, 0.0561, 0.0611, 0.0536,
                      0.0725, 0.0870, 0.0293, 0.0340, 0.0428]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'wag':
            mt.val = [0.0866279, 0.043972, 0.0390894, 0.0570451, 0.0193078,
                      0.0367281, 0.0580589, 0.0832518, 0.0244313, 0.048466,
                      0.086209,  0.0620286, 0.0195027, 0.0384319, 0.0457631,
                      0.0695179, 0.0610127, 0.0143859, 0.0352742, 0.0708956]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'rtRev':
            mt.val = [0.0646, 0.0453, 0.0376, 0.0422, 0.0114, 0.0606,
                      0.0607, 0.0639, 0.0273, 0.0679, 0.1018, 0.0751,
                      0.0150, 0.0287, 0.0681, 0.0488, 0.0622, 0.0251,
                      0.0318, 0.0619]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum

        elif spec == 'tmjtt94':
            mt.val = [0.105068479, 0.015695291, 0.018494452, 0.008897331,
                      0.021893432, 0.014095771, 0.009697091, 0.075777267,
                      0.016794962, 0.118764371, 0.163450965, 0.011196641,
                      0.033290013, 0.077676697, 0.025992202, 0.056782965,
                      0.052284315, 0.022293312, 0.032390283, 0.119464161]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'tmlg99':
            mt.val = [0.100632, 0.014017, 0.014706, 0.010371, 0.030668,
                      0.015152, 0.011343, 0.069235, 0.017501, 0.107722,
                      0.155161, 0.009723, 0.038730, 0.086453, 0.031761,
                      0.064333, 0.044847, 0.028277, 0.036988, 0.112380]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'lg':
            mt.val = [0.079066, 0.055941, 0.041977, 0.053052, 0.012937,
                      0.040767, 0.071586, 0.057337, 0.022355, 0.062157,
                      0.099081, 0.064600, 0.022951, 0.042302, 0.044040,
                      0.061197, 0.053287, 0.012066, 0.034155, 0.069147]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'blosum62':
            mt.val = [0.074, 0.052, 0.045, 0.054, 0.025,
                      0.034, 0.054, 0.074, 0.026, 0.068,
                      0.099, 0.058, 0.025, 0.047, 0.039,
                      0.057, 0.051, 0.013, 0.032, 0.073]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'hivb':
            mt.val = [0.060490222, 0.066039665, 0.044127815, 0.042109048, 0.020075899,
                      0.053606488, 0.071567447, 0.072308239, 0.022293943, 0.069730629,
                      0.098851122, 0.056968211, 0.019768318, 0.028809447, 0.046025282,
                      0.05060433, 0.053636813, 0.033011601, 0.028350243, 0.061625237]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'mtart':
            mt.val = [0.054116, 0.018227, 0.039903, 0.020160, 0.009709,
                      0.018781, 0.024289, 0.068183, 0.024518, 0.092638,
                      0.148658, 0.021718, 0.061453, 0.088668, 0.041826,
                      0.091030, 0.049194, 0.029786, 0.039443, 0.057700]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'mtzoa':
            mt.val = [0.068880,    0.021037,    0.030390,    0.020696,    0.009966,
                      0.018623,    0.024989,    0.071968,    0.026814,    0.085072,
                      0.156717,    0.019276,    0.050652,    0.081712,    0.044803,
                      0.080535,    0.056386,    0.027998,    0.037404,    0.066083]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'gcpREV':
            mt.val = [0.079510, 0.056001, 0.040459, 0.033220, 0.009051,
                      0.037505, 0.049675, 0.080233, 0.021880, 0.080496,
                      0.107512, 0.049324, 0.020776, 0.047731, 0.039916,
                      0.073820, 0.053615, 0.016705, 0.030790, 0.071781]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'stmtREV':
            mt.val = [0.046181, 0.053408, 0.036197, 0.023332, 0.023417,
                      0.039040, 0.034128, 0.038916, 0.016464, 0.089153,
                      0.161731, 0.055134, 0.023326, 0.091125, 0.034471,
                      0.077108, 0.041860, 0.020078, 0.030542, 0.064385]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'vt':  #from iqtree
            mt.val = [0.0770764620135024, 0.0500819370772208, 0.0462377395993731, 0.0537929860758246, 
                      0.0144533387583345, 0.0408923608974345, 0.0633579339160905, 0.0655672355884439, 
                      0.0218802687005936, 0.0591969699027449, 0.0976461276528445, 0.0592079410822730, 
                      0.0220695876653368, 0.0413508521834260, 0.0476871596856874, 0.0707295165111524, 
                      0.0567759161524817, 0.0127019797647213, 0.0323746050281867, 0.0669190817443274]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'pmb':  # from iqtree
            mt.val = [0.076, 0.054, 0.038, 0.045, 0.028, 
                      0.034, 0.053, 0.078, 0.030, 0.060, 
                      0.096, 0.052, 0.022, 0.045, 0.042, 
                      0.068, 0.056, 0.016, 0.036, 0.071]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum


        return mt

    def newRMatrix(self, partNum=0, free=0, spec='ones', val=None, symbol=None):
        """Make, attach, and return a new RMatrix instance.

        spec should be one of:

        -   'ones'         - for JC, poisson, F81
        -   '2p'           - for k2p and hky
        -   'specified'
        -   'cpREV'
        -   'd78'
        -   'jtt'
        -   'mtREV24'
        -   'mtmam'
        -   'wag'
        -   'rtRev'
        -   'tmjtt94'
        -   'tmlg99'
        -   'lg'
        -   'blosum62'
        -   'hivb'
        -   'mtart'
        -   'mtzoa'
        -   'gcpREV'
        -   'stmtREV'
        -   'vt'
        -   'pmb'

        See var.rMatrixProteinSpecs

        You do not set the 'val' arg unless the spec is 'specified' or
        '2p'.  If spec='2p', then you set val to kappa.

        If the spec is 'specified', you specify all the numerical values
        in a list given as the 'val' arg.  The length of that list will be
        (((dim * dim) - dim) / 2), so for DNA, where dim=4, you would
        specify a list containing 6 numbers.  """

        # not implemented:
        # 'blosum62a'
        # 'blosum62b'
        # 'phat70'

        complaintHead = '\nTree.newRMatrix()'
        gm = [complaintHead]
        self._checkModelThing(partNum, symbol, complaintHead)
        if self.model.cModel:
            self.deleteCStuff()
        mt = RMatrix()
        mt.partNum = partNum
        #mt.dim = self.data.parts[partNum].dim
        mt.free = free
        if spec not in var.rMatrixSpecs:
            gm.append("Got unknown rMatrix spec '%s'." % spec)
            gm.append("Should be one of: %s" % var.rMatrixSpecs)
            raise P4Error(gm)
        mt.spec = spec
        mt.num = len(self.model.parts[partNum].rMatrices)
        if symbol:
            mt.symbol = symbol
        else:
            mt.symbol = var.modelSymbols[mt.num]

        self.model.parts[partNum].rMatrices.append(mt)

        # assign val
        dim = self.model.parts[partNum].dim
        if var.rMatrixNormalizeTo1:
            goodLen = int((((dim * dim) - dim) / 2))
        else:
            goodLen = int((((dim * dim) - dim) / 2) - 1)

        v = None
        if spec == 'specified':
            if val:
                # should check that values are all floats
                if len(val) == goodLen:
                    v = numpy.array(val, numpy.float)
                    if var.rMatrixNormalizeTo1:
                        v /= v.sum()
                elif var.rMatrixNormalizeTo1 and len(val) == goodLen - 1:
                    gm.append("var.rMatrixNormalizeTo1 is set, val length should be %i, got %i" % (
                        goodLen, len(val)))
                    raise P4Error(gm)
                else:
                    gm.append("Bad length for arg val.  Length %i, should be %i" % (len(val), goodLen))
                    raise P4Error(gm)
            else:
                gm.append("spec is 'specified', but there are no specified rMatrix values.")
                gm.append("Specify rMatrix values by eg val=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]")
                raise P4Error(gm)
        elif spec == 'ones':
            v = numpy.array([1.0] * goodLen)
            if var.rMatrixNormalizeTo1:
                v /= v.sum()
        elif spec == '2p':
            try:
                v = float(val)
            except (ValueError, TypeError):
                gm.append("Kappa ('val' arg) should be a float.  Setting to 2.0")
                v = 2.0
            if v < var.KAPPA_MIN:
                gm.append("Kappa is too small.   Setting to %f" %
                          var.KAPPA_MIN)
                v = var.KAPPA_MIN
            elif v > var.KAPPA_MAX:
                gm.append("Kappa is too big.  Setting to %f" % var.KAPPA_MAX)
                v = var.KAPPA_MAX
            v = numpy.array([v], numpy.float)
        elif spec in var.rMatrixProteinSpecs:
            if self.data.parts[partNum].dataType != 'protein':
                gm.append("A protein matrix has been specified, but the dataType for part %i is %s." % (
                    partNum, self.data.parts[partNum].dataType))
                raise P4Error(gm)
            if free:
                gm.append('The rMatrix should not be free if it is an empirical protein matrix.')
                raise P4Error(gm)

        mt.val = v  # type numpy.ndarray, or None for specified protein
        return mt

    def newGdasrv(self, partNum=0, free=0, val=None, symbol=None):
        gm = ['Tree.newGdasrv()']

        if not self.model:
            gm.append("Set the data first.  Eg myTree.data = Data()")
            raise P4Error(gm)

        if self.model.cModel:
            self.deleteCStuff()

        # check if there is an nGammaCat > 1:
        if self.model.parts[partNum].nGammaCat == 1:
            gm.append("For this part (%s), the number of nGammaCat has been set to 1." % partNum)
            gm.append("So gdasrv won't work.")
            gm.append("You can set the nGammaCat with yourTree.setNGammaCat(partNum=x, nGammaCat=y)")
            raise P4Error(gm)

        # check val
        if val == None:
            gm.append("Please specify a val, a positive float.")
            raise P4Error(gm)
        try:
            v = float(val)
        except:
            gm.append("Arg val must be a float.  Got '%s'" % val)
            raise P4Error(gm)

        # This week, we have in defines.h
        # define GAMMA_SHAPE_MIN 0.000001
        # define GAMMA_SHAPE_MAX 300.0
        if v <= 0.000001 or v >= 300.0:
            gm.append("Arg val must be between 0.000001 and 300.  Got %f" % v)
            raise P4Error(gm)

        self._checkModelThing(partNum, symbol, gm[0])

        mt = Gdasrv()
        mt.nGammaCat = self.model.parts[partNum].nGammaCat
        mt.partNum = partNum
        mt.free = free
        # no spec or dim
        mt.num = len(self.model.parts[partNum].gdasrvs)
        if symbol:
            mt.symbol = symbol
        else:
            mt.symbol = var.modelSymbols[mt.num]

        self.model.parts[partNum].gdasrvs.append(mt)
        mt.freqs = numpy.zeros(mt.nGammaCat, numpy.float)
        mt.rates = numpy.zeros(mt.nGammaCat, numpy.float)
        mt._val[0] = v
        mt.calcRates()
        return mt

    def setPInvar(self, partNum=0, free=0, val=0.0):
        complaintHead = '\nTree.setPInvar()'
        gm = [complaintHead]

        # check val
        try:
            v = float(val)
        except:
            gm.append("Arg val must be a float.  Got '%s'" % val)
            raise P4Error(gm)

        if v < 0.0 or v >= 1.0:
            gm.append(
                "Arg val must be zero or more, and less than 1.  Got %f" % v)
            raise P4Error(gm)

        self._checkModelThing(partNum, None, complaintHead)
        if self.model.cModel:
            self.deleteCStuff()
        mt = PInvar()
        mt.partNum = partNum
        mt.free = free
        mt.val = v
        self.model.parts[partNum].pInvar = mt

    def setRelRate(self, partNum=0, val=0.0):
        complaintHead = '\nTree.setRelRate()'
        gm = [complaintHead]

        # check val
        try:
            v = float(val)
        except:
            gm.append("Arg val must be a float.  Got '%s'" % val)
            raise P4Error(gm)

        if v <= 0.0 or v >= 1000.0:
            gm.append(
                "Arg val must be more than zero, and less than 1000 (arbitrarily).  Got %f" % v)
            raise P4Error(gm)

        self._checkModelThing(partNum, None, complaintHead)
        if self.model.cModel:
            self.deleteCStuff()
        self.model.parts[partNum].relRate = v

    def setModelComponentOnNode(self, theModelComponent, node=None, clade=1):
        complaintHead = '\nTree.setModelComponentOnNode()'
        gm = [complaintHead]

        if theModelComponent and \
            (isinstance(theModelComponent, Comp) or
                isinstance(theModelComponent, RMatrix) or
                isinstance(theModelComponent, Gdasrv)):
            pass
        else:
            gm.append("Expecting a model component instance of some sort.")
            gm.append("Ie a comp, rMatrix, or gdasrv, instance.")
            gm.append("Got theModelComponent = %s" % theModelComponent)
            raise P4Error(gm)

        if self.model.cModel:
            self.deleteCStuff()

        partNum = theModelComponent.partNum

        if node == None:
            theNode = self.root
        else:
            theNode = self.node(node)

        isBad = 0
        if isinstance(theModelComponent, Comp):
            if theModelComponent != self.model.parts[partNum].comps[theModelComponent.num]:
                isBad = 1
        elif isinstance(theModelComponent, RMatrix):
            if theModelComponent != self.model.parts[partNum].rMatrices[theModelComponent.num]:
                isBad = 1
        elif isinstance(theModelComponent, Gdasrv):
            if theModelComponent != self.model.parts[partNum].gdasrvs[theModelComponent.num]:
                isBad = 1
        else:  # This will never happen-- we checked above.  Overkill.
            gm.append("I don't recognise theModelComponent.")
            raise P4Error(gm)
        if isBad:
            gm.append("The model component can only be set on the tree that made it.")
            raise P4Error(gm)

        # For the root, we set comps and nothing else.  For other nodes we
        # set anything.
        if theNode == self.root:
            if isinstance(theModelComponent, Comp):
                theNode.parts[partNum].compNum = theModelComponent.num
        else:
            if isinstance(theModelComponent, Comp):
                theNode.parts[partNum].compNum = theModelComponent.num
            elif isinstance(theModelComponent, RMatrix):
                theNode.br.parts[partNum].rMatrixNum = theModelComponent.num
            elif isinstance(theModelComponent, Gdasrv):
                theNode.br.parts[partNum].gdasrvNum = theModelComponent.num

        if clade:
            aboves = self.getNodeNumsAbove(theNode, leavesOnly=0)
            for i in aboves:
                if isinstance(theModelComponent, Comp):
                    self.nodes[i].parts[partNum].compNum = theModelComponent.num
                elif isinstance(theModelComponent, RMatrix):
                    self.nodes[i].br.parts[
                        partNum].rMatrixNum = theModelComponent.num
                elif isinstance(theModelComponent, Gdasrv):
                    self.nodes[i].br.parts[
                        partNum].gdasrvNum = theModelComponent.num

    def setModelThingsRandomly(self, forceRepresentation=2):
        """This method has been renamed to setModelComponentsOnNodesRandomly"""

        gm = ["Tree.setModelThingsRandomly() has been renamed Tree.setModelComponentsOnNodesRandomly()"]
        raise P4Error(gm)

    def setModelComponentsOnNodesRandomly(self, forceRepresentation=2):
        """Place model components (semi-)randomly on the tree.

        For example, if there are 2 compositions in model part partNum,
        this method will decorate each node of the tree with zeros and
        ones, randomly. The actual component set is
        node.parts[partNum].compNum.  If the model is homogeneous,
        it will just put zeros in all the nodes.

        We want to have each model component on the tree somewhere, and so it
        is not really randomly set.  If the model component numbers were
        assigned randomly on the tree, it may occur that some model component
        numbers by chance would not be represented.  This is not allowed,
        and you can set forceRepresentation to some positive integer, 1 or
        more.  That number will be the lower limit allowed on the number
        of nodes that get assigned the model component number.  For example,
        if forceRepresentation is set to 2, then each model component must get
        assigned to at least 2 nodes."""

        gm = ['Tree.setModelComponentsOnNodesRandomly()']

        if not self.model or not self.model.nParts:
            gm.append("No model parts?")
            raise P4Error(gm)

        if self.model.cModel:
            self.deleteCStuff()
        # self.model.dump()

        if not isinstance(forceRepresentation, int) or forceRepresentation < 1:
            gm.append("Arg 'forceRepresentation' should be 1 or more.")
            gm.append("Got forceRepresentation = %s" % forceRepresentation)
            raise P4Error(gm)

        for i in self.preOrder:
            if i == var.NO_ORDER:
                gm.append("This method does not work if any nodes are not used in the tree.")
                raise P4Error(gm)

        for pNum in range(self.model.nParts):
            mp = self.model.parts[pNum]

            if mp.ndch2:
                continue

            # First do comps
            if mp.nComps == 1:
                for n in self.nodes:
                    n.parts[pNum].compNum = 0
            elif mp.nComps > 1:
                nNodes = len(self.nodes)
                if (mp.nComps * forceRepresentation) > nNodes:
                    gm.append("Part %i" % pNum)
                    gm.append(
                        "There are not enough nodes (%i) to put %i" % (nNodes, mp.nComps))
                    gm.append(
                        "comps on at least forceRepresentation (%i) nodes." % forceRepresentation)
                    raise P4Error(gm)
                nList = self.nodes[:]
                random.shuffle(nList)
                # get the forceRepresentation out of the way first
                for mtNum in range(mp.nComps):
                    for fr in range(forceRepresentation):
                        n = nList.pop()
                        n.parts[pNum].compNum = mtNum
                # Now do the rest
                for n in nList:
                    n.parts[pNum].compNum = random.randrange(mp.nComps)
            else:
                gm.append("No comps in part %i" % pNum)
                raise P4Error(gm)

            # Second do rMatrices
            if mp.nRMatrices == 1:
                for n in self.nodes:
                    if n != self.root:
                        n.br.parts[pNum].rMatrixNum = 0
            elif mp.nRMatrices > 1:
                nNodes = len(self.nodes) - 1
                if (mp.nRMatrices * forceRepresentation) > nNodes:
                    gm.append("Part %i" % pNum)
                    gm.append(
                        "There are not enough nodes (%i) to put %i" % (nNodes, mp.nRMatrices))
                    gm.append(
                        "rMatrices on at least forceRepresentation (%i) nodes." % forceRepresentation)
                    raise P4Error(gm)
                nList = self.nodes[:]
                nList.remove(self.root)
                random.shuffle(nList)
                # get the forceRepresentation out of the way first
                for mtNum in range(mp.nRMatrices):
                    for fr in range(forceRepresentation):
                        n = nList.pop()
                        n.br.parts[pNum].rMatrixNum = mtNum
                # Now do the rest
                for n in nList:
                    n.br.parts[pNum].rMatrixNum = random.randrange(
                        mp.nRMatrices)

            else:
                gm.append("No rMatrices in part %i" % pNum)
                raise P4Error(gm)

            # Third do gdasrvs
            if mp.nGammaCat > 1:
                if mp.nGdasrvs == 1:
                    for n in self.nodes:
                        if n != self.root:
                            n.br.parts[pNum].gdasrvNum = 0
                elif mp.nGdasrvs > 1:
                    nNodes = len(self.nodes) - 1
                    if (mp.nGdasrvs * forceRepresentation) > nNodes:
                        gm.append("Part %i" % pNum)
                        gm.append(
                            "There are not enough nodes (%i) to put %i" % (nNodes, mp.nGdasrvs))
                        gm.append(
                            "gdasrvs on at least forceRepresentation (%i) nodes." % forceRepresentation)
                        raise P4Error(gm)
                    nList = self.nodes[:]
                    nList.remove(self.root)
                    random.shuffle(nList)
                    # get the forceRepresentation out of the way first
                    for mtNum in range(mp.nGdasrvs):
                        for fr in range(forceRepresentation):
                            n = nList.pop()
                            n.br.parts[pNum].gdasrvNum = mtNum
                    # Now do the rest
                    for n in nList:
                        n.br.parts[pNum].gdasrvNum = random.randrange(
                            mp.nGdasrvs)
                else:
                    gm.append(
                        "No gdasrvs in part %i and yet nGammaCat > 1" % pNum)
                    raise P4Error(gm)

        # self.dump(model=True)

    def setModelComponentsNNodes(self):
        """Set nNodes for all model components"""

        gm = ['Tree.setModelComponentsNNodes()']

        if not self.model or not self.model.nParts:
            gm.append("No model parts?")
            raise P4Error(gm)

        for pNum in range(self.model.nParts):
            mp = self.model.parts[pNum]

            if not mp.nComps:
                gm.append("No comps in model part %i." % pNum)
                raise P4Error(gm)
            elif not mp.nRMatrices:
                gm.append("No rMatrices in model part %i." % pNum)
                raise P4Error(gm)

        for pNum in range(self.model.nParts):
            mp = self.model.parts[pNum]

            # First do comps
            if mp.nComps == 1:
                pass
            elif mp.nComps > 1:
                for mtNum in range(mp.nComps):
                    mp.comps[mtNum].nNodes = 0
                for n in self.iterNodes():
                    mp.comps[n.parts[pNum].compNum].nNodes += 1

            # Second do rMatrices
            if mp.nRMatrices == 1:
                pass
            elif mp.nRMatrices > 1:
                for mtNum in range(mp.nRMatrices):
                    mp.rMatrices[mtNum].nNodes = 0
                for n in self.iterNodesNoRoot():
                    mp.rMatrices[n.br.parts[pNum].rMatrixNum].nNodes += 1

            # Third do gdasrvs
            if mp.nGammaCat > 1:
                if mp.nGdasrvs == 1:
                    pass
                elif mp.nGdasrvs > 1:
                    for mtNum in range(mp.nGdasrvs):
                        mp.gdasrvs[mtNum].nNodes = 0
                    for n in self.iterNodesNoRoot():
                        mp.gdasrvs[n.br.parts[pNum].gdasrvNum].nNodes += 1
                else:
                    gm.append("No gdasrvs in part %i" % pNum)
                    raise P4Error(gm)

    def summarizeModelComponentsNNodes(self):
        """Summarize nNodes for all model components if isHet"""

        gm = ['Tree.summarizeModelComponentsNNodes()']

        if not self.model or not self.model.nParts:
            gm.append("No model parts?")
            raise P4Error(gm)
        if not self.model.isHet:
            gm.append("This method is for hetero models")
            raise P4Error(gm)

        for pNum in range(self.model.nParts):
            mp = self.model.parts[pNum]

            if not mp.nComps:
                gm.append("No comps in model part %i." % pNum)
                raise P4Error(gm)
            elif not mp.nRMatrices:
                gm.append("No rMatrices in model part %i." % pNum)
                raise P4Error(gm)

        for pNum in range(self.model.nParts):
            print("\n%6s %s:" % ("Part", pNum))
            mp = self.model.parts[pNum]

            # First do comps
            if mp.nComps == 1:
                pass
            elif mp.ndch2:
                print("%16s" % "ndch2 is on")
            elif mp.nComps > 1:
                for mtNum in range(mp.nComps):
                    # print "  comp %i nNodes=%i" % (mtNum,
                    # mp.comps[mtNum].nNodes)
                    print("%16s %i %s = %i" % ("composition", mtNum, "nNodes",
                                               mp.comps[mtNum].nNodes))

            # Second do rMatrices
            if mp.nRMatrices == 1:
                pass
            elif mp.ndrh2:
                print("%16s" % "ndrh2 is on")
            elif mp.nRMatrices > 1:
                for mtNum in range(mp.nRMatrices):
                    print("%16s %i %s = %i" % ("rate matrix", mtNum,
                                               "nNodes", mp.rMatrices[mtNum].nNodes))

            # Third do gdasrvs
            if mp.nGammaCat > 1:
                if mp.nGdasrvs == 1:
                    pass
                elif mp.nGdasrvs > 1:
                    for mtNum in range(mp.nGdasrvs):
                        print("  gdasrv %i nNodes =%i" % (mtNum, mp.gdasrvs[mtNum].nNodes))
                else:
                    gm.append("No gdasrvs in part %i" % pNum)
                    raise P4Error(gm)

    def setTextDrawSymbol(self, theSymbol='-', node=None, clade=1):
        gm = ['\nTree.setTextDrawString()']

        if not theSymbol or not isinstance(theSymbol, str)  or len(theSymbol) != 1:
            gm.append("theSymbol should be a single character string.")
            raise P4Error(gm)

        if not node:
            theNode = self.root
        else:
            theNode = self.node(node)

        if theNode == self.root:
            pass
        else:
            theNode.br.textDrawSymbol = theSymbol

        if clade:
            aboves = self.getNodeNumsAbove(theNode, leavesOnly=0)
            for i in aboves:
                self.nodes[i].br.textDrawSymbol = theSymbol

    def setNGammaCat(self, partNum=0, nGammaCat=1):
        gm = ['\nTree.setNGammaCat()']
        if not self.data or not self.model:
            gm.append("No data?")
            raise P4Error(gm)
        if self.model.cModel:
            self.deleteCStuff()
        if partNum < 0 or partNum >= self.model.nParts:
            gm.append("PartNum %s is out of range of %s parts." %
                      (partNum, self.model.nParts))
            raise P4Error(gm)

        try:
            x = int(nGammaCat)
        except ValueError:
            gm.append("'%s' does not appear to be an integer." % i)
            raise P4Error(gm)
        if x < 1:
            gm.append("nGammaCat should not be less than 1.")
            raise P4Error(gm)
        elif x > 16:
            gm.append("nGammaCat '%s' exceeds the arbitrary limit of 16." % x)
            raise P4Error(gm)
        self.model.parts[partNum].nGammaCat = nGammaCat


    def modelSanityCheck(self, resetEmpiricalComps=True):
        """Check that the tree, data, and model specs are good to go.

        Complain and exit if there is anything wrong that might prevent a
        likelihood evaluation from being done.  We are assuming that a
        data object exists and is attached, and that model stuff has been
        set.

        Check that each part has at least 1 each from comps, rMatrices,
        and gdasrvs (if nGammaCat is > 1).

        Check that each node has a comp, rMatrix, and gdasr.  Check that all
        comps, rMatrices, gdasrvs are used on a node somewhere.

        Here relRate, ie the relative rate of each data partition, is
        adjusted based on the size of the data partitions.

        newRelRate_p = oldRelRate_p * (Sum_p[oldRelRate_i * partLen_i] / Sum[partLen_i])

        That ensures that Sum(relRate_i * partLen_i) = totalDataLength, ie
        that the weighted mean of the rates is 1.0.

        This method also tallies up the number of free prams in the whole
        model, and sets self.model.nFreePrams.

        """

        complaintHead = '\nTree.modelSanityCheck()'
        gm = [complaintHead]
        # print("\nTree.modelSanityCheck() here. self.model.nParts=%s" % self.model.nParts)
        # print("\nTree.modelSanityCheck() here.  resetEmpiricalComps=%s" % resetEmpiricalComps)
        isBad = 0
        complaints = []
        if not self.data:
            complaints.append('    No data.')
            isBad = 1
        if not self.model:
            complaints.append('    No model.')
            isBad = 1

        # Set isHet.
        for pNum in range(self.model.nParts):
            mp = self.model.parts[pNum]
            mp.isHet = 0
            if mp.nComps > 1 or mp.nRMatrices > 1:
                mp.isHet = 1
            if mp.nGammaCat > 1 and mp.nGdasrvs > 1:
                mp.isHet = 1

        # This week ndch2 does not play well with other hetero models like ndch,
        # so insist that all partitions are or are not ndch2.
        firstPartIsNdch2 = self.model.parts[0].ndch2
        for pNum in range(self.model.nParts):
            mp = self.model.parts[pNum]
            if mp.ndch2 != firstPartIsNdch2:
                complaints.append("    Can't mix ndch2 with non-ndch2 models.")
                isBad = 1

        # Check that all parts have all the required stuff.  Make a list
        # of errors.  If there is something missing or wrong, don't die
        # right away, but add the problem to the list, and write it all
        # out at the end.  It gives the user a chance to fix more than one
        # error at a time.
        for pNum in range(self.model.nParts):
            complaints.append('  Part %i' % pNum)
            partIsBad = 0
            mp = self.model.parts[pNum]

            # Check if essential things have been set
            if not mp.nComps:
                complaints.append('    No comps in part %s' % pNum)
                partIsBad = 1
            if not mp.nRMatrices:
                complaints.append('    No rMatrices in part %s' % pNum)
                partIsBad = 1
            if mp.nGammaCat > 1:
                if not mp.nGdasrvs:
                    complaints.append('    No gdasrvs in part %s' % pNum)
                    partIsBad = 1
            if mp.nGammaCat == 1:
                if mp.nGdasrvs:
                    complaints.append(
                        '    There should be no gdasrvs in part %s, with nGammaCat=1' % pNum)
                    partIsBad = 1
            if not mp.pInvar:
                complaints.append('    No pInvar in part %s' % pNum)
                partIsBad = 1

            if mp.ndch2:
                if mp.nComps != len(list(self.iterNodes())):
                    complaints.append('Part %i, ndch2 needs a comp for each node' % pNum)
                    partIsBad = 1

            if partIsBad:
                gm.append("  (Indices are zero-based.)")
                gm += complaints
                raise P4Error(gm)

            # Check if comp values have been set.
            for mt in mp.comps:
                if mt.spec != 'empirical' or not resetEmpiricalComps:
                    if mt.val is None:
                        complaints.append(
                            '    No composition val in part %s' % pNum)
                        partIsBad = 1
                    if len(mt.val) != mp.dim:
                        complaints.append('    Composition val is wrong length (%i), but dim is %i' % (
                            len(mt.val), mp.dim))
                        partIsBad = 1

            # We don't want multiple rMatrices or free rMatrices if mp.dim is 2
            if mp.dim == 2:
                if mp.nRMatrices > 1:
                    complaints.append(
                        '    Part %s is dim 2, but we have more than one rMatrix' % pNum)
                    partIsBad = 1
                mt = mp.rMatrices[0]  # hopefully only one
                if mt.free:
                    complaints.append(
                        '    Part %s is dim 2, but rMatrix 0 is free' % pNum)
                    partIsBad = 1

            mp.nCat = mp.nGammaCat
            # If the model part isHet, we need to check that all nodes
            # have something assigned, and that all model components are
            # used.  If the model part is not het, we can skip that,
            # but we need to check that all the
            # node.parts[pNum].compNum are 0, and all the
            # node.br.parts[pNum].rMatrixNum and
            # node.br.parts[pNum].gdasrvNum are set to 0.
            if not mp.isHet:
                # print "model part %i is not het" % pNum
                for n in self.iterNodes():
                    # print("pNum = %i, n.nodeNum=%i, len n.parts = %i" % (pNum, n.nodeNum, len(n.parts)))
                    n.parts[pNum].compNum = 0
                    if n != self.root:
                        n.br.parts[pNum].rMatrixNum = 0
                        if mp.nGammaCat > 1:
                            n.br.parts[pNum].gdasrvNum = 0
            else:  # isHet

                # If there is only one comp, rMatrix, or gdasrv, then
                # simply set it.
                if mp.nComps == 1:
                    for n in self.iterNodes():
                        n.parts[pNum].compNum = 0
                if mp.nRMatrices == 1:
                    for n in self.iterNodes():
                        if n != self.root:
                            n.br.parts[pNum].rMatrixNum = 0
                if mp.nGammaCat > 1 and mp.nGdasrvs == 1:
                    for n in self.iterNodes():
                        if n != self.root:
                            n.br.parts[pNum].gdasrvNum = 0

                # print "model part %i is het" % pNum
                # New ad hoc attribute 'isUsed', to keep track of whether
                # any node uses it.
                for mt in mp.comps:
                    mt.isUsed = 0
                for mt in mp.rMatrices:
                    mt.isUsed = 0
                for mt in mp.gdasrvs:
                    mt.isUsed = 0

                # Does every node have all required things?
                for n in self.iterNodes():
                    mtNum = n.parts[pNum].compNum
                    if mtNum >= 0 and mtNum < mp.nComps:
                        mt = mp.comps[mtNum]
                        mt.isUsed = 1
                    else:
                        complaints.append('    Part %s, node %s has no comp.' % (pNum, n.nodeNum))
                        partIsBad = 1

                    if n != self.root:
                        mtNum = n.br.parts[pNum].rMatrixNum
                        if mtNum >= 0 and mtNum < mp.nRMatrices:
                            mt = mp.rMatrices[n.br.parts[pNum].rMatrixNum]
                            mt.isUsed = 1
                        else:
                            complaints.append('    Part %s, node %s has no rMatrix.' % (pNum, n.nodeNum))
                            partIsBad = 1
                        if mp.nGammaCat > 1:
                            mtNum = n.br.parts[pNum].gdasrvNum
                            if mtNum >= 0 and mtNum < mp.nGdasrvs:
                                mt = mp.gdasrvs[n.br.parts[pNum].gdasrvNum]
                                mt.isUsed = 1
                            else:
                                complaints.append('    Part %s, node %s has no gdasrv. nGammaCat=%s' % (
                                    pNum, n.nodeNum, mp.nGammaCat))
                                partIsBad = 1
                        if mp.nGammaCat == 1:
                            if n.br.parts[pNum].gdasrvNum != -1:
                                complaints.append('    Part %s, node %s has a gdasrv, but nGammaCat is 1.' % (
                                    pNum, n.nodeNum))
                                partIsBad = 1

                # Is every model component used?
                for mt in mp.comps:
                    if not mt.isUsed:
                        complaints.append('    Part %s, comp %s is not used.' % (pNum, mt.num))
                        partIsBad = 1
                for mt in mp.rMatrices:
                    if not mt.isUsed:
                        complaints.append(
                            '    Part %s, rMatrix %s is not used.' % (pNum, mt.num))
                        partIsBad = 1
                for mt in mp.gdasrvs:
                    if not mt.isUsed:
                        complaints.append(
                            '    Part %s, gdasrv %s is not used.' % (pNum, mt.num))
                        partIsBad = 1

                # Clean up ad hoc attr 'isUsed'
                for mt in mp.comps:
                    del(mt.isUsed)
                for mt in mp.rMatrices:
                    del(mt.isUsed)
                for mt in mp.gdasrvs:
                    del(mt.isUsed)

            if partIsBad:
                isBad = 1
            else:
                complaints.append('    ok')

        # ##################################
        if resetEmpiricalComps:
            self.setEmpiricalComps()

        # self.model.isHet if any part isHet
        self.model.isHet = 0
        for pNum in range(self.model.nParts):
            if self.model.parts[pNum].isHet:
                self.model.isHet = 1
                break

        # relativeRates
        self.model.doRelRates = 0
        if self.model.nParts > 1:
            for p in self.model.parts:
                if p.relRate != 1.0:  # This week, the default relRate is 1.0
                    self.model.doRelRates = 1
                    break
        if self.model.relRatesAreFree:
            self.model.doRelRates = 1

        if self.model.doRelRates:
            totDataLen = 0
            for p in self.data.parts:
                totDataLen += p.nChar
            fact = 0.0
            for i in range(self.model.nParts):
                fact += (self.model.parts[i].relRate *
                         self.data.parts[i].nChar)
            fact = float(totDataLen) / fact
            for p in self.model.parts:
                p.relRate *= fact
            if 0:
                print("RelativeRates (adjusted for length)")
                for i in range(self.model.nParts):
                    p = self.model.parts[i]
                    print("  part %s,  nChar %5s, relRate %s" % (p.num, self.data.parts[i].nChar, p.relRate))
            if 1:
                total = 0.0
                for i in range(self.model.nParts):
                    total += (self.model.parts[i].relRate *
                              (float(self.data.parts[i].nChar) / float(totDataLen)))
                if abs(total - 1.0) > 1.0e-12:
                    gm.append(
                        'Error in relativeRate calculation (total=%s).' % total)
                    raise P4Error(gm)

        # print "modelSanityCheck. relRatesAreFree=%s, doRelRates=%s" %
        # (self.model.relRatesAreFree, self.model.doRelRates)

        # model.nFreePrams
        self.model.nFreePrams = 0
        for mp in self.model.parts:
            for mt in mp.comps:
                if mt.free:
                    self.model.nFreePrams += mp.dim - 1
            for mt in mp.rMatrices:
                if mt.free:
                    if mt.spec == '2p':
                        self.model.nFreePrams += 1
                    else:
                        self.model.nFreePrams += (
                            ((mp.dim * mp.dim) - mp.dim) / 2) - 1
            for mt in mp.gdasrvs:
                if mt.free:
                    self.model.nFreePrams += 1
            if mp.pInvar.free:
                self.model.nFreePrams += 1
        # print "Tree.modelSanityCheck().  Counted %i free params." %
        # self.model.nFreePrams

        if self.model.doRelRates and self.model.relRatesAreFree:
            self.model.nFreePrams += self.model.nParts - 1

        if isBad:
            gm.append("(Indices are zero-based.)")
            gm += complaints
            raise P4Error(gm)

    def setEmpiricalComps(self):
        """Set any empirical model comps to the comp of the data.

        This is done by self.modelSanityCheck(), but sometimes you may
        want to do it at other times.  For example, do this after
        exchanging Data objects, or after simulating.  In those cases
        there does not seem to be a reasonable way to do it
        automatically."""
        complaintHead = '\nTree.setEmpiricalComps()'
        gm = [complaintHead]
        if not self.model:
            gm.append("This tree has no model.")
            raise P4Error(gm)
        if not self.data:
            gm.append("This tree has no data.")
            raise P4Error(gm)

        for mp in self.model.parts:
            for c in mp.comps:
                if c.spec == 'empirical':
                    # print "got empirical comp, comp %s in part %s. (nComps=%i, isHet=%s)" % (
                    #    c.num, mp.num, mp.nComps, mp.isHet)
                    if not mp.isHet:
                        seqNums = None
                    elif mp.nComps == 1:
                        seqNums = None
                    else:
                        seqNums = []
                        # for n in self.nodes:
                        #    print "node %2i seqNum=%3i n.parts[%i].compNum=%3i" % (
                        # n.nodeNum, n.seqNum, mp.num, n.parts[mp.num].compNum)

                        for n in self.iterNodes():
                            # Is the comp used by the node?
                            if n.parts[mp.num].compNum == c.num:
                                # print "comp %s is used by node %s" % (c.num,
                                # n.nodeNum)
                                if n.isLeaf:
                                    nodeNums = [n.nodeNum]
                                else:
                                    nodeNums = self.getNodeNumsAbove(
                                        n, leavesOnly=1)
                                # gm.append("nodeNums for %s = %s" %
                                # (n.nodeNum, nodeNums)
                                for i in nodeNums:
                                    seqNum = self.nodes[i].seqNum
                                    if seqNum not in seqNums:
                                        seqNums.append(seqNum)

                        # print "setEmpiricalComps() got seqNums = %s" %
                        # seqNums
                        if not seqNums:
                            gm.append(
                                "Something is wrong here.  part %i, comp %i." % (mp.num, c.num))
                            gm.append(
                                "This comp object has no sequences from which to get the empirical comp.")
                            gm.append(
                                "Maybe you need to yourTree.setModelModelComponentOnNode() or ")
                            gm.append("yourTree.setModelComponentsOnNodesRandomly()")
                            #gm.append(
                            #    "Or maybe its an extra comp in an RJ MCMC? -- If so, fix")
                            #gm.append("the comp val to eg 'equal'.")
                            raise P4Error(gm)

                    # dim long, not dim - 1
                    c.val = self.data.parts[mp.num].composition(seqNums)  
                    # print "  seqNums=%s, c.val=%s" % (seqNums, c.val)

                    needsNormalizing = 0
                    for i in range(len(c.val)):
                        if c.val[i] < var.PIVEC_MIN:
                            c.val[i] = var.PIVEC_MIN + (0.2 * var.PIVEC_MIN) + (var.PIVEC_MIN * random.random())
                            needsNormalizing = 1
                    theSum = numpy.sum(c.val)
                    # print "setEmpiricalComps().  Got theSum = %i" % theSum

                    # We may have asked for the comp of an empty sequence,
                    # in which case val is all zeros.  Check for that.
                    if math.fabs(1.0 - theSum) > 0.1:
                        gm.append(
                            "Something is very wrong here.  Empirical comp vals should sum to 1.0")
                        gm.append("The sum of the comp vals for part %s, comp %s, is %s" % (
                            mp.num, c.num, theSum))
                        gm.append(
                            "Probably the sequences from which the composition was taken were blank.")
                        raise P4Error(gm)

                    if needsNormalizing or abs(theSum - 1.0) > 1e-16:
                        c.val /= theSum

    # def __del__(self, freeTree=pf.p4_freeTree, freeNode=pf.p4_freeNode):
    # def __del__(self, freeTree=pf.p4_freeTree, dp_freeTree = pf.dp_freeTree, mysys=sys):
    # def __del__(self, freeTree=pf.p4_freeTree, dp_freeTree = pf.dp_freeTree):
class RMatrix(object):

    def __init__(self):
        self.num = -1
        self.partNum = None
        self.free = None
        self.spec = None
        self.symbol = None
        self.val = None
        self.nNodes = 0


class Gdasrv(object):

    def __init__(self):
        self.num = -1
        self.partNum = None
        self.free = None
        self.symbol = None
        # self.val=None
        self._val = numpy.zeros(1, numpy.float)
        self.freqs = None
        self.rates = None
        self.nGammaCat = None
        self.c = None    # a p4_gdasrvStruct, if it exists.  It is used in calcRates, below.
        self.nNodes = 0

    def _getVal(self):
        return self._val

    def _setVal(self, theVal):
        if theVal < 1.e-16:
            gm = ["Gdasrv._setVal()"]
            gm.append("Attempt to set Gdasrv.val (ie alpha) to %g" % theVal)
            gm.append(
                "However, we cannot calculate the discrete categories with a value so low.")
            raise P4Error(gm)
        self._val[0] = theVal
        self.calcRates()
    
    def _delVal(self):
        gm = ["Don't/Can't delete this Gdasrv property."]
        raise P4Error(gm)


    val = property(_getVal, _setVal, _delVal)

    def calcRates(self):
        # Use either the p4_gdasrvStruct, or just use the NumPy
        # array vals (np = NumPy).
        # print("self.c = %s" % self.c)
        if self.c:
            pf.gdasrvCalcRates(self.c)
        else:
            pf.gdasrvCalcRates_np(self.nGammaCat, self._val[0], self.freqs, self.rates)
        # print('xxx self.rates = %s, val=%s' % (self.rates, self._val[0]))


class Comp(object):

    def __init__(self):
        self.num = -1
        self.partNum = None
        self.free = None
        self.spec = None
        self.symbol = None
        self._val = None
        self.nNodes = 0

    def _getVal(self):
        return self._val

    def _setVal(self, theVal):
        if self._val is None:
            self._val = numpy.array(theVal)
        else:
            #print("Resetting comp val.")
            assert len(self._val) == len(theVal)
            for i in range(len(theVal)):
                self._val[i] = theVal[i]

    def _delVal(self):
        gm = ["Don't/Can't delete this Comp property."]
        raise P4Error(gm)


    val = property(_getVal, _setVal, _delVal)


class PInvar(object):

    def __init__(self):
        self.num = -1
        self.partNum = None
        self.free = None
        self.val = None

