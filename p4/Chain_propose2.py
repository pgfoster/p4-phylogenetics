"""More Chain 'propose' methods."""

import func,pf
from Var import var
from Glitch import Glitch
import math,random
import numpy


def proposeCompWithSlider(self, theProposal):
    gm = ['Chain.proposeCompWithSlider()']

    mt = self.propTree.model.parts[theProposal.pNum].comps[theProposal.mtNum]
    dim = self.propTree.model.parts[theProposal.pNum].dim

    # mt.val is a list, not a numpy array
    #assert type(mt.val) == numpy.ndarray

    indxs = random.sample(range(dim), 2)
    currentAplusB = mt.val[indxs[0]] + mt.val[indxs[1]]
    thisMin = var.PIVEC_MIN / currentAplusB
    thisMax = 1. - thisMin

    minToMaxDiff = thisMax - thisMin
    thisTuning = theProposal.tuning
    
    # It is possible if A and B char states are missing that both A
    # and B values are very close to var.PIVEC_MIN, in which case
    # thisMin and thisMax will both be close to 0.5, and so the tuning
    # will be too much, requiring too many reflections.  In that case,
    # just change the tuning temporarily.
    if thisTuning > minToMaxDiff:
        thisTuning = minToMaxDiff
        #print "temporarily changing the tuning for comp proposal, to", thisTuning
    
    x = mt.val[indxs[0]] / currentAplusB
    y = x + (thisTuning * (random.random() - 0.5))

    # reflection
    safety = -1
    while 1:
        safety += 1
        if safety > 100:
            gm.append("Did more than 100 reflections -- something is wrong.")
            raise Glitch, gm
        if y < thisMin:
            y = thisMin + (thisMin - y)
        elif y > thisMax:
            y = thisMax - (y - thisMax)
        else:
            break
    #if safety > 1:
    #    print "comp reflections: ", safety
    mt.val[indxs[0]] = y * currentAplusB
    mt.val[indxs[1]] = currentAplusB - mt.val[indxs[0]]

    # The following normalization to a sum of 1 may not be needed.
    mySum = 0.0
    for stNum in range(dim):
        mySum += mt.val[stNum]
    for stNum in range(dim):
        mt.val[stNum] /= mySum

    self.logProposalRatio = 0.0
    # The prior here is a flat Dirichlet, ie Dirichlet(1, 1, 1, ...,
    # 1).  If it is informative, then the prior is affected.
    self.logPriorRatio = 0.0

def proposeCompWithDirichlet(self, theProposal):
    gm = ['Chain.proposeCompWithDirichlet()']

    mt = self.propTree.model.parts[theProposal.pNum].comps[theProposal.mtNum]
    dim = self.propTree.model.parts[theProposal.pNum].dim

    # mt.val is a list of floats, not a numpy.ndarray
    #print type(mt.val), type(mt.val[0])

    # The tuning is the Dirichlet alpha.
    #print theProposal.tuning

    # This method uses func.dirichlet1, which is for lists not numpy
    # arrays.  A copy of inSeq is made, and the copy is modified and
    # returned.
    #dirichlet1(inSeq, alpha, theMin, theMax)
    newVal = func.dirichlet1(mt.val, theProposal.tuning, var.PIVEC_MIN, 1 - var.PIVEC_MIN)
    
    self.logProposalRatio = 0.0

    rangeDim = range(dim)
    mySum = 0.0
    for stNum in rangeDim:
        mySum += newVal[stNum] * theProposal.tuning
    x = pf.gsl_sf_lngamma(mySum)
    for stNum in rangeDim:
        x -= pf.gsl_sf_lngamma(newVal[stNum] * theProposal.tuning)
    for stNum in rangeDim:
        x += ((newVal[stNum] * theProposal.tuning) - 1.) * math.log(mt.val[stNum])

    mySum = 0.0
    for stNum in rangeDim:
        mySum += mt.val[stNum] * theProposal.tuning
    y = pf.gsl_sf_lngamma(mySum)
    for stNum in rangeDim:
        y -= pf.gsl_sf_lngamma(mt.val[stNum] * theProposal.tuning)
    for stNum in rangeDim:
        y += ((mt.val[stNum] * theProposal.tuning) - 1.) * math.log(newVal[stNum])
    self.logProposalRatio = x - y
    mt.val = newVal
                
    # The prior here is a flat Dirichlet, ie Dirichlet(1, 1, 1, ...,
    # 1).  If it is informative, then the prior is affected.
    self.logPriorRatio = 0.0

def proposeRjComp(self, theProposal):
    gm = ['Chain.proposeRjComp()']
    theProposal.doAbort = False
    mp = self.propTree.model.parts[theProposal.pNum]
    assert mp.rjComp
    assert mp.rjComp_k >= 1
    assert mp.rjComp_k <= mp.nComps

    # If the pool size k, mp.rjComp_k, is only 1, then we can only split
    if mp.rjComp_k == 1:
        self.proposeSplitComp(theProposal)

    # If k, mp.rjComp_k is mp.nComps, then we can only merge
    if mp.rjComp_k == mp.nComps:
        self.proposeMergeComp(theProposal)

    # Otherwise, choose randomly.
    if random.random() < 0.5:
        self.proposeSplitComp(theProposal)
    else:
        self.proposeMergeComp(theProposal)

def proposeSplitComp(self, theProposal):
    gm = ['Chain.proposeSplitComp()']

    # var.rjCompUniformAllocationPrior  True by default
    # theProposal.tuning 200.  becomes p0 below

    mp = self.propTree.model.parts[theProposal.pNum]
    assert mp.rjComp
    dim = mp.dim
    
    # Check that k is less than k_max, which is nComps.  This should have been checked before, but check again.
    #print gm[0], "rjComp_k is currently %i, with %i comps" % (mp.rjComp_k, mp.nComps)
    assert mp.rjComp_k < mp.nComps

    # Select an existing comp vector from the pool
    pool = [c for c in mp.comps if c.rj_isInPool]
    notInPool = [c for c in mp.comps if not c.rj_isInPool]
    assert notInPool # or else we can't split
    pi0  = random.choice(pool)

    # The nodes currently associated with pi0
    beta0 = [n for n in self.propTree.iterNodes() if n.parts[theProposal.pNum].compNum == pi0.num]
    b0 = float(len(beta0))
    #print gm[0], "comp %i is chosen, f=%f, currently on nodes" % (pi0.num, pi0.rj_f), [n.nodeNum for n in beta0]

    # Divvy up the contents of beta0 into (new) beta1 and beta2, based on probability u
    if var.rjCompUniformAllocationPrior:
        u = 0.5
    else:
        u = random.random()

    beta1 = []
    beta2 = []
    for it in beta0:
        r = random.random()
        if r < u:
            beta1.append(it)
        else:
            beta2.append(it)
    b1 = float(len(beta1))
    b2 = float(len(beta2))
    bPrime0 = b0 + 2.
    bPrime1 = b1 + 1.
    bPrime2 = b2 + 1.

    # Calculation of f1 and f2 depends on u
    f0 = pi0.rj_f
    f1 = u * f0
    f2 = (1.0 - u) * f0

    
    uu = [random.normalvariate(0., 1.) for i in range(dim)]
    p0 = theProposal.tuning
    s0 = random.gammavariate(p0, 1.)
    m0 = [s0 * it for it in pi0.val]
    #print m0

    # I get a math range error here -- needs debugging.
    #m1 = [m0[j] * math.exp((bPrime0 * uu[j])/(bPrime1 * math.sqrt(m0[j])))
    #      for j in range(dim)]
    #m2 = [m0[j] * math.exp((-bPrime0 * uu[j])/(bPrime2 * math.sqrt(m0[j])))
    #      for j in range(dim)]

    safety = 0
    while 1:
        try:
            m1 = [m0[j] * math.exp((bPrime0 * uu[j])/(bPrime1 * math.sqrt(m0[j])))
                  for j in range(dim)]
            m2 = [m0[j] * math.exp((-bPrime0 * uu[j])/(bPrime2 * math.sqrt(m0[j])))
                  for j in range(dim)]
            break
        except OverflowError:
            print "Overflow error in splitComp() (%2i)" % safety
            safety += 1
            if safety >= 100:
                theProposal.doAbort = True
                #print "Too many overflows in splitComp.  Aborting!"
                return
            uu = [random.normalvariate(0., 1.) for i in range(dim)]
                

    if 0:
        # testing ...
        m1 = [m0[j] * math.exp((bPrime0 * uu[j])/bPrime1) for j in range(dim)]
        m2 = [m0[j] * math.exp((-bPrime0 * uu[j])/bPrime2) for j in range(dim)]

    # Long form of the above for debugging --
    if 0:
        m1 = []
        for j in range(dim):
            top = (bPrime0 * uu[j])
            bottom = (bPrime1 * math.sqrt(m0[j]))
            quot = top/bottom
            try:
                myexp = math.exp(quot)
            except OverflowError:
                gm.append("Got overflow error for m1 exp(%f) at j=%i" % (quot, j))
                gm.append("s0 is %f" % s0)
                gm.append("bPrime0 = %f" % bPrime0)
                gm.append("uu[j] = %f" % uu[j])
                gm.append("bPrime1 = %f" % bPrime1)
                gm.append("m0[j] = %f, sqrt=%f" % (m0[j], math.sqrt(m0[j])))
                gm.append("m0 is %s" % m0)
                gm.append("top = %f" % top)
                gm.append("bottom = %f" % bottom)
                raise Glitch, gm
            m1.append(m0[j] * myexp)

        m2 = []
        for j in range(dim):
            top = (-bPrime0 * uu[j])
            bottom = (bPrime2 * math.sqrt(m0[j]))
            quot = top/bottom
            try:
                myexp = math.exp(quot)
            except OverflowError:
                gm.append("Got overflow error for m2 exp(%f) at j=%i" % (quot, j))
                gm.append("s0 is %f" % s0)
                gm.append("-bPrime0 = %f" % -bPrime0)
                gm.append("uu[j] = %f" % uu[j])
                gm.append("bPrime2 = %f" % bPrime2)
                gm.append("m0[j] = %f, sqrt=%f" % (m0[j], math.sqrt(m0[j])))
                gm.append("m0 is %s" % m0)
                gm.append("top = %f" % top)
                gm.append("bottom = %f" % bottom)
                raise Glitch, gm
            m2.append(m0[j] * myexp)

    if 0:
        # Loggified version, as in Gowri-Shankar and Rattray, eqn 7.
        log_m1 = [math.log(m0[j]) + ((bPrime0 * uu[j])/(bPrime1 * math.sqrt(m0[j])))
              for j in range(dim)]
        log_m2 = [math.log(m0[j]) - ((bPrime0 * uu[j])/(bPrime2 * math.sqrt(m0[j])))
              for j in range(dim)]

        try:
            m1 = [math.exp(it) for it in log_m1]
        except OverflowError:
            gm.append("m0 = %s" % m0)
            gm.append("log_m1 = %s" % log_m1)
            gm.append('overflow m1')
            raise Glitch, gm
        try:
            m2 = [math.exp(it) for it in log_m2]
        except OverflowError:
            gm.append("m0 = %s" % m0)
            gm.append("log_m2 = %s" % log_m2)
            gm.append("overflow m2")
            raise Glitch, gm
    
    if 0:
        print m0
        print m1
        print m2
        print 
        
    s1 = sum(m1)
    s2 = sum(m2)
    newVal1 = [it / s1 for it in m1]
    newVal2 = [it / s2 for it in m2]

    #print newVal1
    #print newVal2

    if 1:
        # Peter adds, the following few lines to make sure the vals are more than var.PIVEC_MIN
        isChanged = False
        for vNum in range(len(newVal1)):
            isGood = False
            while not isGood:
                #print "gen %i" % self.mcmc.gen
                if newVal1[vNum] < var.PIVEC_MIN:
                    newVal1[vNum] = (var.PIVEC_MIN - newVal1[vNum]) + var.PIVEC_MIN
                    isChanged = True
                else:
                    isGood = True
        if isChanged:
            s1 = sum(newVal1)
            newVal1 = [it / s1 for it in newVal1]

        isChanged = False
        for vNum in range(len(newVal2)):
            isGood = False
            while not isGood:
                #print "y gen %i" % self.mcmc.gen
                if newVal2[vNum] < var.PIVEC_MIN:
                    newVal2[vNum] = (var.PIVEC_MIN - newVal2[vNum]) + var.PIVEC_MIN
                    isChanged = True
                else:
                    isGood = True
        if isChanged:
            s2 = sum(newVal2)
            newVal2 = [it / s2 for it in newVal2]    
        
    
    #print newVal1
    #print newVal2

    # Log prior ratio
    # We could have a prior on the pool size, reflected in t1.  If all pool sizes are equally probable, then t1 = 0
    t1 = 0.

    if var.rjCompUniformAllocationPrior:
        b = len([n for n in self.propTree.iterNodes()])
        t2 = b * (math.log(mp.rjComp_k) - math.log(mp.rjComp_k + 1))
    else:
        t2 = (b1 * math.log(f1)) + (b2 * math.log(f2)) - (b0 * math.log(f0))

    # t3 is for the prior on comp vectors.  With the Dirichlet prior alpha values all 1, t3 is log Gamma dim
    t3 = pf.gsl_sf_lngamma(dim)
    
    # t4 is for the f values.  
    if var.rjCompUniformAllocationPrior:
        t4 = 0.0 
    else:
        # If its a uniform Dirichlet, then t4 = log k, where k is from before the split
        t4 = math.log(mp.rjComp_k)

    self.logPriorRatio = t1 + t2 + t3 + t4
    
    # Log proposal ratio
    if mp.rjComp_k == 1:
        t1 = math.log(0.5)
    else:
        t1 = 0.
    if var.rjCompUniformAllocationPrior:
        t2 = b0 * math.log(2.) 
    else:
        t2 = (b0 * math.log(f0)) - (b1 * math.log(f1)) - (b2 * math.log(f2))          # this was changed 26 sept
        
    # for t3, below, do some pre-calculations
    sum_uu2 = sum([u * u for u in uu])
    lastTerm = -pf.gsl_sf_lngamma(p0) + (0.5 * sum_uu2) + \
                 ((dim/2.) * math.log(2 * math.pi))
    t3 = (s0 - s1 - s2) + ((p0 - 1.) * (math.log(s1) + math.log(s2) - math.log(s0))) + lastTerm
    self.logProposalRatio = t1 + t2 + t3
    #print t1,t2,t3,s0,s1,s2
    
    #self.logProposalRatio = 0.

    # The Jacobian
    lastTerm = 0.5 * sum([math.log(v) for v in pi0.val])    # added 26 sept
    t1 = ((((3. * dim) - 2.)/2.) * math.log(s0)) - ((dim - 1.) * (math.log(s1) + math.log(s2))) + lastTerm
    t2 = (2. * dim * math.log(bPrime0)) - (dim * (math.log(bPrime1) + math.log(bPrime2)))
    t3 = sum([uu[j]/(math.sqrt(s0 * pi0.val[j])) for j in range(dim)])
    t3 = ((bPrime0 * (bPrime2 - bPrime1))/(bPrime1 * bPrime2)) * t3

    if var.rjCompUniformAllocationPrior:
        self.logJacobian = t1 + t2 + t3 
    else:
        self.logJacobian = t1 + t2 + t3 + math.log(f0)
        
    # We will now make pi1 and pi2.  The pi1 will be made from pi0,
    # and pi2 will be popped from the notInPool list.  We have newVal1
    # and newVal2 which will be their vals, and we assign them to
    # nodes in beta1 and beta2, and give them rj_f values of f1 and
    # f2.
    pi1 = pi0
    pi1.val = newVal1
    pi1.rj_f = f1
    for n in beta1:
        n.parts[theProposal.pNum].compNum = pi1.num  # not needed, its already that.
        pf.p4_setCompNum(n.cNode, theProposal.pNum, pi1.num)
    pi1.nNodes = b1

    pi2 = notInPool.pop()
    pi2.val = newVal2
    pi2.rj_f = f2
    for n in beta2:
        n.parts[theProposal.pNum].compNum = pi2.num
        pf.p4_setCompNum(n.cNode, theProposal.pNum, pi2.num)
    pi2.nNodes = b2
    pi2.rj_isInPool = True

    self.propTree.model.parts[theProposal.pNum].rjComp_k += 1
    #print "...finished proposeSplitComp()"


def proposeMergeComp(self, theProposal):
    gm = ['Chain.proposeMergeComp()']

    mp = self.propTree.model.parts[theProposal.pNum]
    assert mp.rjComp
    dim = mp.dim
    p0 = theProposal.tuning

    # Check that k is more than 1.  This should have been checked before, but check again.
    #print "rjComp_k is currently %i, with %i comps" % (mp.rjComp_k, mp.nComps)
    if mp.rjComp_k <= 1:
        gm.append("part %i, rjComp_k = %i" % (theProposal.pNum, mp.rjComp_k))
        pool = [c for c in mp.comps if c.rj_isInPool]
        gm.append('len of pool = %i (should be the same as rjComp_k)' % len(pool))
        gm.append("rjComp_k, the pool size, should be more than 1 for a merge.  This isn't.")
        raise Glitch, gm

    # Choose two comps (to make into one).  They must be in the pool.
    pool = [c for c in mp.comps if c.rj_isInPool]
    assert len(pool) == mp.rjComp_k
    pi1, pi2 = random.sample(pool, 2)
    #print "proposing to merge comps %i and %i" % (pi1.num, pi2.num)

    beta1 = []
    beta2 = []
    for n in self.propTree.iterNodes():
        theCompNum = n.parts[theProposal.pNum].compNum
        if theCompNum == pi1.num:
            beta1.append(n)
        elif theCompNum == pi2.num:
            beta2.append(n)
    beta0 = beta1 + beta2
    b1 = float(len(beta1))
    b2 = float(len(beta2))
    b0 = float(b1 + b2)
    assert len(beta0) == b0
    bPrime0 = b0 + 2.
    bPrime1 = b1 + 1.
    bPrime2 = b2 + 1.
    f1 = pi1.rj_f
    f2 = pi2.rj_f
    f0 = f1 + f2

    # Obtain composition vector proposal
    s1 = random.gammavariate(p0, 1.)
    s2 = random.gammavariate(p0, 1.)
    m1 = [v * s1 for v in pi1.val]
    m2 = [v * s2 for v in pi2.val]
    #print "m1 = ", m1
    #print "m2 = ", m2
    #print b0, b1, b2, bPrime0, bPrime1, bPrime2
    m0 = [math.exp(((bPrime1/bPrime0) * math.log(m1[j])) + ((bPrime2/bPrime0) * math.log(m2[j])))
          for j in range(dim)]
    
    #print "m0 = ", m0
    s0 = sum(m0)

    newVal0 = [m0k / s0 for m0k in m0]

    # Log prior ratio
    # We could have a prior on the pool size, reflected in t1.  If all pool sizes are equally probable, then t1 = 0
    t1 = 0.

    if var.rjCompUniformAllocationPrior:
        b = len([n for n in self.propTree.iterNodes()])
        t2 = b * (math.log(mp.rjComp_k) - math.log(mp.rjComp_k - 1))
    else:
        t2 = (b0 * math.log(f0)) - (b1 * math.log(f1)) - (b2 * math.log(f2))
        
    # t3 is for the prior on comp vectors.  With the Dirichlet prior alpha values all 1, t3 is log Gamma dim
    t3 = -pf.gsl_sf_lngamma(dim)

    # t4 is for the f values. 
    if var.rjCompUniformAllocationPrior:
        t4 = 0.0 
    else:
        # If its a uniform Dirichlet, then t4 = - log (k - 1), where k is from before the merge
        t4 = -math.log(mp.rjComp_k - 1)
        
    self.logPriorRatio = t1 + t2 + t3 + t4
    
    # Log proposal ratio
    if mp.rjComp_k == mp.nComps:  # nComps is k_max
        t1 = math.log(0.5)
    else:
        t1 = 0.

    if var.rjCompUniformAllocationPrior:
        t2 = - (b0 * math.log(2.))
    else:
        t2 = (b1 * math.log(f1)) + (b2 * math.log(f2)) - (b0 * math.log(f0))
    
    # for t3, below, do some pre-calculations
    uu = [(bPrime1/bPrime0) * math.sqrt(m0[j]) * (math.log(m1[j]) - math.log(m0[j])) for j in range(dim)]
    sum_uu2 = sum([u * u for u in uu])
    lastTerm = pf.gsl_sf_lngamma(p0) - (0.5 * sum_uu2) - \
                 ((dim/2.) * math.log(2 * math.pi))
    t3 = (s1 + s2 - s0) - ((p0 - 1.) * (math.log(s1) + math.log(s2) - math.log(s0))) + lastTerm

    #logSterm = ((1. - p0) * (math.log(s1) + math.log(s2) - math.log(s0)))
    #print "s1=%.1f s2=%.1f s0=%.1f    sum_uu2=%.1f  logGamma(p0)=%.1f, logSterm=%.1f" % (
    #    s1, s2, s0, sum_uu2, pf.gsl_sf_lngamma(p0), logSterm)
    self.logProposalRatio = t1 + t2 + t3

    #self.logProposalRatio = 20.

    # The Jacobian
    lastTerm = 0.5 * sum([math.log(v) for v in newVal0])  # new 26 sept
    t1 = ((((3. * dim) - 2.)/2.) * math.log(s0)) - ((dim - 1.) * (math.log(s1) + math.log(s2))) + lastTerm
    t2 = (2. * dim * math.log(bPrime0)) - (dim * (math.log(bPrime1) + math.log(bPrime2)))
    t3 = sum([uu[j]/(math.sqrt(s0 * newVal0[j])) for j in range(dim)])
    t3 = ((bPrime0 * (bPrime2 - bPrime1))/(bPrime1 * bPrime2)) * t3
    if var.rjCompUniformAllocationPrior:
        self.logJacobian = -(t1 + t2 + t3)
    else:
        self.logJacobian = -(t1 + t2 + t3 + math.log(f0))


    # Merge pi1 and pi2 => pi0, where pi0 is actually pi1, re-used, by
    # giving "0" values to pi1 = pi0
    pi1.rj_f = f0
    pi1.val = newVal0
    for n in beta0:
        n.parts[theProposal.pNum].compNum = pi1.num
        pf.p4_setCompNum(n.cNode, theProposal.pNum, pi1.num)
    pi1.nNodes = b0
    mp.rjComp_k -= 1
    pi2.rj_isInPool = False
    pi2.nNodes = 0
    



def proposeRMatrixWithSlider(self, theProposal):

    #print "rMatrix proposal. the tuning is %s" % theProposal.tuning

    assert var.rMatrixNormalizeTo1
    mtCur = self.curTree.model.parts[theProposal.pNum].rMatrices[theProposal.mtNum]
    mtProp = self.propTree.model.parts[theProposal.pNum].rMatrices[theProposal.mtNum]
    if mtProp.spec == '2p':
        # For 2p, its actually a Dirichlet, not a slider.  All this is
        # stolen from MrBayes, where the default tuning is 50.  In
        # MrBayes, the "alphaDir" is a 2-item list of Dirichlet
        # parameters (not the multiplier) but they are both by default
        # 1, which makes the prior ratio 1.0 and the logPriorRatio
        # zero.

        
        old = [0.0, 0.0]
        old[0] = mtCur.val / (mtCur.val + 1.0)
        old[1] = 1.0 - old[0]
        new = func.dirichlet1(old, theProposal.tuning, var.KAPPA_MIN, var.KAPPA_MAX)
        mtProp.val[0] = new[0] / new[1]

        theSum = 0.0
        for i in range(2):
            theSum += new[i] * theProposal.tuning
        x = pf.gsl_sf_lngamma(theSum)
        for i in range(2):
            x -= pf.gsl_sf_lngamma(new[i] * theProposal.tuning)
        for i in range(2):
            x += ((new[i] * theProposal.tuning) - 1.0) * math.log(old[i])
        theSum = 0.0
        for i in range(2):
            theSum += old[i] * theProposal.tuning
        y = pf.gsl_sf_lngamma(theSum)
        for i in range(2):
            y -= pf.gsl_sf_lngamma(old[i] * theProposal.tuning)
        for i in range(2):
            y += ((old[i] * theProposal.tuning) -1.0) * math.log(new[i])
        self.logProposalRatio = x - y

    
    else: # specified, ones, eg gtr
        mt = self.propTree.model.parts[theProposal.pNum].rMatrices[theProposal.mtNum]

        # mt.val is a numpy array
        assert type(mt.val) == numpy.ndarray

        nRates = len(mt.val)  # eg 6 for dna gtr, not 5
        indxs = random.sample(range(nRates), 2)
        currentAplusB = mt.val[indxs[0]] + mt.val[indxs[1]]
        thisMin = var.RATE_MIN / currentAplusB
        thisMax = 1. - thisMin

        minToMaxDiff = thisMax - thisMin
        thisTuning = theProposal.tuning

        # It is possible that both A
        # and B values are very close to var.RATE_MIN, in which case
        # thisMin and thisMax will both be close to 0.5, and so the tuning
        # will be too much, requiring too many reflections.  In that case,
        # just change the tuning temporarily.
        if thisTuning > minToMaxDiff:
            thisTuning = minToMaxDiff
            #print "temporarily changing the tuning for rMatrix proposal, to", thisTuning

        x = mt.val[indxs[0]] / currentAplusB
        y = x + (thisTuning * (random.random() - 0.5))

        # reflect
        safety = -1
        while 1:
            safety += 1
            if safety > 20:
                gm.append("Did more than 20 reflections -- something is wrong.")
                raise Glitch, gm
            if y < thisMin:
                y = thisMin + (thisMin - y)
            elif y > thisMax:
                y = thisMax - (y - thisMax)
            else:
                break
        #if safety > 1:
        #    print "rMatrix reflections: ", safety
        mt.val[indxs[0]] = y * currentAplusB
        mt.val[indxs[1]] = currentAplusB - mt.val[indxs[0]]

        mySum = 0.0
        for stNum in range(nRates):
            mySum += mt.val[stNum]
        for stNum in range(nRates):
            mt.val[stNum] /= mySum

        self.logProposalRatio = 0.0
        
    self.logPriorRatio = 0.0

def proposeRjRMatrix(self, theProposal):
    gm = ['Chain.proposeRjRMatrix()']
    mp = self.propTree.model.parts[theProposal.pNum]
    assert mp.rjRMatrix
    assert mp.rjRMatrix_k >= 1
    assert mp.rjRMatrix_k <= mp.nRMatrices

    # If the pool size k, mp.rjRMatrix_k, is only 1, then we can only split
    if mp.rjRMatrix_k == 1:
        self.proposeSplitRMatrix(theProposal)

    # If k, mp.rjRMatrix_k is mp.nRMatrices, then we can only merge
    if mp.rjRMatrix_k == mp.nRMatrices:
        self.proposeMergeRMatrix(theProposal)

    # Otherwise, choose randomly.
    if random.random() < 0.5:
        self.proposeSplitRMatrix(theProposal)
    else:
        self.proposeMergeRMatrix(theProposal)

def proposeSplitRMatrix(self, theProposal):
    gm = ['Chain.proposeSplitRMatrix()']
    # var.rjRMatrixUniformAllocationPrior  True by default
    # theProposal.tuning 300.  becomes p0 below

    mp = self.propTree.model.parts[theProposal.pNum]
    assert mp.rjRMatrix
    rDim = ((mp.dim * mp.dim) - mp.dim) / 2
    
    # Check that k is less than k_max, which is nRMatrices.  This should have been checked before, but check again.
    #print gm[0], "rjRMatrix_k is currently %i, with %i rMatrices" % (mp.rjRMatrix_k, mp.nRMatrices)
    assert mp.rjRMatrix_k < mp.nRMatrices

    # Select an existing rMatrix from the pool
    pool = [c for c in mp.rMatrices if c.rj_isInPool]
    assert mp.rjRMatrix_k == len(pool)
    notInPool = [c for c in mp.rMatrices if not c.rj_isInPool]
    assert notInPool # or else we can't split
    assert mp.nRMatrices == len(pool) + len(notInPool)
    rm0  = random.choice(pool)

    # The nodes currently associated with rm0
    beta0 = [n for n in self.propTree.iterNodesNoRoot() if n.br.parts[theProposal.pNum].rMatrixNum == rm0.num]
    b0 = float(len(beta0))
    #print gm[0], "rMatrix %i is chosen, f=%f, currently on nodes" % (rm0.num, rm0.rj_f), [n.nodeNum for n in beta0]

    # Divvy up the contents of beta0 into (new) beta1 and beta2, based on probability u
    if var.rjRMatrixUniformAllocationPrior:
        u = 0.5
    else:
        u = random.random()

    beta1 = []
    beta2 = []
    for it in beta0:
        r = random.random()
        if r < u:
            beta1.append(it)
        else:
            beta2.append(it)
    b1 = float(len(beta1))
    b2 = float(len(beta2))
    bPrime0 = b0 + 2.
    bPrime1 = b1 + 1.
    bPrime2 = b2 + 1.

    # Calculation of f1 and f2 depends on u
    f0 = rm0.rj_f
    f1 = u * f0
    f2 = (1.0 - u) * f0

    
    uu = [random.normalvariate(0., 1.) for i in range(rDim)]
    p0 = theProposal.tuning
    s0 = random.gammavariate(p0, 1.)
    m0 = [s0 * it for it in rm0.val]
    #print m0

    # I get a math range error here -- needs debugging.
    #m1 = [m0[j] * math.exp((bPrime0 * uu[j])/(bPrime1 * math.sqrt(m0[j])))
    #      for j in range(rDim)]
    #m2 = [m0[j] * math.exp((-bPrime0 * uu[j])/(bPrime2 * math.sqrt(m0[j])))
    #      for j in range(rDim)]
    safety = 0
    while 1:
        try:
            m1 = [m0[j] * math.exp((bPrime0 * uu[j])/(bPrime1 * math.sqrt(m0[j])))
                  for j in range(rDim)]
            m2 = [m0[j] * math.exp((-bPrime0 * uu[j])/(bPrime2 * math.sqrt(m0[j])))
                  for j in range(rDim)]
            break
        except OverflowError:
            print "Overflow error in splitRMatrix() (%2i)" % safety
            safety += 1
            if safety >= 100:
                theProposal.doAbort = True
                print "Too many overflows in splitComp.  Aborting!"
                return
            uu = [random.normalvariate(0., 1.) for i in range(rDim)]

    if 0:
        # Long form of the above for debugging --
        m1 = []
        for j in range(rDim):
            top = (bPrime0 * uu[j])
            bottom = (bPrime1 * math.sqrt(m0[j]))
            quot = top/bottom
            try:
                myexp = math.exp(quot)
            except OverflowError:
                gm.append("Got overflow error for m1 exp(%f) at j=%i" % (quot, j))
                gm.append("bPrime0 = %f" % bPrime0)
                gm.append("uu[j] = %f" % uu[j])
                gm.append("bPrime1 = %f" % bPrime1)
                gm.append("m0[j] = %f, sqrt=%f" % (m0[j], math.sqrt(m0[j])))
                gm.append("m0 is %s" % m0)
                gm.append("top = %f" % top)
                gm.append("bottom = %f" % bottom)
                raise Glitch, gm
            m1.append(m0[j] * myexp)

        m2 = []
        for j in range(rDim):
            top = (-bPrime0 * uu[j])
            bottom = (bPrime2 * math.sqrt(m0[j]))
            quot = top/bottom
            try:
                myexp = math.exp(quot)
            except OverflowError:
                gm.append("Got overflow error for m2 exp(%f) at j=%i" % (quot, j))
                gm.append("-bPrime0 = %f" % -bPrime0)
                gm.append("uu[j] = %f" % uu[j])
                gm.append("bPrime2 = %f" % bPrime2)
                gm.append("m0[j] = %f, sqrt=%f" % (m0[j], math.sqrt(m0[j])))
                gm.append("m0 is %s" % m0)
                gm.append("top = %f" % top)
                gm.append("bottom = %f" % bottom)
                raise Glitch, gm
            m2.append(m0[j] * myexp)
            
        
        
    s1 = sum(m1)
    s2 = sum(m2)
    newVal1 = [it / s1 for it in m1]
    newVal2 = [it / s2 for it in m2]

    #print newVal1
    #print newVal2

    if 1:
        # Peter adds, the following few lines to get the vals more than var.RATE_MIN
        isChanged = False
        for vNum in range(len(newVal1)):
            isGood = False
            while not isGood:
                #print "gen %i" % self.mcmc.gen
                if newVal1[vNum] < var.RATE_MIN:
                    newVal1[vNum] = (var.RATE_MIN - newVal1[vNum]) + var.RATE_MIN
                    isChanged = True
                else:
                    isGood = True
        if isChanged:
            s1 = sum(newVal1)
            newVal1 = [it / s1 for it in newVal1]

        isChanged = False
        for vNum in range(len(newVal2)):
            isGood = False
            while not isGood:
                #print "y gen %i" % self.mcmc.gen
                if newVal2[vNum] < var.RATE_MIN:
                    newVal2[vNum] = (var.RATE_MIN - newVal2[vNum]) + var.RATE_MIN
                    isChanged = True
                else:
                    isGood = True
        if isChanged:
            s2 = sum(newVal2)
            newVal2 = [it / s2 for it in newVal2]    
        
    
    #print newVal1
    #print newVal2

    # Log prior ratio
    # We could have a prior on the pool size, reflected in t1.  If all pool sizes are equally probable, then t1 = 0
    t1 = 0.

    if var.rjRMatrixUniformAllocationPrior:
        b = len([n for n in self.propTree.iterNodesNoRoot()])
        t2 = b * (math.log(mp.rjRMatrix_k) - math.log(mp.rjRMatrix_k + 1))
    else:
        t2 = (b1 * math.log(f1)) + (b2 * math.log(f2)) - (b0 * math.log(f0))

    # t3 is for the prior on rMatrices.  With the Dirichlet prior alpha values all 1, t3 is log Gamma rDim
    t3 = pf.gsl_sf_lngamma(rDim)
    
    # t4 is for the f values.  
    if var.rjRMatrixUniformAllocationPrior:
        t4 = 0.0 
    else:
        # If its a uniform Dirichlet, then t4 = log k, where k is from before the split
        t4 = math.log(mp.rjRMatrix_k)

    self.logPriorRatio = t1 + t2 + t3 + t4
    
    # Log proposal ratio
    if mp.rjRMatrix_k == 1:
        t1 = math.log(0.5)
    else:
        t1 = 0.
    if var.rjRMatrixUniformAllocationPrior:
        t2 = b0 * math.log(2.) 
    else:
        t2 = (b0 * math.log(f0)) - (b1 * math.log(f1)) - (b2 * math.log(f2))          # this was changed 26 sept
        
    # for t3, below, do some pre-calculations
    sum_uu2 = sum([u * u for u in uu])
    lastTerm = -pf.gsl_sf_lngamma(p0) + (0.5 * sum_uu2) + \
                 ((rDim/2.) * math.log(2 * math.pi))
    t3 = (s0 - s1 - s2) + ((p0 - 1.) * (math.log(s1) + math.log(s2) - math.log(s0))) + lastTerm
    self.logProposalRatio = t1 + t2 + t3
    #print t1,t2,t3,s0,s1,s2
    
    #self.logProposalRatio = 0.

    # The Jacobian
    lastTerm = 0.5 * sum([math.log(v) for v in rm0.val])    # added 26 sept
    t1 = ((((3. * rDim) - 2.)/2.) * math.log(s0)) - ((rDim - 1.) * (math.log(s1) + math.log(s2))) + lastTerm
    t2 = (2. * rDim * math.log(bPrime0)) - (rDim * (math.log(bPrime1) + math.log(bPrime2)))
    t3 = sum([uu[j]/(math.sqrt(s0 * rm0.val[j])) for j in range(rDim)])
    t3 = ((bPrime0 * (bPrime2 - bPrime1))/(bPrime1 * bPrime2)) * t3

    if var.rjRMatrixUniformAllocationPrior:
        self.logJacobian = t1 + t2 + t3 
    else:
        self.logJacobian = t1 + t2 + t3 + math.log(f0)
        
    # We will now make rm1 and rm2.  The rm1 will be made from rm0,
    # and rm2 will be popped from the notInPool list.  We have newVal1
    # and newVal2 which will be their vals, and we assign them to
    # nodes in beta1 and beta2, and give them rj_f values of f1 and
    # f2.
    rm1 = rm0
    for rNum in range(rDim):
        rm1.val[rNum] = newVal1[rNum]
    rm1.rj_f = f1
    for n in beta1:
        n.br.parts[theProposal.pNum].rMatrixNum = rm1.num  # not needed, its already that.
        pf.p4_setRMatrixNum(n.cNode, theProposal.pNum, rm1.num)
    rm1.nNodes = b1

    rm2 = notInPool.pop()
    for rNum in range(rDim):
        rm2.val[rNum] = newVal2[rNum]
    rm2.rj_f = f2
    for n in beta2:
        n.br.parts[theProposal.pNum].rMatrixNum = rm2.num
        pf.p4_setRMatrixNum(n.cNode, theProposal.pNum, rm2.num)
    rm2.nNodes = b2
    rm2.rj_isInPool = True

    self.propTree.model.parts[theProposal.pNum].rjRMatrix_k += 1
    #print "...finished proposeSplitRMatrix()"


def proposeMergeRMatrix(self, theProposal):
    gm = ['Chain.proposeMergeRMatrix()']

    mp = self.propTree.model.parts[theProposal.pNum]
    assert mp.rjRMatrix
    rDim = ((mp.dim * mp.dim) - mp.dim) / 2
    p0 = theProposal.tuning

    # Check that k is more than 1.  This should have been checked before, but check again.
    #print "rjRMatrix_k is currently %i, with %i rMatrices" % (mp.rjRMatrix_k, mp.nRMatrices)
    if mp.rjRMatrix_k <= 1:
        gm.append("part %i, rjRMatrix_k = %i" % (theProposal.pNum, mp.rjRMatrix_k))
        pool = [c for c in mp.rMatrices if c.rj_isInPool]
        gm.append('len of pool = %i (should be the same as rjRMatrix_k)' % len(pool))
        gm.append("rjRMatrix_k, the pool size, should be more than 1 for a merge.  This isn't.")
        raise Glitch, gm

    # Choose two rMatrices (to make into one).  They must be in the pool.
    pool = [c for c in mp.rMatrices if c.rj_isInPool]
    assert len(pool) == mp.rjRMatrix_k
    rm1, rm2 = random.sample(pool, 2)
    #print "proposing to merge rMatrices %i and %i" % (rm1.num, rm2.num)

    beta1 = []
    beta2 = []
    for n in self.propTree.iterNodesNoRoot():
        theRMatrixNum = n.br.parts[theProposal.pNum].rMatrixNum
        if theRMatrixNum == rm1.num:
            beta1.append(n)
        elif theRMatrixNum == rm2.num:
            beta2.append(n)
    beta0 = beta1 + beta2
    b1 = float(len(beta1))
    b2 = float(len(beta2))
    b0 = float(b1 + b2)
    assert len(beta0) == b0
    bPrime0 = b0 + 2.
    bPrime1 = b1 + 1.
    bPrime2 = b2 + 1.
    f1 = rm1.rj_f
    f2 = rm2.rj_f
    f0 = f1 + f2

    # Obtain rMatrix proposal
    s1 = random.gammavariate(p0, 1.)
    s2 = random.gammavariate(p0, 1.)
    m1 = [v * s1 for v in rm1.val]
    m2 = [v * s2 for v in rm2.val]
    #print "m1 = ", m1
    #print "m2 = ", m2
    #print b0, b1, b2, bPrime0, bPrime1, bPrime2
    m0 = [math.exp(((bPrime1/bPrime0) * math.log(m1[j])) + ((bPrime2/bPrime0) * math.log(m2[j])))
          for j in range(rDim)]
    
    #print "m0 = ", m0
    s0 = sum(m0)

    newVal0 = [m0k / s0 for m0k in m0]

    # Log prior ratio
    # We could have a prior on the pool size, reflected in t1.  If all pool sizes are equally probable, then t1 = 0
    t1 = 0.

    if var.rjRMatrixUniformAllocationPrior:
        b = len([n for n in self.propTree.iterNodes()])
        t2 = b * (math.log(mp.rjRMatrix_k) - math.log(mp.rjRMatrix_k - 1))
    else:
        t2 = (b0 * math.log(f0)) - (b1 * math.log(f1)) - (b2 * math.log(f2))
        
    # t3 is for the prior on rMatrices.  With the Dirichlet prior alpha values all 1, t3 is log Gamma rDim
    t3 = -pf.gsl_sf_lngamma(rDim)

    # t4 is for the f values. 
    if var.rjRMatrixUniformAllocationPrior:
        t4 = 0.0 
    else:
        # If its a uniform Dirichlet, then t4 = - log (k - 1), where k is from before the merge
        t4 = -math.log(mp.rjRMatrix_k - 1)
        
    self.logPriorRatio = t1 + t2 + t3 + t4
    
    # Log proposal ratio
    if mp.rjRMatrix_k == mp.nRMatrices:  # nRMatrices is k_max
        t1 = math.log(0.5)
    else:
        t1 = 0.

    if var.rjRMatrixUniformAllocationPrior:
        t2 = - (b0 * math.log(2.))
    else:
        t2 = (b1 * math.log(f1)) + (b2 * math.log(f2)) - (b0 * math.log(f0))
    
    # for t3, below, do some pre-calculations
    uu = [(bPrime1/bPrime0) * math.sqrt(m0[j]) * (math.log(m1[j]) - math.log(m0[j])) for j in range(rDim)]
    sum_uu2 = sum([u * u for u in uu])
    lastTerm = pf.gsl_sf_lngamma(p0) - (0.5 * sum_uu2) - \
                 ((rDim/2.) * math.log(2 * math.pi))
    t3 = (s1 + s2 - s0) - ((p0 - 1.) * (math.log(s1) + math.log(s2) - math.log(s0))) + lastTerm

    #logSterm = ((1. - p0) * (math.log(s1) + math.log(s2) - math.log(s0)))
    #print "s1=%.1f s2=%.1f s0=%.1f    sum_uu2=%.1f  logGamma(p0)=%.1f, logSterm=%.1f" % (
    #    s1, s2, s0, sum_uu2, pf.gsl_sf_lngamma(p0), logSterm)
    self.logProposalRatio = t1 + t2 + t3

    #self.logProposalRatio = 20.

    # The Jacobian
    lastTerm = 0.5 * sum([math.log(v) for v in newVal0])  # new 26 sept
    t1 = ((((3. * rDim) - 2.)/2.) * math.log(s0)) - ((rDim - 1.) * (math.log(s1) + math.log(s2))) + lastTerm
    t2 = (2. * rDim * math.log(bPrime0)) - (rDim * (math.log(bPrime1) + math.log(bPrime2)))
    t3 = sum([uu[j]/(math.sqrt(s0 * newVal0[j])) for j in range(rDim)])
    t3 = ((bPrime0 * (bPrime2 - bPrime1))/(bPrime1 * bPrime2)) * t3
    if var.rjRMatrixUniformAllocationPrior:
        self.logJacobian = -(t1 + t2 + t3)
    else:
        self.logJacobian = -(t1 + t2 + t3 + math.log(f0))


    # Merge rm1 and rm2 => rm0, where rm0 is actually rm1, re-used, by
    # giving "0" values to rm1 = rm0
    rm1.rj_f = f0
    for rNum in range(rDim):
        rm1.val[rNum] = newVal0[rNum]
    for n in beta0:
        n.br.parts[theProposal.pNum].rMatrixNum = rm1.num
        pf.p4_setRMatrixNum(n.cNode, theProposal.pNum, rm1.num)
    rm1.nNodes = b0
    mp.rjRMatrix_k -= 1
    rm2.rj_isInPool = False
    rm2.nNodes = 0
    


def proposeGdasrv(self, theProposal):

    # This is a multiplier proposal.

    gm = ["Chain.proposeGdasrv()"]
    mt = self.propTree.model.parts[theProposal.pNum].gdasrvs[theProposal.mtNum]

    # We can't have alpha less than about 1.e-16, or DiscreteGamma hangs.
    # But that is moot, as var.GAMMA_SHAPE_MIN is much bigger

    # Dont' do something like the following, cuz mt.val is a property that invokes a function.
    #mt.val = newVal
    #mt.val /= theProposal.tuning

    # mt.val is a numpy.ndarray type, an array with 1 element.
    assert type(mt.val) == numpy.ndarray
    oldVal = mt.val
    newVal = oldVal * math.exp(theProposal.tuning * (random.random() - 0.5))

    isGood = False
    while not isGood:
        if newVal < var.GAMMA_SHAPE_MIN:
            newVal = var.GAMMA_SHAPE_MIN * var.GAMMA_SHAPE_MIN / newVal
        elif newVal > var.GAMMA_SHAPE_MAX:
            newVal = var.GAMMA_SHAPE_MAX * var.GAMMA_SHAPE_MAX / newVal
        else:
            isGood = True

    #print type(self.logProposalRatio), type(self.logPriorRatio),
    self.logProposalRatio = math.log(newVal / oldVal)

    self.logPriorRatio = 0.0
    # as in proposeBrLen()
    #self.logPriorRatio = self.mcmc.tunings.parts[theProposal.pNum].gdasrvPriorLambda * float(oldVal - newVal)
    mt.val = newVal
    assert type(mt.val) == numpy.ndarray
    #print type(self.logProposalRatio), type(self.logPriorRatio),
    #print self.logProposalRatio, self.logPriorRatio

def proposePInvar(self, theProposal):
    mt = self.propTree.model.parts[theProposal.pNum].pInvar

    # Slider proposal
    mt.val += (random.random() - 0.5) * theProposal.tuning

    # Linear reflect
    isGood = False
    #while (mt.val < var.PINVAR_MIN) or (mt.val > var.PINVAR_MAX):
    while not isGood:
        if mt.val < var.PINVAR_MIN:
            mt.val = (var.PINVAR_MIN - mt.val) + var.PINVAR_MIN
        elif mt.val > var.PINVAR_MAX:
            mt.val = var.PINVAR_MAX - (mt.val - var.PINVAR_MAX)
        else:
            isGood = True

    self.logProposalRatio = 0.0
    self.logPriorRatio = 0.0

def proposeRelRate(self, theProposal):
    for pNum in range(self.propTree.model.nParts):
        mp = self.propTree.model.parts[pNum]
        ran = (random.random() - 0.5) * theProposal.tuning
        mp.relRate += ran
        isGood = False
        #while (mp.relRate < var.RELRATE_MIN) or (mp.relRate > var.RELRATE_MAX):
        while not isGood:
            if mp.relRate < var.RELRATE_MIN:
                mp.relRate = (var.RELRATE_MIN - mp.relRate) + var.RELRATE_MIN
            elif mp.relRate > var.RELRATE_MAX:
                mp.relRate = var.RELRATE_MAX - (mp.relRate - var.RELRATE_MAX)
            else:
                isGood = True
                
    totDataLen = 0
    for p in self.propTree.data.parts:
        totDataLen += p.nChar
    fact = 0.0
    for pNum in range(self.propTree.model.nParts):
        fact += (self.propTree.model.parts[pNum].relRate * self.propTree.data.parts[pNum].nChar)
    fact = float(totDataLen) / fact
    for p in self.propTree.model.parts:
        p.relRate *= fact
    
    self.logProposalRatio = 0.0
    self.logPriorRatio = 0.0

def proposeCompLocation(self, theProposal):
    gm = ["Chain.proposeCompLocation()"]
    mp = self.propTree.model.parts[theProposal.pNum]
    #nMT = self.propTree.model.parts[theProposal.pNum].nComps
    if mp.rjComp:
        pool = [c for c in mp.comps if c.rj_isInPool]
        # We need at least 2 comps, because one of them will be the
        # current comp for the node chosen below, and we need at least
        # one other to change to.
        if len(pool) < 2:
            theProposal.doAbort = True
            return True
    validNodeNums = [n for n in self.propTree.preOrder if n != var.NO_ORDER]
    validNodes = [self.propTree.nodes[n] for n in validNodeNums]
    if mp.rjComp:
        validNodes =  [n for n in validNodes if (
        mp.comps[n.parts[theProposal.pNum].compNum].rj_isInPool)]
    else:
        validNodes = [n for n in validNodes if (
            mp.comps[n.parts[theProposal.pNum].compNum].nNodes > 1)]
    if not validNodes:
        theProposal.doAbort = True
        return True
    theNode = random.choice(validNodes)
    currentNum = theNode.parts[theProposal.pNum].compNum
    if mp.rjComp:
        validCompNumbers = [c.num for c in pool if c.num is not currentNum]
        if not validCompNumbers:
            theProposal.doAbort = True
            return True
    else:
        validCompNumbers = [c.num for c in mp.comps if c.num is not currentNum]
    #proposedNum = currentNum
    #while proposedNum == currentNum:
    #    proposedNum = random.randrange(nMT)
    proposedNum = random.choice(validCompNumbers)
    if 0 and self.mcmc.gen == 399:
        self.propTree.draw()
        print "proposeCompLocation().  node %i, before=%i, new=%s" % (theNode.nodeNum, currentNum, proposedNum)
    self.propTree.model.parts[theProposal.pNum].comps[currentNum].nNodes -= 1
    self.propTree.model.parts[theProposal.pNum].comps[proposedNum].nNodes += 1
    theNode.parts[theProposal.pNum].compNum = proposedNum
    self.logProposalRatio = 0.0
    #self.logPriorRatio = 0.0
    self.logPriorRatio = theProposal.tuning


def proposeRMatrixLocation(self, theProposal):
    #gm = ["proposeRMatrixLocation()"]
    mp = self.propTree.model.parts[theProposal.pNum]
    #nMT = mp.nRMatrices
    if mp.rjRMatrix:
        pool = [c for c in mp.rMatrices if c.rj_isInPool]
        # We need at least 2 rMatrices, because one of them will be the
        # current rMatrix for the node chosen below, and we need at least
        # one other to change to.
        if len(pool) < 2:
            theProposal.doAbort = True
            return True
    validNodeNums = [n for n in self.propTree.preOrder if n != var.NO_ORDER]
    validNodeNums.remove(self.propTree.root.nodeNum)
    validNodes = [self.propTree.nodes[n] for n in validNodeNums]
    
    if mp.rjRMatrix:
        validNodes =  [n for n in validNodes if (
        mp.rMatrices[n.br.parts[theProposal.pNum].rMatrixNum].rj_isInPool)]
    else:
        validNodes = [n for n in validNodes if (
            mp.rMatrices[n.br.parts[theProposal.pNum].rMatrixNum].nNodes > 1)]
    
    if not validNodes:
        theProposal.doAbort = True
        return True
    theNode = random.choice(validNodes)
    currentNum = theNode.br.parts[theProposal.pNum].rMatrixNum
    
    if mp.rjRMatrix:
        validRMatrixNumbers = [c.num for c in pool if c.num is not currentNum]
        if not validRMatrixNumbers:
            theProposal.doAbort = True
            return True
    else:
        validRMatrixNumbers = [c.num for c in mp.rMatrices if c.num is not currentNum]
        
    #proposedNum = currentNum
    #while proposedNum == currentNum:
    #    proposedNum = random.randrange(nMT)
    #print "proposeRMatrixLocation().  node %i, before=%i, new=%s" % (theNode.nodeNum, currentNum, proposedNum)
    proposedNum = random.choice(validRMatrixNumbers)
    
    self.propTree.model.parts[theProposal.pNum].rMatrices[currentNum].nNodes -= 1
    self.propTree.model.parts[theProposal.pNum].rMatrices[proposedNum].nNodes += 1
    theNode.br.parts[theProposal.pNum].rMatrixNum = proposedNum
    self.logProposalRatio = 0.0
    #self.logPriorRatio = 0.0
    self.logPriorRatio = theProposal.tuning

def proposeGdasrvLocation(self, theProposal):
    #gm = ["proposeGdasrvLocation()"]
    nMT = self.propTree.model.parts[theProposal.pNum].nGdasrvs
    validNodeNums = [n for n in self.propTree.preOrder if n != var.NO_ORDER]
    validNodeNums.remove(self.propTree.root.nodeNum)
    validNodes = [self.propTree.nodes[n] for n in validNodeNums]
    validNodes = [n for n in validNodes if (
        self.propTree.model.parts[theProposal.pNum].gdasrvs[n.br.parts[theProposal.pNum].gdasrvNum].nNodes > 1)]
    if not validNodes:
        theProposal.doAbort = True
        return True
    theNode = random.choice(validNodes)
    currentNum = theNode.br.parts[theProposal.pNum].gdasrvNum
    proposedNum = currentNum
    while proposedNum == currentNum:
        proposedNum = random.randrange(nMT)
    #print "proposeGdasrvLocation().  node %i, before=%i, new=%s" % (theNode.nodeNum, currentNum, proposedNum)
    self.propTree.model.parts[theProposal.pNum].gdasrvs[currentNum].nNodes -= 1
    self.propTree.model.parts[theProposal.pNum].gdasrvs[proposedNum].nNodes += 1
    theNode.br.parts[theProposal.pNum].gdasrvNum = proposedNum
    self.logProposalRatio = 0.0
    self.logPriorRatio = 0.0

def proposeCmd1CompDir(self, theProposal):
    gm = ['Chain.proposeCmd1CompDir()']

    #print gm[0], theProposal.pNum, theProposal.mtNum

    mp = self.propTree.model.parts[theProposal.pNum]

    # The proposal mtNum is -1, meaning do all, or any
    assert theProposal.mtNum == -1
    mt = random.choice(mp.comps)

    # mt.val is a list of floats, not a numpy.ndarray
    #print type(mt.val), type(mt.val[0]), mt.val, mt.num

    # This method uses func.dirichlet1, which is for lists not numpy
    # arrays.  A copy of inSeq is made, and the copy is modified and
    # returned.
    #dirichlet1(inSeq, alpha, theMin, theMax)
    newVal = func.dirichlet1(mt.val, mp.cmd1_p, var.PIVEC_MIN, 1 - var.PIVEC_MIN)

    # proposal ratio
    dirPrams = [mp.cmd1_p * v for v in newVal]
    logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(mt.val))
    dirPrams = [mp.cmd1_p * v for v in mt.val]
    logPdfProps = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(newVal))
    self.logProposalRatio = logPdfProps - logPdfCurrs

    # prior ratio
    dirPrams = [mp.cmd1_alpha * v for v in mp.cmd1_pi0]
    logPdfProps = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(newVal))
    logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(mt.val))
    self.logPriorRatio = logPdfProps - logPdfCurrs

    mt.val = newVal
    
def proposeCmd1Comp0Dir(self, theProposal):
    gm = ['Chain.proposeCmd1Comp0Dir()']

    #print gm[0], theProposal.pNum, theProposal.mtNum

    mp = self.propTree.model.parts[theProposal.pNum]

    # The proposal mtNum is -1, meaning do all, or any
    assert theProposal.mtNum == -1
    curVal = mp.cmd1_pi0

    # mt.val is a list of floats, not a numpy.ndarray
    #print type(mt.val), type(mt.val[0]), mt.val, mt.num

    # This method uses func.dirichlet1, which is for lists not numpy
    # arrays.  A copy of inSeq is made, and the copy is modified and
    # returned.
    #dirichlet1(inSeq, alpha, theMin, theMax)
    newVal = func.dirichlet1(curVal, mp.cmd1_q, var.PIVEC_MIN, 1 - var.PIVEC_MIN)

    # proposal ratio
    dirPrams = [mp.cmd1_q * v for v in newVal]
    logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(curVal))
    dirPrams = [mp.cmd1_q * v for v in curVal]
    logPdfProps = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(newVal))
    self.logProposalRatio = logPdfProps - logPdfCurrs

    # prior ratio
    dirPrams = [mp.cmd1_s * v for v in [1.0] * mp.dim]
    logPdfProps = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(newVal))
    logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(curVal))
    self.logPriorRatio = logPdfProps - logPdfCurrs

    for pi_i in mp.comps:
        dirPrams = [mp.cmd1_alpha * v for v in newVal]
        logPdfProps = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(pi_i.val))
        dirPrams = [mp.cmd1_alpha * v for v in curVal]
        logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(pi_i.val))
        diff = logPdfProps - logPdfCurrs
        self.logPriorRatio += diff

    mp.cmd1_pi0 = newVal


def proposeCmd1AllCompDir(self, theProposal):
    gm = ['Chain.proposeCmd1AllCompDir()']
    # all the comps, including pi0, in one go.
    
    mp = self.propTree.model.parts[theProposal.pNum]

    # The proposal mtNum is -1, meaning do all, or any
    assert theProposal.mtNum == -1

    # First do proposal for pi0

    # mt.val is a list of floats, not a numpy.ndarray
    #print type(mt.val), type(mt.val[0]), mt.val, mt.num
    # This method uses func.dirichlet1, which is for lists not numpy
    # arrays.  A copy of inSeq is made, and the copy is modified and
    # returned.
    #dirichlet1(inSeq, alpha, theMin, theMax)
    myU = 0.0
    pi0_newVal = func.dirichlet1(mp.cmd1_pi0, mp.cmd1_q, var.PIVEC_MIN, 1 - var.PIVEC_MIN, u=myU)

    # Now do proposals for all the comps in mp.comps, now using u
    # added to the dirichlet prams within func.dirichlet1().  This of
    # course needs to be taken into account when calculating the
    # proposal ratio below.
    piNewVals = [func.dirichlet1(mt.val, mp.cmd1_p, var.PIVEC_MIN, 1 - var.PIVEC_MIN, u=myU) for mt in mp.comps]

    # proposal ratio for pi0
    #myU = 0.0
    dirPrams = [(mp.cmd1_q * v) + myU for v in pi0_newVal]
    logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(mp.cmd1_pi0))
    dirPrams = [(mp.cmd1_q * v) + myU for v in mp.cmd1_pi0]
    logPdfProps = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(pi0_newVal))
    self.logProposalRatio = logPdfProps - logPdfCurrs

    # proposal ratios for all the comps in mp.comps
    for mtNum in range(len(mp.comps)):
        mt = mp.comps[mtNum]
        newVal = piNewVals[mtNum]
        dirPrams = [(mp.cmd1_p * v) + myU for v in newVal]
        logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(mt.val))
        dirPrams = [(mp.cmd1_p * v) + myU for v in mt.val]
        logPdfProps = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(newVal))
        self.logProposalRatio += (logPdfProps - logPdfCurrs)

    # prior ratio for pi0 component
    dirPrams = [mp.cmd1_s * v for v in [1.0] * mp.dim]
    logPdfProps = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(pi0_newVal))
    logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(mp.cmd1_pi0))
    self.logPriorRatio = logPdfProps - logPdfCurrs

    # prior ratios for all the comps in mp.comps
    for mtNum in range(len(mp.comps)):
        mt = mp.comps[mtNum]
        newVal = piNewVals[mtNum]
        dirPrams = [mp.cmd1_alpha * v for v in pi0_newVal]
        logPdfProps = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(newVal))
        dirPrams = [mp.cmd1_alpha * v for v in mp.cmd1_pi0]
        logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(mt.val))
        diff = logPdfProps - logPdfCurrs
        self.logPriorRatio += diff

    # assign the proposals to the prop tree
    mp.cmd1_pi0 = pi0_newVal
    for mtNum in range(len(mp.comps)):
        mt = mp.comps[mtNum]
        newVal = piNewVals[mtNum]
        mt.val = newVal

def proposeCmd1Alpha(self, theProposal):
    gm = ['Chain.proposeCmd1Alpha()']
    MIN = 1.
    MAX = 1000.

    #print gm[0], theProposal.pNum, theProposal.mtNum

    mp = self.propTree.model.parts[theProposal.pNum]

    # The proposal mtNum is -1, meaning do all, or any
    assert theProposal.mtNum == -1
    curVal = mp.cmd1_alpha

    assert curVal <= MAX
    assert curVal >= MIN

    # mt.val is a list of floats, not a numpy.ndarray
    #print type(mt.val), type(mt.val[0]), mt.val, mt.num

    if 0:
        # Do a log scale proposal
        if 0:
            # Make proposals, and if it is outside MIN, MAX, then try again.
            while 1:
                newVal = curVal * math.exp(mp.cmd1_alphaLogScaleProposalTuning * (random.random() - 0.5))
                if (newVal < MIN) or (newVal > MAX):
                    continue
                else:
                    break
        else:
            # If it is outside MIN, MAX, then do logarithmic reflect
            newVal = curVal * math.exp(mp.cmd1_alphaLogScaleProposalTuning * (random.random() - 0.5))
            if 1:
                # Logarithmic reflect if needed
                while (newVal < MIN) or (newVal > MAX):
                    if newVal < MIN:
                        newVal = MIN * MIN / newVal
                    elif newVal > MAX:
                        newVal = MAX * MAX / newVal
    else:
        # Do linear proposal
        newVal = curVal + (mp.cmd1_alphaLinearScaleProposalTuning * (random.random() - 0.5))

        # Linear reflect if needed
        while (newVal < MIN) or (newVal > MAX):
            if newVal < MIN:
                newVal = (MIN - newVal) + MIN
            elif newVal > MAX:
                newVal = MAX - (newVal - MAX)

    #print curVal, newVal

    # proposal ratio
    #self.logProposalRatio = math.log(newVal/curVal)
    self.logProposalRatio = 0.0

    # prior ratio
    if 0:
        # As in Tom's blurb
        # logNormal for self
        zeta = mp.cmd1_LN_a
        sigma = mp.cmd1_LN_t
        pdfProp = pf.gsl_ran_lognormal_pdf(newVal, zeta, sigma)
        pdfCurr = pf.gsl_ran_lognormal_pdf(curVal, zeta, sigma)
        priorRatio = pdfProp/pdfCurr
        self.logPriorRatio = math.log(priorRatio)
    else:
        self.logPriorRatio = 0.0 # flat!

    for pi_i in mp.comps:
        dirPrams = [newVal * v for v in mp.cmd1_pi0]
        logPdfProps = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(pi_i.val))
        dirPrams = [curVal * v for v in mp.cmd1_pi0]
        logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(mp.dim, numpy.array(dirPrams), numpy.array(pi_i.val))
        diff = logPdfProps - logPdfCurrs
        self.logPriorRatio += diff

    mp.cmd1_alpha = newVal

    
