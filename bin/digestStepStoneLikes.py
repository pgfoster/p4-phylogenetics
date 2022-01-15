#!/usr/bin/env python3

if 1:
    # This is needed for the SE calculation
    import rpy2
    from rpy2.robjects.packages import importr
    from rpy2.robjects import FloatVector
    from rpy2.robjects import r
    importr('coda')

import statistics
import math
import numpy
# from scipy.stats import beta
import sys
import argparse

# https://docs.python.org/3/howto/argparse.html
# https://docs.python.org/3/library/argparse.html
# https://docs.python.org/3/library/argparse.html#the-add-argument-method

myDescription = """Digests log likelihoods for stepping stone method.

Wangang Xie, Paul O. Lewis, Yu Fan, Lynn Kuo, Ming-Hui Chen 2011.
Improving Marginal Likelihood Estimation for Bayesian Phylogenetic Model Selection.
Systematic Biology, 60 150-–160 https://doi.org/10.1093/sysbio/syq085

Thanks to Mario dos Reis for explanations and code on his website, https://dosreislab.github.io/2018/08/31/bppr-stepstones.html 

For the variance calculation, this implementation uses the effectiveSize function in the coda package in R, so that needs to be available on your system, as well as the rpy2 Python module (https://rpy2.github.io/doc.html) to enable Python to use R.

The log likelihoods need to be all in a single column in one file, and should go from high (all Posterior) to low beta (all Prior) values.

As in Xie et al, the beta values are given as β_k = (k/K)^{1/α}, where alpha is less than one.  That is a different parameterization to the Beta(alpha,1.0) as given in Scipy.

If you provide an alpha value, this script will use the Xie et al β_k formula to to calculate beta values.  Alternatively, if your beta values are in a file, provide the file name.  The latter option also provides K. 

"""

# See https://stackoverflow.com/questions/35857610/how-to-preserve-newlines-in-argparse-version-output-while-letting-argparse-auto
# for info about RawDescriptionHelpFormatter
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description=myDescription)

parser.add_argument("logLikesFile", type=str, help="Name of the file containing the log likelihoods")
parser.add_argument("-sps", "--samplesPerStep", type=int, help="The number of power posteriors per step")
parser.add_argument("-K", "--K", type=int, help="K, the number of stepping stones ")
parser.add_argument("-alpha", "--alpha", type=float, help="The alpha for beta_k values, as in Xie et al pg 154, eg 0.25 or 0.3")
parser.add_argument("-bt", "--betaValsFile", type=str, help="Name of the file containing the beta values")
parser.add_argument("-doSE", "--doSE", action="store_true", help="Calculate the SE (standard error).  Optional --- default is False")
parser.add_argument("-gr", "--gram", action="store_true", help="Plot using gram, default is False. Gram.baseName is meanLogLikes")
# parser.add_argument("-v", "--verbosity", type=int, choices=[0, 1, 2], default=0, help="requires an argument, must be one of the choices")

# If there is no point in proceeding ...
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()
# print(args)

assert args.samplesPerStep, "Need to specify the --samplesPerStep, the number of samples in each stepping stone"

# Get the beta values, either from a file, or regenerate them.
if args.betaValsFile:
    bb = [float(l) for l in open("beta.txt")]
    K = len(bb)
else:
    assert args.alpha
    assert args.K

    alpha = args.alpha
    K = args.K
    bb = [math.pow((k/K),(1/alpha)) for k in range(K)]
    #print("bb ", bb)

bb.append(1.0)
bdiffs = [bb[i] - bb[i-1] for i in range(1, K + 1)]

maxima = numpy.zeros(K)
zr = numpy.zeros(K)           # mean scaled likelihoods (not log likelihoods)
meanLogLikes = numpy.zeros(K)
ess = numpy.zeros(K)
vzr = numpy.zeros(K)
likes = []

#print(f"bb {bb}")
#print(f"bdiffs {bdiffs}")

# All likelihoods from all beta values
allLikes = [float(l) for l in open(args.logLikesFile)]
myMsg = f"Got {len(allLikes)} log likelihood values, should be samplesPerStep {args.samplesPerStep} * K {K} = {args.samplesPerStep * K}"
assert len(allLikes) == args.samplesPerStep * K, myMsg
allLikes.reverse()


for k in range(K):
    # Get the log likes into a numpy array, logLikes
    start = args.samplesPerStep * k
    logLikes = numpy.array(allLikes[start: start + args.samplesPerStep])

    meanLogLikes[k] = statistics.mean(logLikes)
    #print(f"meanLogLikes {meanLogLikes}")

    logLikes *= bdiffs[k]
    maxima[k] = max(logLikes)
    #print(f"maxima {maxima}")

    likes.append(numpy.exp(logLikes - maxima[k]))

    zr[k] = numpy.mean(likes[k])

for k in range(K):
    # From Mario
    #ess[k] <- coda::effectiveSize(likes[[k]])
    fvR = FloatVector(likes[k])
    ret = r['effectiveSize'](fvR)
    myFloat = float(ret[0])
    ess[k] = myFloat

    vzr[k] = numpy.var(likes[k]) / ess[k]

    # the delta approximation does not work well if vzr/zr^2 > 0.1
    if (vzr[k] / (zr[k] * zr[k])) > 0.1:
       print("unreliable se: var(r_k)/r_k^2 = ", vzr[k] / (zr[k] * zr[k]), f" > 0.1 for k={k}, b=", bb[k])

if 0:
    print("%3s %12s %12s %12s %12s %12s %14s" % ("k", "b", "var(likes)", "ess", "vzr", "zr", "var(r_k)/r_k^2"))
    for k in range(K):
        print(f"{k:3}", end=' ')
        print(f"{bb[k]:12.6}", end=' ')
        print(f"{numpy.var(likes[k]):12.4f}", end=' ')
        print(f"{ess[k]:12.4f}", end=' ')
        print(f"{vzr[k]:12.8f}", end=' ')
        print(f"{zr[k]:12.6f}", end=' ')
        print(f"{vzr[k] / (zr[k] * zr[k]):14.6f}", end=' ')
        print()
    

vmlnl = sum(vzr / (zr * zr))


lnml = sum(numpy.log(zr) + maxima)

if args.gram:
    from gram import Plot
    gp = Plot()
    xx = list(bb[:-1])
    yy = list(meanLogLikes)
    #gp.svgPxForCm = 100
    gp.baseName = 'meanLogLikes'
    #gp.scatter(xx, yy)
    gp.line(xx, yy)
    gp.xAxis.title = r'$\beta_k$'
    gp.yAxis.title = 'mean log likelihood'
    gp.minXToShow = 0
    gp.maxXToShow = 1
    gp.pdf()



#print(meanLogLikes)
#print(maxima)
#print(zr)
#print(ess)
#print(vzr)
# print(vmlnl)
print("log marginal likelihood:", lnml)

