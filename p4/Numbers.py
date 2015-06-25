import os,sys,math,string
import func
from Glitch import Glitch
from Var import var
if var.usePfAndNumpy:
    import numpy
    
class Numbers(object):
    """Simple 1-dimensional data handling.  Emphasis on 'simple'.

    Feed this a list of floats, or the name of a file containing
    floats.  If it is a file, lines are read in only if the first
    character is a digit.  So lines starting with words, '#', or blank
    lines, are ignored (and those ignored lines do not contribute to
    the 'skip', below).  So this should be able to read in MrBayes ``*.p``
    files, sensibly skipping the first lines.

    Options:
    
        col
                  If inThing is a file, it can have multiple columns.
                  This says which one to use.  Zero-based!
        skip
                  If inThing is a file, you can skip lines at the
                  beginning.  This says how many.  It skips 'skip'
                  data lines in addition to ignored lines, 'cuz ignored
                  lines are skipped anyway.

    The ``histo()`` method does a simple 1-D histogram, putting the data
    into bins.  It would be usual to not supply a ``binSize``, in which
    case a suitable binSize is chosen for you.  But if you don't like
    the result, you can optionally specify a binSize,
    (eg myNumbers.binSize=10.0)

    If your numbers are in more than one list or file, you can read in
    more stuff with the ``read()`` method.

    """
    
    def __init__(self, inThing, col=0, skip=0):
        self.data = []
        self.bins = None
        # self.binSize is a property
        self._binSize = None
        self.nBins = None
        #self.min = None
        #self.max = None
        #self.range = None
        #self.col
        #self.skip
        if inThing:
            self.read(inThing, col, skip)

    def _setBinSize(self, binSize):
        try:
            theBinSize = float(binSize)
        except (ValueError,TypeError):
            raise Glitch, "Arg binSize, if set, should be a float."
        if theBinSize <= 0.0:
            raise Glitch, "Arg binSize, if set, should be a positive float."
        self._binSize = theBinSize
            
    def _delBinSize(self):
        self._binSize = None
        
    binSize = property(lambda self: self._binSize, _setBinSize, _delBinSize)

    def read(self, inThing, col=0, skip=0):
        """Slurp in some more numbers.

        You can use this repeatedly, eg if your numbers are in more
        than one list or file.  """

        gm = ['Numbers.read()']
        if type(inThing) == type('string') and os.path.isfile(inThing):
            # col and skip only come into play if its a file.
            try:
                self.col = int(col)
                self.skip = int(skip)
            except (ValueError, TypeError):
                gm.append("Args col and skip must be ints")
                raise Glitch, gm

            flob = file(inThing)
            theLines = flob.readlines()
            flob.close()
            #if skip:
            #    if len(theLines) <= self.skip:
            #        gm.append("File '%s' has %i lines, " % (inThing, len(theLines)))
            #        gm.append("but skip is set to %i." % self.skip)
            #        raise Glitch, gm
            skipsDone = 0
            digitsPlusMinus = '0123456789+-'
            for aLine in theLines:
                ll = aLine.lstrip()
                if ll.startswith("#"):
                    pass
                elif not ll:
                    pass
                # Read MrBayes *.p files
                #elif ll.startswith("[ID: "):
                #    pass
                #elif ll.startswith("Gen"):
                #    pass
                elif ll[0] not in digitsPlusMinus:
                    pass
                else:
                    if skipsDone < self.skip:
                        skipsDone += 1
                    else:
                        splitLine = aLine.split()
                        try:
                            theOne = splitLine[col]
                        except IndexError:
                            gm.append("Line '%s'.  " % string.rstrip(aLine))
                            gm.append("Can't get the item at (zero-based) index %i  " % col) 
                            raise Glitch, gm
                        try:
                            aFloat = float(theOne)
                            self.data.append(aFloat)
                        except (ValueError, TypeError):
                            gm.append("Line '%s'.  " % string.rstrip(aLine))
                            gm.append("Can't make sense of '%s'" % theOne)
                            raise Glitch, gm
        elif type(inThing) == type([]):
            for thing in inThing:
                try:
                    aFloat = float(thing)
                    self.data.append(aFloat)
                except (ValueError,TypeError):
                    gm.append("Can't make sense of '%s'" % thing)
                    gm.append("I was expecting a float.")
                    raise Glitch, gm
        else:
            gm.append("Can't understand inThing.  Should be a file or a list.")
            raise Glitch, gm
        #print "got %i data points" % len(self.data)
        #print self.data
        
        self.min = min(self.data)
        self.max = max(self.data)
        self.range = self.max - self.min

    def mean(self):
        if len(self.data):
            return sum(self.data) / len(self.data)
        else:
            return None

    def _makeBins(self, padMin, padMax):

        if not self.binSize:
            # Try to figure out a good binSize
            orderOfMag = 0
            rng = self.range
            orderOfMag = int(math.floor(math.log(rng, 10)))
            rng *= math.pow(10, -orderOfMag)
            #print "rng = %f" % rng
            #print "orderOfMag = %i" % orderOfMag

            # This gives 10-20 bins
            #if rng < 2.0:
            #    step = 0.1
            #elif rng < 4.0:
            #    step = 0.2
            #else:
            #    step = 0.5

            # This gives 9-18 bins (I think)
            if rng < 1.5:
                step = 0.1
            elif rng < 3.6:
                step = 0.2
            elif rng < 9.0:
                step = 0.5
            else:
                step = 1.0
            
            step = step * math.pow(10, orderOfMag)
            #print "step = %f" % step
            self.binSize = step

        # At this point we have a binSize.  Now calculate a
        # niceMin, so that if the bins included zero, then
        # zero would be a bin border.

        if padMin != None:
            theMin = padMin
            #print "theMin = padMin = %f" % padMin
        else:
            theMin = self.min
            #print "theMin = self.min = %f" % self.min
        
        if padMax != None:
            theMax = padMax
            #print "theMax = padMax = %f" % padMax
        else:
            theMax = self.max
            #print "theMax = self.max = %f" % self.max
        
        nStepsToMin = int(math.floor(theMin / self.binSize))
        niceMin = nStepsToMin * self.binSize

        #print "nStepsToMin = %i, niceMin = %f" % (nStepsToMin, niceMin)
        #sys.exit()

        self.nBins = int(math.ceil((theMax - niceMin) / float(self.binSize)))
        #print "self.max = %f, niceMin + (self.nBins * self.binSize) = %f" % (
        #    self.max, niceMin + (self.nBins * self.binSize))
        assert niceMin + (self.nBins * self.binSize) >= theMax
        # If the maximum point lies on a bin border, we need another bin.
        if niceMin + (self.nBins * self.binSize) == theMax:
            #print "...adding another bin."
            self.nBins += 1

        if 0:
            print "binSize= %f" % self.binSize
            print "self.max = %f, self.min = %f, niceMin=%f" % (self.max, self.min, niceMin)
            print "padMin=%s, theMin=%s, padMax=%s, theMax=%s" % (padMin, theMin, padMax, theMax)
            print "self.range = %f" % self.range
            print "theMax - niceMin = %f" % (theMax - niceMin)
            print "self.nBins = %i" % self.nBins
            #sys.exit()

        # Make the bins
        self.bins = []
        binBorders = []
        for i in range(self.nBins):
            bb = i * self.binSize
            self.bins.append([niceMin + bb, 0])
            binBorders.append(niceMin + bb)
        #niceRange = self.nBins * self.binSize
        for i in self.data:
            # This next line did not work for 0.7 on a 0.7 border!
            #bNum = int(math.floor(((i - niceMin) / niceRange) * self.nBins))
            # Nor did this next bit get it right.  Seems to be insoluble.
            bNum = 0
            for j in range(1, self.nBins):
                if i >= binBorders[j]:
                    bNum = j
                else:
                    break
            assert bNum >= 0 and bNum < self.nBins
            self.bins[bNum][1] += 1
        
    def histo(self, verbose=True, binSize=None, padMin=None, padMax=None):
        """Put the data nicely into bins.

        After you do this, the bins are available in self.bins, a list
        of pairs.  

        Args *padMin* and *padMax* extend the range up or down, where
        the extended bins would have zero content, so they would be
        placeholders.  It would be good for making two different
        histos have exactly the same range, eg so you can plot them in
        the same plot.
        """

        
        gm = ['Numbers.histo()']
        if padMin != None:
            assert padMin <= self.min
        if padMax != None:
            assert padMax >= self.max
        if binSize:
            self.binSize = binSize  # a property, so it checks to make sure it is a float
            
        # It is possible that the data are all the same, so its really
        # not clear how to make bins.  Unless there is a binSize defined.
        if self.binSize:
            pass
        else:
            if not self.range:
                self.dump()
                gm.append("The data are all the same.  max=min.  That will not work.")
                raise Glitch, gm

        self._makeBins(padMin, padMax)
        if not self.bins:
            gm.append("No bins.")
            raise Glitch, gm
        if verbose:
            print "%i data points, min=%s, max=%s, binSize=%s, nBins=%i" % (
                len(self.data), self.min, self.max, self.binSize, self.nBins)
            if padMin != None or padMax != None:
                print "padMin=%s, padMax=%s" % (padMin, padMax)
            print "%i points at min, %i points at max" % (self.data.count(self.min), self.data.count(self.max))
            for bin in self.bins:
                print " %-8s  %i" % (bin[0], bin[1]) 


    def dump(self):
        print "%i data points, " % len(self.data),
        print "min=%s, " % self.min,
        print "max=%s, " % self.max,
        print "mean=%s, " % self.mean(),
        print "binSize=%s, " % self.binSize,
        print "nBins=%s" % self.nBins
     

    def plot(self, term='x11'):
        """A horrible hack to plot stuff with GnuPlot.

        A file, gnuplot_instructs, is written.  Then gnuplot is
        called.  This of course needs X-windows.

        Couldn't this be done in a clever way with (existing) proper
        python hooks into GnuPlot?  -yes, but it would require installing
        other modules."""

        weirdName = 'tmPFoRGnUPloT'
        f1 = file(weirdName, 'w')
        for i in self.data:
            f1.write('%s\n' % i)
        f1.close()

        instructionsFileName = 'gNupLot_inStruCts'
        f1 = open(instructionsFileName, 'w')
        f1.write('set term %s\n' % term)   # my new gnuplot on my mac, from macports, has aqua as the default.
        f1.write('plot "%s" notitle\n' % weirdName)
        f1.close()
        
        os.system('gnuplot -persist %s' % instructionsFileName)
        os.system('rm %s' % weirdName)
        os.system('rm %s' % instructionsFileName)

    def tailAreaProbability(self, theStat, verbose=1):
        """Access to :func:`func.tailAreaProbability`

        Supply *theStat*, ie the test quantity, and it uses
        ``self.data`` to make the distribution.
        """
        return func.tailAreaProbability(theStat, self.data, verbose)

    # def mplPlot(self):
    #     """Plot with mathplotlib.  Does this work?"""
    #     from pylab import plot,show
    #     plot(self.data)
    #     show()

    # def harmonicMeanOfLogs2(self):
    #     theMax = -min(self.data)
    #     diff = 600. - theMax
    #     theSum = 0.0
    #     for i in range(len(self.data)):
    #         val = -self.data[i] + diff
    #         theSum += math.exp(val)
    #     return -(math.log(theSum) - diff)

if var.usePfAndNumpy:
    def gsl_meanVariance(self):
        """Uses gsl.  Returns a tuple of the mean and the variance."""
        
        if len(self.data):
            a = numpy.array(self.data, numpy.float)
            m = numpy.zeros([1], numpy.float)
            v = numpy.zeros([1], numpy.float)
            func.gsl_meanVariance(a, m, v)
            return (m, v)
        else:
            return None
    Numbers.gsl_meanVariance = gsl_meanVariance
    del(gsl_meanVariance)

    def arithmeticMeanOfLogs(self):
        ar = numpy.array(self.data, numpy.float)
        #print ar

        if 0:
            numpy.exp(ar, ar)
            print ar
            m = sum(ar) / len(ar)
            return math.log(m)

        # From MrBayes.  Thanks again, guys!
        scaler = max(ar) - 100.0
        a = aOld = n = 0.0
        for x in ar:
            y = x
            y -= scaler
            if y < -100.0:
                print "Numbers.arithmeticMeanOfLogs()  Ignoring outlier  %f" % x
                continue
            else:
                x = numpy.exp(y)
            if n < 0.5:
                a = x
            else:
                aOld = a
                a = aOld + (x -aOld) / (n + 1.0)
            n += 1.0
        mean = numpy.log(a) + scaler
        #print "mean = %f" % mean
        return mean
    Numbers.arithmeticMeanOfLogs = arithmeticMeanOfLogs
    del(arithmeticMeanOfLogs)
    
        
    def harmonicMeanOfLogs(self):
        """Returns log of the harmonic mean of logs.

        Self is in log form.  The harmonic mean is calculated after
        exponentiating, and the result is returned as a log.

        The method was stolen uncritically from MrBayes.  This gives
        the same result as MrBayes.
        """

        # The numbers are scaled relative to the lowest log value.  If
        # the highest log values are more than 200 log units more than
        # the lowest log value, then they are ignored.
        
        ar = numpy.array(self.data, numpy.float)
        #print ar

        # From MrBayes.
        scaler = (0.0 - min(ar)) - 100.0
        #print "scaler = %f" % scaler
        
        a = aOld = n = 0.0
        for x in ar:
            x *= -1.0
            y = x
            y -= scaler
            #print "x=%f, y = %f" % (x, y)
            if y < -100.0:
                print "Numbers.harmonicMeanOfLogs()  Ignoring outlier %f" % -x
                continue
            else:
                x = numpy.exp(y)
            if n < 0.5:
                a = x
            else:
                aOld = a
                a = aOld + (x - aOld) / (n + 1.0)
            n += 1.0
        harm_mean = - numpy.log(a) - scaler
        return harm_mean

    Numbers.harmonicMeanOfLogs = harmonicMeanOfLogs
    del(harmonicMeanOfLogs)




