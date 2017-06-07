# Read a pickled mcmc tunings (from autoTune()) and dump the info.
try:
    fName = var.argvAfterDoubleDash[0]
except IndexError:
    print("Tell me the tunings file name after a double dash, eg")
    print("sTuningsDump.py -- tunings_0.pickle")
    sys.exit()

tf = open(fName, 'rb')
tunings = pickle.load(tf)
tf.close()

tunings.dump()
