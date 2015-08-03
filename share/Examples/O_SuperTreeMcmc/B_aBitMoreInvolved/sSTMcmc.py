#os.system("rm -f mcmc*")
# (self, inTrees, modelName='SR2008_rf_aZ', beta=1.0, stRFCalc='purePython1', runNum=0, sampleInterval=100, checkPointInterval=None)

# Choose one of these, the fastest available
myCalc='purePython1'
myCalc='bitarray'
#myCalc='fastReducedRF'

read('inTrees.phy')

stm = STMcmc(var.trees, modelName='SR2008_rf_aZ', beta=1.2, stRFCalc=myCalc, runNum=0, sampleInterval=20, checkPointInterval=10000)
stm.run(20000)
stm = STMcmc(var.trees, modelName='SR2008_rf_aZ', beta=1.2, stRFCalc=myCalc, runNum=1, sampleInterval=20, checkPointInterval=10000)
stm.run(20000)
