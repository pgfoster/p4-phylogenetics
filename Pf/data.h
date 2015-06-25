data *newData(int nTax, int nParts);
void freeData(data *theData);
void dumpData(data *theData);
PyObject *rell(int nBoots, rellStuff *rStuff);
void bootstrapData(data *reference, data *toFill, const gsl_rng *gsl_rng);
