from p4.supertreesupport import SuperTreeSupport

sts = SuperTreeSupport('supertree.nex', 'input.nex')
sts.doSaveDecoratedTree = False
sts.decoratedFilename='mytree.nex'
sts.verbose=2
sts.doDrawTree=True
sts.superTreeSupport()
