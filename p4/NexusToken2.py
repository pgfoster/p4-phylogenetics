import pf
from Var import var
import numpy,string
from Glitch import Glitch

"""A faster version of nextTok(), using memory allocated (once only)
using numpy, and using functions written in C.  The slow, pure
python module is NexusToken.py.  This version is about twice as fast.

Which one is used is under the control of var.nexus_doFastNextTok.
This one does not work for CStrings, so we need to revert to the old
way whenever CStrings are encountered."""


class NexusToken(object):
    def __init__(self, max):
        self.max = numpy.array([max], numpy.int32)
	self.tokLen = numpy.array([0], numpy.int32)
	self.tok = numpy.array(['x'] * int(self.max), 'c')
	self.embeddedCommentLen = numpy.array([0], numpy.int32)
	self.embeddedComment = numpy.array(['x'] * int(self.max), 'c')
	self.savedCommentLen = numpy.array([0], numpy.int32)
	self.filePtr = None
        self.nexusToken = pf.newNexusToken(var._nexus_writeVisibleComments,
                                           var._nexus_getP4CommandComments,
                                           var._nexus_getWeightCommandComments,
                                           var._nexus_getAllCommandComments,
                                           var._nexus_getLineEndingsAsTokens,
                                           self.max,
                                           self.tokLen,
                                           self.tok,
                                           self.embeddedCommentLen,
                                           self.embeddedComment,
                                           self.savedCommentLen)
        #self.previousTok = None
        #self.previousEmbeddedComment = None


nt = NexusToken(300)

def checkLineLengths(flob):
    global nt
    #print 'NexusToken2.checkLineLengths here.'
    flob.seek(0,0)
    longest = pf.nexusTokenCheckLineLengths(nt.nexusToken, flob)
    flob.seek(0,0)
    #print 'The longest line length is %i' % longest
    if longest > nt.max:
        nt = NexusToken(longest)

def nextTok(flob):
    #print 'NexusToken2.nextTok() here.  nt.nexusToken = %i, max=%s, tokLen=%s, type(tokLen)=%s' % (nt.nexusToken, nt.max, nt.tokLen[0], type(nt.tokLen))
    #assert type(nt.tokLen) == type(numpy.array([0], numpy.int32))
    #print "NexusToken2.nextTok().  nt.wordIsFinished[0]=%i, nt.tokLen=%i, previousTok=%s, previousComment=%s" % (nt.wordIsFinished[0], nt.tokLen[0], nt.previousTok, nt.previousEmbeddedComment) 
    #if nt.wordIsFinished[0]:
    #    assert nt.tokLen[0]
    #    ret = nt.tok[:int(nt.tokLen[0])].tostring()
    #    nt.tokLen[0] = 0
    #    nt.wordIsFinished[0] = 0
    #    #nt.previousTok = ret
    #    return ret
    #print '    x1 NexusToken2.nextTok() here. savedCommentLen=%i' % nt.savedCommentLen[0]
    if nt.savedCommentLen[0]:
        ret = nt.embeddedComment[:int(nt.savedCommentLen[0])].tostring()
        nt.savedCommentLen[0] = 0
        return ret
    pf.nextToken(nt.nexusToken, flob)
    #print '    x2 tokLen = %i, embeddedCommentLen[0] = %i' % (nt.tokLen[0], nt.embeddedCommentLen[0])
    if nt.embeddedCommentLen[0]:
        ret = nt.embeddedComment[:int(nt.embeddedCommentLen[0])].tostring()
        nt.embeddedCommentLen[0] = 0
        #nt.previousEmbeddedComment = ret
        return ret
    else:
        if nt.tokLen[0]:
            ret = nt.tok[:int(nt.tokLen[0])].tostring()
            nt.tokLen[0] = 0
            #nt.previousTok = ret
            return ret
        else:
            return None
    
def safeNextTok(flob, caller=None):
    t = nextTok(flob)
    if not t:
        if caller:
            gm = ["safeNextTok(), called from %s" % caller]
        else:
            gm = ["safeNextTok()"]
        gm.append("Premature Death.")
        gm.append("Ran out of understandable things to read in nexus file.")
        raise Glitch, gm
    else:
        return t

def nexusSkipPastNextSemiColon(flob):
    pf.nexusSkipPastNextSemiColon(nt.nexusToken, flob)


def nexusSkipPastBlockEnd(flob):
    """Read up to and including a block 'end' or 'endblock'."""
    # This should only ever be issued after a semi-colon

    complaintHead = '\nNexus: nexusSkipPastBlockEnd()'
    if hasattr(flob, 'name'):
        complaintHead += " file: %s" % flob.name
    while 1:
        tok = nextTok(flob)
        if tok:
            lowTok = string.lower(tok)
            if lowTok == 'end' or lowTok == 'endblock':
                tok2 = nextTok(flob)
                if not tok2 or tok2 != ';':
                    gm = [complaintHead]
                    gm.append("    Expecting a semicolon after %s" % tok)
                    if not tok2:
                        gm.append("Got nothing.")
                    else:
                        gm.append("Got '%s'" % tok2)
                    raise Glitch, gm
                return
            elif lowTok == ';':  # for pathological cases where the last command is a ';' by itself.
                continue
            else:
                pf.nexusSkipPastNextSemiColon(nt.nexusToken, flob)
        else:
            break
    gm = [complaintHead]
    gm.append("Failed to find either 'end' or 'endblock'")
    gm.append("Premature end of file?")
    raise Glitch, gm
