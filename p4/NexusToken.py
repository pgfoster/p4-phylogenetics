import os,string,sys,types
from Var import var
from Glitch import Glitch


############################################################################
#
# Token generation stuff.
#
#       safeNextTok()  # checks for None, and dies
#       nextTok()      # may return None
#
# Handling comments:
#   - Behaviour is under the control of some variables in the var module.
#   - Visible comments are never returned.  They may be printed, depending
#     on variable settings.
#   - Command comments may be returned (depending on the variable settings).
#     While it is true that comments are not tokens, the Nexus format
#     specifies informative command comments, and so they are returned by
#     nextTok() and safeNextTok().
#
# Some differences with the old nextToken() stuff:
#   - comments within quoted words are not skipped,
#   - comments are returned in one piece
#
# Below are some possibly useful functions:
#    nexusSkipPastNextSemiColon()
#    nexusSkipPastBlockEnd()
#    nexusNextCommand()        <-- does not appear to be used anywhere
#
#
############################################################################

pieces = []  # This is used in _getWord(), where it is global.
             # It needs to survive multiple calls to nextTok() and _getWord()
#wordIsFinished = False
comment = None

##Ignore
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

def nextTok(flob):
    """Returns Nexus tokens or (some) comments.

    You feed this function a flob, an open file or file-like object.

    Comments are handled under the control of a few P4 variables,
    namely

        p4.var.nexus_writeVisibleComments      # write, but do not get, all [!...]
        p4.var.nexus_getP4CommandComments      # all [&&p4 ...]
        p4.var.nexus_getWeightCommandComments  # all [&w ...]
        p4.var.nexus_getAllCommandComments     # all [&...]

    """

    #global pieces
    #global wordIsFinished
    #if not flob:
    #    if 1:
    #        print '\nNexusToken.nextTok()'
    #        print "    No flob?"
    #        import sys; sys.exit()
    #    else:
    #        return None

    # It may be that the last tok returned was a comment, but we
    # know that there is a word in pieces that is finished.  So
    # just return it.
    #if wordIsFinished:
    #    assert pieces
    #    ret = ''.join(pieces)
    #    wordIsFinished = False
    #    pieces = []
    #    return ret

    global comment

    if comment:
        ret = comment
        comment = None
        return ret

    while 1:

        
        c = flob.read(1)
        if 0:
            if c == '\n' or c == '\r':
                print 'nt %3i  c = line ending' % flob.tell()
            else:
                print 'nt %3i  c = %s' % (flob.tell(), c)
        if not c:
            return None
        #if 0:
        #    if c in string.whitespace:
        #        print " c %3i  whitespace" % flob.tell()
        #    else:
        #        print " c %3i  %s" % (flob.tell(), c)
        if var.nexus_getLineEndingsAsTokens:
            if c == '\n' or c == '\r':
                return c
        if c in string.whitespace:
            continue
        elif c == '[':
            ret = _handleComment(flob)
            if ret:
                return ret
            else:
                continue
        elif c == "'":
            return _getQuotedStuff(flob)
        #elif c in var.nexus_punctuation:
        elif c in var.punctuation:
            return c
        else:
            return _getWord(flob)




def _handleComment(flob):
    """Do the right thing with [a nexus comment].

    The only comments we are interested in are the tree weights, eg
    [&W 1/2] and command comments, eg [& foo].  Visible comments are
    not returned, but if var.nexus_writeVisibleComments is turned on
    they are written to stdout."""

    #   nexus_writeVisibleComments      # write, but do not get, all [!...]
    #   nexus_getP4CommandComments      # all [&&p4 ...]
    #   nexus_getWeightCommandComments  # all [&w ...]
    #   nexus_getAllCommandComments     # all [&...]


    # We are here having read the opening '[' character.
    commentStartPos = flob.tell() - 1

    c = flob.read(1)
    if not c: # empty string, ie at the end of the flob.  level can be assumed to be more than zero.
        gm = ["NexusToken._skipComment()"]
        gm.append("Reached the end while still in a comment.")
        raise Glitch, gm
    if 0:
        if c in string.whitespace:
            print "hc %3i whitespace"  % flob.tell()
        else:
            print "hc %3i  %s" % (flob.tell(), c)
    if c == ']':
        return None
    if c == '!' and var.nexus_writeVisibleComments:
        flob.seek(commentStartPos, 0)
        theComment = _getComment(flob)
        print theComment
        return None
    elif c in ['&']: #, '\\']:
        if var.nexus_getAllCommandComments:
            flob.seek(commentStartPos, 0)
            return  _getComment(flob)
        c2 = flob.read(1)
        if not c2:
            gm = ["NexusToken._skipComment()"]
            gm.append("Reached the end while still in a comment.")
            raise Glitch, gm
        if c == '&' and string.lower(c2) == 'w' and var.nexus_getWeightCommandComments:
            flob.seek(commentStartPos, 0)
            return _getComment(flob)
        elif c2 == '&' and var.nexus_getP4CommandComments:
            # ask whether the next 3 characters are 'p4 ' (ie p, 4, space)
            c2 = flob.read(3)
            if len(c2) != 3 or string.lower(c2) != 'p4 ':
                flob.seek(commentStartPos + 1, 0)
                _skipComment(flob)
                return
            flob.seek(commentStartPos, 0)
            return _getComment(flob)
        else:
            flob.seek(commentStartPos + 1, 0)
            return _skipComment(flob)
    else:
        _skipComment(flob)
        return


def _skipComment(flob):
    level = 1 # ie level of nested comments.  This assumes that we have already read one '['.
    while 1:
        c = flob.read(1)
        if 0:
            if c in string.whitespace:
                print "sc %3i  whitespace" % flob.tell()
            else:
                print "sc %3i  %s" % (flob.tell(), c)
        if not c: # empty string, ie at the end of the flob.  level can be assumed to be more than zero.
            gm = ["NexusToken._skipComment()"]
            gm.append("Reached the end while still in a comment.")
            raise Glitch, gm
        if c == '[':
            level = level + 1
        elif c == ']':
            level = level - 1
        if level == 0:
            return


def _getComment(flob):
    startPos = flob.tell()
    level = 0 # ie level of nested comments.  This assumes that we have not already read one '['.
    while 1:
        c = flob.read(1)
        if 0:
            if c in string.whitespace:
                print "gc %3i  whitespace" % flob.tell()
            else:
                print "gc %3i  %s" % (flob.tell(), c)
        if not c:
            gm = ["NexusToken._getComment()"]
            gm.append("Reached the end while still in a comment.")
            raise Glitch, gm
        if c == '[':
            level = level + 1
        elif c == ']':
            level = level - 1
        if level == 0:
            endPos = flob.tell()
            flob.seek(startPos)
            theLen = endPos - startPos
            return flob.read(theLen)


def _getWord(flob):
    """Get a nexus token that isnt punctuation or a comment.

    The word might be broken by comments (the authors of the Nexus
    format were mad).  Internal command comments (eg mid[&a
    comment]dle, or aWord[&a comment at the end]) are returned before
    the enclosing word, because the comment ends before the enclosing
    word."""

    #pieces = [] # It might be broken by comments.
    global pieces
    #global wordIsFinished
    global comment

    startPos = flob.tell() - 1
    while 1:
        c = flob.read(1)
        if 0:
            if not c:
                print "gw %3i  empty (position given is that of the last char)" % flob.tell()
            if c in string.whitespace:
                print "gw %3i  whitespace" % flob.tell()
            else:
                print "gw %3i  %s" % (flob.tell(), c)
        if c == '[':
            endPos = flob.tell() - 1
            theLen = endPos - startPos
            flob.seek(startPos)
            theWordPiece = flob.read(theLen)
            #print "gw got theWordPiece = %s" % theWordPiece
            pieces.append(theWordPiece)
            commentStartPos = flob.tell()
            flob.seek(1, 1) # skip the '['
            #_skipComment(flob)
            ret = _handleComment(flob)
            startPos = flob.tell()
            #print 'gw %3i  comment=%s' % (flob.tell(), ret)

            if ret:
                # Now we need to ask if the word is finished, or whether it continues immediately after the comment.
                wordIsFinished = False
                c2 = flob.read(1)
                if c2 in string.whitespace or c2 in var.punctuation or not c2:
                    wordIsFinished = True
                if c2:
                    flob.seek(-1, 1) # back up one space only if we went forward one space

                if wordIsFinished:
                    comment = ret
                else:
                    return ret
            else:
                continue

            #if 0:
            #    # If it is the end of a word, like this:
            #    # aWord[aComment], then maybe we do not want to skip
            #    # it.  That is why we saved the commentStartPos.  Test
            #    # for whether it was at the end of a word.
            #    c2 = flob.read(1)
            #    if c2 in string.whitespace or c2 in var.nexus_punctuation or not c2:
            #        # It was at the end of a word, so back up and break
            #        flob.seek(commentStartPos, 0)
            #        break
            #    else:
            #        flob.seek(-1,1) # un-do the move of the previous read(1)
            #startPos = flob.tell()
        #elif c in string.whitespace or c in var.nexus_punctuation or not c:
        elif c in string.whitespace or c in var.punctuation or not c:
            if c:
                flob.seek(-1, 1) # back up one space only if we went forward one space
            endPos = flob.tell()
            theLen = endPos - startPos
            if theLen:
                flob.seek(startPos)
                theWordPiece = flob.read(theLen)
                #print "gw got theWordPiece = %s" % theWordPiece
                pieces.append(theWordPiece)
            break
    theWord = string.join(pieces, '')
    pieces = []
    return theWord


def _getQuotedStuff(flob):
    complaintHead = '\nNexusToken._getQuotedStuff()'
    local_pieces = []
    startPos = flob.tell() # we have passed the opening single quote

    # Check that we do not have a single quote immediately following.
    # 2 single quotes is ok, tho.
    c = flob.read(1)
    if c:
        if c == '\'':
            c2 = flob.read(1)
            if c2:
                if c2 == '\'':
                    # The opening single quote has been directly
                    # followed by 2 single quotes. Thats ok.
                    flob.seek(-2, 1) # Back up 2 spaces.
                else:
                    # The opening single quote was followed by a
                    # single quote and then something else.  Thats
                    # bad.
                    gm = [complaintHead]
                    gm.append("Got 2 single quotes in a row, not properly within single quotes.")
                    raise Glitch, gm
            else:
                gm = [complaintHead]
                gm.append("Got 2 single quotes in a row, not properly within single quotes.")
                gm.append("And then the file ended.  Bad.")
                raise Glitch, gm
        else:
            flob.seek(-1, 1) # it wasn't a single quote, so back up one space
    else:
        gm = [complaintHead]
        gm.append("    File ended with an un-matched single quote.")
        raise Glitch, gm

    #print 'At the start of the while(1) loop, the file position is %i' % flob.tell()
    while 1:
        c = flob.read(1)
        #print '! c= %s' % c
        if 0:
            if not c:
                print "gq %3i  empty (position is eof)" % flob.tell()
            elif c in string.whitespace:
                print "gq %3i  whitespace" % flob.tell()
            else:
                print "gq %3i  %s" % (flob.tell(), c)
        if c == '\'':
            c2 = flob.read(1)
            if c2:
                flob.seek(-1, 1)
            endPos = flob.tell() - 1
            theLen = endPos - startPos
            flob.seek(startPos)
            if theLen:
                thePiece = flob.read(theLen)
                #print "gq got thePiece = %s" % thePiece
                local_pieces.append(thePiece)

            if c2 and c2 == '\'': # ie 2 single quotes in a row
                local_pieces.append('\'\'')
                flob.seek(2, 1)
                startPos = flob.tell()
                continue
            else:
                flob.seek(1,1)
                break

        elif not c:
            gm = [complaintHead]
            gm.append("File ended while still in a quoted word.")
            raise Glitch, gm
    return string.join(['\''] + local_pieces + ['\''], '')

# (Thats it for the token generation stuff)


##Ignore
def nexusSkipPastNextSemiColon(flob):
    complaintHead = '\nNexus.nexusSkipPastNextSemiColon()'
    while 1:
        c = flob.read(1)
        #print "nexusSkipPastNextSemiColon() c=%s" % c
        if c:
            if c == '[':
                _skipComment(flob)
            elif c == ';':
                return
        else:
            break
    gm = [complaintHead]
    gm.append("Every Nexus command must end in a semicolon.")
    gm.append("No semicolon was found.  Premature end of file?")
    raise Glitch, gm

##Ignore
def nexusSkipPastBlockEnd(flob):
    """Read up to and including a block 'end' or 'endblock'."""
    # This should only ever be issued after a semi-colon
    gm = ['NexusToken.nexusSkipPastBlockEnd()']
    while 1:
        tok = nextTok(flob)
        #print "nexusSkipPastBlockEnd() tok=%s" % tok
        if tok:
            lowTok = string.lower(tok)
            if lowTok == 'end' or lowTok == 'endblock':
                tok2 = nextTok(flob)
                if not tok2 or tok2 != ';':
                    print complaintHead
                    gm.append("Expecting a semicolon after %s" % tok)
                    if not tok2:
                        gm.append("Got nothing.")
                    else:
                        gm.append("Got '%s'" % tok2)
                    raise Glitch, gm
                return
            elif lowTok == ';':  # for pathological cases where the last command is a ';' by itself.
                continue
            else:
                nexusSkipPastNextSemiColon(flob)
        else:
            break
    gm.append("Failed to find either 'end' or 'endblock'")
    gm.append("Premature end of file?")
    raise Glitch, gm


def nexusNextCommand(flob):
    """Return a nexus command as a list of tokens.

    Does not return the final semi-colon."""

    toks = []
    tok = nextTok(flob)
    if not tok: # end of file, not an error
        return None
    if string.lower(tok) == '#nexus':
        tok = nextTok(flob)
    while tok:
        if tok == ';':
            return toks
        toks.append(tok)
        tok = nextTok(flob)
    gm = ["NexusToken.nexusNextCommand()"]
    gm.append("Could not find a semi-colon to end the nexus command.")
    gm.append("Premature end of file?")
    raise Glitch, gm


