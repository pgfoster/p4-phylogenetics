from __future__ import print_function
import os
import string
import sys
from p4.var import var
from p4.p4exceptions import P4Error


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


def safeNextTok(flob, caller=None):
    t = nextTok(flob)
    if not t:
        if caller:
            gm = ["safeNextTok(), called from %s" % caller]
        else:
            gm = ["safeNextTok()"]
        gm.append("Premature Death.")
        gm.append("Ran out of understandable things to read in nexus file.")
        raise P4Error(gm)
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


    global comment

    if comment:
        ret = comment
        comment = None
        return ret

    while 1:
        c = flob.read(1)
        if 0:
            if c == '\n' or c == '\r':
                print('nt %3i  c = line ending' % flob.tell())
            else:
                print('nt %3i  c = %s' % (flob.tell(), c))
        if not c:
            return None
        # if 0:
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
            ret = _handleComment(flob, c)
            if ret:
                return ret
            else:
                continue
        elif c == "'":
            return _getQuotedStuff(flob)
        elif c in var.punctuation:
            return c
        else:
            # print("nt %3i  %s" % (flob.tell(), c))
            return _getWord(flob, c)


def _handleComment(flob, c0):
    """Do the right thing with [a nexus comment].

    The only comments we are interested in are the tree weights, eg
    [&W 1/2] and command comments, eg [& foo].  Visible comments are
    not returned, but if var.nexus_writeVisibleComments is turned on
    they are written to stdout."""

    #   nexus_writeVisibleComments      # write, but do not get, all [!...]
    #   nexus_getP4CommandComments      # all [&&p4 ...]
    #   nexus_getWeightCommandComments  # all [&w ...]
    #   nexus_getAllCommandComments     # all [&...]

    gm = ["nexustoken._handleComment()"]
    localPieces = [c0]
    # presumably c0 is '['

    c = flob.read(1)
    localPieces.append(c)
    if not c:
        # empty string, ie at the end of the flob.  level can be assumed to be
        # more than zero.
        gm.append("Reached the end while still in a comment. No comment close.")
        raise P4Error(gm)
    if c == ']':
        return None   # appears to be '[]'
    if c == '!' and var.nexus_writeVisibleComments:
        theComment = _getComment(flob, localPieces)
        print(theComment)
        return None
    elif c in ['&']:  # , '\\']:
        if var.nexus_getAllCommandComments:
            return _getComment(flob, localPieces)
        c2 = flob.read(1)
        if not c2:
            gm.append("Reached the end while still in a comment.")
            raise P4Error(gm)
        localPieces.append(c2)
        if c == '&' and c2.lower() == 'w' and var.nexus_getWeightCommandComments:
            return _getComment(flob, localPieces)
        elif c2 == '&' and var.nexus_getP4CommandComments:
            # We may need to get back here, and relative seek will not work.
            startSpot = flob.tell()
            # ask whether the next 3 characters are 'p4 ' (ie p, 4, space)
            c3 = flob.read(3)
            if len(c3) != 3 or c3.lower() != 'p4 ':
                flob.seek(startSpot)
                _skipComment(flob)
                return
            # So it is a p4 comment.  Get it.
            localPieces.append(c3)
            return _getComment(flob, localPieces)
        else:
            return _skipComment(flob)  # None, normally
    else:
        _skipComment(flob)
        return


def _skipComment(flob):
    # ie level of nested comments.  This assumes that we have already read one
    # '['.
    level = 1
    while 1:
        c = flob.read(1)
        if 0:
            if c in string.whitespace:
                print("sc %3i  whitespace" % flob.tell())
            else:
                print("sc %3i  %s" % (flob.tell(), c))
        # empty string, ie at the end of the flob.  level can be assumed to be
        # more than zero.
        if not c:
            gm = ["nexustoken._skipComment()"]
            gm.append("Reached the end while still in a comment.")
            raise P4Error(gm)
        if c == '[':
            level = level + 1
        elif c == ']':
            level = level - 1
        if level == 0:
            return


def _getComment(flob, localPieces):
    # level, ie level of nested comments.  This assumes that we have already read one '['.
    level = 1
    while 1:
        c = flob.read(1)
        localPieces.append(c)
        if 0:
            if c in string.whitespace:
                print("gc %3i  whitespace" % flob.tell())
            else:
                print("gc %3i  %s" % (flob.tell(), c))
        if not c:
            gm = ["nexustoken._getComment()"]
            gm.append("Reached the end while still in a comment.")
            raise P4Error(gm)
        if c == '[':
            level = level + 1
        elif c == ']':
            level = level - 1
        if level == 0:
            return ''.join(localPieces)


def _getWord(flob, c0):
    """Get a nexus token that isnt punctuation or a comment.

    The word might be broken by comments (the authors of the Nexus
    format were mad).  Internal command comments (eg mid[&a
    comment]dle, or aWord[&a comment at the end]) are returned before
    the enclosing word, because the comment ends before the enclosing
    word."""

    global pieces
    #global wordIsFinished
    global comment

    localPieces = [c0]
    while 1:
        prevPosn = flob.tell()
        c = flob.read(1)
        if 0:
            if not c:
                print("gw %3i  empty (position given is that of the last char)" % flob.tell())
            if c in string.whitespace:
                print("gw %3i  whitespace" % flob.tell())
            else:
                print("gw %3i  %s" % (flob.tell(), c))
        if c == '[':
            theWordPiece = ''.join(localPieces)
            # print("gw comm got theWordPiece '%s'" % theWordPiece) 
            pieces.append(theWordPiece)
            localPieces = []
            ret = _handleComment(flob, c)
            # print('gw %3i  comment=%s' % (flob.tell(), ret))

            # Here it might be useful to have wordIsFinished as a global, but here it is local.

            if ret:
                # Now we need to ask if the word is finished, or whether it
                # continues immediately after the comment.
                beforeC2Posn = flob.tell()
                # print("beforC2Posn is %i" % beforeC2Posn)
                wordIsFinished = False
                c2 = flob.read(1)
                if 0:
                    if c2 in string.whitespace:
                        print("gw %3i  c2=whitespace" % flob.tell())
                    else:
                        print("gw %3i  c2=%s" % (flob.tell(), c2))

                if c2 in string.whitespace or c2 in var.punctuation or not c2:
                    wordIsFinished = True
                if c2:
                    # back up one space only if we need to
                    if c2 in string.whitespace:
                        if var.nexus_getLineEndingsAsTokens:
                            if c2 == '\n' or c2 == '\r':
                                #print("seeking %i" % beforeC2Posn)
                                flob.seek(beforeC2Posn)
                    elif c in var.punctuation:
                        #print("seeking %i" % beforeC2Posn)
                        flob.seek(beforeC2Posn)

                if wordIsFinished:
                    # The word (ie pieces) will be returned first, then comment will be returned on the following nextTok()
                    comment = ret
                else:
                    # There is more word to come, so return the comment now.
                    return ret
            else:
                # We threw the comment away.  Look for more word.
                continue

        elif c in string.whitespace or c in var.punctuation or not c:
            if c:
                # We need to back up one (unicode) char, sometimes
                if c in string.whitespace:
                    if var.nexus_getLineEndingsAsTokens:
                        if c == '\n' or c == '\r':
                            #print("seeking %i" % prevPosn)
                            flob.seek(prevPosn)
                elif c in var.punctuation:
                    #print("seeking %i" % prevPosn)
                    flob.seek(prevPosn)
            theWordPiece = ''.join(localPieces)
            # print("gw got theWordPiece '%s'" % theWordPiece) 
            pieces.append(theWordPiece)
            break
        else:
            localPieces.append(c)
    theWord = ''.join(pieces)
    pieces = []
    return theWord


def _getQuotedStuff(flob):
    complaintHead = 'nexustoken._getQuotedStuff()'
    local_pieces = ["'"]

    # Check that we do not have a single quote immediately following.
    # 2 single quotes is ok, tho.
    c = flob.read(1)
    if c:
        local_pieces.append(c)
        if c == "'":
            c2 = flob.read(1)
            if c2:
                local_pieces.append(c2)
                if c2 == "'":
                    # The opening single quote has been directly
                    # followed by 2 single quotes. Thats ok.
                    pass
                else:
                    # The opening single quote was followed by a
                    # single quote and then something else.  Thats
                    # bad.
                    gm = [complaintHead]
                    gm.append("Got 2 single quotes in a row, not properly within single quotes.")
                    raise P4Error(gm)
            else:
                gm = [complaintHead]
                gm.append(
                    "Got 2 single quotes in a row, not properly within single quotes.")
                gm.append("And then the file ended.  Bad.")
                raise P4Error(gm)
        else:
            pass
    else:
        gm = [complaintHead]
        gm.append("File ended with an un-matched single quote.")
        raise P4Error(gm)

    # print 'At the start of the while(1) loop, the file position is %i' % flob.tell()
    while 1:
        c = flob.read(1)
        # print '! c= %s' % c
        if 0:
            if not c:
                print("gq %3i  empty (position is eof)" % flob.tell())
            elif c in string.whitespace:
                print("gq %3i  whitespace" % flob.tell())
            else:
                print("gq %3i  %s" % (flob.tell(), c))
        if c:
            local_pieces.append(c)
            if c == "'":
                posnBeforeC2 = flob.tell()
                c2 = flob.read(1)
                if c2:
                    if c2 == "'":  # ie 2 single quotes in a row
                        local_pieces.append(c2)
                        continue
                    else:
                        # It was a single quote, not followed by a single quote,
                        # so we are finished.  But in finding out that fact, we
                        # have gone too far.  So back up.
                        flob.seek(posnBeforeC2)
                        break
                    
                else:
                    # It was a closing single quote, followed by nothing.
                    break

        elif not c:
            gm = [complaintHead]
            gm.append("File ended while still in a quoted word.")
            raise P4Error(gm)
    return ''.join(local_pieces)

# (Thats it for the token generation stuff)


def nexusSkipPastNextSemiColon(flob):
    complaintHead = 'nexustoken.nexusSkipPastNextSemiColon()'
    while 1:
        c = flob.read(1)
        # print "nexusSkipPastNextSemiColon() c=%s" % c
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
    raise P4Error(gm)


def nexusSkipPastBlockEnd(flob):
    """Read up to and including a block 'end' or 'endblock'."""
    # This should only ever be issued after a semi-colon
    gm = ['nexustoken.nexusSkipPastBlockEnd()']
    while 1:
        tok = nextTok(flob)
        # print "nexusSkipPastBlockEnd() tok=%s" % tok
        if tok:
            lowTok = tok.lower()
            if lowTok == 'end' or lowTok == 'endblock':
                tok2 = nextTok(flob)
                if not tok2 or tok2 != ';':
                    print(complaintHead)
                    gm.append("Expecting a semicolon after %s" % tok)
                    if not tok2:
                        gm.append("Got nothing.")
                    else:
                        gm.append("Got '%s'" % tok2)
                    raise P4Error(gm)
                return
            # for pathological cases where the last command is a ';' by itself.
            elif lowTok == ';':
                continue
            else:
                nexusSkipPastNextSemiColon(flob)
        else:
            break
    gm.append("Failed to find either 'end' or 'endblock'")
    gm.append("Premature end of file?")
    raise P4Error(gm)


def nexusNextCommand(flob):                           # does not appear to be used anywhere
    """Return a nexus command as a list of tokens.

    Does not return the final semi-colon."""

    toks = []
    tok = nextTok(flob)
    if not tok:  # end of file, not an error
        return None
    if tok.lower() == '#nexus':
        tok = nextTok(flob)
    while tok:
        if tok == ';':
            return toks
        toks.append(tok)
        tok = nextTok(flob)
    gm = ["nexustoken.nexusNextCommand()"]
    gm.append("Could not find a semi-colon to end the nexus command.")
    gm.append("Premature end of file?")
    raise P4Error(gm)
