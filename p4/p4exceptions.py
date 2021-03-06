import textwrap

class P4Error(Exception):

    """A class for exceptions in p4.

    You can raise this with a string, or a list of strings.  If its
    a single string, it gets wrapped.  If its a list of 2 strings, the
    first one is output flush and unwrapped, and the second is
    indented and wrapped."""

    def __init__(self, msg='', flavour=''):
        myIndent = ' ' * 4
        if isinstance(msg, str):
            try:
                if msg.startswith('\n\n'):
                    firstLine = '%s' % msg
                elif msg.startswith('\n'):
                    firstLine = '\n%s' % msg
                else:
                    firstLine = '\n\n%s' % msg
            except:
                firstLine = ''
            self.msg = firstLine

        elif isinstance(msg, list):
            try:
                if msg[0].startswith('\n\n'):
                    firstLine = '%s' % msg[0]
                elif msg[0].startswith('\n'):
                    firstLine = '\n%s' % msg[0]
                else:
                    firstLine = '\n\n%s' % msg[0]
            except:
                firstLine = ''
            niceMsgList = [firstLine]
            for i in range(len(msg))[1:]:
                if isinstance(msg[i], str):
                    #  If it is short, use it as is.  If it is long, wrap it.
                    if len(msg[i]) < 66:
                        niceMsgList.append(myIndent + msg[i])
                    else:
                        wLine = textwrap.fill(
                            msg[i], 70, initial_indent=myIndent, subsequent_indent=myIndent)
                        niceMsgList.append(wLine)
                else:
                    pass

            self.msg = '\n'.join(niceMsgList)

        else:
            self.msg = ''
        Exception.__init__(self)
        self.flavour = flavour

    def __str__(self):
        return self.msg
