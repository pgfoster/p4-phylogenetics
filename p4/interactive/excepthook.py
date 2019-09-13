import sys,traceback
import os
from p4.var import var
import p4.func

def invoke_editor(extractedTraceback):
    p4.func.writeInColour("\nWhere do you want to go? ...\n", colour="blue")

    te = extractedTraceback
    for teItemNum in range(len(te)):  # te is a traceback.extract_tb() result
        theTeItem = te[teItemNum]
        p4.func.writeInColour("%2i" % teItemNum, colour='red')
        print("  line %4i,  %s" %  (theTeItem[1], theTeItem[0]))
    p4.func.setTerminalColour("blue")
    ret = input('Tell me a number (or nothing to do nothing): ')
    p4.func.unsetTerminalColour()
    #print("Got %s" % ret)
    retNum = None
    if ret == '':
        pass
    else:
        try:
            retNum = int(ret)
            if retNum < 0 or retNum >= len(te):
                retNum = None
        except ValueError:
            pass
    if retNum != None:
        theTeItem = te[retNum]
        theFileName = theTeItem[0]
        if os.path.isfile(theFileName):
            try:
                theLineNum = int(theTeItem[1])
                theCommand = "%s +%i %s" % (var.excepthookEditor, theLineNum, theFileName)
                #print(theCommand)
                # It was sometimes not working well, so do it twice
                os.system(theCommand)
                os.system(theCommand)
            except:
                print("...could not make an int from theLineNum '%s'" % theLineNum)
                pass
        else:
            print("-> '%s' is not a regular file" % theFileName)
    #sys.exit()
       

if var.interactiveHelper == 'bpython':
    # modified from bpython
    def my_showtraceback(self):
        """This needs to override the default traceback thing
        so it can put it into a pretty colour and maybe other
        stuff, I don't know"""
        if 1:
            t, v, tb = sys.exc_info()
            sys.last_type = t
            sys.last_value = v
            sys.last_traceback = tb
            tblist = traceback.extract_tb(tb)
            del tblist[:1]
            # Set the right lineno (encoding header adds an extra line)
            #if not py3:
            for i, (fname, lineno, module, something) in enumerate(tblist):
                if fname == '<input>':
                    tblist[i] = (fname, lineno - 1, module, something)

            l = traceback.format_list(tblist)
            if l:
                l.insert(0, "Traceback (most recent call last peter):\n")
            l[len(l):] = traceback.format_exception_only(t, v)
        # finally:
        #     tblist = tb = None

        #self.writetb(l)
        #self.write("Over-riding bpython's showtraceback with the p4 version to invoke emacsclient...\n")
        for line in l:
            self.write(line)
        if var.excepthookEditor:
            invoke_editor(tblist)

    from bpython.repl import Interpreter
    Interpreter.showtraceback = my_showtraceback
    del(my_showtraceback)

def myExceptHook(exctype, excvalue, tb):
    try:
        te = traceback.extract_tb(tb)
    except:
        print("tb is %s" % tb)
        te = traceback.extract_tb(tb)

    # This next line does a regular traceback.
    sys.__excepthook__(exctype, excvalue, tb)

    if var.excepthookEditor:
        invoke_editor(te)

sys.excepthook = myExceptHook
del(myExceptHook)
    
    
