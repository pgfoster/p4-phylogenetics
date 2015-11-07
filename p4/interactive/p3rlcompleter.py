"""Peter's re-write of word completion for GNU readline

The foundation for this is the wonderful rlcompleter module that comes
with Python.  It has been re-written to make it more informative.

When you import this module, it does this for you:
    readline.parse_and_bind("tab: complete")
which makes the tab key the complete key, which you probably want.

As you would expect: When you hit the tab you get completion.  When
you hit the tab twice, you get choices.

An informative peculiarity: Completions for Classes and routines end
in parens, eg foobar() or MyClass().  Completions for modules and
variables do not.  Lists are shown as aList[], and NumPy arrays are
shown as someArray[N]

More: When you hit the tab twice after something ending in a dot, you
get documentation, if it exists.  If you hit the tab after some
routine or class name that ends in a single open paren ie
    foobar(<tab>
then you get the argspec.  This doesn't work for builtin or non-python
functions and such, because inspect can't get the argspec.  And it
does not work on the Python that comes with Mac OS 10.5, because that
python is built with editline, not readline.

For example: import this module, by typing
    import p3rlcompleter
Now type
    p3rl<tab>
It gets completed to p3rlcompleter.  Type
    p3rlcompleter.<tab><tab>
and you get this documentation, and the names in the module.  Type
    p3rlcompleter.C<tab>
and then
    p3rlcompleter.Completer.<tab><tab>
and you get the documentation (there isn't much!).  Then type
    p3rlcompleter.Completer(<tab>
and you get the argspec for the Completer.__init__ method.

(Not a particularly useful example, but hopefully you get the point.)

The names (choices) that you get when you hit <tab> twice do not
usually include those names starting with '_'; this is a feature, not a bug.
If you want those, you can give it a '_' after the dot to complete.  """

import readline,inspect,re
import __builtin__
import __main__

__all__ = ["Completer"]

try:
    from completionIgnores import *
except ImportError:
    #print 'Could not find completionIgnores.py'
    completionIgnores = []
    funcIgnores = []
#print completionIgnores

try:
    import numpy
    aNumPyArray = numpy.array([1], numpy.int32)
except ImportError:
    aNumPyArray = None
    

class Completer:
    """This is the class doc string for the Completer class.

    This doc string was made so that the tiny tutorial in the module
    doc string on how to use p3rlcompleter with itself as an example
    would have something to read for the class doc for Completer."""


    def __init__(self, namespace = None):
        """Create a new completer for the command line.

        Completer([namespace]) -> completer instance.

        If unspecified, the default namespace where completions are performed
        is __main__ (technically, __main__.__dict__). Namespaces should be
        given as dictionaries.

        Completer instances should be used as the completion mechanism of
        readline via the set_completer() call:

        readline.set_completer(Completer(my_namespace).complete)
        """

        if namespace and not isinstance(namespace, dict):
            raise TypeError,'namespace must be a dictionary'

        # Don't bind to namespace quite yet, but flag whether the user wants a
        # specific namespace or to use __main__.__dict__. This will allow us
        # to bind to __main__.__dict__ at completion time, not now.
        if namespace is None:
            self.use_main_ns = 1
        else:
            self.use_main_ns = 0
            self.namespace = namespace
        self.toggle = 1

    def complete(self, text, state):
        """Return the next possible completion for 'text'.

        This is called successively with state == 0, 1, 2, ... until it
        returns None.  The completion should begin with 'text'.

        """
        if self.use_main_ns:
            self.namespace = __main__.__dict__

        #print "complete()  text = |%s|" % text
        if state == 0:
            if "." in text:
                #print "\ncomplete(), calling attr_matches() with text '%s'" % text
                self.matches = self.attr_matches(text)
            else:
                #print "complete(), calling global_matches()"
                self.matches = self.global_matches(text)
            #print "complete(). state=0 self.matches = %s" % self.matches
        try:
            return self.matches[state]
        except IndexError:
            return None

    def global_matches(self, text):
        """Compute matches when text is a simple name.

        Return a list of all keywords, built-in functions and names currently
        defined in self.namespace that match.

        """
        #import keyword
        matches = []
        n = len(text)
        # keyword.kwlist is ['and', 'assert', 'break', 'class',
        # 'continue', 'def', 'del', 'elif', 'else', 'except', 'exec',
        # 'finally', 'for', 'from', 'global', 'if', 'import', 'in',
        # 'is', 'lambda', 'not', 'or', 'pass', 'print', 'raise',
        # 'return', 'try', 'while', 'yield']

        #for list in [keyword.kwlist,
        #             __builtin__.__dict__,
        #             self.namespace]:
##        for list in [self.namespace]:
##            for word in list:
##                #print 'considering word %s, isfunction=%s' % (word, inspect.isfunction(word))
##                if word[:n] == text and word != "__builtins__":
##                    print 'got word %s.  type = %s, isfunction=%s' % (word, type(word), inspect.isfunction(word))
##                    if inspect.isfunction(word):
##                        matches.append('%s()' % word)
##                    else:
##                        matches.append(word)

        #print 'gm text = |%s|' % text
        ks = self.namespace.keys()
        #ignores = ['__builtins__', '__doc__', '__file__', '__name__']
        for k in ks:
            #print '%20s %s' % (k, self.namespace[k])
            if not text and k.startswith('_'):
                pass
            elif k in completionIgnores:
                pass
            elif k[:n] == text: #and k not in ignores:
                if inspect.isfunction(self.namespace[k]) or inspect.isclass(self.namespace[k]):
                    theComments = inspect.getcomments(self.namespace[k])
                    if theComments and theComments.find('##Ignore') != -1:
                        pass
                    else:
                        matches.append('%s()' % k)
                elif inspect.ismodule(self.namespace[k]):
                    # Skip all modules.  Maybe I only want to skip
                    # ##Ignore'd modules?  If that is the case, I can
                    # inspect.getcomments(self.namespace[k]) as above.
                    # Comments in the module are at the beginning of
                    # the module, as you would expect.
                    #pass

                    # Changed my mind.  Now modules get added
                    matches.append(k)
                else:
                    matches.append(k)


            # If the user is asking for the argspec of the function or a class
            elif text.endswith('(') and k == text[:-1] and inspect.isfunction(self.namespace[k]):
                try:
                    args, varargs, varkw, defaults = inspect.getargspec(self.namespace[k])
                    argspec = inspect.formatargspec(args, varargs, varkw, defaults)
                    matches.append('%s%s' % (k,argspec))
                except TypeError:
                    matches.append('%s(?)' % k)
            elif text.endswith('(') and k == text[:-1] and inspect.isclass(self.namespace[k]):
                try:
                    args, varargs, varkw, defaults = inspect.getargspec(self.namespace[k].__init__)
                    if args[0] == 'self':
                        args = args[1:]
                    argspec = inspect.formatargspec(args, varargs, varkw, defaults)
                    matches.append('%s%s' % (k,argspec))
                except TypeError:
                    matches.append('%s(?)' % k)

            #print 'matches = %s' % matches

        return matches

    def attr_matches(self, text):
        """Compute matches when text contains a dot.

        Assuming the text is of the form NAME.NAME....[NAME], and is
        evaluatable in self.namespace, it will be evaluated and its attributes
        (as revealed by dir()) are used as possible completions.  (For class
        instances, class members are are also considered.)

        WARNING: this can still invoke arbitrary C code, if an object
        with a __getattr__ hook is evaluated.

        """
        
        #print "attr_matches() here, with text '%s'" % text
        #import re
        #m = re.match(r"(\w+(\.\w+)*)\.(\w*)", text)
        # peter changed the above so that parens would be part of the afterDot
        #m = re.match(r"(\w+(\.\w+)*)\.([\w\(]*)", text)
        # peter changed the above so that square brackets would be part of the beforeDot, so that it would complete lists
        m = re.match(r"([\w\[\]]+(\.[\w\[\]]+)*)\.([\w\(]*)", text)
        if not m:
            #print "\nattr_matches().  Failed to re.match() the text '%s'" % text
            return
        beforeDot, afterDot = m.group(1, 3)  # everything before the final dot, and everything after the final dot.
        #print '\nattr_matches here.  beforeDot = \'%s\', afterDot = \'%s\'' % (beforeDot, afterDot)
        theThing = eval(beforeDot, __main__.__dict__)
        thingDir = dir(theThing)
        import types
        #print 'theThing is \'%s\', type %s' % (theThing, type(theThing))
        #print 'completionIgnores = %s' % completionIgnores
        
        words = []
        #print "thingDir = %s" % thingDir
        for w in thingDir:
            #if w not in completionIgnores:
            #    print 'considering w=%s' % w
            #    print 'type %s' % type(getattr(theThing,w))
            #else:
            #    print "ignoring %s" % w
            #print 'ismodule=%s' % inspect.ismodule(getattr(theThing,w))
            if w[0] == '_' and not(afterDot and len(afterDot) and afterDot[0] == '_'):
                #print 'q(w=%s)' % w
                pass
            elif w in completionIgnores:
                #print "completionIgnores"
                pass
            elif beforeDot == 'func' and w in funcIgnores:
                pass

            # These next few lines put square brackets after any lists,
            # to tell you that they are lists.  This might be a good
            # idea, or maybe not.  If not, comment out these next few
            # lines ...
            elif type(getattr(theThing,w)) == type([]): # if its a list...
                #print "                its a list"
                words.append('%s[]' % w)

            # same idea with NumPy arrays
            elif aNumPyArray and type(getattr(theThing,w)) == type(aNumPyArray): # if its a NumPy array...
                #print "                 its an array"
                words.append('%s[N]' % w)

            elif inspect.isroutine(getattr(theThing,w)) or inspect.isclass(getattr(theThing,w)):
                # Look for a comment line that includes the string '##Ignore'.
                # Comment lines preceed the routine or Class def, at the same indentation level.
                # See for example this method-- I have put in a sample comment.
                #print "                  its a routine or a Class"
                theComments = inspect.getcomments(getattr(theThing,w))
                #print "theComments: ", theComments
                if theComments and theComments.find('##Ignore') != -1:
                    #print "        ##Ignore %s" % w
                    pass
                elif afterDot and len(afterDot) and afterDot.endswith('('):
                    if inspect.isroutine(getattr(theThing,w)):
                        #print 'x(w=%s)' % w,
                        try:
                            args, varargs, varkw, defaults = inspect.getargspec(getattr(theThing, w))
                            #print '(args=%s)' % args,
                            if args and len(args) and args[0] == 'self':
                                args = args[1:]
                            argspec = inspect.formatargspec(args, varargs, varkw, defaults)
                            #print '(argspec=%s)' % argspec,
                            words.append('%s%s' % (w, argspec))
                        except TypeError:
                            words.append('%s(?)' % w)
                    elif inspect.isclass(getattr(theThing,w)):
                        #print 'y(w=%s)' % w,
                        try:
                            if hasattr(getattr(theThing,w), '__init__'):
                                args, varargs, varkw, defaults = \
                                      inspect.getargspec(getattr(getattr(theThing,w), '__init__'))
                            else:
                                args, varargs, varkw, defaults = [], None, None, None
                            if args and len(args) and args[0] == 'self':
                                args = args[1:]
                            argspec = inspect.formatargspec(args, varargs, varkw, defaults)
                            words.append('%s%s' % (w, argspec))
                        except TypeError:
                            words.append('%s(?)' % w)
                    else: # when would this happen?  -- never!
                        words.append(w)
                else:
                    words.append('%s()' % w)
            else:
                #print "                its 'else'"
                words.append(w)
            #if len(words):
            #    print "x words=%s" % words
            #print "next word ..."

        #print 'y words = %s' % words

        if type(theThing) == types.InstanceType:
            if self.toggle:
                self.toggle = 0
            else:
                self.toggle = 1
                if len(afterDot) > 1:
                    pass
                else:
                    print '\n\nInstance of class %s' % theThing.__class__

        elif inspect.isroutine(theThing):
            #print 'Its a routine.  toggle=%s' % self.toggle
            if self.toggle == 0:
                self.toggle = 1
                typeOfTheThing = type(theThing)
                thingName = theThing.__name__
                thingModule = inspect.getmodule(theThing)
                if not thingModule:
                    if hasattr(theThing, '__module__'):
                        thingModule = theThing.__module__
                    else:
                        thingModule = 'unknown module'
                if typeOfTheThing == types.MethodType:
                    print '\n\n%s method \'%s\', in %s' % (theThing.im_class, thingName, thingModule)
                elif typeOfTheThing == types.FunctionType:
                    print '\n\nfunction \'%s\', in %s' % (thingName, thingModule)
                elif typeOfTheThing == types.BuiltinFunctionType:
                    print '\n\nbuiltin function \'%s\', in module %s' % (thingName, thingModule)
                else:
                    print '\n\nroutine \'%s\', in %s' % (thingName, thingModule)
                #try:
                #    print 'defined in file %s' % inspect.getfile(theThing)
                #except TypeError:
                #    print 'boing!',
                #    pass
                argspec = None
                try:
                    args, varargs, varkw, defaults = inspect.getargspec(theThing)
                    #if args[0] == 'self':
                    #    args = args[1:]
                    argspec = inspect.formatargspec(args, varargs, varkw, defaults)
                except TypeError:
                    pass
                if argspec:
                    print ''
                    print theThing.__name__ + argspec
                    print ''
                else:
                    print ''
                theDocString = inspect.getdoc(theThing)
                if theDocString:
                    print theDocString
                else:
                    print 'No documentation available'
                if afterDot and len(afterDot):
                    pass
                else:
                    return ['', '']
            else:
                self.toggle = 0
                return ['', '']

        elif inspect.isclass(theThing) or inspect.ismodule(theThing):
            #print 'Its a class.'
            if self.toggle:
                self.toggle = 0
                #print '\n\nClass %s, from module %s' % (theThing.__name__, inspect.getmodule(theThing))
                # if there is nothing after the dot, give the documentation
                if not afterDot or len(afterDot) == 0:
                    #print 'words = %s' % words
                    #if len(words) <= 1:
                    #    return ['', '']
                    return ['', '']
            else:
                self.toggle = 1
                if not afterDot or len(afterDot) == 0:
                    if inspect.isclass(theThing):
                        print '\n\nClass %s, from module %s' % (theThing.__name__, inspect.getmodule(theThing))
                    elif inspect.ismodule(theThing):
                        print '\n\nModule %s' % theThing.__name__
                    print ''
                    theDocString = inspect.getdoc(theThing)
                    if theDocString:
                        print theDocString
                    else:
                        print 'No documentation available'
                    if not len(words):
                        return ['', '']


        matches = []
        #print 'words = %s' % words

        n = len(afterDot)  # afterDot is only the stuff after the last dot, so it might be nothing
        #print 'afterDot = \'%s\'' % afterDot
        for word in words:
            if word[:n] == afterDot and word != "__builtins__":
                matches.append("%s.%s" % (beforeDot, word))
        #print 'c attr_matches returning %s' % matches
        return matches


##def get_class_members(klass):
##    ret = dir(klass)
##    if hasattr(klass,'__bases__'):
##        for base in klass.__bases__:
##            ret = ret + get_class_members(base)
##    return ret

readline.set_completer(Completer().complete)

# peter added, so that it would not split words in the input line on a "("
import string
delims = readline.get_completer_delims()
#print "delims = %s" % delims
delimList = list(delims)
if '(' in delimList:
    delimList.remove('(')
if 1:
    delimList.remove('[')
    delimList.remove(']')
#print delimList
readline.set_completer_delims(string.join(delimList, ''))

# make the tab work for completion
try:
    from Var import var
    #print "var.readlineUsesEditline is %s" % var.readlineUsesEditline
    if var.readlineUsesEditline:
        # For Mac OS 10.5, which uses editline, not readline, and has a different syntax.
        readline.parse_and_bind("bind ^I rl_complete")
    else:
        readline.parse_and_bind("tab: complete")
        
except ImportError:
    readline.parse_and_bind("tab: complete")


