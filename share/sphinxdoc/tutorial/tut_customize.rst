============
Customizing
============



P4 is meant to be customized.

   * You can change variables found in :class:`p4.var.Var` 

   * You can semi-permanently add your own code into files that p4
     reads on startup.  Doing that you can extend or over-ride the
     behaviour of p4.  

     - You want p4 Alignments to do something else? - Add it!  

     - You don't like the way that p4 Tree objects are written?  - Change it!

Setting var variables
---------------------

You can change the default behaviour by setting or un-setting some
variables in ``var`` (which is an instance of the :class:`p4.var.Var` class, but you don't
really need to know that).  For example, when p4 reads files it is by
default not verbose, not describing what it is doing as it goes.  This
is good everyday behaviour, but not very good for debugging when
something goes wrong.  You can make it verbose by setting::

    var.verboseRead=True

Generally you will get your phylogenetic stuff into p4 from files, using
the read() function.  The read() function also reads strings, which is
occasionally handy.  If you hand the read() function a string, it will
first ask if there is a file by that name, and if so it will attempt to
read it.  If there is no file, it will try it as a string.  The problem
is that if you mis-spell a file name, you are liable to get a bunch of
confusing error messages as it attempts to read the file name as some
sort of phylogenetic data.  So the default behaviour when it can't find
a file is to give a little warning, and you can turn that warning off by
setting:: 

    var.warnReadNoFile=False

There are several other things in the :class:`p4.var.Var` class that you might want to know about.
Take a :ref:`look <var-class-ref-label>`.

.. note::
   You can make these settings to ``var``, either 

   * temporarily

     - in interactive p4   
     - at the top of individual scripts
   * or in startup scripts

     - possibly in your ``~/.p4`` directory
 

Add your own code
-----------------

You can of course add your own code in your p4 script to extend p4 -
perhaps to add your own function, or to extend a p4 Class with an
additional method or two.  You might put that part of the code in a
separate file in your working directory and ask your script to load
that module.  But if you find your new code generally useful you might
want to semi-permanently add it to your p4, so that it is read in for
any script in any directory, without you needing to explicitly ask it
to do so.  To do that, you can put your code in a directory known to
p4, and p4 will read it on startup.  One such directory is ``~/.p4``.  So
if you put python files in there, p4 will find them and read them on
startup.  Other directories can be specified by the environment
variable ``P4_STARTUP_DIRS``.  That variable is a colon-separated list of
directory names.

So for example when p4 starts up, by default it prints out a **splash
screen**.  Whether that is done or not is controlled by ``var.doSplash``.
When you get tired of seeing the splash screen, you can turn it off by
setting ``var.doSplash=False``.  But wait! - By the time you have the
opportunity to set that variable, its already too late-- the splash
screen has already displayed itself.  The solution is to put
``var.doSplash=False`` in a file ending in ``.py`` in a ``.p4`` directory in
your home directory, or in another of the ``P4_STARTUP_DIRS``.  Any file
ending in ``.py`` in one of those directories is read on startup.  Its a
good place to put ``var.warnReadNoFile=False``, for example.

Since it is Python, you can add methods to existing classes.  For
example, if you have this in your script::

     def greet(self):
         print "Howdy!"
     Alignment.greet = greet
     del(greet)

and you had an Alignment instance ``a``, you could then do::

     p4> a.greet()
     Howdy!
     p4>

You can over-ride existing methods that way also.  If you don't like the
way p4 works, change it to suit yourself.

A good place to put these new methods and functions is in files in the
``~/.p4`` directory, or in one of your specified ``P4_STARTUP_DIRS``.
You can easily turn off files, but allow them to remain there for
later, by changing the file name so that it no longer ends in ``.py``, or
just put them in a sub-directory.

One special file that is read in after the others, but only if you are
interactive, is the ``~/.p4/interactive`` file, and note that while it is
a Python file, it does not end in ``.py``.  I like to instantiate a few
blank objects in there, for quick access to the doc strings via
completion (see below).  Here is my current ``interactive`` file, for what it's worth::

     #print('Here is the interactive startup file.')

     # Make some blank objects.
     if 1:
         a0 = Alignment()
         d0 = Data([])
         t0 = Tree()
         n0 = Node()
         #tt0 = Trees()

     # Give single Alignments, Trees, or SequenceLists names.
     if pyFileCount == 0:
         if len(var.alignments) == 1:
             print("There is one alignment in var.alignments-- I am naming it 'a'")
             a = var.alignments[0]
         if len(var.trees) == 1:
             print("There is one tree in var.trees-- I am naming it 't'")
             t = var.trees[0]
         if len(var.sequenceLists) == 1:
             print("There is one sequenceList in var.sequenceLists-- I am naming it 'sl'")
             sl = var.sequenceLists[0]

     #import sys
     #sys.ps1 = '>>> '


.. _InteractiveHelpers:

Using interactive helpers
-------------------------

If you are using p4 interactively you can set it up so that you get

* completion
* previous commands
* signatures of functions and methods
* doc strings

This week, there are three ways to do it ---

* p3rlcompleter (comes with p4, based on rlcompleter that comes with Python)
* `bpython <http://bpython-interpreter.org>`_
* `ipython <http://ipython.org>`_

You can turn it on by setting :attr:`p4.var.Var.interactiveHelper`, which is by default ``None``, for example::

    var.interactiveHelper = 'p3rlcompleter'


In p3rlcompleter, 

* you get signatures (argspecs) inserted in place
* in the attribute listing, you get an indication of whether the attribute is a function or method (with a ``()``) or a list, a numpy array, or a plain variable

However, ipython and bpython look a lot nicer!  With colour!



Using p3rlcompleter
-------------------

Completion is a memory aid and can save you a lot of keystrokes. To use
this, you need to have the ``readline`` library linked to your Python. I
have modified the wonderful ``rlcompleter`` module that comes with
Python so that it is slightly more helpful.
Completion will remind you of methods and functions, and their
arguments, classes, variables, and documentation.

If you type ``<tab>``, a partially typed name is completed up to the point
of ambiguity. If it is unambiguous, the whole thing is completed. If
what you typed is ambiguous, then the function or method is competed
only up to the point of ambiguity. At a point of ambiguity, typing
``<tab><tab>`` tells you your options.

.. (If you are using the mac, see :ref:`completion_on_the_mac`.)

For example, if you type::

	func.chi<tab>

then, since this is unambiguous (there is only one function in the func
module that starts with chi) p4 will complete it, resulting in ::

	func.chiSquaredProb()

If you just type::

	func.<tab>

then nothing happens, because it is ambiguous at that point. Typing a
second <tab> tells you the available possibilities. Then you can type
more of the name to resolve the ambiguity and finish it up with a <tab>.

You can also get completion starting from nothing. If you just type ::

	<tab><tab>

then you get all the top level names that p4 knows about.

All this works for method names also. For example, type the name of the
class or object, and the dot. If you then type <tab><tab> all the
available methods and instance variable names appropriate to that
object are given.

.. 
   (Actually, only a subset of names are given. That subset is the "user
   interface". There may be more names available that don't show up with
   completion, but those invisible names are considered to be the
   programmer's interface, not the user's interface.  Everything shows up
   with help() and dir())

If you type a function or method name up to the argument list, but not
including the opening parenthesis, and then a dot, and then <tab><tab>,
then the documentation for the function or method is given, if it
exists. For example, do this to get the documentation for the read
function::

	read.<tab><tab>

To get the documentation for the Alignment method translate()::

	Alignment.translate.<tab><tab>

Generally, functions and methods tell you that they are functions and
methods by printing out a pair of parens after the name (eg foobar()).
You can get the argspec (the stuff that goes inside the parens) by
completing the name followed by a single open paren. So, for example,
using an Alignment instance a, to get the argspec for the method
translate(), you might say::

	a.tr<tab>

which gets completed to::

	a.translate()

You can back over the unparen and get the argspec (signature) by::

	a.translate(<tab>

which then tells you::

	a.translate(transl_table=1)

In this example, I asked for completion using an Alignment instance ``a``,
but it works using the Alignment class as well. Note that in methods,
the first argument is ``self``; I skip that in argspecs via completion, but
retain it in the documentation output via completion.

.. tip::
   You can delete the entire line backwards with control-u.
   Delete everything to the end of the line with control-k

Interfacing with your editor
----------------------------

In my world, I use a terminal to run p4 and emacs to edit my scripts and code.  When you use emacs to run your Python code and you get a traceback, you can use that traceback to easily go to those locations in the source.  That is brilliant functionality, but I have never got the hang of running Python in emacs  --- so I use a terminal.  Using the terminal is simpler, but I still want that traceback functionality.

Using the terminal, when I get a traceback from an exception I want to go to where the error is, or perhaps one of the previous places in the traceback stack, and I want to go quickly and easily and start editing.  For that I use an exceptionhook that invokes the editor and goes to the right line in the right file.  It is under the control of :attr:`p4.var.Var.excepthookEditor`. It is by default ``None``, but you can turn it on by::

    var.excepthookEditor = 'emacsclient -n'

It appears to work for vim, as well, although I have not tested it well.


 
