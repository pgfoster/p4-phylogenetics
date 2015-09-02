============
Introduction
============


**Help** is available in various places--

   * This tutorial (that you are reading now) is a tutorial introduction,
     but it is not a comprehensive reference.  It skips a lot of detail.

   * The examples.  The examples should have been installed along with
     the rest of p4, in a share directory.  See the splash screen to
     find out where they are.  Hopefully they work as advertised.  Copy
     them to a working directory, and shake 'em 'til you break 'em.

   * The class and module descriptions :ref:`here
     <classes_and_modules>`.  It should be fairly comprehensive and
     definitive.

   * Those descriptions mentioned in the previous point are also
     available in the Python help() and dir() functions. For example,
     try help(Alignment) from within p4 to read all the docstrings in
     the Alignment class.

   * Completion, for quick reminders.  See :ref:`Completion`.


P4 is used in Python -- it is a Python package.  You can import p4
into your Python as usual, ::

  import p4

or (and I prefer this, but you might not) you can use the p4
script.  The p4 script loads the package as::

  from p4 import *

so all the p4 stuff goes to the 'top level' (which you may not like).  The p4 script will
read in files from the command line, and if you are interactive it will
set up completion.

This week, p4 will recognize and read data files in the form of Nexus
files, fasta files, Phylip (sequential or interleaved), GDE files
(proper gde files only, not gde flat files), and clustalw aln files. P4
will recognize trees in Nexus and Phylip or raw Newick format. And of
course p4 recognizes Python files. P4 understands Nexus data, taxa,
characters, trees, and sets blocks, but it does not understand all the
commands within those blocks.

| Q: So why is it called p4?
| A: Its an acronym- pppp
| Q: So what's it short for?
| A: So its easy to type.


You need to know this
=====================

When you import p4, or start the p4 script, you import ``var``, which
is an object that holds things like lists and flags that are generally
useful to you and p4.

With the p4 script you can read things from the command line, via::

  p4 aFile anotherFile

If those files contain phylogenetic data, they get put in

- ``var.sequenceLists``  (a list)
- ``var.alignments`` (a list)
- ``var.trees`` (a list)
- ``var.nexusSets`` (not a list; if it exists it is a NexusSets object)


For the impatient
=================

Here we provide a taster for some of the toolkit aspects of p4.  These
examples are very easy, and (except for the likelihood calculation) can
be done using p4 interactively.  You will probably find that for
anything bigger writing a script would be easier.

These quickstart examples can be found in the Examples, in the
A_quickStart directory.  I suggest that you make a copy of the Examples,
then change directory to A_quickStart in your copy, and then cd to
directories therein.  Be sure to read the files.


Draw a tree on the screen
-------------------------

Lets say you have a file that has a tree description in Newick or nexus
format, and you want to draw that tree to the screen to see what it is.
Say, to your shell::

     p4 -d t.nex

(or whatever your file name is)


Convert between file formats
----------------------------

Lets say that you have an alignment in Phylip format and you want to
convert it to nexus format.  This is a small job, so we can use p4
interactively.  First, read in the data, then give the alignment a name,
and then tell the alignment to write itself in nexus format to a file ::

     $ p4 data.phy                # Read the data at the command line
     p4> a = var.alignments[0]    # Give the alignment a name
     p4> a.writeNexus('data.nex') # Tell it to write itself to a file
     p4>                          # Control-d to quit.


Make an eps picture of a tree
-----------------------------

Lets say that you have a file with a tree description, and you want to
make a nice picture of it (rather than just a text screen picture).  The
following makes an encapsulated postscript (eps) file ::

     $ p4 t.nex                 # Read the tree from the command line
     p4> t = var.trees[0]       # Give the tree a name
     p4> t.eps()                # Tell it to make an eps file.
     p4>                        # Quit with control-d

You can view the resulting file with for example ``gv``, ie ``ghostview``,
a previewer for ghostscript.


Extract a charset from an alignment
-----------------------------------

Lets say you have a multi-gene alignment, and you want to extract one of
the genes and put it in an alignment file by itself.  Here, the
placement of the various genes in the alignment is defined in a nexus
sets block, in this example in the file ``sets.nex``.  Do something like
this::

     $ p4 d.phy sets.nex                   # Read everything
     p4> a=var.alignments[0]               # Name the alignment
     p4> b = a.subsetUsingCharSet('gene1') # b is a new alignment
     p4> b.writeNexus('gene1.nex')         # Tell b to write itself
     p4>                                   # Quit with control-d


Tabulate the character frequencies in the data
----------------------------------------------

::

     $ p4 d.nex             # Read the data
     p4> d=Data()           # Make a Data object
     p4> d.compoSummary()   # Ask for a table of composition


Compare the topology of two trees
---------------------------------

Lets say that you have two best trees from two different phylogenetic
analyses using the same data.  The trees might have the same topology,
but the trees might not look much like each other, because they are
rooted on different nodes, and many branches are rotated on their stems.
However, ignoring branch lengths, the trees might still have the same
topology, regardless of the various permutations.  You can quickly find
out how different they are, by doing something like the following::

     $ p4 tree1.nex tree2.nex          # Read in both trees
     p4> t1 = var.trees[0]             # Name the first one
     p4> t2 = var.trees[1]             # Name the second one
     p4> print t1.topologyDistance(t2) # Symmetric difference
     0                                 # Zero means they are the same

The default metric for the ``topologyDistance()`` method is the symmetric
difference, aka the unweighted Robinson-Foulds distance, which is the
number of splits in one tree that are not in the other, plus the number
of splits in the other tree that are not in the one.  In this example,
the trees are the same, and so the difference is zero.  If the two
trees had only one difference, the symmetric difference would be 2.

See also :meth:`Tree.Tree.tvTopologyCompare`


A very simple likelihood calculation
------------------------------------

This example is a bit more involved, and is not well suited to
interactive use.  The usual way to use p4 would be to make a script, and
that is what we do here.  Make a file with the following, and save it as
``s.py``::

     read(""" 2 2
     one
     ac
     two
     gt
     """)
     read('(one,two);')
     t = var.trees[0]
     t.data = Data()
     t.newComp()
     t.newRMatrix()
     t.setPInvar()
     t.calcLogLike()

Usually p4 scripts refer to other files for the data and the tree, but
here it is all in the one script.  Sequence data, a tree, and a model
are described and then the likelihood is calculated without
optimization.  To make the script happen, say, to your command line::

     p4 s.py


Using p4
========


Using the p4 script
-------------------

As usual for a Python package, you can import p4 from within Python, by
``import p4``.  However, my favourite way of using it is with the ``p4``
script (which should have been installed, hopefully somewhere in your
path).  The p4 script does a ``from p4 import *``, so your top-level
namespace becomes instantly cluttered, but then you do not need to
precede everything with ``p4.``.  The p4 script allows you to read in
phylogenetic data and tree files at the command line.  To see if it is
installed and working, say::

     p4 --help


Some peculiarities, bugs, and "features"
-----------------------------------------

   * P4 doesn't do searches with ML.  It only evaluates explicit trees.
     It will do tree searches with Bayesian analysis.

   * Arrays and such in Nexus files are 1-based, while in Python and p4
     they are zero-based.  Nexus files are case-insensitive in p4, as
     they should be; but elsewhere in p4, as in Python, they are
     case-sensitive.

   * Nexus trees do not require a Taxa or Data block.  Counter to the
     Nexus standard, having different trees with different taxa in the
     same trees block is acceptable in p4.  However, if you do supply a
     taxa block, p4 will use it to enforce those taxon names in the
     trees.

   * Often the results are not written automatically: you have to tell
     p4 to write the result.  However, you can tell it how you would
     like the output formatted.  Since the interface is a programming
     language you can do whatever you want with the numbers.

   * You give starting values for optimization.  If the params are
     already optimized, then it shouldn't take long for re-optimization.

   * The default branch length in trees is 0.1

   * When reading Nexus files, p4 recognizes datatype DNA, but not RNA.
     If you give it datatype nucleotide, p4 will assume you mean DNA.

   * The Nexus spec does not allow names that are all numerals, and by
     default p4 follows that rule.  You can (eg temporarily) override
     that behaviour and allow all-digit names by True-ing the variable
     ``var.nexus_allowAllDigitNames``

   * The symbols ``-``, ``?``, and ``.`` are the only ones allowed for gaps,
     unknown, and matchchars, respectively.


.. _Completion:

Using completion
----------------

Who needs pull-down menus or dialog boxes when you have ... completion!

Completion is a memory aid and can save you a lot of keystrokes. To use
this, you need to have the ``readline`` library linked to your Python. I
have re-written the wonderful ``rlcompleter`` module that comes with
Python so that it is slightly more helpful. Completion of course only
works if you are using p4 interactively.

Completion will remind you of methods and functions, and their
arguments, classes, variables, and documentation.

If you type ``<tab>``, a partially typed name is completed up to the point
of ambiguity. If it is unambiguous, the whole thing is completed. If
what you typed is ambiguous, then the function or method is competed
only up to the point of ambiguity. At a point of ambiguity, typing
``<tab><tab>`` tells you your options.

(If you are using the mac, see :ref:`completion_on_the_mac`.)

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

(Actually, only a subset of names are given. That subset is the "user
interface". There may be more names available that don't show up with
completion, but those invisible names are considered to be the
programmer's interface, not the user's interface.  Everything shows up
with help() and dir())

If you type a function or method name up to the argument list, but not
including the opening parenthesis, and then a dot, and then <tab><tab>,
then the documentation for the function or method is given, if it
exists. For example, do this to get the documentation for the read
function:

     read.<tab><tab>

To get the documentation for the Alignment method translate(),

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

You can back over the unparen and get the argspec by::

     a.translate(<tab>

which then tells you::

     a.translate(transl_table=1)

In this example, I asked for completion using an Alignment instance a,
but it works using the Alignment class as well. Note that in methods,
the first argument is ``self``; I skip that in argspecs via completion, but
retain it in the documentation output via completion.

(Tip: you can delete the entire line backwards with control-u.)



The ``read()`` function
-----------------------

The usual way to get phylogenetic data and trees into p4 is with the
:func:`func.read` function.  This function is used in scripts that you might
write; it is also used within the ``p4`` script to read in files from the
command line.  

Nexus files or files containing trees often contain multiple things
which would be turned into multiple Python objects, so a command
like::

     myAlignment = readNexus('myData.nex') # Doesn't work

would not always work.  The ``read()`` function does not return anything;
rather, when you read in alignments from files, Alignment objects are
made and stuffed into var.alignments, a list.  P4 puts alignments in
``var.alignments``, trees in ``var.trees``, unaligned sequences in
``var.sequenceLists``, and Nexus sets in ``var.nexusSets``.  So a typical script
might start out by reading in data and trees, like this::

     read('myAlignment.nex')
     a = var.alignments[0]
     read('myTreeFile.nex')
     t = var.trees[0]

There might be several trees in ``myTreeFile.nex``, and when the file is
read in they are all converted to Tree objects and put in the ``var.trees``
list.  So if you want one of those trees, you have to say which one.  To
get the first one, say something like:: 

     firstTree = var.trees[0]

When you read in files from the command line using the p4 script, you
can use file name globbing, ie use wildcards, as::

     p4 *.phy

In that case the shell does the file name expansion.  You can also do
the same sort of thing with ``read()``, as in::

     read('*.phy')

If you want to make p4 forget the trees that it has read, say::

     var.trees = []  

You can do the same for alignments and sequenceLists.

If you are sure that the file that you are trying to read has only one
thing (Tree, Alignment, SequenceList), then you can use
:func:`func.readAndPop`.  This is handy if for example you are reading
in a bunch of sequence list files and one of them happens to have all
its sequences the same length -- so it gets promoted to an Alignment
object and gets put in ``var.alignments``.  If you use
:func:`func.readAndPop` then you would not need to check. ::

    # awkward using read()
    for fName in myFileList:
        var.alignments = []
        var.sequenceLists = []
        read(fName)
        try:
            sl = var.sequenceLists[0]
        except IndexError:
            sl = var.alignments[0]
       
    # easier using readAndPop()
    for fName in myFileList:
        sl = func.readAndPop(fName)   # SequenceList or Alignment


The ``dump()`` function and methods
-----------------------------------

There is a function, :func:`func.dump` that gives a quick summary of files
that you have read, and objects that have been made and placed in
``var.trees``, ``var.alignments``, and so on.  It does not know about alignments
and such that are not in var.\*.

Several classes have ``dump()`` methods as well.  For example, to see inside
trees in fine detail, you can use the :meth:`Tree.Tree.dump` method, for example::

     t.dump()

or::

     t.dump(all=True)

To see details about models, use Model.dump(), for example::

     t.model.dump()


Customizing p4
--------------

   * You can change p4 variables found in p4.var

   * You can semi-permanently add your own code into files that p4
     reads on startup.  Doing that you can extend or over-ride the
     behaviour of p4.  

     - You want p4 Alignments to do something else? - Add it!  

     - You don't like the way that p4 Tree objects are written?  - Change it!


You can change the default behaviour by setting or un-setting some
variables in var (which is an instance of the Var class, but you don't
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

There are several other things in Var.py that you might want to know about.
Take a :ref:`look <var-class-ref-label>`.

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
``~/.p4`` directory, or in another of your specified ``P4_STARTUP_DIRS``.
You can easily turn off files, but allow them to remain there for
later, by changing the file name so that it no longer ends in ``.py``, or
just put them in a sub-directory.

One special file that is read in after the others, but only if you are
interactive, is the ``~/.p4/interactive`` file, and note that while it is
a Python file, it does not end in ``.py``.  I like to instantiate a few
blank objects in there, for quick access to the doc strings via
completion.  Here is my current ``interactive`` file, for what its worth::

     #print 'Here is the interactive startup file.'

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
             print "There is one alignment in var.alignments-- I am naming it 'a'"
             a = var.alignments[0]
         if len(var.trees) == 1:
             print "There is one tree in var.trees-- I am naming it 't'"
             t = var.trees[0]
         if len(var.sequenceLists) == 1:
             print "There is one sequenceList in var.sequenceLists-- I am naming it 'sl'"
             sl = var.sequenceLists[0]

     #import sys
     #sys.ps1 = '>>> '

