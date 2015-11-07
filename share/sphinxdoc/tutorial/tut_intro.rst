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

   * Completion, for quick reminders.  See :ref:`InteractiveHelpers`.


P4 is used in Python -- it is a Python package.  You can import p4
into your Python as usual, ::

  import p4

or (and I prefer this, but you might not) you can use the p4
script.  The p4 script loads the package as::

  from p4 import *

so all the p4 stuff goes to the 'top level' (which you may not like).  The p4 script will
read in files from the command line.

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
|
|


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


..
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


