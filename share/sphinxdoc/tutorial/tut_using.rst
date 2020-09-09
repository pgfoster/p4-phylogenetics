========
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
     ``var.nexus_allowAllDigitNames``::
 
       var.nexus_allowAllDigitNames = True

   * The symbols ``-``, ``?``, and ``.`` are the only ones allowed for gaps,
     unknown, and matchchars, respectively.


.. _tut_using_var_bucket_label:
   
The var object, the global bucket
---------------------------------



.. note::
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
   

.. _tut_using_read_function_label:

The ``read()`` function
-----------------------

The usual way to get phylogenetic data and trees into p4 is with the
:func:`~p4.func.read` function.  This function is used in scripts that you might
write; it is also used within the ``p4`` script to read in files from the
command line.  

Nexus files or files containing trees often contain multiple things
which would be turned into multiple Python objects, so a command
like::

     myAlignment = readNexus('myData.nex') # Doesn't work

would not always work.  The :func:`~p4.func.read` function does not return anything;
rather, when you read in alignments from files, Alignment objects are
made and stuffed into ``var.alignments``, a list, as explained above 
in the section :ref:`tut_using_var_bucket_label`.  So a typical script
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

When you read in files from the command line using the ``p4`` script, you
can use file name **globbing**, *ie* use wildcards, as::

     p4 *.phy

In that case the shell does the file name expansion.  You can also do
the same sort of thing with ``read()``, as in::

     read('*.phy')

If you want to make p4 forget the trees that it has read, say::

     var.trees = []  

You can do the same for alignments and sequenceLists.  To make p4
forget any Nexus sets info that it has read, you can do::

     var.nexusSets = None

If you are sure that the file that you are trying to read has only one
thing (Tree, Alignment, SequenceList), then you can use
:func:`p4.func.readAndPop`.  This is handy if for example you are reading
in a bunch of sequence list files and one of them happens to have all
its sequences the same length -- so it gets promoted to an Alignment
object and gets put in ``var.alignments``.  If you use
:func:`p4.func.readAndPop` then you would not need to check. ::

    # awkward using read()
    for fName in myFileList:
        var.alignments = []
        var.sequenceLists = []
        read(fName)
        try:
            sl = var.sequenceLists[0]
        except IndexError:
            sl = var.alignments[0]
       
    # easier using func.readAndPop()
    for fName in myFileList:
        sl = func.readAndPop(fName)   # SequenceList or Alignment


The ``dump()`` function and methods
-----------------------------------

There is a function, :func:`p4.func.dump` that gives a quick summary of files
that you have read, and objects that have been made and placed in
``var.trees``, ``var.alignments``, and so on.  It does not know about alignments
and such that are not in ``var.*``.

Several classes have ``dump()`` methods as well.  For example, to see inside
trees in fine detail, you can use the :meth:`p4.tree.Tree.dump` method, for example::

     t.dump()

or::

     t.dump(all=True)

To see details about models, use :meth:`p4.model.Model.dump`, for example::

     t.model.dump()

