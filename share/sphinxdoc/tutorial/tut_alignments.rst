===================
Alignments and data
===================

You can read in an alignment from a file either with the p4 script with
the file name as a command-line argument, or using the read() function
within p4.  Any alignments that are made are put in var.alignments, a
list.  (If its a fasta file and the sequences are not all the same
length, it is not an alignment, it is a SequenceList, and it is put in
var.sequenceLists, another list.)  So a lot of p4 interactive sessions
might start out like this::

     yourShellPrompt$ p4 myData.nex
     p4> a = var.alignments[0]

and a lot of p4 scripts might start out like this::

     read('myData.nex')
     a = var.alignments[0]

Saying ``read()`` here reads the file, makes an Alignment object from
it, and put that object in ``var.alignments``.  Then you give it the
name ``a``.  Then you can do useful things with ``a``.



Alignment, Part, and Data classes
=================================

There are 3 classes.

   - **Alignment**.  A SequenceList subclass, with a lot of functionality
     in Python.

   - **Part**.  A data partition, with the sequence data in the C language, for
     likelihood calculations.  Each alignment has at least one part,
     maybe more.

   - **Data**.  All the partitions that you want to use.  It is at
     least one Alignment (but perhaps more), with at least one Part
     (at least one Part for each Alignment, but perhaps more).  A Data object is
     needed for likelihood calculations.  It has a C-language
     representation of the sequences as a list of partitions.


An alignment of phylogenetic data in a file can be read in to make an
Alignment object, which is good for simple tasks that you might want to
do with alignments.  However, if you want to calculate the likelihood of
the data, then things get more complicated because then we transfer the
alignment into the C-language.  A Part object encapsulates an alignment
(or an alignment partition) in C.  A Data object contains all the Part
objects that you want to use, together with the alignments that the
parts came from.

When you do likelihood calculations, your data might be simple-- they
might be in just one alignment, with only one data partition.  In that
case your Data object has only one Alignment object, and one Part
object.  Or you might have 2 alignments, each of which has a single part
each, and so your Data object would have 2 alignments and 2 parts.  Or
perhaps you have a single alignment which is partitioned into 5 genes,
in which case your Data object would have one alignment and 5 parts.

The way to get phylogenetic data into p4 is with the 'read()' function.
When you read in alignments from files, Alignment
objects are made and stuffed into var.alignments, a list.  So a typical
script to do something simple (ie not likelihood) with an alignment
might start out by reading in a file, like this::

     read('myAlignment.nex')
     a = var.alignments[0]

For likelihood calculations you need a Data object, so a typical script
for that might start out like this::

     read('myAlignment.nex')
     d = Data()

In either case when you read in the file, an Alignment object is made
and deposited in var.alignments.  It is there whether or not you give it
a name, and whether or not you use it to make a Data object.

If you need to partition your alignment, then you might do something
like::

    read("myAlignment.nex")
    read("mySets.nex")  # with the partition info
    a = var.alignments[0]
    a.setCharPartition('myCharPartitionName')
    d = Data()  # now with more Parts that if the Alignment had not been partitioned

Part objects are made when you make a Data object.

**The data that you want to use might be in more than one alignment, and
that is perfectly acceptable to p4.**  In fact, it is essential if you
are using more than one datatype (eg the DNA sequence for one gene,
analysed simultaneously with the protein sequence for another gene).  A
MrBayes-style data matrix with more than one datatype is not allowed by
p4 (nor is it Nexus compliant).  So you might do something like the
following to read in two alignments::

    read('myFirstAlignment.nex')
    read('mySecondAlignment.nex')
    d = Data()

P4 analyses data partition-by-partition.  Separate alignments
necessarily make separate data partitions, but you can also split up
alignments into separate partitions--that would be done using Nexus
sets blocks.  To do this you might do something like::

    read('myAlignment.nex')
    a = var.alignments[0]
    read('myNexusSetsBlock.nex') # Defines cp1
    a.setCharPartition('cp1')    # Applies cp1 to the alignment
    d = Data()  # now with more Parts that if the Alignment had not been partitioned

In order to get a Data object from one alignment that you want to
partition into 2 parts together with one unpartitioned alignment, you
might do something like the following::

    read('myFirstAlignment.nex')
    a = var.alignments[0]
    read('myNexusSetsBlock.nex') # defines cp1
    a.setCharPartition('cp1')
    read('mySecondAlignment.nex')
    d = Data()

The ``d`` Data instance now contains 3 Part objects (in ``d.parts``) and 2
Alignment objects (in ``d.alignments``).

If you are going to use a Data object together with a tree, eg for ML
evaluations, Bayesian analysis, or simulation, you need to attach the
Data object to the tree, as::

     t.data = d


Nexus charpartitions
--------------------

This week, p4 understands Nexus charset and charpartition commands
within Nexus sets blocks.  When a Nexus sets block is read in, its
contents are added to var.nexusSets, a NexusSets object which can hold
charsets and charpartitions.  When you want to use a charpartition on an
alignment, you would say something like::

     a = var.alignments[0]
     a.setCharPartition(theCharPartitionName)

That attaches a copy of var.nexusSets to the alignment.  To un-set a
charpartition from an alignment, you can say ``a.setCharPartition(None)``.

Setting a charPartition does not actually do anything to the sequences
in an Alignment, rather a charPartition only comes into play when you
make a Data object.  When you make a data object from your partitioned
alignments, you can check that you are indeed dealing with partitioned
data by doing a ``dump()`` method on your data object.


Subsetting alignments
---------------------

By *subsetting* I mean extracting a part of an alignment into its own
new alignment.  You can subset alignments based on charsets or mask
strings.

If you do not have a Nexus sets block, you can still use the Nexus
character list format to make a mask on the fly using the
``func.maskFromNexusCharacterList()`` function, for example::

     p4> print func.maskFromNexusCharacterList("1 3-5 7 11", 12, invert=1)
     010001011101


Nexus standard datatype
-----------------------

P4 supports Nexus standard datatype.

Something I have been looking into is **recoding protein alignments into
amino acid groups**, for example the 6 Dayhoff groups, for analysis.  That uses the
:meth:`Alignment.Alignment.recodeDayhoff()` method or the
:meth:`Alignment.Alignment.recodeProteinIntoGroups` method.  Here the
amino acid symbols are recoded to the numerals.  For example in
Dayhoff recoding, since C (Cysteine) is in a group by itself, it is
simply recoded as 1, but the amino acids S, T, P, A, and G are all
recoded as 2, as they are all in the same Dayhoff group.  It is rather
like recoding DNA into transversions to decrease the effects of
saturation.  It is based on the notion that changes within groups
might be saturated, biased, and difficult to model, while changes
between groups are less problematic.  Certainly there is loss of
information, but (hopefully) the information that is retained is of
higher quality.  See :ref:`kosiol_ais`.


(Non-parametric) bootstrapping of data
--------------------------------------

Bootstrapping an alignment is trivially easy, and lots of programs do
it, p4 included.  It is a method in the Alignment class.

However, **bootstrapping of partitioned data** is also supported by p4.  It
is a method of the Data class.  The alignment and partition structure is
maintained in the bootstrapping process.


Dealing with duplicate sequences
--------------------------------

There is no point in doing a deep phylogenetic analysis with
duplicated sequences. (There may be a reason for popgen level work.)
If you have any duplicated sequences, p4 will point them out when you
read in your data file.  P4 will remove the sequences if you tell it
to.  But often you want those dupes in your final result - generally a
tree.  So you will want to remove the duplicated sequences before you
do your analysis, and then restore the taxa to the resulting tree
after.  For example, lets say that you read in your data like this::

     yourShellPrompt$ p4 d3.nex

and then you get a long message complaining about duplicated sequences,
starting with::

      Alignment from file 'd3.nex'
      This alignment has duplicate sequences!
      Sequence numbers below are 1-based.
         sequence 1 (A) is the same as sequence 3 (C).
         sequence 1 (A) is the same as sequence 8 (H).
         sequence 2 (B) is the same as sequence 5 (E).
         sequence 2 (B) is the same as sequence 7 (G).
     ...

To remove the duplicated sequences and make a Python dictionary file
that will be useful later for restoring those removed taxa, you can do
something like this::

     read('d3.nex')
     a = var.alignments[0]
     a.checkForDuplicateSequences(removeDupes=True, makeDict=True)
     a.writePhylip(fName='d3_noDupes.phy')

(Above I wrote the subsetted alignment as a phylip file, but that need
not be- you could write it in Nexus or fasta format as you choose.)

That makes a dictionary file, by default named ``p4DupeSeqRenameDict.py``.
Don't lose it, because it is needed when you restore the deleted duped
taxa to the resulting tree.

It removes all but one of the duplicate sequences, and re-names the
remaining sequence with a name like p4Dupe1 or p4Dupe2.

Then you can do your analysis, which will go faster because you don't
have duped sequences.  At the end, you will generally have a tree.  You
can restore the removed taxa to the tree using the dictionary file made
above.  You would do something like::

     read('myTree.nex')
     t = var.trees[0]
     t.restoreDupeTaxa()
     t.writeNexus(fName='myRestoredTree.nex')

The restored taxa will by default have compound names.
