"""Stuff to ignore when using p3rlcompleter.

The p3rlcompleter is meant to be for the user, to expose the user
interface.  As such, it generally does not want to display all that
dir() can provide.  This file is one way to tell the completer to
ignore.  The other way is with ##Ignore comments in the method doc."""

completionIgnores = [
    'Alignment_logDet',
    'Alignment_muck',
    'Alignment_readWrite',
    'Glitch',
    'Nexus',
    'NodeBranch',
    'Chain_propose1',
    'Chain_propose2',
    'Model',    # Model has no user interface
    'NexusToken',
    'NexusSets',
    'Tree_biRootStuff.py',
    'Tree_fit',
    'Tree_model',
    'Tree_muck',
    'Tree_optSim',
    'Tree_write',
    'Trees_converge',
    'ScaleBar',
    'completionIgnores',
    'func2',
    'p3rlcompleter',
    'Var',
    ]

# The func module imports a lot of stuff for its own purposes, and is
# itself generally imported as "import func", which drags all that
# rubbish along with it.  We don't want to look at it when we use the
# completer.
funcIgnores = [
    'os','sys','re','string','math','cStringIO', 'sets', 'random',
    'var','pf',
    'Sequence',
    'SequenceList',
    'Alignment',
    'Nexus',
    'nextTok',
    'Tree',
    'Node'
    ]
