========================
Scripts for common tasks
========================

This is a cook book of reminders, suggestions, and points of departure
for customization.  You can copy and paste these into your own files.
Or you can use the p4 recipes function, ``func.recipes()``, to write new
files.

In these scripts I suggest variations, often by enclosing the variation
in an ``if 1:`` block to turn it on, or an ``if 0:`` block to turn it off.
Hopefully the meaning is clear and the (sometimes mutually exclusive)
possibilities are obvious.



Make a consensus tree, with variations
--------------------------------------

.. literalinclude:: ../Examples/W_recipes/sCon.py

Make a consensus tree, uncomplicated
------------------------------------

.. literalinclude:: ../Examples/W_recipes/sConSimple.py

Calculate a likelihood
----------------------

.. literalinclude:: ../Examples/W_recipes/sLike.py

Calculate likelihood with more than one data partition
------------------------------------------------------

.. literalinclude:: ../Examples/W_recipes/sLikeMultiPart.py

Do an MCMC
----------

.. literalinclude:: ../Examples/W_recipes/sMcmc.py

Do an MCMC with more than one data partition
--------------------------------------------

.. literalinclude:: ../Examples/W_recipes/sMcmcMultiPart.py

Read checkPoints from an MCMC
-----------------------------

.. literalinclude:: ../Examples/W_recipes/sReadMcmcCheckPoints.py

Restart an MCMC
---------------

.. literalinclude:: ../Examples/W_recipes/sRestartMcmc.py

Simulate data
-------------

.. literalinclude:: ../Examples/W_recipes/sSim.py

Simulate data with more than one data partition
-----------------------------------------------

.. literalinclude:: ../Examples/W_recipes/sSimMultiPart.py

Simulate very hetero data
-------------------------

.. literalinclude:: ../Examples/W_recipes/sSimVeryHetero.py

