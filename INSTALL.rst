============
Installation
============

The source code is hosted at `<https://github.com/pgfoster/p4-phylogenetics>`_

P4 needs Python 3; it no longer works with Python 2.

I've installed it on Linux and Mac OS X.  In either case, you need to
have the basic C-language programming tools, including a compiler,
libraries, headers, and so on.   

**Preparations for the full install on the Mac**

You probably want to install `Homebrew <http://brew.sh>`_ and install Python3 from there.
The command-line tools for development need to be installed for Homebrew to work, but 
the last time I installed Homebrew on my mac it did the installation of those command-line tools itself.  
If that does not work, you may be able to install them via::

    xcode-select --install

You don't need the full Xcode for this (but if you want it, it is available from the Apple App Store).

You will need the ``gsl`` library (Gnu
Scientific Library), and the ``nlopt`` library.  For this, `Homebrew <http://brew.sh>`_ is recommended.

You can install ``scipy`` with ``pip3``; it will install ``numpy`` as well.  You will also need ``bitarray``.


 
**Preparations for the full install on Ubuntu Linux**

If the installation process complains about lack of Python.h, then you
need what on Ubuntu would be called 'python-dev'. 

I have recently installed p4 on Ubuntu 16.04, and had to::

    sudo apt-get install libgsl-dev
    sudo apt-get install libnlopt-dev
    sudo apt-get install python3-dev

You can use ``pip3`` to install ``scipy`` and ``bitarray``.

And if you want to use the GUI tree-drawing::

    sudo apt-get install python-tk

Presumably other Ubuntu versions will be similar or identical.


Installing p4 using setup.py
============================

This is the usual way that Python packages are installed; however if you plan on
keeping up with the git repo you it may be easier to install p4 in-place as
described in another section below.

If you are upgrading and you want a clean install, you can un-install the
previous version with the p4 func.uninstall() function --- see below.

Presumably you are starting with the file ``p4-x.y.z.tar.gz``.  Unpack it in
your favourite source directory.  In the newly-created directory note
the file ``setup.py``.  That file controls the build and installation.  

**Simple install**

First you can build it, without installing it, by saying::

    python3 setup.py build

(no need to be root or use sudo for the above step)

After building it, you then install it.  The default location for
installation is where python libraries are installed, and you as
JoeUser may not have file-writing permission to put files there, so
you may need to be root or use sudo for the next step.  Eg if you sudo
it, you can say::

    sudo python3 setup.py install


**Installation variations**

To get a description of the various installation options, do::
    
    python3 setup.py install --help. 

To install it in your home directory, say::

    python3 setup.py install --prefix=~

Or some unusual place in your home directory::

    python3 setup.py install --prefix=~/Unusual

If you install it in your home directory, 
    
- there is no need to be root or to use sudo

- if you do this you may need to setenv your ``PYTHONPATH`` to eg
  ``~/lib/python``.  Eg in your ``~/.bash_profile`` or ``~/.profile`` you can put the
  line::

      export PYTHONPATH=$HOME/lib/python

  (depending on where your P4 lib directory is, and what it is called), or
  you can add ::

      export PYTHONPATH=$PYTHONPATH:$HOME/lib/python

  if you already have a ``PYTHONPATH`` defined.

- you may also need to set your ``PATH`` environment variable to
  include ``~/bin``.  In many cases this will already be done, but if it is
  not, and you are using the bash shell, you can do something like::

      PATH=$PATH:$HOME/bin

  and then as usual after all your paths have been set, you should have a line that says::

      export PATH



**Where things go**


The default installation location has a "root", which might be ``/usr`` or
``/usr/local``, or your home directory.  These three things are installed:

    1.  **The p4 package.**          Goes where 3rd party packages go
        Eg ``/usr/local/lib/python3.6/site-packages``

    2.  **The p4 script.**           Goes somewhere in your path
        Eg ``/usr/local/bin``

    3.  **The examples.**             Goes in a share/doc directory
        Eg ``/usr/local/share/doc``

The default location for installation of the modules is something like::

    /usr/lib/python3.6/site-packages

or::

    ~/lib/python3.6

depending on the "root" of the installation, of course.

The default location for the script p4 is something like::

    /usr/local/bin

The default location for the examples is something like::

    /usr/local/share/doc/p4-1.x.y/Examples


Installing it in-place
======================

My fave way of using the git version of p4 is to install it in-place rather than
installing it with ``setup.py``.  The advantage is that it makes it easier to
keep up with the changes made to the git repo.  The first thing would be to
clone it from GitHub.  After that, you need to make it usable.


To make it usable in-place, you need to do three things, which in overview are

1. Add the p4 git directory, eg ``/usr/local/src/P4Git`` to your ``PYTHONPATH``

2. Add the p4 git bin directory, eg ``/usr/local/src/P4Git/bin`` to your ``PATH``

3. Build the ``pf`` module, installing it in-place

Now look at those three steps in detail.
For example if you install it in your home directory, to add the p4
git directory to your ``PYTHONPATH``, you might add something like the
following line to your ``~/.profile`` or ``~/.bash_profile``::

  export PYTHONPATH=$HOME/src/P4Git

(depending on where your P4 lib directory is, and what it is called), or
you can add ::

  export PYTHONPATH=$PYTHONPATH:$HOME/src/P4Git

if you already have a ``PYTHONPATH`` defined.

The second thing you will want to do is to add the location of the p4
script to your ``PATH``.  Similar to adjusting the ``PYTHONPATH``
above, you can add a line like this to your  ``~/.profile`` or ``~/.bash_profile``::

  export PATH=$PATH:$HOME/src/P4Git/bin

depending on where your P4 git directory is, and what it is called.

To build the ``pf`` module, say::

   python3 setup.py build_ext -i

It might actually work.  If it doesn't, note the error messages that
flew by.  The earliest error message is usually a clue.


**Updating from git**


The motivation for installing it in-place is that it makes it easy to
update.  Generally all you need to do is to go to the p4 git directory
and say::

  git pull

That is usually sufficient.  

Occasionally there may have been changes to the C-language code in the ``pf``
module.  If that is the case (would you be able to see those files as they are
updated?), and you use the ``pf`` module then you would need to do::

  python3 setup.py build_ext -i

You would also need to do that when you install it in-place for the
first time, or if you make any changes to the C-language code
yourself.  If you are not sure rebuilding is needed, it's OK to do it anyway.


Installing scqdist, the sub-cubic quartet distance module
=========================================================

See the directory Qdist in the source, with its own instructions.


To see if it works
==================

If, in your shell, you are still in the same directory that you built it from,
go to some other directory, or the following test will not work.  Even better,
use a new shell.

To see if you can load the package, start up python3 and then::

    import p4

To see if the p4 script works, say (perhaps from a new terminal) to
your shell (not in interactive python)::

    p4 --help

(Once it gets installed, if everything went perfectly and it still
does not work, try it in a new shell, or maybe even restart your
terminal program to refresh your PATH and PYTHONPATH.)



Deinstallation
==============

There is a func.uninstall() function, which may work.  You may need to
run it as root, or use sudo.

If that does not work, then recall that things get installed in 3
places.  Search out the Python package, the p4 script, and the
examples.



 
If you want to statically link your gsl libs
============================================

For those who may not want to do the usual dynamic linking of gsl
libs, it is possible to statically link the gsl libs to the pf.so
module when you build it.  See the ``setup.py``
file, and uncomment and adjust the ``extra_link_args`` line.



