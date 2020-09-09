from distutils.core import setup
from distutils.extension import Extension

myFiles = """pyquartet_dist.cpp
  HDT.cpp
  HDTCountingCCToC.cpp
  HDTCountingCTransform.cpp
  HDTCountingG.cpp
  HDTCountingIGToC.cpp
  HDTCountingLeaf.cpp
  HDTFactory.cpp
  HDTListUtils.cpp
  NewickParser.cpp
  RootedTree.cpp
  RootedTreeFactory.cpp
  AbstractDistanceCalculator.cpp
  QuartetDistanceCalculator.cpp
  int_stuff.cpp""".split()



setup(ext_modules=[
        Extension("pytqdist", myFiles,
                  # Adjust the following to be able to find boost-python and, boost-numpy.
                  # On my mac, using Homebrew, I did not need include_dirs or library_dirs, as it was installed
                  # in a "usual" location.
                  # On my old Ubuntu box, I installed boost "in-place" (or rather I did not "install" it as such) 
                  # and so I had to specify both dirs.
                  #include_dirs = ["/usr/local/src/Boost/boost_1_66_0"],
                  #library_dirs = ["/usr/local/src/Boost/boost_1_66_0/stage/lib"],
                  define_macros=[('quartetsToo', None)],

                  # Adjust and use this next line if you want to build with
                  # statically linked libs.  You might want to do this if your
                  # libs are in a non-standard location, and (although it
                  # compiles well) when you run it it complains that it cannot
                  # find the boost libs.  A related reason might be if you are
                  # using more than one boost.
                  #extra_link_args = ['-Wl,-rpath=/usr/local/src/Boost/boost_1_66_0/stage/lib' ],  # static build
                  libraries = ["boost_python38"]
        )
    ])
