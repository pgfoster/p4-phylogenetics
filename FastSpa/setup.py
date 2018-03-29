from distutils.core import setup
from distutils.extension import Extension

setup(name="fastSpa",
    ext_modules=[
        Extension("fastspa", ['fastSpa.cpp'],
                  # Adjust the following to be able to find boost-python and, boost-numpy.
                  # On my mac, using Homebrew, I did not need include_dirs or library_dirs, as it was installed
                  # in a usual location.
                  # On my old Ubuntu box, I installed boost "in-place" (or rather I did not "install" it as such) 
                  # and so I had to specify both dirs.
                  include_dirs = ["/usr/local/src/Boost/boost_1_66_0"],
                  library_dirs = ["/usr/local/src/Boost/boost_1_66_0/stage/lib"],
                  libraries = ["boost_python3", "boost_numpy3"]   # needed and worked on both linux and mac
        )
    ])