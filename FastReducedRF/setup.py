from distutils.core import setup
from distutils.extension import Extension

setup(name="p4",
    ext_modules=[
        Extension("fastReducedRF", ['fastReducedRF.cpp'],
            # Adjust the following to be able to find pyublas, numpy, and boost.  Maybe need library_dirs as well?
            include_dirs = ["/usr/local/lib/python2.7/site-packages/pyublas/include",
              "/usr/local/lib/python2.7/site-packages/numpy/core/include",
                            "/usr/local/include"],
                  libraries = ["boost_python"])   # or "boost_python-mt" if multi-threaded, on ubuntu
    ])
