from distutils.core import setup
from distutils.extension import Extension

setup(name="p4",
    ext_modules=[
        Extension("fastReducedRF", ['fastReducedRF.cpp'],
            # Adjust the following to be able to find pyublas, numpy, and boost.  Maybe need library_dirs as well?
            #include_dirs = ["/usr/local/lib/python3.8/site-packages/pyublas/include",
            #                "/usr/local/lib/python3.8/site-packages/numpy/core/include",
            #                "/usr/local/include"],    # worked on my mac
            include_dirs = ["/usr/local/lib/python3.7/site-packages/pyublas/include",
                            "/usr/local/lib/python3.7/site-packages/numpy/core/include",
                            "/usr/local/include"],    # worked on my ubuntu 16.04
                  libraries = ["boost_python37"]   # worked on ubuntu 16.04.  No mt; mt not needed?
                  #libraries = ["boost_python38"]   # worked on my mac
        )
    ])
