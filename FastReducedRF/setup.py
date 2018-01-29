from distutils.core import setup
from distutils.extension import Extension

setup(name="p4",
    ext_modules=[
        Extension("fastReducedRF", ['fastReducedRF.cpp'],
            # Adjust the following to be able to find pyublas, numpy, and boost.  Maybe need library_dirs as well?
            include_dirs = ["/usr/local/lib/python3.6/site-packages/pyublas/include",
                            "/usr/local/lib/python3.6/site-packages/numpy/core/include",
                            "/usr/local/include"],    # worked on my mac
            # include_dirs = ["/usr/local/lib/python3.4/dist-packages/pyublas/include",
            #                 "/usr/local/lib/python3.4/dist-packages/numpy/core/include",
            #                 ],    # worked on my ubuntu 14.04
                  # libraries = ["boost_python"]   # or "boost_python-mt" if multi-threaded, on ubuntu
                  # libraries = ["boost_python-py34"]   # worked on ubuntu 14.04.  No mt; mt not needed?
                  libraries = ["boost_python3"]   # worked on my mac
        )
    ])
