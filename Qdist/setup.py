from distutils.core import setup
from distutils.extension import Extension
import glob

mySources = glob.glob('*.cpp')
toRemove = [
    'main.cpp',
    'testMatrix.cpp',
    'testQDist.cpp'
    ]
for fN in toRemove:   
    mySources.remove(fN)

setup(name="PackageName",
    ext_modules=[
          Extension("scqdist", mySources,
                    # Adjust the following to be able to find boost root and cblas/blas/atlas.
                    #include_dirs = ["/opt/local/include"],
                    #library_dirs = ["/opt/local/lib"],
                    # Might not need atlas?  Might be blas rather than cblas?
                    # libraries = ["cblas", "atlas", "boost_python"]) 
                    # libraries = ["cblas", "boost_python"]             # this worked on ubuntu 16.04 with Py2
                    # libraries = ["cblas", "boost_python-py34"]         # this worked on ubuntu 16.04 with Py3
                    # libraries = ["cblas", "boost_python"]             # this worked on my mac with Py2
                    libraries = ["cblas", "boost_python36"]             # this worked on my mac with Py3.6
          )
    ])
