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
                    libraries = ["cblas", "boost_python"])             # this worked on my mac
    ])
