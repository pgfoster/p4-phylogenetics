from distutils.core import setup
from distutils.extension import Extension
import glob

mySources = glob.glob('*.cpp')
toRemove = [
    'all_pairs_quartet_distance.cpp',
    'all_pairs_triplet_distance.cpp',
    'pairs_quartet_distance.cpp',
    'pairs_triplet_distance.cpp',
    'quartet_dist.cpp',
    'triplet_dist.cpp',
    'rQuartetDist.cpp',
    'rTripletDist.cpp',
    'pyQuartetDist.cpp',
    'pyTripletDist.cpp',
    'test_quartet.cpp',
    'test_triplet.cpp',
    ]

for fN in toRemove:   
    mySources.remove(fN)

setup(name="PackageName",
    ext_modules=[
          Extension("tqdist", mySources,
                    # Adjust the following to be able to find boost_python and, if needed, cblas/blas/atlas.
                    #include_dirs = ["/opt/local/include"],
                    #library_dirs = ["/opt/local/lib"], 
                    # Might not need atlas?  Might be blas rather than cblas?
                    # libraries = ["cblas", "atlas", "boost_pythonXX"]) 
                    libraries = ["boost_python37"]              # this worked on ubuntu 16.04 with Py3.7
                    #libraries = ["boost_python38"]             # this worked on my mac with Py3.8
          )
    ])
