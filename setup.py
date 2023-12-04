theLongDescription = """P4 does Bayesian and maximum likelihood phylogenetic analyses on
molecular sequences.  It's specialty is that you can use heterogeneous
models, where the model parameters can differ in different parts of
the tree, or over different parts of the data.
"""

# distutils is deprecated and removal is planned for Python 3.12
from distutils.core import setup, Extension
import sys
import os

if (sys.version_info < (3, 0)):
    print('P4 uses Python3')
    sys.exit(1)

sys.path.insert(0, "p4")

try:
    import numpy
except ImportError:
    print("P4 needs numpy")
    sys.exit()
try:
    import scipy
except ImportError:
    print("P4 needs scipy")
    sys.exit()

############################################################

# Find the location of the gsl library, and its header files, and the
# location of the numpy header files.  If this attempt fails, find
# them yourself, and add them by hand to the following two lists.  Eg
#
#my_include_dirs = ["/my/weird/include"]
#my_lib_dirs = ["/my/weird/lib"]

my_include_dirs = ["/share/apps/include"]
my_lib_dirs = ["/usr/local/lib64", "/share/apps/lib64"]

likelyDirs = [ "/usr",
               "/usr/local",
               "/sw",
               "/opt/local",
               "/opt/homebrew",
               os.path.expanduser('~'),
               os.path.join(os.path.expanduser('~'), ".linuxbrew"),
	       "/home/linuxbrew/.linuxbrew"
               ]

import glob

# Get GSL (GNU Scientific Library) lib and header locations
found_libgsl = False
found_gsl_headers = False

for lpth in my_lib_dirs:
    if os.path.exists(lpth):
        fList = glob.glob(os.path.join(lpth, "libgsl*"))
        #print(fList)
        if fList:
            found_libgsl = True
for ipth in my_include_dirs:
    if os.path.exists(os.path.join(ipth, 'gsl')):
        found_gsl_headers = True
        break
#print("found_libgsl = %s" % found_libgsl)
#print("found_gsl_headers = %s" % found_gsl_headers)

if not found_libgsl or not found_gsl_headers:
    for lD in likelyDirs:
        for libdir in ['lib', 'lib64', 'lib/x86_64-linux-gnu', 'lib/i386-linux-gnu']:
            lpth = os.path.join(lD, libdir)
            if os.path.exists(lpth):
                fList = glob.glob(os.path.join(lpth, "libgsl*"))
                #print(fList)
                if fList:
                    found_libgsl = True
                    ipth = os.path.join(lD, 'include')   # We assume the header is near the lib, but check it next line
                    if os.path.exists(os.path.join(ipth, 'gsl')):
                        found_gsl_headers = True
                        my_include_dirs.append(ipth)
                        my_lib_dirs.append(lpth)
                        break

if not found_libgsl or not found_gsl_headers:
    print("The setup.py script could not find the libgsl or gsl headers.")
    print("So that is not going to work.")
    sys.exit()

hasNumpyHeaders = True
# Get numpy header location.  We have already successfully imported numpy, above.
ipth = numpy.__file__
ipth = os.path.dirname(ipth)
ipth = os.path.join(ipth, "core")
ipth = os.path.join(ipth, "include")
if os.path.exists(ipth):
    if ipth not in my_include_dirs:
        my_include_dirs.append(ipth)
else:
    print("The setup.py script was not able to find the numpy headers.")
    print("So that is not going to work.")
    sys.exit()

pfSources = []
sourceDir = 'Pf'
allFiles = os.listdir(os.path.join(os.curdir, sourceDir))
for f in allFiles:
    if f.endswith('.c'):
        if f.endswith('nexusToken.c'):  # hack, this week
            continue
        pfSources.append(f)
for i in range(len(pfSources)):
    pfSources[i] = os.path.join(sourceDir, pfSources[i])
	

setup(name="p4",
      description="Phylogenetic analysis with heterogeneous models",
      author="Peter Foster",
      author_email="p.foster _at_ nhm ac uk",
      url="https://p4.nhm.ac.uk",
      license="GPL",
      long_description=theLongDescription,
      scripts = ["bin/p4"],
      #data_files= ['share/Examples', 'share/sphinxdoc'],
      packages=["p4", "p4.interactive"],
      ext_package="p4",
      ext_modules=[Extension("pf",
                             pfSources,
                             include_dirs = my_include_dirs,
                             library_dirs = my_lib_dirs,
                             # If you need to link the gsl stuff statically, uncomment the extra_link_args line below,
                             # and adjust the rpath to the location of your gsl libs.
                             #extra_link_args = ['-Wl,-rpath=/home/peter/Secret/lib' ],
                             libraries=["gsl", "gslcblas", "nlopt"])])
    
