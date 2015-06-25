theLongDescription = """P4 does Bayesian and maximum likelihood phylogenetic analyses on
molecular sequences.  It's specialty is that you can use heterogeneous
models, where the model parameters can differ in different parts of
the tree, or over different parts of the data.
"""

from distutils.core import setup, Extension
from distutils.command.install_data import install_data
from distutils.command.install_lib import install_lib
from distutils.command.install_scripts import install_scripts
import sys,os,shutil
sys.path.append(os.path.join(os.getcwd(),'p4'))
from version import versionString

# Check for version 2, then for version 2.4.  Using sys.version_info
# is more convenient than sys.version, but it is only in v2.
majorVersion = int(sys.version[0])
if majorVersion < 2:
    print 'p4 wants Python version 2.4 or better.'
    sys.exit()
if sys.version_info[1] < 4:
    print 'p4 wants Python version 2.4 or better.'
    sys.exit()

usePfAndNumpy = True
hasNumpy = False
hasNumpyHeaders = False
try:
    import numpy
    hasNumpy = True
except ImportError:
    usePfAndNumpy = False
try:
    import readline
except ImportError:
    print "P4 wants a Python built with 'readline', and your Python does not appear to have it."
    sys.exit()

############################################################

# Find the location of the gsl library, and its header files, and the
# location of the numpy header files.  If this attempt fails, find
# them yourself, and add them by hand to the following two lists.  Eg
#
#my_include_dirs = ["/my/weird/include"]
#my_lib_dirs = ["/my/weird/lib"]

my_include_dirs = []
my_lib_dirs = []

likelyDirs = [ "/usr",
               "/usr/local",
               "/sw",
               "/opt/local",
               os.path.expanduser('~')
               ]

import glob

# Get GSL (GNU Scientific Library) lib and header locations
found_libgsl = False
found_gsl_headers = False

for lpth in my_lib_dirs:
    if os.path.exists(lpth):
        fList = glob.glob(os.path.join(lpth, "libgsl*"))
        #print fList
        if fList:
            found_libgsl = True
for ipth in my_include_dirs:
    if os.path.exists(os.path.join(ipth, 'gsl')):
        found_gsl_headers = True
        break
#print "found_libgsl = %s" % found_libgsl
#print "found_gsl_headers = %s" % found_gsl_headers

if not found_libgsl or not found_gsl_headers:
    for lD in likelyDirs:
        lpth = os.path.join(lD, 'lib')
        if os.path.exists(lpth):
            fList = glob.glob(os.path.join(lpth, "libgsl*"))
            #print fList
            if fList:
                found_libgsl = True
                ipth = os.path.join(lD, 'include')   # We assume the header is near the lib, but check it next line
                if os.path.exists(os.path.join(ipth, 'gsl')):
                    found_gsl_headers = True
                    my_include_dirs.append(ipth)
                    my_lib_dirs.append(lpth)
                    break

if not found_libgsl or not found_gsl_headers:
    usePfAndNumpy = False

if hasNumpy:
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
        hasNumpyHeaders = False
        usePfAndNumpy = False

#print "found_libgsl = %s" % found_libgsl
#print "found_gsl_headers = %s" % found_gsl_headers
#print "my_include_dirs", my_include_dirs
#print "my_lib_dirs", my_lib_dirs
#sys.exit()

########################################################################################
#usePfAndNumpy=False
if usePfAndNumpy:
    pfSources = []
    sourceDir = 'Pf'
    allFiles = os.listdir(os.path.join(os.curdir, sourceDir))
    for f in allFiles:
        if f.endswith('.c'):
            pfSources.append(f)
    for i in range(len(pfSources)):
        pfSources[i] = os.path.join(sourceDir, pfSources[i])
	
instFileName = "installation.py"
if os.path.exists(instFileName):
    print "The file '%s' exists, but it should not. Remove or re-name it." % instFileName
    sys.exit()

# For use in shutil.copytree(), below
def my_ignore(adir, filenames):
    return [fN for fN in filenames if fN.startswith(".svn")]

class P4_install_data(install_data):
    def run(self):
        # This is only called when installing, I think.
        print 'P4_install_data.run() here, sitting in for install_data.run()'
        # data_files is a list of one string, defined below in setup(), ['share/doc/p4-0.xx']
        print 'data_files = %s' % self.data_files   
        print 'install_dir = %s' % self.install_dir # eg /usr
        #fList = self.copy_tree('share', os.path.join(self.install_dir, self.data_files[0]))
        # Don't copy .svn directories
        shutil.copytree('share', 
                        os.path.join(self.install_dir, self.data_files[0]), 
                        ignore=my_ignore)

        
        mySharePath = os.path.join(self.install_dir, self.data_files[0])
        mySphinxIndexPath = os.path.join(mySharePath, 'sphinxdoc/_build/html/index.html')
        myExamplesPath = os.path.join(mySharePath, 'Examples')

        try:
            loc = {}
            execfile("%s" % instFileName, {}, loc)  # get the p4_lib_dir
            p4_lib_dir = loc['p4LibDir']
            instFile = file(instFileName, 'a')
            instFile.write("p4DocDir = '%s'\n" % mySharePath)
            instFile.write("p4SphinxIndexPath = '%s'\n" % mySphinxIndexPath)
            instFile.write("p4ExamplesDir = '%s'\n" % myExamplesPath)
            instFile.close()
            os.system("cp %s %s" % (instFileName, p4_lib_dir))
            from py_compile import compile
            compile(os.path.join(p4_lib_dir, instFileName))
            os.system("rm -f %s" % instFileName)
        except IOError:
            print "The file '%s' cannot be found." % instFileName


class P4_install_lib(install_lib):
    def run(self):
        print "P4_install_lib()"
        print "self.install_dir = %s" % self.install_dir
        instFile = file(instFileName, 'w')
        instFile.write("p4LibDir = '%s'\n" % os.path.join(self.install_dir, 'p4'))
        instFile.close()
        install_lib.run(self)


class P4_install_scripts(install_scripts):
    def run(self):
        print "P4_install_scripts()"
        print "self.install_dir = %s" % self.install_dir
        instFile = file(instFileName, 'a')
        instFile.write("p4ScriptPath = '%s'\n" % os.path.join(self.install_dir, 'p4'))
        instFile.close()
        install_scripts.run(self)
        
            
        
if usePfAndNumpy:
    setup(name="p4",
          #version=P4_VERSION,
          version=versionString,
          description="Phylogenetic analysis with heterogeneous models",
          author="Peter Foster",
          author_email="p.foster _at_ nhm ac uk",
          url="http://www.bmnh.org/~pf/p4.html",
          license="GPL",
          long_description=theLongDescription,
          cmdclass = {'install_data':P4_install_data,
                      'install_lib':P4_install_lib,
                      'install_scripts':P4_install_scripts},
          scripts = ["bin/p4"],
          #data_files= ['share/doc/p4-%s' % versionString],
          data_files= ['share/doc/p4'],
          packages=["p4"],
          ext_package="p4",
          ext_modules=[Extension("pf",
                                 pfSources,
                                 include_dirs = my_include_dirs,
                                 library_dirs = my_lib_dirs,
                                 # If you need to link the gsl stuff statically, uncomment the extra_link_args line below,
                                 # and adjust the rpath to the location of your gsl libs.
                                 #extra_link_args = ['-Wl,-rpath=/home/peter/Secret/lib' ],
                                 libraries=["gsl", "gslcblas"])])
else:
    setup(name="p4",
          #version=P4_VERSION,
          version=versionString,
          description="Phylogenetic analysis with heterogeneous models",
          author="Peter Foster",
          author_email="p.foster _at_ nhm ac uk",
          url="http://www.bmnh.org/~pf/p4.html",
          license="GPL",
          long_description=theLongDescription,
          cmdclass = {'install_data':P4_install_data,
                      'install_lib':P4_install_lib,
                      'install_scripts':P4_install_scripts},
          scripts = ["bin/p4"],
          #data_files= ['share/doc/p4-%s' % versionString],
          data_files= ['share/doc/p4'],
          packages=["p4"],
          ext_package="p4")
    
if not usePfAndNumpy:
    print "\n\n=============================================================\n"
    print "Some parts of p4 involving likelihood calculations, mcmc, "
    print "etc, cannot be installed.  You may not really need those parts, "
    print "in which case fine.\n"

    print "However, if you do need those parts, they have pre-requisites "
    print "and it seems that those pre-requisites were not found.\n"
    
    if not hasNumpy:
        print "The 'numpy' module does not appear to be installed."
    if not hasNumpyHeaders:
        print "The numpy headers could not be found."
        
    if not found_libgsl:
        print "The GSL (Gnu Scientific Library) libgsl could not be found."
        print "(If you want to install it, "
        print "if you use a package manager, you may need to install both gsl"
        print " and gsl-dev, for the header files.)"
        sys.exit()
    if found_libgsl and not found_gsl_headers:
        print "The GSL (Gnu Scientific Library) libgsl was found."
        print "However, the gsl header files could not be found."
        print "If you use a package manager, maybe you need to install 'gsl-dev'?"

    print "\n(At least that is what setup.py thinks.  But it is not very smart."
    print "If it is wrong, see and modify 'setup.py' accordingly.)\n"

################################################################################

