# This might be useful for those who install "in-place".  It builds the pf.so
# module, and then puts it in ./p4

# The c-language stuff does not change much over time, so this need
# not be done every update.  But occasionally something in the pf
# module changes, and so this build_ext will need to be done again.
# If you are not sure, it is ok to do it anyway.

#rm -f p4/pf.so
rm -f p4/pf.cpython-36m-darwin.so

# Python2 or Python3?
python3 setup.py build_ext -i

# Default -Wunreachable-code is noisy, because I have such a lot of it.
# CFLAGS="-Wno-unreachable-code" python3 setup.py build_ext -i
