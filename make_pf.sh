# This might be useful for those who install "in-place".  It builds the pf.so
# module, and then puts it in ./p4

# The c-language stuff does not change much over time, so this need
# not be done every update.  But occasionally something in the pf
# module changes, and so this build_ext will need to be done again.
# If you are not sure, it is ok to do it anyway.

# Might be python 36 or python 37
rm -f p4/pf.cpython-3?m-darwin.so

python3 setup.py build_ext -i

