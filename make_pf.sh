# This might be useful for those who install "in-place" It builds the
# pf.so module, and then puts it in ./p4

# The c-language stuff does not change much over time, so this need
# not be done every hg update.  But occasionally something in the pf
# module changes, and so this build_ext will need to be done again.

rm -f p4/pf.so

python setup.py build_ext -i
