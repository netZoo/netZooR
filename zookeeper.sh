pip3 install --user numpy scipy pandas

R CMD build .
R CMD check *tar.gz --no-manual
