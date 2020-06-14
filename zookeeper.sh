pip3 install --user numpy scipy pandas

Rscript travis_script.r
R CMD build .
R CMD check *tar.gz --no-manual
