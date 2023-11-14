library(testthat)
library(sars)

test_check("sars")
data("galap")


#devtools::test("sars")
#https://github.com/r-lib/testthat/issues/257

#don't forget to manually un-hash and run some of the GDM, threshold and
#test_crit_method (just the last one for beta p) tests. GDM as issues with 
#having BAT as a dependency, and threshold and betap due
#to speed and running time on Circle CI. Also one sar_habitat test hashed
#out due to needing external dataset.

#also, if making lots of commits in short amount of time, hash out
#code coverage from the circle CI yml file to avoid it do every time.
