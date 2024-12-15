# TODO: enable C++ tests. Currently they only work with devtools::test()
# the problem seems to be the test-mqs and test-genotype-masking files import 
# test_utils. change this pattern
test_that("C++ tests", {
  run_cpp_tests("mmer")
})
