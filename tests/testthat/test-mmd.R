# X = matrix(rnorm(300*5), nrow=300)
X = amazon[1:100, 1:5]
Y = amazon[50:60, 1:5]
# Y = matrix(rnorm(400*5), nrow=400)
# R.matlab::writeMat("amazon_data.mat",Xt=amazon)

test_that("mmd works", {
  expect_equal(round(fast_mmd(X, Y, 2, 64, seed=1)$unbiased, 3), .283)
})
