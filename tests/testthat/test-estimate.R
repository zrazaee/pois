
Xt = order_cols_by_freq(amazon, 10)
cat(sprintf('sparsity = %3.1f%%',100*(1-Matrix::nnzero(Xt)/prod(dim(Xt)))))

out <- fit_pois(Xt, solver="global", method="glmnet")

test_that("pois_fit works", {
  expect_equal(round(norm(out$gam)), 1)
})

