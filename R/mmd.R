
#' @export
pair_complement_mmd = function(X, Y,
                               sigvec = 10^seq(-2, 2, by = 0.2),
                               agg_func = mean,
                               max_npairs = Inf,
                               nBasis = 2^6,
                               seed = NULL,
                               ncores = 4) {
  p = ncol(X)
  all_indices = 1:p
  if (choose(p,2) > max_npairs) {
    combs = matrix_to_list_of_rows( unique(matrix(sample(p, 2*max_npairs, T), ncol=2)) )
  } else {
    combs = combn(all_indices, 2, simplify = F)
  }

  curr_time = tic()
  printf('  [Computing pair-comp. MMD (nBasis = %d, nPairs = %d) ... ', nBasis, length(combs))
  out = parallel::mclapply(combs,
  # lapply(combn(all_indices, 2, simplify = F),
                     function(ind) {
                        ind =  setdiff(all_indices, ind)
                        agg_func(fast_mmd(X[, ind], Y[, ind], sigvec, nBasis, seed)$unbiased)
                     },
                     mc.cores = ncores) # %>% unlist()
  values = unlist(out)

  printf('%2.1f (s).  value = %2.3f]\n', toc(curr_time), mean(values))

  return(values)
  # doParallel::registerDoParallel(ncores)
  # p = ncol(X)
  # C = combn(1:p, 2, simplify = F)
  # nC = length(C)
  #
  # `%dopar%` <- foreach::`%dopar%`
  # foreach::foreach(i=1:nC, .combine = "c") %dopar% {
  #    ind = setdiff(1:p, C[[i]])
  #    mean( fast_mmd(X[, ind], Y[, ind], sigvec, nBasis, seed)$unbiased )
  # }
}


pair_complement_mmd2 = function(X, Y,
                               sigvec = 10^seq(-2, 2, by = 0.2),
                               nBasis = 2^6,
                               seed = NULL,
                               nWorkers = 4) {

  all_indices = 1:ncol(X)
  out = lapply(combn(all_indices, 2, simplify = F),
                     function(ind) {
                       ind =  setdiff(all_indices, ind)
                       mean(fast_mmd(X[, ind], Y[, ind], sigvec, nBasis, seed)$unbiased)
                     })
  unlist(out)
  # p = ncol(X)
  # C = combn(1:p, 2, simplify = F)
  # nC = length(C)
  # mean_mmd = vector("double", nC)
  # for (ii in 1:nC) {
  #   ind = setdiff(1:p, C[[ii]])
  #   mean_mmd[ii] = mean( fast_mmd(X[, ind], Y[, ind], sigvec, nBasis, seed)$unbiased )
  # }
  # mean_mmd
}



l2norm_squared = function(X) Re(sum(Conj(t(X)) %*% X))

fast_mmd = function(X, Y,
                    sigvec = 10^seq(-2, 2, by = 0.2),
                    nBasis = 2^6, seed = NULL) {
  # The function use the Random Kitchen Sink idea (i.e., random Fourier
  # features) to compute the MMD based on a Gaussian kernel.
  # nBasis: Number of Fourier basis vectors

  if (!is.null(seed)) set.seed(seed)
  # rng('default');
  k0 = 1 # k0 = K(0,0)
  d = ncol(X)
  if (ncol(Y) != d) stop('# of columns of X and Y should match.')

  n = nrow(X)
  m = nrow(Y)
  N = nBasis

  Z = matrix(rnorm(N*d), nrow=N) # randn(N,d);
  # ZXt = Z %*% Matrix::t(X) # N x n
  # ZYt = Z %*% Matrix::t(Y) # N x m
  # ZXt = Matrix::t( X %*% t(Z) ) # N x n
  # ZYt = Matrix::t( Y %*% t(Z) ) # N x m
  ZXt = Z %*% t(as.matrix(X))
  ZYt = Z %*% t(as.matrix(Y))
  phi = function(sigma) rowMeans( exp(1i*ZXt/sigma) / sqrt(N))
  psi = function(sigma) rowMeans( exp(1i*ZYt/sigma) / sqrt(N))

  K = length(sigvec)
  biased = unbiased = rep(0, K)
  for (k in 1:K) {
    sigma = sigvec[k]
    u = phi(sigma)  # N x 1
    v = psi(sigma)  #  N x 1
    temp = l2norm_squared(u - v)
    biased[k] = sqrt(temp); #  norm(phi(sigma) - psi(sigma));
    unbiased[k] = temp + l2norm_squared(u)/(n-1) + l2norm_squared(v)/(m-1) - k0*(n+m-2)/(n-1)/(m-1)
    unbiased[k] = sqrt(max(unbiased[k],0))
  }
  list(biased = biased, unbiased = unbiased)
}

# l2norm = function(x) sqrt(sum(x*x))
# # The port to R does not work as in MATLAB; something is not right
# MMDFourierFeature = function(xPos, xNeg,
#                              allSgm = 10^seq(-2,2,by = 0.2),
#                              nBasis = 2^6, MAX_SIZE = 1e7) {
#   # Approximate the MMD by Random Fourier Features (Random Kitchen Sinks)
#   #  with block process
#   # The kernel is Gaussian: k(x,y) = exp(-||x-y||_2^2 / (2*sigma^2)).
#   # Reference:
#   #  [1] Ji Zhao, Deyu Meng. Ensemble of Circular Discrepancy for Efficient Two-Sample Test.
#   #      NIPS Workshop on Randomized Methods for Machine Learning (RMML2013), 2013.
#   # [2] Ji Zhao, Deyu Meng. FastMMD: Ensemble of Circular Discrepancy for Efficient Two-Sample Test.
#   #      Neural Computation, 2015.
#   # Input:
#   #  xPos, xNeg: two sample sets
#   #  allSgm: bandwidth parameter for Gaussian kernel, scale or vector
#   #  nBasis: number of basis for approximating p(w), see our paper.
#   # Output
#   #  d1: estimate of biased MMD
#   #  d2: estimate of unbiased MMD
#
#   # Ji Zhao@CMU
#   # zhaoji84@gmail.com
#   # 02/18/2014
#
#
#   k0 = 1 # K(0,0)=1
#   nDim = ncol(xPos)
#   Nsig = length(allSgm)
#   d1 = d2 = rep(0, Nsig)
#
#   nPos = nrow(xPos)
#   nNeg = nrow(xNeg)
#   xPos = as.matrix(xPos)
#   xNeg = as.matrix(xNeg)
#
#   bsz = max(ceiling(MAX_SIZE/nDim), 1)
#
#   nBlock1 = ceiling(nBasis/bsz); # for W
#   nBlock2 = ceiling(nPos/bsz); # for dataset 1
#   nBlock3 = ceiling(nNeg/bsz); # for dataset 2
#   cat(nBlock1)
#   cat(nBlock2)
#   cat(nBlock3)
#
#   phiPos = phiNeg = matrix(0, nrow=nBasis*2, ncol = Nsig)
#   #phiPos = zeros(nBasis*2, numel(allSgm));
#   #phiNeg = zeros(nBasis*2, numel(allSgm));
#
#   for (ii in 1:nBlock1) {
#     i1 = (ii-1)*bsz + 1
#     i2 = min(ii*bsz, nBasis)
#     i1 = i1*2-1 # need double space to store [cos sin]
#     i2 = i2*2
#     tmp = min(bsz, nBasis-(ii-1)*bsz)
#     W = matrix(rnorm(tmp*nDim), nrow = tmp) # randn(tmp, nDim);
#     for (jj in 1:nBlock2) {
#       j1 = (jj-1)*bsz + 1
#       j2 = min(jj*bsz, nPos)
#       tt = W %*% t(xPos[j1:j2, ]) # W*xPos(j1:j2, :)';
#       for (kk in 1:Nsig) {
#         sgm = allSgm[kk];
#         t1 = tt/sgm;
#         t2 = cbind(cos(t1), sin(t1)) # [cos(t1); sin(t1)];
#         phiPos[i1:i2, kk] = phiPos[i1:i2, kk] + rowSums(t2)
#       }
#     }
#     for (jj in 1:nBlock3) {
#       j1 = (jj-1)*bsz + 1
#       j2 = min(jj*bsz, nNeg)
#       tt = W %*% t(xNeg[j1:j2, ])
#       for (kk in 1:Nsig) {
#         sgm = allSgm[kk]
#         t1 = tt/sgm
#         t2 = cbind(cos(t1), sin(t1))
#         phiNeg[i1:i2, kk] = phiNeg[i1:i2, kk] + rowSums(t2)
#       }
#     }
#   }
#
#   phiPos = phiPos/nPos * nBasis^(-1/2)
#   phiNeg = phiNeg/nNeg * nBasis^(-1/2)
#
#   for (ii in 1:Nsig) {
#     d1[ii] = l2norm(phiPos[ , ii] - phiNeg[, ii])
#     tt = d1[ii]^2 + l2norm(phiPos[ , ii])^2/(nPos-1) +
#       l2norm(phiNeg[ , ii])^2/(nNeg-1) -
#       k0*(nPos+nNeg-2)/(nPos-1)/(nNeg-1)
#     d2[ii] = (max(tt,0))^(1/2)
#   }
#
#   list(d1 = d1, d2 = d2)
# }
#




