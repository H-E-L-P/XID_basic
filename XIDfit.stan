//Full Bayesian inference fit on XID sims
data {
  int<lower=0> npix;//number of pixels
  int<lower=0> nsrc;//number of sources
  vector[npix] data;//flattened map
  vector[npix] sigma;//flattened uncertianty map (assuming no covariance between pixels)
  matrix[npix,nsrc] A;//image matrix
}
parameters {
  vector[nsrc] src_f;//source vector
}

model {
  data ~ normal(A*src_f,sigma)
    }
