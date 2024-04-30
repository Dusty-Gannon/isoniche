

data{

  int N;                        // number of individuals
  int P[2];                     // number of explanatory variables for mean
  int K[2];                     // number of explanatory variables for sds
  int J;                        // number of explanatory variables for correlation
  matrix[N, P[1] + P[2]] X;     // design matrix for means
  matrix[N, K[1] + K[2]] Z;     // design matrix for sds
  matrix[N, J] G;               // design matrix for correlation
  matrix[N, 2] y;               // response vectors

}

parameters{

  vector[P[1]] beta_1;
  vector[P[2]] beta_2;
  vector[K[1]] zeta_1;
  vector[K[2]] zeta_2;
  vector[J] gamma;

}

model{

  // Some additional object definitions
  matrix[P[1] + P[2], 2] B;
  matrix[K[1] + K[2], 2] Zeta;
  matrix[N, 2] mu;
  matrix[N, 2] sigma;
  vector[N] rho;
  matrix[2, 2] Omega[N];
  matrix[2, 2] Sigma[N];

  // construct B matrix
  B[1:P[1], 1] = beta_1;
  B[(P[1] + 1):(P[1] + P[2]), 1] = rep_vector(0, P[2]);
  B[1:P[1], 2] = rep_vector(0, P[1]);
  B[(P[1] + 1):(P[1] + P[2]), 2] = beta_2;

  // construct Zeta matrix
  Zeta[1:K[1], 1] = zeta_1;
  Zeta[(K[1] + 1):(K[1] + K[2]), 1] = rep_vector(0, K[2]);
  Zeta[1:K[1], 2] = rep_vector(0, K[1]);
  Zeta[(K[1] + 1):(K[1] + K[2]), 2] = zeta_2;

  mu = X * B;
  sigma = exp(Z * Zeta);
  for(i in 1:N){
    rho[i] = 2 * inv_logit(G[i, ] * gamma) - 1;
    Omega[i] = diag_matrix(rep_vector(1, 2));
    Omega[i][1, 2] = rho[i];
    Omega[i][2, 1] = rho[i];
  }

  // construct covariance matrices
  for(i in 1:N){
    Sigma[i] = diag_matrix(sigma[i, ]') * Omega[i] * diag_matrix(sigma[i, ]');
  }

  // priors
  beta_1 ~ normal(0, 5);
  beta_2 ~ normal(0, 5);
  zeta_1 ~ std_normal();
  zeta_2 ~ std_normal();
  gamma ~ std_normal();

  // likelihood
  for(i in 1:N){
    y[i, ]' ~ multi_normal(mu[i, ]', Sigma[i]);
  }

}
