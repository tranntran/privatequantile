functions {
  real n_lpdf(real y, real mu, real sigma2) {
     real ans;
     real ans1;
     
     ans1 = 2*pi()*sigma2;
     ans = log(1/sqrt(ans1))*(-pow((y-mu),2))/(2*sigma2);
     
     return ans;
  }
}
data {
  real mu;         	
  real<lower=0> sigma2;

}
parameters {
  real y;
}

model {
  y ~ n(mu, sigma2);
}
