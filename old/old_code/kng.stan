functions {
  real kng_lpdf(matrix beta, int N, int dim, real eps, real tau,
 		matrix y, matrix x, matrix sumX) {
     matrix[N, 1] left;
     matrix[N, 1] less;
     real ans1;
     matrix[1, 1] ans2;
     real ans;

     left = y - x * beta;
     for (i in 1:N){
	less[i,1] = left[i,1] <= 0;
     }
     ans1 = max(fabs(-tau*sumX + (less)' * x));
     ans2 = (beta)' * beta;
     ans = -(eps/2) * ans1 / ((1-tau)*2*max(fabs(x))) - 0.005*ans2[1, 1];
     
     return ans;
  }
}
data {
  int<lower=1> N;         	// sample size
  int<lower=1> dim;
  real<lower=0> eps;
  real<lower=0, upper=1> tau;
  matrix[N, 1] y;     	// add prior
  matrix[N, dim] x; 		// 
  matrix[1, dim] sumX;
}
parameters {
   matrix<lower=0>[dim, 1] beta;      // coefficient matrix
}
model {
  beta~ kng(N, dim, eps, tau, y, x, sumX);
}
