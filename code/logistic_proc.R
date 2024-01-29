rm(list=ls())
library("rstan")

## Data preparation ------------------------------------------------------------

data.in = read.csv("data_example_logistic.csv", row.names=1)

head(data.in)
n = ncol(data.in)
m = nrow(data.in)-1

times = as.numeric(data.in[1, ])
N = unname(as.matrix(data.in[2:(m+1), ]))

data = list(n = n,
            m = m, 
            t = times,
            N = N,
            f = 0.1)

str(data)

## Stan model ------------------------------------------------------------------

stanmodelcode = "
functions{
  vector odemodel(real t, vector N, vector p){
    // p[1]=r, p[2]=K
    vector[1] dNdt;
    dNdt[1] = p[1]*(1-(N[1]/p[2]))*N[1];
    return dNdt;
  }
}

data{
  int n; // observations
  int m; // replicates
  array[n] real t; 
  array[m,n] int N;
  real f;
}

parameters{
  real<lower=0> r;
  real<lower=0> K;
  real<lower=0> tau;
}

model{
  vector[2] p;
  array[1] vector[1] Nsim; // simulated values, matrix. dim1 = 1 point in time, dim2 = dim_ODE = 1

  // priors (uninformative)
  r ~ lognormal(-2,1);
  K ~ lognormal(9,1);
  tau ~ gamma(2,0.1);

  // parameters for integrator
  p[1] = r;
  p[2] = K;
  
  for(j in 1:m){
    for(i in 1:(n-1)){
      if(N[j,i]>0){
        // integrate ODE
        Nsim = ode_rk45(odemodel, [N[j,i]/f]', t[i], {t[i+1]}, p);
        // likelihood
        N[j,i+1] ~ neg_binomial_2(Nsim[1,1]*f, tau);
      }  
    }
  }
}
"

## Model fitting ---------------------------------------------------------------

# stan options
chains = 3
rstan_options(auto_write = TRUE)
options(mc.cores = chains)
iter   =  4000
warmup =  2000

# initial values for sampling 
init=rep(list(list(r=0.1,
                   K=1e4,
                   tau=10))
         ,chains)

stanmodel = stan_model(model_code=stanmodelcode)

fit = sampling(stanmodel,
               data=data,
               iter=iter,
               warmup=warmup,
               chains=chains,
               init=init
)

## Model diagnostics -----------------------------------------------------------

print(fit, digits=3, pars=c("r","K","tau"), probs=c(0.025, 0.5, 0.975))

library("coda")
samples=As.mcmc.list(fit)
plot(samples[, c("r","K")])

pairs(fit, pars=c("r","K"))

## Posterior predictions -------------------------------------------------------

library("deSolve")
ode.model = function(t,N,p){
  dNdt =  p[1]*N[1]*(1.0-N[1]/p[2])
  return(list(dNdt))
}

post = as.matrix(fit)
n.post = nrow(post)
n.post = 1000

par(mfrow=c(3,2), mar=c(3,3,0,0), oma=c(3,3,1,1))
for(i in 1:m){
  
  N.pred = matrix(NA, nrow=n.post, ncol=length(data$t)-1)
  for(j in 1:n.post){
    for(k in 1:(length(data$t)-1)){
      N.pred[j,k] = as.data.frame(lsoda(y     = c( N = data$N[i,k]/data$f ),
                                        times = data$t[c(k,k+1)],
                                        func  = ode.model,
                                        parms = c(post[j, "r"],
                                                  post[j, "K"])
      ))$N[2]
    }
  }
  N.pred.qs = apply(N.pred, 2, function(x) quantile(x, probs=c(0.025,0.500,0.975), na.rm=TRUE))
  
  plot((data$t[2:n]), (data$N[i,2:n]/data$f),ylim=c(0,1.2e4), type="n", xlab="", ylab="")
  for(k in 1:(n-1)){
    lines(c(data$t[k+1],data$t[k+1]), c(N.pred.qs[1,k],N.pred.qs[3,k]), col="red")
  }
  points((data$t[2:n]), N.pred.qs[2, ], col="red", pch="+", cex=2)
  points((data$t[2:n]), (data$N[i,2:n]/data$f), type="p")
}
mtext("Time [h]", side=1, outer=TRUE, line=1)
mtext("Abundance [ind]", side=2, outer=TRUE, line=1)



