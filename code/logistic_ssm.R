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
  real<lower=0> sigma;
  real<lower=0> tau;
  matrix<lower=0>[m,n] Ntrue; // states in times t
}

model{
  vector[2] p;
  array[1] vector[1] Nsim; // simulated values, matrix. dim1 = 1 point in time, dim2 = dim_ODE = 1
  
  // priors 
  r ~ lognormal(-2,1);
  K ~ lognormal(9,1);
  sigma ~ gamma(2,1);
  tau ~ gamma(2,0.1);
  
  for(i in 1:m){
    Ntrue[i,1] ~ normal(0,100);
  }

  // parameters for integrator
  p[1] = r;
  p[2] = K;
  
  // process equations
  for(j in 1:m){
    for(i in 1:(n-1)){
      // integrate ODE
      Nsim = ode_rk45(odemodel, [Ntrue[j,i]]', t[i], {t[i+1]}, p);
      // process error
      Ntrue[j,i+1] ~ lognormal(log(Nsim[1,1]), sigma);
    }
  }
  
  // observation equations
  for(j in 1:m){
    for(i in 1:n){
      N[j,i] ~ neg_binomial_2(Ntrue[j,i]*f, tau); 
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
                   sigma=0.1,
                   tau=10,
                   Ntrue=data$N/data$f
))
,chains)

stanmodel = stan_model(model_code=stanmodelcode)

fit = sampling(stanmodel,
               data=data,
               iter=iter,
               warmup=warmup,
               chains=chains,
               init=init,
               control=list(adapt_delta=0.8, max_treedepth=10)
)

## Model diagnostics -----------------------------------------------------------

print(fit, digits=3, pars=c("r","K","sigma","tau"), probs=c(0.025, 0.5, 0.975))

library("coda")
samples=As.mcmc.list(fit)
plot(samples[, c("r","K")])

pairs(fit, pars=c("r","K"))

## Posterior predictions -------------------------------------------------------

summary = summary(fit)$summary
N.pred.qs = matrix(NA, nrow=3, ncol=length(data$t))

par(mfrow=c(3,2), mar=c(3,3,0,0), oma=c(3,3,1,1))

for(i in 1:m){
  
  plot(data$t, data$N[i, ]/data$f, type="n", 
       ylim=c(0,1.1e4), xlab="Time [h]", ylab="Abundance [ind]", log="")
  
  for(j in 1:length(times)){
    N.pred.qs[, j] = summary[paste0("Ntrue[",i,",",j,"]"), c("2.5%", "50%", "97.5%") ]
  }
  
  polygon(c(data$t, rev(data$t)), c(N.pred.qs[1, ], rev(N.pred.qs[3, ])), col = adjustcolor("red",alpha.f=0.25),  border = NA)
  lines(data$t, N.pred.qs[2, ], col="red", lwd=1.5)

  points(data$t, data$N[i, ]/data$f) 
}
mtext("Time [h]", side=1, outer=TRUE, line=1)
mtext("Abundance [ind]", side=2, outer=TRUE, line=1)

