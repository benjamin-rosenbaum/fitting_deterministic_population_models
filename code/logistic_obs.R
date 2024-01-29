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
  vector<lower=0>[m] N0sim; // individual initial values for replicates
}

model{
  vector[2] p;
  array[n-1] vector[1] Nsim; // simulated values, matrix. dim1 = time without t0, dim2 = dim_ODE = 1
  
  // priors 
  r ~ lognormal(-2,1);
  K ~ lognormal(9,1);
  N0sim ~ normal(0,100);
  tau ~ gamma(2,0.1);
  
  // parameters for integrator
  p[1] = r;
  p[2] = K;
  
  for (j in 1:m){
    // integrate ODE
    Nsim = ode_rk45(odemodel, [N0sim[j]]', t[1], t[2:n], p);
    // likelihood
    N[j,1] ~ neg_binomial_2(N0sim[j]*f, tau);
    for (i in 2:n){
      N[j,i] ~ neg_binomial_2(Nsim[i-1,1]*f, tau);
    }  
  }
}

generated quantities{
  real alpha = r/K;
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
                   N0sim=data$N[, 1]/data$f,
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

print(fit, digits=7, pars=c("alpha"), probs=c(0.025, 0.5, 0.975))

plot(samples[, "alpha"])

## Posterior predictions -------------------------------------------------------

library("deSolve")

ode.model = function(t,N,p){
  dNdt =  p[1]*N[1]*(1.0-N[1]/p[2])
  return(list(dNdt))
}

post = as.matrix(fit)
head(post)

n.post = 1000 # nrow(post)
t.pred = seq(from=min(data$t), to=max(data$t), by=1)
N.pred = matrix(NA, nrow=n.post, ncol=length(t.pred))

par(mfrow=c(3,2), mar=c(3,3,0,0), oma=c(3,3,1,1))

for(i in 1:6){
  for(j in 1:n.post){
    N.pred[j, ] = as.data.frame(lsoda(y     = c( N = post[j, paste0("N0sim[",i,"]")] ),
                                      times = t.pred,
                                      func  = ode.model,
                                      parms = c(post[j, "r"],
                                                post[j, "K"])
    ))$N
  }
  
  N.pred.qs = apply(N.pred, 2, function(x) quantile(x, probs=c(0.025,0.500,0.975), na.rm=TRUE))
  
  plot(data$t, data$N[i, ]/data$f, type="n")
  points(data$t, data$N[i, ]/data$f) 
  polygon(c(t.pred, rev(t.pred)), c(N.pred.qs[1, ], rev(N.pred.qs[3, ])), 
          col = adjustcolor("red",alpha.f=0.25),  border = NA)
  lines(t.pred, N.pred.qs[2, ], col="red", lwd=1.5)
}
mtext("Time [h]", side=1, outer=TRUE, line=1)
mtext("Abundance [ind]", side=2, outer=TRUE, line=1)




