rm(list=ls())
library("rstan")

## Data preparation ------------------------------------------------------------

data.in = read.csv("data_example_consumer_resource.csv", row.names=1)

n = 18
m = 6
mc = 3

times = as.numeric(data.in[1, ])
N1 = data.in[2:(m+1), ]
N2 = data.in[(m+2):(2*m+1), ]
N1C = data.in[(2*m+2):(2*m+1+mc), ]
N2C = data.in[(2*m+2+mc):((2*m+1+2*mc)), ]

data = list(n = n,
            m = m,
            mc = mc,
            t = times,
            N1 = N1,
            N2 = N2,
            N1C  = N1C,
            N2C  = N2C,
            f = 0.1)

str(data)

## Stan model ------------------------------------------------------------------

stanmodelcode="
functions{
  vector odemodel(real t, vector N, vector p){
    // p[1]=r, p[2]=K, p[3]=b, p[4]=h, p[5]=e, p[6]=d
    vector[2] dNdt;
    dNdt[1] = p[1]*N[1]*(1.0-N[1]/p[2])-p[3]*N[1]*N[2]/(1.0+p[3]*p[4]*N[1]);
    dNdt[2] = p[5]*p[3]*N[1]*N[2]/(1.0+p[3]*p[4]*N[1])-p[6]*N[2];
    return dNdt;
  }
}

data{
  int n; // observations
  int m; // replicates
  int mc; // replicates control
  real t[n];
  array[m,n] int N1;
  array[m,n] int N2;
  array[mc,n] int N1C;
  array[mc,n] int N2C;
  real f;
}

parameters{
  real<lower=0> r;
  real<lower=0> K;
  real<lower=0> b;
  real<lower=0> h;
  real<lower=0, upper=1> e;
  real<lower=0> d;
  vector<lower=0>[2] tau;
  vector<lower=0>[m] N10sim; // individual initial values for replicates
  vector<lower=0>[m] N20sim; // individual initial values for replicates
  vector<lower=0>[mc] N1C0sim; // individual initial values for replicates
  vector<lower=0>[mc] N2C0sim; // individual initial values for replicates
}

model{
  vector[6] p;
  array[n-1] vector[2] Nsim; // simulated values, matrix. dim1 = time without t0, dim2 = dim_ODE = 2
  array[n-1] vector[1] NCsim; // simulated values, matrix. dim1 = time without t0, dim2 = dim_ODE = 1
  
  // priors 
  r ~ lognormal(-2,1);
  K ~ lognormal(9,1);
  b ~ gamma(2,1);
  h ~ gamma(2,1);
  // e ~ beta(2,2);
  d ~ gamma(2,1);
  tau ~ gamma(2,0.1);
  
  N10sim ~ normal(0,2000);
  N20sim ~ normal(0,2000);
  N1C0sim ~ normal(0,2000);
  N2C0sim ~ normal(0,2000);
  
  //---------------------------------------------------------------------------
  // 2-species data
  //---------------------------------------------------------------------------
    
  // parameters for integrator
  p[1] = r;
  p[2] = K;
  p[3] = b;
  p[4] = h;
  p[5] = e;
  p[6] = d;
  
  for (j in 1:m){
    // integrate ODE
    Nsim = ode_rk45(odemodel, [N10sim[j],N20sim[j]]', t[1], t[2:n], p);
    // likelihood
    N1[j,1] ~ neg_binomial_2(N10sim[j]*f,tau[1]);
    N2[j,1] ~ neg_binomial_2(N20sim[j]*f,tau[2]);
    for (i in 2:n){
      N1[j,i] ~ neg_binomial_2(Nsim[i-1,1]*f,tau[1]);
      N2[j,i] ~ neg_binomial_2(Nsim[i-1,2]*f,tau[2]);
    }  
  }
  
  //---------------------------------------------------------------------------
  // resource: control data
  //---------------------------------------------------------------------------
    
  for (j in 1:mc){
    // integrate ODE
    for(i in 2:n){
      NCsim[i-1,1] = K / ( 1+(K-N1C0sim[j])/N1C0sim[j] * exp(-r*(t[i]-t[1]))); // analytical solution
    }  
    // likelihood
    N1C[j,1] ~ neg_binomial_2(N1C0sim[j]*f,tau[1]);
    for (i in 2:n){
      N1C[j,i] ~ neg_binomial_2(NCsim[i-1,1]*f,tau[1]);
    }  
  }

  //---------------------------------------------------------------------------
  // consumer: control data
  //---------------------------------------------------------------------------
  
  for (j in 1:mc){
    // integrate ODE
    for(i in 2:n){
      NCsim[i-1,1] = N2C0sim[j] * exp(-d*(t[i]-t[1])); // analytical solution
    }  
    // likelihood
    N2C[j,1] ~ neg_binomial_2(N2C0sim[j]*f,tau[2]);
    for (i in 2:n){
      N2C[j,i] ~ neg_binomial_2(NCsim[i-1,1]*f,tau[2]);
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
                   K=2e4,
                   b=0.0003,
                   h=0.5,
                   e=0.1,
                   d=0.1,
                   N10sim=data$N1[, 1]/data$f, 
                   N20sim=data$N2[, 1]/data$f,
                   N1C0sim=data$N1C[, 1]/data$f,
                   N2C0sim=data$N2C[, 1]/data$f,
                   tau=c(10,10)
))
,chains)

stanmodel = stan_model(model_code=stanmodelcode)

fit = sampling(stanmodel,
               data=data,
               iter=iter,
               warmup=warmup,
               chains=chains,
               init=init,
               refresh=10
)

## Model diagnostics -----------------------------------------------------------

print(fit, digits=3, pars=c("r","K","b","h","e","d","tau"), probs=c(0.025, 0.5, 0.975))

library("coda")
samples=As.mcmc.list(fit)
plot(samples[, c("r","K","b","h","e","d")])

pairs(fit, pars=c("r","K","b","h","e","d"))

## Posterior predictions -------------------------------------------------------

library("deSolve")

ode.model = function(t,N,p){
  dN1dt = p[1]*N[1]*(1.0-N[1]/p[2])-p[3]*N[1]*N[2]/(1.0+p[3]*p[4]*N[1])
  dN2dt = p[5]*p[3]*N[1]*N[2]/(1.0+p[3]*p[4]*N[1])-p[6]*N[2]
  return(list(c(dN1dt,dN2dt)))
}

post = as.matrix(fit)
n.post = nrow(post)
n.post = 1000
t.pred = seq(from=min(data$t), to=max(data$t), by=1)
N1.pred = matrix(NA, nrow=n.post, ncol=length(t.pred))
N2.pred = matrix(NA, nrow=n.post, ncol=length(t.pred))
p = rep(NA, 6)

par(mfrow=c(3,2), mar=c(3,3,0,0), oma=c(3,3,1,1))

for(i in 1:m){
  for(j in 1:n.post){
    Ns = as.data.frame(lsoda(y     = c(N=as.numeric(post[j, c(paste0("N10sim[",i,"]"),
                                                              paste0("N20sim[",i,"]"))])),
                             times = t.pred,
                             func  = ode.model,
                             parms = c(post[j, c("r","K","b","h","e","d")])
    ))[, 2:3]
    N1.pred[j, ]=Ns[, 1]
    N2.pred[j, ]=Ns[, 2]
  }

  plot(c(data$t,data$t), c(data$N1[i, ]/data$f,data$N2[i, ]/data$f), type="n", xlab="", ylab="", log="y")

  N1.pred.qs = apply(N1.pred, 2, function(x) quantile(x, probs=c(0.025,0.500,0.975), na.rm=TRUE))
  polygon(c(t.pred, rev(t.pred)), c(N1.pred.qs[1, ], rev(N1.pred.qs[3, ])), col = adjustcolor("red",alpha.f=0.25),  border = NA)
  lines(t.pred, N1.pred.qs[2, ], col="red", lwd=1.5)
  points(data$t, data$N1[i, ]/data$f) 
  
  N2.pred.qs = apply(N2.pred, 2, function(x) quantile(x, probs=c(0.025,0.500,0.975), na.rm=TRUE))
  polygon(c(t.pred, rev(t.pred)), c(N2.pred.qs[1, ], rev(N2.pred.qs[3, ])), col = adjustcolor("blue",alpha.f=0.25),  border = NA)
  lines(t.pred, N2.pred.qs[2, ], col="blue", lwd=1.5)
  points(data$t, data$N2[i, ]/data$f) 
}
mtext("Time [h]", side=1, outer=TRUE, line=1)
mtext("Abundance [ind]", side=2, outer=TRUE, line=1)
