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
  vector<lower=0>[2] sigma;
  vector<lower=0>[2] tau;
  matrix<lower=0>[m,n] N1true; // states in times t
  matrix<lower=0>[m,n] N2true; // states in times t
  matrix<lower=0>[mc,n] N1Ctrue; // states in times t
  matrix<lower=0>[mc,10] N2Ctrue; // states in times t, only 10 datapoints
}

model{
  vector[6] p;
  array[1] vector[2] Nsim; // simulated values, matrix. dim1 = 1 point in time, dim2 = dim_ODE = 2
  array[1] vector[1] NCsim; // simulated values, matrix. dim1 = 1 point in time, dim2 = dim_ODE = 1

  // priors 
  r ~ lognormal(-2,1);
  K ~ lognormal(9,1);
  b ~ gamma(2,1);
  h ~ gamma(2,1);
  // e ~ beta(2,2);
  d ~ gamma(2,1);
  // sigma ~ gamma(2,1);
  sigma ~ lognormal(-2,1);
  tau ~ gamma(2,0.1);

  for(i in 1:m){
    N1true[i,1] ~ normal(0,2000);
    N2true[i,1] ~ normal(0,2000);
  }
  for(i in 1:mc){
    N1Ctrue[i,1] ~ normal(0,2000);
    N2Ctrue[i,1] ~ normal(0,2000);
  }

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

  // process equation
  for(j in 1:m){
    for(i in 1:(n-1)){
      Nsim = ode_rk45(odemodel, [N1true[j,i],N2true[j,i]]', t[i], {t[i+1]}, p);
      N1true[j,i+1] ~ lognormal(log(Nsim[1,1]), sigma[1]);
      N2true[j,i+1] ~ lognormal(log(Nsim[1,2]), sigma[2]);
    } 
  }
  // observation equation
  for(j in 1:m){
    for(i in 1:n){
      N1[j,i] ~ neg_binomial_2(N1true[j,i]*f, tau[1]); 
      N2[j,i] ~ neg_binomial_2(N2true[j,i]*f, tau[2]); 
    }
  }

  //---------------------------------------------------------------------------
  // resource: control data
  //---------------------------------------------------------------------------
  
  // process equation
  for(j in 1:mc){
    for(i in 1:(n-1)){
      NCsim[1,1] = K/( 1.0+(K-N1Ctrue[j,i])/N1Ctrue[j,i] * exp(-r*(t[i+1]-t[i]))); // analytical solution
      N1Ctrue[j,i+1] ~ lognormal(log(NCsim[1,1]), sigma[1]);
    }
  }
  // observation equation
  for(j in 1:mc){
    // observation error
    for(i in 1:n){
      N1C[j,i] ~ neg_binomial_2(N1Ctrue[j,i]*f, tau[1]); 
    }
  }

  //---------------------------------------------------------------------------
  // consumer: control data
  //---------------------------------------------------------------------------
  
  // process equation
  for(j in 1:mc){
    for(i in 1:(10-1)){ // only 10 datapoints
      NCsim[1,1] = N2Ctrue[j,i] * exp(-d*(t[i+1]-t[i])); // analytical solution
      N2Ctrue[j,i+1] ~ lognormal(log(NCsim[1,1]), sigma[2]);
    } 
  }
  // observation equation
  for(j in 1:mc){
    for(i in 1:10){ // only 10 datapoints
      N2C[j,i] ~ neg_binomial_2(N2Ctrue[j,i]*f, tau[2]); 
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
                   sigma=c(0.1,0.1),
                   tau=c(10,10),
                   N1true=data$N1/data$f,
                   N2true=data$N2/data$f,
                   N1Ctrue=data$N1C/data$f,
                   N2Ctrue=data$N2C[, 1:10]/data$f
))
,chains)

stanmodel = stan_model(model_code=stanmodelcode)

fit = sampling(stanmodel,
               data=data,
               iter=iter,
               warmup=warmup,
               chains=chains,
               init=init,
               refresh=10,
               control=list(adapt_delta=0.95, max_treedepth=12)
)

## Model diagnostics -----------------------------------------------------------

print(fit, digits=3, pars=c("r","K","b","h","e","d","sigma","tau"))

library("coda")
samples=As.mcmc.list(fit)
plot(samples[, c("r","K","b","h","e","d")])

pairs(fit, pars=c("r","K","b","h","e","d"))

## Posterior predictions -------------------------------------------------------

N1.pred.qs = matrix(NA, nrow=3, ncol=length(data$t))
N2.pred.qs = matrix(NA, nrow=3, ncol=length(data$t))

par(mfrow=c(3,2), mar=c(3,3,0,0), oma=c(3,3,1,1))

for(i in 1:m){
  for(j in 1:length(times)){
    N1.pred.qs[, j] = summary[paste0("N1true[",i,",",j,"]"), c("2.5%", "50%", "97.5%") ]
    N2.pred.qs[, j] = summary[paste0("N2true[",i,",",j,"]"), c("2.5%", "50%", "97.5%") ]
  }
  
  plot(data$t, data$N1[i, ]/data$f,
       ylim=c(1e2,2e4), type="n", xlab="", ylab="", log="y")
  
  polygon(c(data$t, rev(data$t)), c(N1.pred.qs[1, ], rev(N1.pred.qs[3, ])), 
          col = adjustcolor("red",alpha.f=0.25),  border = NA)
  lines(data$t, N1.pred.qs[2, ], col="red", lwd=1.5)
  points(data$t, data$N1[i, ]/data$f) 
  
  polygon(c(data$t, rev(data$t)), c(N2.pred.qs[1, ], rev(N2.pred.qs[3, ])), 
          col = adjustcolor("blue",alpha.f=0.25),  border = NA)
  lines(data$t, N2.pred.qs[2, ], col="blue", lwd=1.5)
  points(data$t, data$N2[i, ]/data$f) 
  
  if(i==1) legend("topleft", c("res","con"), lwd=c(2,2), col=c("red","blue"))
  
}
mtext("Time [h]", side=1, outer=TRUE, line=1)
mtext("Abundance [ind], logscale", side=2, outer=TRUE, line=1)
