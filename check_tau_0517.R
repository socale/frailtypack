
# Codes for checking the data generation for the Clayton copula #

G=500;N=1000

theta_true=5 ## parameter for the Clayton copula ##
Tau_true=theta_true/(theta_true+2)

beta1_true=1
beta2_true=1
r1_true=1
r2_true=1

eta_true=1
alpha_true=0

##### Monte carlo #####
S_vec=T_vec=Z=group=NULL
ij=0

for(i in 1:G){
  set.seed(i)
  u_i=rgamma(1,shape=1/eta_true,scale=eta_true)

  for(j in 1:N){
    
    ij=ij+1
    
    group[ij]=i
    Z[ij]=rbinom(n=1,size=1,prob=0.5)
    r1_ij=r1_true*u_i*exp(beta1_true*Z[ij])
    r2_ij=r2_true*(u_i^alpha_true)*exp(beta2_true*Z[ij])
    V1=runif(1)
    V2=runif(1)
    S_ij=-1/r1_ij*log(1-V1);W=(1-V1)^(-theta_true)
    T_ij=1/theta_true/r2_ij*log(  1-W+W*(1-V2)^(-theta_true/(theta_true+1))  )
    S_vec[ij]=S_ij
    T_vec[ij]=T_ij
  }
  
}

## Compute Kendall's tau at each study ##
Tau=numeric(G)
for(i in 1:G){
  G_i=c(group==i)
  Tau[i]=cor(S_vec[G_i&Z],T_vec[G_i&Z],method="kendall")
}

## True vs. empirical ##
Tau_true
mean(Tau)
## difference at most 0.01##

