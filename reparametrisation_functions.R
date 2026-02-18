
##### two-parameter GEV distribution (upper/lower-bound fixed to zero => mu = sig/xi) #####

GEV_density_known_mu <- function(x,par) {
  
  # classic two-parameter GEV density
  
  sig <- par[1]
  
  xi <- par[2]
  
  if (abs(xi) <= 10e-3){return("xi must be non null")}
  
  else{
    
    t <- (xi*x/sig)
    
    return((1/sig)*t^(-1-1/xi)*exp(-t^(-1/xi)))}
}

GEV_nll_known_mu <- function(x,par) {
  
  # negative log-likelihood of the classic two-parameter GEV distribution
  
  sig <- par[1]
  
  xi <- par[2]
  
  t <- xi*x/sig
  
  if (sig <= 10e-3){l <- 10^6}
  
  if (xi > 0 & min(t) <= 0){l <- 10^6}
  
  if (xi < 0 & min(t) >= 0){l <- 10^6}
  
  else{l <- -sum(log(GEV_density_known_mu(x,par)))}
  
  l
}

GEV_density_known_mu_fixed_xi <- function(x,par) {
  
  # reparametrised, with fixed shape parameter, two-parameter GEV density
  
  rho <- par[1] # scale parameter in the new parametrisation
  
  xi <- par[2]
  
  gamma <- -digamma(1)
  
  if (xi > 10e-3){
    
    t <- x/(rho*exp((1-gamma)*xi))
    
    return((1/(rho*xi*exp((1-gamma)*xi)))*t^(-1-1/xi)*exp(-t^(-1/xi)))
  }
  
  if (xi < -10e-3){
    
    t <- -x/(rho*exp((1-gamma)*xi))
    
    return((-1/(rho*xi*exp((1-gamma)*xi)))*t^(-1-1/xi)*exp(-t^(-1/xi)))
  }
  
  if (abs(xi) <= 10e-3){return("xi must be non null")}
}

GEV_nll_known_mu_fixed_xi <- function(x,par) {
  
  # negative log-likelihood of the reparametrised, with fixed shape parameter, two-parameter GEV distribution
  
  gamma <- -digamma(1)
  
  rho <- par[1] # scale parameter in the new parametrisation
  
  xi <- par[2]
  
  t <- sign(xi)*x/(rho*exp((1-gamma)*xi))
  
  if ((rho*exp((1-gamma)*xi)) <= 10e-3){l <- 10^6}
  
  if (xi > 0 & min(t) <= 0){l <- 10^6}
  
  if (xi < 0 & min(t) >= 0){l <- 10^6}
  
  if (abs(xi) <= 10e-3){l <- 10^6}
  
  else{l <- -sum(log(GEV_density_known_mu_fixed_xi(x,par)))}
  
  l
}

### Score functions, i.e. partial first derivatives of the log-likelihood, of the classic two-parameter GEV distribution

score_sig_GEV  <- function(x,sig,xi) {
  
  # partial first derivative, according to the scale parameter, of the log-likelihood of the two-parameter GEV distibution
  
  if (abs(xi) <= 10e-3){return("xi must be non null")}
  
  else{
    t <- x*xi/sig
    
    return((1/(xi*sig))*(1-t^(-1/xi)))
  }
}

score_xi_GEV  <- function(x,sig,xi) {
  
  # partial first derivative, according to the shape parameter, of the log-likelihood of the two-parameter GEV distibution
  
  if (abs(xi) <= 10e-3){return("xi must be non null")}
  
  else{
    t <- x*xi/sig
    
    return((1/xi^2)*(log(t)-1)*(1-t^(-1/xi))-1/xi)
  }
}

### Information functions, i.e. partial second derivatives of the log-likelihood, of the classic two-parameter GEV distribution

IF_sig_GEV  <- function(x,sig,xi) {
  
  # partial second derivative, twice according to the scale parameter, of the log-likelihood of the GEV distibution
  
  if (abs(xi) <= 10e-3){return("xi must be non null")}
  
  else{
    t <- x*xi/sig
    
    return(-(1/(xi*sig^2))*(1-(t^(-1/xi))+(1/xi)*t^(-1/xi)))
  }
}

IF_xi_GEV  <- function(x,sig,xi) {
  
  # partial second derivative, twice according to the shape parameter, of the log-likelihood of the GEV distibution
  
  if (abs(xi) <= 10e-3){return("xi must be non null")}
  
  else{
    t <- x*xi/sig
    
    return((-2/xi^3)*log(t)+(1/xi^2)+(3/xi^3)-(t^(-1/xi))*((-2/xi^3)*log(t)+(3/xi^3)+((1/xi^2)*log(t)-1/xi^2)^2))
  }
}

IF_xi_sig_GEV  <- function(x,sig,xi) {
  
  # partial second derivative, once according to both parameters, of the log-likelihood of the GEV distibution
  
  if (abs(xi) <= 10e-3){return("xi must be non null")}
  
  else{
    
    t <- x*xi/sig
    
    return((1/(sig*xi^2))*(-1+t^(-1/xi)-(1/xi)*(t^(-1/xi))*(log(t)-1)))
  }
}

### Fisher information of the classic two-parameter GEV distribution

Fisher_sig_GEV  <- function(sig,xi) {
  
  # (1,1)-component of the Fisher information matrix i_{\sigma\sigma}
  
  return(1/(xi*sig)^2)
}

Fisher_xi_GEV  <- function(sig,xi) {
  
  # (2,2)-component of the Fisher information matrix i_{\xi\xi}
  
  gamma<- -digamma(1)
  
  return((-1/xi^2)*((-pi^2/6)-1+2*gamma-gamma^2+(2/xi)*(gamma-1)-1/xi^2))
}

Fisher_xi_sig_GEV  <- function(sig,xi) {
  
  # (1,2)-component of the Fisher information matrix i_{\sigma\xi}
  
  gamma<- -digamma(1)
  
  return((gamma-1)/(sig*xi^2)-1/(sig*xi^3))
}





##### two-parameter GP distribution (mu = zero) ######


GPD_density <- function(x,par) {
  
  # classic two-parameter GP density
  
  # the shape parameter is always assumed non null, since the problem of reparametrisation is not interesting in the case of a one-parameter distribution
  
  sig <- par[1]
  
  xi <- par[2]
  
  if (abs(xi) <= 10e-3){return((1/sig)*exp(-x/sig))}
  
  if (xi > 10e-3){
    
    t <- (1+xi*x/sig)*(x > 0)
    
    return((1/sig)*t^(-1-1/xi))
  }
  
  if (xi  < -10e-3){
    
    t <- (1+xi*x/sig)*(x >= 0 & x <= -sig/xi)
    
    return((1/sig)*t^(-1-1/xi))
  }
}

GPD_nll <- function(x,par) {
  
  # negative log-likelihood of the classic two-parameter GP distribution
  
  sig <- par[1]
  
  xi <- par[2]
  
  t <- 1 + xi*x/sig
  
  if (min(t) <= 0){l <- 10^6}
  
  else{l <- -sum(log(GPD_density(x,par)))}
  
  l
  
}

GPD_density_fixed_xi <- function(x,par) {
  
  # reparametrised, with fixed shape parameter, two-parameter GP density
  
  rho <- par[1] # scale parameter in the reparametrised settings
  
  xi <- par[2]
  
  if (rho <= 10e-4){return(10e9)}
  
  if (abs(xi) <= 10e-3){return((1/rho)*exp(-x/rho))}
  
  else{
    
    t <- (1+xi*(xi+1)*x/rho)
    
    return(((xi+1)/rho)*t^(-1-1/xi))
  }
}

GPD_nll_fixed_xi <- function(x,par) {
  
  # negative log-likelihood of the reparametrised, with fixed shape parameter, two-parameter GP distribution
  
  rho <- par[1]
  
  xi <- par[2]
  
  t <- 1+xi*(xi+1)*x/rho
  
  if (min(t) <= 0){l <- 10^6}
  
  else{l <- -sum(log(GPD_density_fixed_xi(x,par)))}
  
  l
}

GPD_density_fixed_sig <- function(x,par) {
  
  # reparametrised, with fixed scale parameter, two-parameter GP density
  
  sig <- par[1] 
  
  zeta <- par[2] # shape parameter in the reparametrised settings
  
  if (sig <= 10e-4){return(10e9)}
  
  if (abs(zeta-log(sig)/2) <= 10e-3){return((1/sig)*exp(-x/sig))}
  
  else{
    
    t <- (1+(zeta-log(sig)/2)*x/sig)
    
    return((1/sig)*t^(-1-1/(zeta-log(sig)/2)))
  }
}

GPD_nll_fixed_sig <- function(x,par) {
  
  # negative log-likelihood of the reparametrised, with fixed scale parameter, two-parameter GP distribution
  sig <- par[1]
  
  zeta <- par[2]
  
  y <- 1+(zeta-log(sig)/2)*x/sig
  
  if (min(y) <= 0){l <- 10^6}
  
  else{return(-sum(log(GPD_density_fixed_sig(x,par))))}
}



### Score functions of the classic two-parameter GPD distribution


score_sig_GPD  <- function(x,sig,xi) {
  
  # partial first derivative, according to the scale parameter, of the log-likelihood of the GPD distribution
  
  t <- x*xi/sig
  
  return(-1/sig+(xi+1)*x/(sig^2*(1+t)))
}

score_xi_GPD  <- function(x,sig,xi) {
  
  # partial first derivative, according to the shape parameter, of the log-likelihood of the GPD distribution
  
  t <- x*xi/sig
  
  return(log(1+t)/xi^2-(1+1/xi)*x/(sig*(1+t)))
}


### Information functions of the classic two-parameter GPD distribution


IF_sig_GPD  <- function(x,sig,xi) {
  
  # partial second derivative, twice according to the scale parameter, of the log-likelihood of the GPD distribution
  
  t <- x*xi/sig
  
  return((1/sig^2)*(
    1+(1+1/xi)*t^2*(1/(1+t))^2-2*(1+1/xi)*t*(1/(1+t))
  )
  )
}

IF_xi_GPD  <- function(x,sig,xi) {
  
  # partial second derivative, twice according to the shape parameter, of the log-likelihood of the GPD distribution
  
  t <- x*xi/sig
  
  return((-2/xi^3)*log(1+t)+2*x/(sig*xi^2*(1+t))+(1+1/xi)*(x/sig)^2*(1/(1+t)^2))
}

IF_xi_sig_GPD  <- function(x,sig,xi) {
  
  # partial second derivative, once according to both parameters, of the log-likelihood of the GPD distribution
  
  return((x*(sig-x))/(sig*(sig+xi*x)^2))
}

### Fisher matrix of the classic two-parameter GP distribution


Fisher_sig_GPD  <- function(sig,xi) {
  
  # (1,1)-component of the Fisher information matrix i_{\sigma\sigma}
  
  return(1/(sig^2*(1+2*xi)))
}

Fisher_xi_GPD  <- function(sig,xi) {
  
  # (2,2)-component of the Fisher information matrix i_{\xi\xi}
  
  return(2/((1+xi)*(1+2*xi)))
}

Fisher_xi_sig_GPD  <- function(sig,xi) {

  # (1,2)-component of the Fisher information matrix i_{\xi\sigma}
  
  return(1/(sig*(1+xi)*(1+2*xi)))
}


### Gumbel distribution


Gumbel_density <- function(x,par) {
  
  # classic Gumbel density
  
  mu <- par[1]
  
  sig <- par[2]
  
  t <- (x-mu)/sig
  
  return((1/sig)*exp(-t)*exp(-exp(-t)))
}

Gumbel_nll <- function(x,par) {
  
  # negative log-likelihood of the classic Gumbel distribution
  
  return(-sum(log(Gumbel_density(x,par))))
}

Gumbel_density_fixed_sig <- function(x,par) {
  
  # reparametrised, with fixed scale parameter, Gumbel density
  
  nu <- par[1]
  
  sig <- par[2]
  
  gamma<- -digamma(1)
  
  t <- (x-nu)/sig-(1-gamma)
  
  return((1/sig)*exp(-t)*exp(-exp(-t)))
}

Gumbel_nll_fixed_sig <- function(x,par) {
  
  # negative log-likelihood of the reparametrised, with fixed scale parameter, Gumbel distribution
  
  return(-sum(log(Gumbel_density_fixed_sig(x,par))))
}

Gumbel_density_fixed_mu <- function(x,par) {
  
  # reparametrised, with fixed location parameter, Gumbel density
  
  mu <- par[1]
  
  rho <- par[2]
  
  gamma<- -digamma(1)
  
  K <- (1-gamma)/(pi^2/6+gamma^2-2*gamma+1)
  
  t <- (x-mu)/(K*mu+rho)
  
  return((1/(K*mu+rho))*exp(-t)*exp(-exp(-t)))
}


Gumbel_nll_fixed_mu <- function(x,par) {
  
  # negative log-likelihood of the reparametrised, with fixed location parameter, Gumbel distribution
  
  return(-sum(log(Gumbel_density_fixed_mu(x,par))))
}

### Score functions of the classic Gumbel distribution

score_mu_Gumbel  <- function(x,mu,sig) {
  
  # partial first derivative, according to the location parameter, of the log-likelihood of the Gumbel distribution
  
  t <- (x-mu)/sig
  
  return((1/sig)*(1-exp(-t)))
}

score_sig_Gumbel  <- function(x,mu,sig) {
  
  # partial first derivative, according to the scale parameter, of the log-likelihood of the Gumbel distribution
  
  t <- (x-mu)/sig
  
  return((t/sig)*(1-exp(-t))-1/sig)
}



### Information functions of the classic Gumbel distribution


IF_mu_Gumbel  <- function(x,mu,sig) {
  
  # partial second derivative, twice according to the location parameter, of the log-likelihood of the Gumbel distribution
  
  t <- (x-mu)/sig
  
  return((-1/sig^2)*exp(-t))
}

IF_sig_Gumbel  <- function(x,mu,sig) {
  
  # partial second derivative, twice according to the scale parameter, of the log-likelihood of the Gumbel distribution
  
  t <- (x-mu)/sig
  
  return((1/sig^2)*(1-2*t*(1-exp(-t))-t^2*exp(-t)))
}

IF_mu_sig_Gumbel  <- function(x,mu,sig) {
  
  # partial second derivative, once according to both parameters, of the log-likelihood of the Gumbel distribution
  
  t <- (x-mu)/sig
  
  return((-1/sig^2)*(1-exp(-t)+t*exp(-t)))
}

### Fisher matrix Gumbel


Fisher_mu_Gumbel  <- function(mu,sig) {
  
  # (1,1)-component of the Fisher information matrix i_{\mu\mu}
  
  return(1/sig^2)
}

Fisher_sig_Gumbel  <- function(mu,sig) {
  
  # (2,2)-component of the Fisher information matrix i_{\sigma\sigma}
  
  gamma<- -digamma(1)
  
  return((pi^2/6+gamma^2-2*gamma+1)/sig^2)
}

Fisher_mu_sig_Gumbel  <- function(mu,sig) {
  
  # (1,2)-component of the Fisher information matrix i_{\mu\sigma}
  
  gamma<- -digamma(1)
  
  return((gamma-1)/sig^2)
}


##### Fisher information of the classic three-parameter GEV distribution


Fisher_mu_GEV3  <- function(mu,sig,xi) {
  
  # (1,1)-component of the Fisher information matrix i_{\mu\mu}
  
  p_xi <- (1+xi^2)*gamma(1+2*xi)
  
  return(p_xi/sig^2)
}


Fisher_sig_GEV3  <- function(mu,sig,xi) {
  
  # (2,2)-component of the Fisher information matrix i_{\sigma\sigma}
  
  p_xi <- (1+xi^2)*gamma(1+2*xi)
  
  return((1-2*gamma(2+xi)+p_xi)/(sig*xi)^2)
}


Fisher_xi_GEV3  <- function(mu,sig,xi) {
  
  # (3,3)-component of the Fisher information matrix i_{\xi\xi}
  
  gamma<- -digamma(1)
  
  p_xi <- (1+xi^2)*gamma(1+2*xi)
  
  q_xi <- gamma(2+xi)*(digamma(1+xi)+1/xi+1)
  
  return((1/xi^2)*(pi^2/6 +(1-gamma+1/xi)^2-2*q_xi/xi+p_xi/xi^2))
}

Fisher_mu_sig_GEV3  <- function(mu,sig,xi) {
  
  # (1,2)-component of the Fisher information matrix i_{\mu\sigma}
  
  p_xi <- (1+xi^2)*gamma(1+2*xi)
  
  return(-(p_xi-gamma(2+xi))/(sig^2*xi))
}

Fisher_mu_xi_GEV3  <- function(mu,sig,xi) {
  
  # (1,3)-component of the Fisher information matrix i_{\mu\xi}
  
  p_xi <- (1+xi^2)*gamma(1+2*xi)
  
  q_xi <- gamma(2+xi)*(digamma(1+xi)+1/xi+1)
  
  return(-(q_xi-p_xi/xi)/(sig*xi))
}

Fisher_sig_xi_GEV3  <- function(mu,sig,xi) {
  
  # (2,3)-component of the Fisher information matrix i_{\sigma\xi}
  
  gamma <- -digamma(1)
  
  p_xi <- (1+xi^2)*gamma(1+2*xi)
  
  q_xi <- gamma(2+xi)*(digamma(1+xi)+1/xi+1)
  
  return(-(1-gamma+(1-gamma(2+xi))/xi-q_xi+p_xi/xi)/(sig*xi^2))
}



### Fisher information of the classic three-parameter GP distribution

Fisher_mu_GPD3  <- function(mu,sig,xi) {
  
  # (2,2)-component of the Fisher information matrix i_{\sigma\sigma}
  
  return((xi+1)^2/(sig^2*(1+2*xi)))
}

Fisher_sig_GPD3  <- function(mu,sig,xi) {
  
  # (2,2)-component of the Fisher information matrix i_{\sigma\sigma}
  
  return(1/(sig^2*(1+2*xi)))
}

Fisher_xi_GPD3  <- function(mu,sig,xi) {
  
  # (3,3)-component of the Fisher information matrix i_{\xi\xi}
  
  return(2/((1+xi)*(1+2*xi)))
}

Fisher_mu_sig_GPD3  <- function(mu,sig,xi) {
  
  # (1,2)-component of the Fisher information matrix i_{\mu\sigma}
  
  return(-xi/(sig^2*(2*xi+1)))
}

Fisher_mu_sig_GPD3  <- function(mu,sig,xi) {
  
  # (1,3)-component of the Fisher information matrix i_{\mu\xi}
  
  return(-xi/(sig*(xi+1)*(2*xi+1)))
}

Fisher_xi_sig_GPD3  <- function(mu,sig,xi) {
  
  # (2,3)-component of the Fisher information matrix i_{\xi\sigma}
  
  return(1/(sig*(1+xi)*(1+2*xi)))
}





### Gamma distribution

Gamma_density <- function(x,par) {
  
  # classic Gumbel density
  
  alpha <- par[1]
  
  theta <- par[2]
  
  return(x^(alpha-1)*exp(-x/theta)/(gamma(alpha)*theta^alpha))
}

Gamma_nll <- function(x,par) {
  
  # negative log-likelihood of the classic Gumbel distribution
  
  return(-sum(log(Gamma_density(x,par))))
}

Gamma_density_fixed_alpha <- function(x,par) {
  
  # reparametrised, with fixed scale parameter, Gumbel density
  
  alpha <- par[1]
  
  omega <- par[2]
  
  return(x^(alpha-1)*alpha^alpha*exp(-x*alpha/omega)/(gamma(alpha)*omega^alpha))
}


Gamma_nll_fixed_alpha <- function(x,par) {
  
  # negative log-likelihood of the classic Gumbel distribution
  
  return(-sum(log(Gamma_density_fixed_alpha(x,par))))
}
