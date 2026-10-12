functions {
  real cgp_gaussian_log_density(vector state, data matrix distance2,
    data matrix group_cov, data vector observed, data int kernel);
  #include censor-transform.staninc
}
data {
  int<lower=1> N;
  int<lower=1> G;
  int<lower=0> No;
  int<lower=0> Nl;
  int<lower=0> Nr;
  int<lower=0> Ni;
  int<lower=0,upper=1> parameterization;
  int<lower=0,upper=1> custom_gradient;
  vector[No] yo;
  vector[Nl] left_upper;
  vector[Nr] right_lower;
  vector[Ni] interval_lower;
  vector[Ni] interval_upper;
  array[N] real t;
  array[N] int<lower=1,upper=G> group;
}
transformed data {
  int Nc=Nl+Nr+Ni;
  matrix[N,N] distance2;
  matrix[N,N] group_cov;
  vector[parameterization==0?Nl:0] bound_left=head(left_upper,parameterization==0?Nl:0);
  vector[parameterization==0?Nr:0] bound_right=head(right_lower,parameterization==0?Nr:0);
  vector[parameterization==0?Ni:0] bound_interval_lower=head(interval_lower,parameterization==0?Ni:0);
  vector[parameterization==0?Ni:0] bound_interval_upper=head(interval_upper,parameterization==0?Ni:0);
  vector[Nc] b_lower=append_row(rep_vector(negative_infinity(),Nl),append_row(right_lower,interval_lower));
  vector[Nc] b_upper=append_row(left_upper,append_row(rep_vector(positive_infinity(),Nr),interval_upper));
  array[Nc] int kind;
  if(No + Nc != N) reject("Observation and censoring counts must sum to N");
  if(Ni>0) for(i in 1:Ni) {
    if(interval_lower[i]>=interval_upper[i]) reject("Interval lower bound must be below upper bound");
  }
  for(i in 1:N) for(j in 1:N) {
    distance2[i,j]=square(t[i]-t[j]);
    group_cov[i,j]=G==1 ? 0 : (group[i]==group[j]?1.0:-1.0/(G-1.0));
  }
  if(Nl>0) for(i in 1:Nl) kind[i]=1;
  if(Nr>0) for(i in 1:Nr) kind[Nl+i]=2;
  if(Ni>0) for(i in 1:Ni) kind[Nl+Nr+i]=3;
}
parameters {
  real<lower=0> alpha_mu;
  real<lower=0> rho_mu;
  real<lower=0> alpha_eta;
  real<lower=0> rho_eta;
  real<lower=0> sigma;
  vector<upper=bound_left>[parameterization==0?Nl:0] yc_left;
  vector<lower=bound_right>[parameterization==0?Nr:0] yc_right;
  vector<lower=bound_interval_lower,upper=bound_interval_upper>[parameterization==0?Ni:0] yc_interval;
  vector<lower=0,upper=1>[parameterization==1?Nc:0] u;
}
transformed parameters {
  vector[Nc] y_cens;
  real log_lik;
  if(custom_gradient==1 && (parameterization==0 || Nc==0)) {
    vector[5+Nc] state;
    y_cens=append_row(yc_left,append_row(yc_right,yc_interval));
    state[:5]=[alpha_mu,rho_mu,alpha_eta,rho_eta,sigma]';
    if(Nc>0) state[6:]=y_cens;
    log_lik=cgp_gaussian_log_density(state,distance2,group_cov,yo,1);
  } else {
    matrix[N,N] K=square(alpha_mu)*exp(-0.5*distance2/square(rho_mu));
    if(G>1) K+=group_cov.*(square(alpha_eta)*exp(-0.5*distance2/square(rho_eta)));
    {
      matrix[N,N] L=cholesky_decompose(add_diag(K,square(sigma)));
      if(parameterization==0 || Nc==0) {
        y_cens=append_row(yc_left,append_row(yc_right,yc_interval));
        log_lik=multi_normal_cholesky_lpdf(append_row(yo,y_cens)|rep_vector(0,N),L);
      } else {
        tuple(vector[Nc],real) result=censor_transform(L,yo,u,b_lower,b_upper,kind);
        y_cens=result.1;
        log_lik=result.2;
        if(No>0) log_lik+=multi_normal_cholesky_lpdf(yo|rep_vector(0,No),L[:No,:No]);
      }
    }
  }
}
model {
  alpha_mu~normal(0,1.5);
  alpha_eta~normal(0,1);
  rho_mu~lognormal(log(0.3),0.7);
  rho_eta~lognormal(log(0.2),0.7);
  sigma~normal(0,0.5);
  target+=log_lik;
}
