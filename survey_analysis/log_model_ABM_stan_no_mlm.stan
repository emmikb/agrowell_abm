data {
  int<lower=0> n;             //number of farmers
vector[n] aw;
vector[n] major;
vector[n] female;
vector[n] sinhala;
vector[n] ses;
vector[n] land_owner;
vector[n] he;
vector[n] fo;
vector[n] y;

}

parameters{
real beta_aw;
real beta_major;
real beta_female;
real beta_sinhala;
real beta_ses;
real beta_land_owner;
real beta_he;
real beta_fo;

real alpha;

real<lower=0> sigma;
}

model {

vector[n] mu;

beta_aw ~ cauchy(0,2.5);
beta_major ~ cauchy(0,2.5);
beta_female ~ cauchy(0,2.5);
beta_sinhala ~ cauchy(0,2.5);
beta_ses ~ cauchy(0,2.5);
beta_land_owner ~ cauchy(0,2.5);
beta_he ~ cauchy(0,2.5);
beta_fo ~ cauchy(0,2.5);

alpha ~ cauchy(20,20);
sigma ~ cauchy(0, 10);

mu <- alpha + beta_aw*aw + beta_major*major + beta_female*female + beta_sinhala*sinhala +
beta_ses*ses + beta_land_owner*land_owner + beta_he*he + beta_fo*fo;

y ~ normal(mu, sigma);

}
