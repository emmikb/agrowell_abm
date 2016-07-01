require(rstan)
require(shinystan)

rstan_options(auto_write = TRUE)

#source('C:/Users/Emily Burchfield/Box Sync/WF/Survey/load_seads_data.R')
#source('C:/Users/Emily Burchfield/Box Sync/WF/Survey/variable_creation.R')

###############################################
#DATA PREP
##############################################
gn.names <- as.vector(c1$HI4)
uq <- unique(gn.names)
ngn<- length(uq)
gn <- rep(NA, ngn)
for (i in 1:ngn){
  gn[gn.names == uq[i]] <- i
  sample.size <- as.vector(table(gn))
}
c1$gn <- gn

#outcome variable
y <- c1$FAR1B_20  #Do you regularly cultivate OFC in your paddy field?
y <-ifelse(y==1, 1, 0)
n <- length(y)
k <- 8


#############################################
#STAN MODEL
#############################################

ofc_model <- 
'data {
  int<lower=0> n;             //number of farmers
  int<lower=0> ngn;           //number of gns
  int k;                      //number of predictors

  int<lower=0, upper=1> y[n];  //outcome
  row_vector[k] x[n];         //predictors
  int gn[n];                  //mapping gn to group

}

parameters{
  real alpha;               //constant
  real a[ngn];              //group-specific random effects
  vector[k] beta;           //predictor coefficiencts
  real sigma;               //standard deviation of random effects
}

model {
  alpha ~ normal(0,100);
  a ~ normal(0, sigma);
  beta ~ normal(0, 100);
  for (i in 1:n) {
    y[i] ~ bernoulli( inv_logit(alpha + a[gn[i]] + x[i]*beta));
  }

}'
  
ofc_data <- list(aw = c1$agrowell_user, major = c1$major_flag,
                 female = c1$female, sinhala = c1$sinhalese,
                 ses = c1$Standardized_SES, land_owner = c1$owner, 
                 he = c1$head_end, fo = c1$fo)

ofc_data <- as.data.frame(ofc_data)

model_data <-  list(n=n, ngn = ngn, k=8, y = y, x= ofc_data[,], gn = gn)

lfit <- stan(model_code = ofc_model, model_name = "OFC", data = model_data, iter = 2000, chains = 2)
my_sso <- launch_shinystan(lfit)
save(lfit, file = 'lfit.Rda')



ofc_model <- 
  'data {
int<lower=0> n; 
int<lower=0> ngn;
vector<lower=0, upper=1>[n] he;
vector<lower=0, upper=1>[n] aw;
int<lower=0, upper=ngn> gn[n];
int<lower=0, upper=1>y[n];
}

parameters{
vector[ngn] a;
vector[2] beta;
real<lower=0, upper=100> sigma_a;
real mu_a;
}

transformed parameters{
vector[n] y_hat;

for (i in 1:n)
y_hat[i] <- beta[1] * aw[i] + beta[2] * he[i] + a[gn[i]];
}

model {
mu_a ~ normal(0,1);
beta ~ normal(0,1);
a ~ normal(mu_a, sigma_a);
y ~ bernoulli_logit(y_hat);
}'