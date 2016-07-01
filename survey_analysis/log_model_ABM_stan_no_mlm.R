library(rstan)
library(shinystan)
library(ggmcmc)
library(stringr)
library(dplyr)
library(purrr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source('load_seads_data.R')
source('variable_creation.R')

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

#############################################
#STAN MODEL
#############################################

ofc_data <- list(n=n, y = y, aw = c1$agrowell_user, major = c1$major,
                 female = c1$female, sinhala = c1$sinhalese,
                 ses = c1$Standardized_SES, land_owner = c1$owner,
                 he = c1$prop.he, fo = c1$fo)

#ofc_data <- as.data.frame(ofc_data)

lfit_no_mlm <- stan('log_model_ABM_stan_no_mlm.stan', model_name = "OFC no MLM", data = ofc_data, iter = 2000, chains = 2)
save(lfit_no_mlm, file = "lfit_no_mlm.Rda")
# my_sso <- launch_shinystan(lfit_no_mlm)

par_names <- names(lfit_no_mlm) %>% purrr::keep(~str_detect(.,"^beta_"))
new_par_names <- str_replace_all(par_names,
                                 c('aw$'='AW','ses$'='status', 'he$'='HE','fo$'='FO',
                                   'sinhala$'='Sinhala','land_owner$' = 'landowner',
                                   'beta_(.+)$'='beta[\\1]' ))
par_labels = data.frame(Parameter = par_names, Label = new_par_names)


g <- ggs(lfit_no_mlm, family = 'beta_', par_labels = par_labels)

parsed_labels <- g %>% group_by(Parameter) %>% summarize(m = median(value)) %>%
  ungroup() %>% arrange(m) %>%
  mutate(Parameter = as.character(Parameter)) %>% dplyr::select(Parameter) %>%
  unlist() %>% unname() %>% parse(text = .)

pdf('ofc_regression.pdf', height = 3, width = 4.5)
ggs_caterpillar(g, greek = F, thick_ci = c(0.17, 0.83)) +
  geom_vline(xintercept = 0, size=0.1) +
  geom_point(aes(x = median), size = 1) +
  scale_y_discrete(labels = parsed_labels) +
  labs(x = "Value", y = "Coefficient") + theme_bw(base_size = 10)
dev.off()
