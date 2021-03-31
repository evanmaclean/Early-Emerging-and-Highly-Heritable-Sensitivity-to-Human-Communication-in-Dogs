library(tidyverse)
library(here)
library(rstan)
library(RNOmni)
library(kinship2)
library(bayesplot)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load('h2_data.RData')

####
## Build A matrix
###
A <- 2*kinship(NPED$idanim,NPED$isire,NPED$idam)
###

# Create the Stan Model
unimodel = "

data {
int<lower=1> Q; // number of independent variables
int<lower=0> nanim; // number of animals
int<lower=0> nrec;  // number of records in 1 trait
real Y[nrec]; // response variable
matrix[nrec,Q] X; // design matrix
int pid[nrec];  // animal number - consecutive
matrix[nanim, nanim]    A; // relationship matrix
}

transformed data{
matrix[nanim, nanim] LA;
LA = cholesky_decompose(A);
}

parameters {
vector[Q] beta;
vector[nanim] ebv; // breeding values
real<lower=0> sg;  // genetic sd
real<lower=0> se;  // residual sd
}

transformed parameters {
 vector[nanim] a;
 real<lower=0,upper=1> h2;
 h2 = (sg*sg)/(sg*sg+se*se);
 a = sg * (LA * ebv);
}

model{
real mu[nrec];

for(i in 1:nrec){
mu[i] = (X[i,]*beta + a[pid[i]] ); 
}
Y ~ normal(mu,se);  

beta ~ normal(0,1000);
ebv   ~ normal(0, 1);
sg ~ cauchy(0, 5);
se ~ cauchy(0, 5);
}

generated quantities {
vector[nrec] log_lik;
vector[nanim] bv;
real pred[nrec];
bv = sg * (LA * ebv);
for(i in 1:nrec){
pred[i] = normal_rng(X[i,]*beta + bv[pid[i]],se); 
log_lik[i] = normal_lpdf(Y[i] | X[i,]*beta + bv[pid[i]],se);
}
}

"
comp_mod = stan_model(model_code=unimodel,model_name="unimodel")
######
######
# pick trait
my.traits <- c('Pointing_PerCorrect','Marker_PerCorrect','HumInt_avglookingtime','HumInt_avginteracttime','Unsolvable_avglooktime')

results.list <- list()
beta.list <- list()
j = 0
for (trait in my.traits) {
print(trait)
j = j+1
dat2 <- dat[!is.na(dat[,trait]),]

##
# Build the data list for stan
##

Y <- dat2[,trait]
Y <- RankNorm(Y)
#######
pid <- dat2$idanim
X <- model.matrix(~ sex + breed + puptestage + raise_location ,data=dat2)
beta.list[[j]] <- colnames(X) # saving the beta names here

Q <- length(X[1,])
nrec <- length(Y)
nanim <- max(dat$idanim)
data_list <- list(Q=Q,nrec=nrec,nanim=nanim,A=A,
                  X=X,Y=Y,pid=pid)
#########

#### Run the model

my.iter = 75000
my.warm = 5000
my.thin = 20
my.chains = 4
my.cores = 4

my.fit = sampling(comp_mod,
               data=data_list, pars = c("beta","h2","sg","se", "log_lik","pred"), iter=my.iter,warmup=my.warm, thin=my.thin, chains=my.chains, cores = my.cores)
results.list[[j]] <- my.fit
}
names(results.list) <- my.traits

#######
## Summarize results
#######

my.posteriors <- matrix(nrow = 14000, ncol = length(results.list))# 1400 draws
my.posteriors <- as.data.frame(my.posteriors)
colnames(my.posteriors) <- names(results.list)

for (i in 1:length(results.list)) {
  my.var <- names(results.list[i])
  tmp <- results.list[[i]]
  my.smry <- summary(tmp, pars = 'h2', probs = c(0.05, 0.5, 0.95))$summary
  beta.names <- beta.list[[i]]
  my.betas <- summary(tmp, pars = c('beta[2]','beta[3]','beta[4]','beta[5]','beta[6]'), probs = c(0.05, 0.5, 0.95))$summary
  # replace uninformative beta names with actual names
  row.names(my.betas) <- beta.names[2:6]
  my.smry <- as.data.frame(my.smry)
  my.betas <- as.data.frame(my.betas)
  my.smry <- bind_rows(my.smry, my.betas)
  
  my.smry$divergences = sum(get_divergent_iterations(tmp),inc_warmup = F)
  list_of_draws <- extract(tmp)
  my.smry$trait <- names(results.list[i])
  my.smry <- dplyr::select(my.smry, trait, everything())
  names(list_of_draws)
  length(list_of_draws$h2)
  my.posteriors[,my.var] <- list_of_draws$h2
  my.smry$iterations <- dim(list_of_draws$h2)
  my.smry <-  rownames_to_column(my.smry, var = 'param')
  
  if (i == 1) {
    results.out <- my.smry
  } else {
    results.out <- bind_rows(results.out, my.smry)
  }
}

my.posteriors <- dplyr::select(my.posteriors, 'Pointing_PerCorrect','HumInt_avglookingtime','HumInt_avginteracttime','Marker_PerCorrect','Unsolvable_avglooktime')

my.intervals <- mcmc_intervals_data(my.posteriors)

theme_set(theme_bw())
mycolors <- c('#4CAECD','#ED9E09','#878787')

p <- ggplot(my.intervals, aes(x = parameter, y = m)) +
  geom_linerange(aes(ymin = ll, ymax = hh), size = 0.75, color = 'black') +
  geom_linerange(aes(ymin = l, ymax = h), size = 5, color = '#4CAECD' ) +
  geom_point(color = 'black', size = 5) +
  geom_point(color = 'white', size = 4) +
  #geom_hline(yintercept = 0, color = 'black', size = 0.5) +
  coord_flip() +
  labs(x = "", y = "h2") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(axis.text.x = element_text(size = 22), axis.title.x = element_text(size = 22)) + 
  theme(axis.text.y = element_text(size = 16)) + 
  labs(y = 'heritability (posterior distribution)')
p
