## ---- echo = FALSE--------------------------------------------------------------------
NCORES = 36L
CACHE_DIR = "pomp_cache/"
run_level = 1
NP = switch(run_level, 50, 1e3, 3e3) # number of particles
NMIF_S = switch(run_level, 5, 50, 100) # number of filtering iterations - small
NMIF_L = switch(run_level, 10, 100, 200) # - large
NREPS_EVAL = switch(run_level, 5, 20, 40) # number of replications in likelihood evaluation
NREPS_LOCAL = switch(run_level, 10, 20, 30) # number of replications in local search
NSTART = switch(run_level, 50, 500, 800) # number of starting points in the global search
NSIM = switch(run_level, 50, 100, 500) # number of simulations

PARAMS_FILE = sprintf("seir_params_runlevel=%i.csv", run_level)

suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
  library(doRNG)
  library(tidyverse)
  library(pomp)
})
cl = makeCluster(NCORES)
registerDoParallel(cl)
registerDoRNG(625904618)


## -------------------------------------------------------------------------------------
cases = read.csv("cases_deaths_by_county_date.csv",
                 colClasses = list(Date = "Date")) %>%
  # select records of confirmed cases in Washtenaw County
  filter(!is.na(Date), CASE_STATUS == "Confirmed", COUNTY == "Washtenaw") %>%
  select(-Updated, -CASE_STATUS, -COUNTY, -ends_with("Cumulative")) %>%
  # select records in 2020 only
  filter(Date < as.Date("2021-01-01")) %>%
  # transform `Date` to numeric form
  mutate(Time = 1:n())
covid_t0 = min(cases$Time)


## -------------------------------------------------------------------------------------
seir_step = Csnippet("
  double Beta;
  if(intervention == 1) Beta = b1;
  else if(intervention == 2) Beta = b2;
  else if(intervention == 3) Beta = b3;
  else if(intervention == 4) Beta = b4;
  else Beta = b5;
  
  double dN_SE = rbinom(S, 1 - exp(-Beta * I / N * dt));
  double dN_EI = rbinom(E, 1 - exp(-mu_EI * dt));
  double dN_IR = rbinom(I, 1 - exp(-mu_IR * dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  H += dN_IR;
")

seir_rinit = Csnippet("
  S = nearbyint(eta * N);
  E = 100;
  I = 200;
  H = 0;
")

dmeas <- Csnippet("
  double tol=1.0e-25;
  double mean =rho*H;
  double sd =sqrt(pow(tau*H,2)+rho*H);
  if(Cases>0.0){
    lik=pnorm(Cases+0.5,mean,sd,1,0)-pnorm(Cases-0.5,mean,sd,1,0)+tol;
  } else {
    lik=pnorm(Cases+0.5,mean,sd,1,0)+tol;
  }
  if(give_log) lik=log(lik);
")

rmeas <- Csnippet("
  Cases = rnorm(rho*H, sqrt(pow(tau*H,2)+rho*H));
  if(Cases>0.0){
    Cases=nearbyint(Cases);
  } else {
    Cases=0.0;
  }
")

seir_covar <- covariate_table(
  t = cases$Time,
  intervention = c(rep(1, 23),
                   rep(2, 77),
                   rep(3, 20),
                   rep(4, 75),
                   rep(5, 111)),
  times = "t")

covidSEIR = cases %>% select(Time, Cases) %>%
  pomp(
    times = "Time", t0 = covid_t0,
    rprocess = euler(seir_step, delta.t = 1), # delta.t set to 1 day
    rinit = seir_rinit,
    rmeasure = rmeas,
    dmeasure = dmeas,
    accumvars = "H",
    partrans=parameter_trans(
      log = c("mu_EI", "mu_IR", "tau", "b1", "b2", "b3", "b4", "b5"),
      logit = c("rho", "eta")
    ),
    statenames = c("S", "E", "I", "H"),
    paramnames = c("b1", "b2", "b3", "b4", "b5", "mu_EI", "mu_IR", 
                   "eta", "rho", "N", "tau"),
    covar = seir_covar
  )


## -------------------------------------------------------------------------------------
pop_washtenaw = 367601
params = c(b1 = 3, b2 = 0.5, b3 = 4, b4 = 1.5, b5 = 3.2, 
            mu_EI = 0.1, mu_IR = 0.1, rho = 0.5, eta = 0.09, 
            tau = 0.001, N = pop_washtenaw)
fixed_params = params[c("N", "mu_EI", "mu_IR")]
params_rw.sd = rw.sd(b1 = 0.02, b2 = 0.02, b3 = 0.02, b4 = 0.02, b5 = 0.02, 
                     rho = 0.02, tau = 0.0001, eta = ivp(0.02))

## ---- fig.width = 6, fig.height = 3---------------------------------------------------
if (FALSE) {
  covidSEIR %>%
    simulate(
      params = params,
      nsim = 10, format='data.frame', include.data=TRUE
    ) -> dat_simulated
  
  dat_simulated %>% 
    ggplot(aes(x=Time, y=Cases, group=.id, color=.id=='data')) +
    theme_bw() +
    geom_line() +
    guides(color = FALSE)
}

run_id = 0
registerDoRNG(1235252)
bake(file = sprintf("%srunlevel=%i_%s", CACHE_DIR, run_level, "lik_starting_values.rds"), {
  foreach(i=1:10, .combine = c) %dopar% {
    library(pomp)
    covidSEIR %>% pfilter(params=params,  Np = NP)
  }
}) -> pf
L_pf = pf %>% logLik() %>% logmeanexp(se = TRUE)
print(L_pf)

coef(pf[[1]]) %>% bind_rows() %>%
  bind_cols(loglik = L_pf[1], loglik.se = L_pf[2], id = run_id) %>%
  write_csv(PARAMS_FILE)


## ---- fig.height = 5, fig.width = 8---------------------------------------------------
run_id = 1
registerDoRNG(482947940)
bake(file = sprintf("%srunlevel=%i_%s", CACHE_DIR, run_level, "local_search.rds"),{
  foreach(i = 1:NREPS_LOCAL, .combine = c) %dopar% {
    suppressPackageStartupMessages({
      library(tidyverse)
      library(pomp)
    })
    covidSEIR %>%
      mif2(
        params = params,
        Np = NP, Nmif = NMIF_S,
        cooling.fraction.50 = 0.5,
        rw.sd = params_rw.sd
      )
  } -> mifs_local
  attr(mifs_local,"ncpu") <- getDoParWorkers()
  mifs_local
}) -> mifs_local
t_loc <- attr(mifs_local,"system.time")
ncpu_loc <- attr(mifs_local,"ncpu")
cat(sprintf("Local search (iterative filtering) finished in %4.3f seconds\n", t_loc["elapsed"]))

# likelihood est. for local search results
registerDoRNG(900242057)
bake(file = sprintf("%srunlevel=%i_%s", CACHE_DIR, run_level, "lik_local.rds"),{
  foreach(mf = mifs_local, .combine = rbind) %dopar% {
    suppressPackageStartupMessages({
      library(tidyverse)
      library(pomp)
    })
    ll = replicate(NREPS_EVAL, logLik(pfilter(mf, Np = NP))) %>% 
         logmeanexp(se = TRUE)
    coef(mf) %>% bind_rows() %>% bind_cols(loglik = ll[1], loglik.se = ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) -> results
t_local <- attr(results,"system.time")
ncpu_local <- attr(results,"ncpu")
cat(sprintf("Local search (likelihood estimation) finished in %4.3f seconds\n", t_local["elapsed"]))


read.csv(PARAMS_FILE) %>%
  bind_rows(results %>% mutate(id = run_id)) %>%
  arrange(-loglik) %>%
  write_csv(PARAMS_FILE)
# results %>% arrange(-loglik) %>% head %>% knitr::kable(digits = 3)

## -------------------------------------------------------------------------------------
run_id = 2

# create a box of starting values (for parameters)
set.seed(2062379496)
guesses = runif_design(
  lower = c(b1 = 0, b2 = 0, b3 = 0, b4 = 0, b5 = 0, 
            rho = 0, eta = 0, tau = 0),
  upper = c(b1 = 10, b2 = 10, b3 = 10, b4 = 10, b5 = 10, 
            rho = 1, eta = 0.3, tau = 0.1),
  nseq = NSTART
)

mf1 = mifs_local[[1]] # take the output of previous IF process (local search)

bake(file = sprintf("%srunlevel=%i_%s", CACHE_DIR, run_level, "global_search.rds"),{
  registerDoRNG(1270401374)
  foreach(guess=iter(guesses, "row"), .combine = rbind) %dopar% {
    suppressPackageStartupMessages({
      library(tidyverse)
      library(pomp)
    })
    mf = mf1 %>% # cooling.fraction.50 = 0.5
      mif2(params = c(unlist(guess), fixed_params), Nmif = NMIF_L) %>%
      mif2(Nmif = NMIF_L) %>%
      mif2(Nmif = NMIF_L)
    mf = mf %>%
      mif2(Nmif = NMIF_L, cooling.fraction.50 = 0.3) %>%
      mif2(Nmif = NMIF_L, cooling.fraction.50 = 0.3) %>%
      mif2(Nmif = NMIF_L, cooling.fraction.50 = 0.1) %>%
      mif2(Nmif = NMIF_L, cooling.fraction.50 = 0.1)
    ll = replicate(NREPS_EVAL, mf %>% pfilter(Np = NP) %>% logLik()) %>%
      logmeanexp(se = TRUE)
    coef(mf) %>% bind_rows() %>%
      bind_cols(loglik = ll[1],loglik.se = ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) %>%
  filter(is.finite(loglik)) -> results
t_global <- attr(results, "system.time")
ncpu_global <- attr(results, "ncpu")
cat(sprintf("Global search finished in %4.3f seconds\n", t_global["elapsed"]))

read.csv(PARAMS_FILE) %>%
  bind_rows(results %>% mutate(id = run_id)) %>%
  filter(is.finite(loglik)) %>%
  arrange(-loglik) %>%
  write_csv(PARAMS_FILE)

if (FALSE) {
  all = read_csv(PARAMS_FILE) %>%
    filter(id <= 2) %>%
    filter(loglik > max(loglik) - 50) %>%
    bind_rows(guesses) %>%
    mutate(type = if_else(is.na(loglik), "guess", "result")) %>%
    arrange(type)
  
  pairs(~loglik + b1 + b2 + b3 + b4 + b5, data = all,
        col = ifelse(all$type == "guess", grey(0.5), "red"), pch = 16)
  pairs(~loglik + mu_EI + mu_IR + eta + rho + tau, data = all,
        col = ifelse(all$type == "guess", grey(0.5), "red"), pch = 16)
}


# Profile likelihood for rho --------------------------------------------------
run_id = 3

guesses = read.csv(PARAMS_FILE) %>%
  group_by(cut = round(rho, 2)) %>%
  filter(rank(-loglik) <= 10) %>%
  ungroup() %>%
  select(-cut, -loglik, -loglik.se)

rw.sd_rho_fixed = rw.sd(
  b1 = 0.02, b2 = 0.02, b3 = 0.02, b4 = 0.02, b5 = 0.02, 
  rho = 0, tau = 0.0001, eta = ivp(0.02)
)

mf1 = mifs_local[[1]]
registerDoRNG(2105684752)
bake(file = sprintf("%srunlevel=%i_%s", CACHE_DIR, run_level, "profile_rho.rds"), {
  foreach(guess = iter(guesses, "row"), .combine = rbind) %dopar% {
    suppressPackageStartupMessages({
      library(tidyverse)
      library(pomp)
    })
    mf = mf1 %>%
      mif2(params = guess, rw.sd = rw.sd_rho_fixed) %>% 
      mif2(Nmif = NMIF_L, cooling.fraction.50 = 0.3) %>%
      mif2(cooling.fraction.50 = 0.1)
    ll = replicate(NREPS_EVAL, mf %>% pfilter(Np = NP) %>% logLik()) %>%
      logmeanexp(se = TRUE)
    coef(mf) %>% bind_rows() %>% bind_cols(loglik = ll[1],loglik.se = ll[2])
  } -> results
  attr(results, "ncpu") = getDoParWorkers()
  results
}) -> results
t_profile = attr(results, "system.time")
ncpu_profile = attr(results,"ncpu")
cat(sprintf(
  "Profile likelihood for rho finished in %4.3f seconds\n", t_profile["elapsed"]
))

read.csv(PARAMS_FILE) %>%
  bind_rows(results %>% mutate(id = run_id)) %>%
  filter(is.finite(loglik)) %>%
  arrange(-loglik) %>%
  write_csv(PARAMS_FILE)

# -----------------------------------------------------------------------------
stopCluster(cl) # shut the cluster down
