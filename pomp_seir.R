# when altering fixed parameters, update:
# (1) fixed_params (2) params_rw.sd (3) guesses
# -------------------------------------------
NCORES = 25L
CACHE_DIR = "pomp_cache/"
DO_PLOT = FALSE

run_level = 2
NP = switch(run_level, 50, 1e3, 3e3) # number of particles
NMIF_S = switch(run_level, 5, 50, 100) # number of filtering iterations - small
NMIF_L = switch(run_level, 10, 75, 150) # - large
NREPS_EVAL = switch(run_level, 5, 10, 20) # number of replications in likelihood evaluation
NREPS_LOCAL = switch(run_level, 10, 20, 30) # number of replications in local search
NSTART = switch(run_level, 50, 500, 800) # number of starting points in the global search
NSIM = switch(run_level, 50, 100, 500) # number of simulations

PARAMS_FILE = sprintf("seir_params_runlevel=%i.csv", run_level)

# -------------------------------------------
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


# -------------------------------------------
cases = read.csv("cases_deaths_by_county_date.csv",
                 colClasses = list(Date = "Date")) %>%
  filter(!is.na(Date), CASE_STATUS == "Confirmed", COUNTY == "Washtenaw") %>%
  select(-Updated, -CASE_STATUS, -COUNTY, -ends_with("Cumulative"))

cases = cases %>%
  # transform Date to numeric form
  mutate(Time = 2020 + as.numeric(Date - as.Date("2020-01-01"))/366)
t0 = min(cases$Time) # 2020.164, i.e. 2020-03-01


# -------------------------------------------
seir_step = Csnippet("
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
  E = 1;
  I = 0;
  H = 0;
")

dmeas = Csnippet("
  lik = dbinom(Cases, H, rho, give_log);
")

rmeas = Csnippet("
  Cases = rbinom(H, rho);
")

measSEIR = cases %>% select(Time, Cases) %>%
  pomp(
    times = "Time", t0 = t0,
    rprocess = euler(seir_step, delta.t = 1),
    rinit = seir_rinit,
    rmeasure = rmeas,
    dmeasure = dmeas,
    accumvars = "H",
    partrans=parameter_trans(
      log = c("Beta", "mu_EI", "mu_IR"),
      logit = c("rho","eta")
    ),
    statenames=c("S", "E", "I", "H"),
    paramnames=c("Beta", "mu_EI", "mu_IR", "eta", "rho", "N")
  )


# -------------------------------------------
pop_washtenaw = 367601
params = c(Beta = 1e05, mu_EI = 1, mu_IR = 1, rho = 0.12, eta = 0.8, N = pop_washtenaw)
fixed_params = params[c("N")]
params_rw.sd = rw.sd(Beta = 0.02, mu_EI = 0.02, mu_IR = 0.02,
                     rho = 0.02, eta = ivp(0.04)) 
pairs_formula = ~loglik + Beta + mu_EI + mu_IR + eta + rho


## ---- echo = FALSE, eval = FALSE------------
## # run a simulation to check that rprocess and rmeasure work
# y = measSEIR %>% simulate(params = params, nsim = 10, format = "data.frame")
## 
## # run a pfilter and check that (rprocess and) dmeasure worksi
# pf = measSEIR %>% pfilter(Np = 1000, params = params)


# -------------------------------------------
# Likelihood est. for starting values of parameters
run_id = 0
system.time({
  foreach(i = 1:10, .combine = c) %dopar% {
    suppressPackageStartupMessages({
      library(pomp)
      library(tidyverse)
    })
    measSEIR %>% pfilter(params = params, Np = NP)
  } -> pf
})

L_pf = pf %>% logLik() %>% logmeanexp(se=TRUE)
print(L_pf)

coef(pf[[1]]) %>% bind_rows() %>%
  bind_cols(loglik = L_pf[1], loglik.se = L_pf[2], id = run_id) %>%
  write_csv(PARAMS_FILE)


run_id = 1
registerDoRNG(482947940)
bake(file = sprintf("%srunlevel=%i_%s", CACHE_DIR, run_level, "local_search.rds"), {
  foreach(i = 1:NREPS_LOCAL, .combine = c) %dopar% {
    suppressPackageStartupMessages({
      library(tidyverse)
      library(pomp)
    })
    measSEIR %>%
      mif2(
        params = params,
        Np = NP, Nmif = NMIF_S,
        cooling.fraction.50 = 0.5,
        rw.sd = params_rw.sd
      )
  } -> mifs_local
  attr(mifs_local,"ncpu") = getDoParWorkers()
  mifs_local
}) -> mifs_local
t_loc = attr(mifs_local,"system.time")
ncpu_loc = attr(mifs_local,"ncpu")
cat(sprintf("Local search (IF2) finished in %4.3f seconds\n", t_loc["elapsed"]))


registerDoRNG(900242057)
bake(file = sprintf("%srunlevel=%i_%s", CACHE_DIR, run_level, "lik_local.rds"),{
  foreach(mf = mifs_local, .combine = rbind) %dopar% {
    suppressPackageStartupMessages({
      library(tidyverse)
      library(pomp)
    })
    ll = replicate(NREPS_EVAL, logLik(pfilter(mf, Np = NP))) %>% logmeanexp(se = TRUE)
    coef(mf) %>% bind_rows() %>% bind_cols(loglik = ll[1], loglik.se = ll[2])
  } -> results
  attr(results,"ncpu") = getDoParWorkers()
  results
}) -> results
t_local = attr(results,"system.time")
ncpu_local = attr(results,"ncpu")
cat(sprintf("Local search (likelihood est.) finished in %4.3f seconds\n", t_local["elapsed"]))


read.csv(PARAMS_FILE) %>%
  bind_rows(results %>% mutate(id = run_id)) %>%
  arrange(-loglik) %>%
  write_csv(PARAMS_FILE)


if (DO_PLOT) {
  mifs_local %>%
    traces() %>%
    melt() %>%
    filter(!variable %in% names(fixed_params)) %>%
    ggplot(aes(x = iteration, y = value, group = L1, color = factor(L1))) +
    theme_bw() +
    geom_line() +
    guides(color = FALSE) +
    facet_wrap(~variable, scales="free_y")
  
  pairs(pairs_formula, data = results, pch = 16)
}

## -------------------------------------------
run_id = 2

# create a box of starting values (for parameters)
set.seed(2062379496)
guesses = runif_design(
  lower = c(Beta = 1e04, mu_EI = 0, mu_IR = 0 , rho = 0.05, eta = 0.5),
  upper=c(Beta = 1e06, mu_EI = 10, mu_IR = 10, rho = 0.3, eta = 1),
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
      # mif2(Nmif = NMIF_L) %>%
      mif2(Nmif = NMIF_L)
    mf = mf %>%
      mif2(Nmif = NMIF_L, cooling.fraction.50 = 0.3) %>%
      # mif2(Nmif = NMIF_L, cooling.fraction.50 = 0.3) %>%
      # mif2(Nmif = NMIF_L, cooling.fraction.50 = 0.1) %>%
      mif2(Nmif = NMIF_L, cooling.fraction.50 = 0.1)
    ll = replicate(NREPS_EVAL, mf %>% pfilter(Np = NP) %>% logLik()) %>%
      logmeanexp(se = TRUE)
    coef(mf) %>% bind_rows() %>%
      bind_cols(loglik = ll[1],loglik.se = ll[2])
  } -> results
  attr(results,"ncpu") = getDoParWorkers()
  results
}) %>%
  filter(is.finite(loglik)) -> results
t_global = attr(results,"system.time")
ncpu_global = attr(results,"ncpu")
cat(sprintf("Global search finished in %4.3f seconds\n", t_global["elapsed"]))


read.csv(PARAMS_FILE) %>%
  bind_rows(results %>% mutate(id = run_id)) %>%
  filter(is.finite(loglik)) %>%
  arrange(-loglik) %>%
  write_csv(PARAMS_FILE)


if (DO_PLOT) {
  read.csv(PARAMS_FILE) %>%
    filter(id <= 2) %>%
    # filter(loglik > max(loglik) - 50) %>%
    bind_rows(guesses) %>%
    mutate(type = if_else(is.na(loglik), "guess", "result")) %>%
    arrange(type) -> all
  
  pairs(pairs_formula, data = all,
        col = ifelse(all$type == "guess", grey(0.5), "red"), pch = 16)
}

## -------------------------------------------
stopCluster(cl) # shut the cluster down

