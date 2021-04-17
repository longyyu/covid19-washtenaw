# test measurement process
# -------------------------------------------
NP = 100
suppressPackageStartupMessages({
  library(tidyverse)
  library(pomp)
})

plot_simulation = function(sim_dat) {
  sim_dat %>%
    ggplot() +
    theme_bw() +
    geom_line(aes(Time, Cases, group = .id, 
                  color = (.id == "data"), alpha = (.id == "data"), 
                  linetype = (.id == "data"))) +
    scale_color_manual(values = c("#18bc9c", "#c61919")) +
    scale_alpha_manual(values = c(0.5, 1)) +
    scale_linetype_manual(values = c(5, 1)) +
    guides(color = FALSE, linetype = FALSE, alpha = FALSE)
} # plot_simulation()


# -------------------------------------------
cases = read.csv("cases_deaths_by_county_date.csv",
                 colClasses = list(Date = "Date")) %>%
  filter(!is.na(Date), CASE_STATUS == "Confirmed", COUNTY == "Washtenaw") %>%
  select(-Updated, -CASE_STATUS, -COUNTY, -ends_with("Cumulative")) %>% 
  filter(Date >= as.Date("2020-06-01"), Date < as.Date("2021-03-01"))

cases = cases %>%
  # transform Date to numeric form
  mutate(Time = 2020 + as.numeric(Date - as.Date("2020-01-01"))/366)
t0 = min(cases$Time)


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
  E = 100;
  I = 0;
  H = 0;
")

dmeas = Csnippet("
  lik = dbinom(Cases, H, rho, give_log);
  // lik = dpois(Cases, H * rho, give_log);
")

rmeas = Csnippet("
  Cases = rbinom(H, rho);
  // Cases = rpois(H *rho);
")

covidSEIR = cases %>% select(Time, Cases) %>%
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
pop_washtenaw = 367000
# # full sample
# params = c(Beta = 70, mu_EI = 80, mu_IR = 14, rho = 0.12, eta = 0.5, N = pop_washtenaw)
# subsample
params = c(Beta = 70, mu_EI = 80, mu_IR = 10, rho = 0.12, eta = 0.5, N = pop_washtenaw)

dat_simulated = covidSEIR %>%
  simulate(params = params, nsim = 10, format = "data.frame", include.data = TRUE)

dat_simulated %>% plot_simulation()

# test measurement process with simulated data --------------------------------
logliks = sapply(c("data", as.list(1:10)), function(group) {
  covidSEIR_sim = dat_simulated %>% filter(.id == group) %>% select(Time, Cases) %>%
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
  pf = sapply(1:10, function(i) {
    covidSEIR_sim %>% pfilter(params = params, Np = NP)
  })
  sapply(pf, logLik) %>% logmeanexp(se = TRUE)
})

logliks # the first element is the log likelihood estimate for the original data
        # the others are estimates for the simulated data
