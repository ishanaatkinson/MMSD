# Dynamic multi-spcies, multi-patch model that is entangled at each timestep
# to represent malaria transmission in Bolivar, Amazonas and Sucre in 
# Venezuela 


# Load packages ####
library(pacman)
p_load(deSolve, tidyverse, doParallel, manipulate,readxl, scales)

# Input definitions ####
# define the number of days to run the model
times <- seq(0, 365*11 -1, 1)



# Initial conditions ####

# Population of venezuela in 2000 ~ 25 million
# Assuming the same proportion of populations across each state in
# 2020 to 2000: Sucre =  3.446%, Bolivar=6.048%, Amazonas=0.633%. 
# Hence, Sucre ~ 861500, Bolivar&Amazonas ~ 1670250 split for 76% 
# viv population and 24% falc population 

# These have been chosen arbitrarily using naive/basic model fitting or 
# calibration techniques to fit expected case data published by the WHO.
# The population and population ratios between each patch is realistic 
# according to World Bank and Statista population data. More details can 
# be found in the report.

start <- c(# P. falciparum in Bolivar and Amazonas
           Sf1 = 350000, 
           Ef1 = 500, 
           Af1 = 500, 
           Cf1 = 500, 
           Tf1= 400, 
           Rf1 = 50000, 
           E2f1=500,
           
           # P. vivax in Bolivar and Amazonas
           Sv1 = 1000000, 
           Ev1 = 3000, 
           Av1 = 1000, 
           Cv1 = 800, 
           Tv1=2000, 
           Lv1=100000, 
           Rv1 = 100000, 
           E2v1=3000,
           
           Cincf1=0, 
           Cincv1=0, 
           Ctrtf1=0, 
           Ctrtv1=0, 
           
           # P. vivax in Sucre
           Sv2 = 600000, 
           Ev2 = 1000, 
           Av2 = 1000, 
           Cv2 = 1000, 
           Tv2=1000, 
           Lv2=50000, 
           Rv2 = 50000, 
           E2v2=1000,
           
           Cincv2=0, 
           Ctrtv2=0,
           CumInc=0, 
           tests = 0
)

# Parameter table ####

# These parameters have mostly been chosen reasonably through research
# of peer-reviewed published literature. In the case where expected values 
# could not be found via this way, either an educated opinion or basic
# model calibration techniques were used. More details can be found in
# the report.  

parms <- c(a = 1/4, # human feeding rate per mosquito per day
           b = 0.55, # transmission efficacy M->H
           c=0.55, # transmission efficacy H->M
           gamma_m = 1/10, # 1/latent period in mosquitos
           mu_m = 1/13, # natural death rate of mosquitos
           mu_h = 1/(70.06*365), # natural birth/death rate
           gamma_hf = 1/10.5, # latent period, falciparum
           gamma_hf2 = 1/(10.5*1.25), # latent period for reinfection, falc
           gamma_hv = 1/15, # latent period, vivax
           gamma_hv2 = 1/(15*1.25), # latent period for reinfection, vivax
           tau_f = 1/5, # rate of ACT treatment seeking
           tau_vba = 1/10, # rate of getting treated in Bolivar&Amazonas
           tau_vs=1/7, # rate of getting treated in Sucre
           rhof = 1/(365), # loss of immunity to falciparum
           rho_v = 1/(9*30), # loss of immunity to vivax
           m=0.85, # relative number of mosquitos to humans in b&a
           m2 = 0.45, # relative number of mosquitos to humans in Sucre 
           relv=0.9,  # relative transmission parameter V
           pa1 = 0.602, # probability of asymptomatic infection, vivax
           pa3 = 0.1, # probability of asymptomatic infection, falciparum
           delta1 = 1/34.4, # natural rate of recovery, asymptomatic vivax 
           delta2 = 1/130, # natural rate of recovery, asymptomatic falciparum 
           omega_f = 1/20, # rate of symptom loss, falciparum
           omega_v = 1/20, # rate of symptom loss, vivax 
           zeta_a = 12.6/27, # relative infectiousness of A to C infections
           zeta_t = 12.6/(27*3), # relative infectiousness of T to C  
           phi = 3/12, # cases peak in Bolivar&Amazonas in feb/april 
           tau3=1/3, # rate of trt recovery (1/pct) ACT
           d_hyp=1/400, # death rate hypnozoites
           relr= 0.0056, # vivax relapse rate
           p_hyp=0.13, # prob of recovery with hyp after PQ+CQ treatment 
           p_hyp2=0.68, # prob of recovery with hyp w/o treatment 
           tau2=1/14, #recovery rate with 14day dose of primaquine
           t4=1, # Entanglement 4 - cross immunity
           beta12=0.8, # cross-immunity from f to v
           beta21=0.8, # cross-immunity from v to f
           phi1 = 2/12, # dec-april dry season so feb=peak dry month
           phi2 = 8/12, # april to nov wet season so aug=peak of wet month
           vivpercent = 0.76, # percent vivax in Bolivar&Amazonas
           itn_own = 0.76, # prop of population at risk who own a LLIN
           itn_use = 0.432, #prop of those with LLIN that use it effectively
           itn_eff = 0.25, #effectiveness of nets in reducing infectious bites
           pimm1 = 0.8, # partial immunity for reinfection, vivax
           pimm2 = 0.8, # partial immunity for reinfection, falciparum 
           avail_ba = 0.6, # avg. test availability in Bolivar and Amazonas
           avail_s = 0.8, # average test availability in Sucre
           testsens_ms=0.957, # microscopy test sensitivity
           testsens_rdt = 0.923 # rdt test sensitivity
)

# Read data ####

# read data for LLIN net distribution and case data from WHO

data<-read.csv("https://raw.githubusercontent.com/ishanaatkinson/MMSD/main/data2.csv")


# Define entangled model function for Bolivar, Amazonas and Sucre ####
venmodel<-function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         Pv2 <- (Sv2+Ev2+Av2+Cv2+Tv2+Lv2+Rv2+E2v2) # population, Sucre to be used later on 
         
         # Between state migration
         # migration to Bolivar&Amazonas function from Sucre 
         mig1 = 0.02+0.2*cos(pi*(t/365-phi1))^2
         
         # migration to Sucre from Bolivar&Amazonas function (less than 
         # migration to Bolivar&Amazonas)
         mig2 = 0.01+0.01*cos(pi*(t/365-phi2))^2
         
         # drug availability fluctuations
         drugstock = 0.5+0.5*cos(pi*(t/365-0.5))^2

         # cases
         data_t<-data$timestep
         cases<-data$cases
         vencases<-approx(data_t,cases,t)$y
         
         # patch 1 - Bolivar&Amazonas
         
         # define variables
         
         # falciparum population Bolivar&Amazonas
         Pf1 <- (Sf1+Ef1+Af1+Cf1+Tf1+Rf1+E2f1) 
         
         # vivax population Bolivar&Amazonas
         Pv1 <- (Sv1+Ev1+Av1+Cv1+Tv1+Lv1+Rv1+E2v1) 
         Pba <- Pf1+Pv1 # Population in Bolivar and Amazonas
         
  
         # All states use ITNs
         #proportion of the population covered by NEW nets
         itn_cov = (data$nets)*1.8/(Pba+Pv2) 
         itn_effcov = itn_use*itn_eff*pmin(itn_cov, 1)
         
         itn<-approx(data_t, itn_effcov, t)$y
         
         # exponential increase in importation from outside Venezuela for mining purposes 
         expon = 1/7 * (25*exp(t/365)-2.5)/(6e7)
         ratemig = expon # 1/7 days migration rate into Bolivar and Amazonas with exponential growth component due to increase in mining activities
         
         Inf_f1 <- Cf1 + zeta_a*Af1+zeta_t*Tf1 # infectious with falciparum in Bolivar&Amazonas
         Inf_v1 <- Cv1 + zeta_a*Av1+zeta_t*Tv1 # infectious with vivax in Bolivar&Amazonas
         seas <- 0.05+0.5*cos(pi*(t/365 - phi))^2 # seasonality in Bolivar&Amazonas
         
         lambdaf1 <- (1-itn)*seas*(a^2*b*c*m*Inf_f1/Pf1)/(a*c*Inf_f1/Pf1+mu_m)*(gamma_m/(gamma_m+mu_m)) + ratemig
         lambdav1 <- (1-itn)*seas*relv*(a^2*b*c*m*Inf_v1/Pv1)/(a*c*Inf_v1/Pv1+mu_m)*(gamma_m/(gamma_m+mu_m)) + ratemig
         lambdav21 <- lambdav1*pimm1 # partial immunity once recovered
         lambdaf21 <- lambdaf1*pimm2 # partial immunity once recovered
         
         alpha1 = (1-t4*beta12*Inf_f1/Pf1)
         alpha2 = (1-t4*beta21*Inf_v1/Pv1)
         
         # P. falciparum equations 
         
         # susceptible
         dSf1 <- mu_h*Pf1 -
           lambdaf1*Sf1*alpha2 +
           rhof*Rf1 -
           mu_h*Sf1 - 
           mig2*Sf1 + 
           (1-vivpercent)*mig1*Sv2 
         
         # exposed
         dEf1 <- lambdaf1*Sf1*alpha2 - 
           gamma_hf*Ef1 - 
           mu_h*Ef1 
         
         # asymptomatic
         dAf1 <- pa3*gamma_hf*Ef1 + 
           pa3*gamma_hf2*E2f1 + 
           omega_f*Cf1 - 
           delta2*Af1 -
           mu_h*Af1 
         
         # clinical 
         dCf1 <- (1-pa3)*gamma_hf*Ef1 +
           (1-pa3)*gamma_hf2*E2f1 -
           omega_f*Cf1 -
           tau_f*avail_ba*testsens_ms*drugstock*Cf1 - 
           mu_h*Cf1      
         
         # treated
         dTf1 <- tau_f*avail_ba*testsens_ms*drugstock*Cf1 - 
           tau3*Tf1 - 
           mu_h*Tf1                 
         
         # recovered
         dRf1 <- tau3*Tf1 + 
           delta2*Af1 - 
           lambdaf21*Rf1*alpha2 - 
           rhof*Rf1 - 
           mu_h*Rf1 
         
         # re-exposed
         dE2f1 <- lambdaf21*Rf1*alpha2 - 
           gamma_hf2*E2f1 - 
           mu_h*E2f1 
         
         #P. vivax equations
         
         # susceptible
         dSv1 <- mu_h*Pv1 - 
           lambdav1*Sv1*alpha1 +
           rho_v*Rv1 + 
           d_hyp*Lv1 - 
           mu_h*Sv1 - 
           mig2*Sv1 + 
           vivpercent*mig1*Sv2 
         
         # exposed
         dEv1 <- lambdav1*Sv1*alpha1 - 
           gamma_hv*Ev1 - 
           mu_h*Ev1 
         
         # asymptomatic
         dAv1 <- pa1*gamma_hv*Ev1 + 
           pa1*gamma_hv2*E2v1 + 
           omega_v*Cv1  - 
           delta1*Av1 - 
           mu_h*Av1 
         
         # clinical
         dCv1 <- (1-pa1)*gamma_hv*Ev1 + 
           (1-pa1)*gamma_hv2*E2v1 -
           omega_v*Cv1 - 
           tau_vba*avail_ba*testsens_ms*drugstock*Cv1 - 
           mu_h*Cv1  
         
         # treated
         dTv1 <- tau_vba*avail_ba*testsens_ms*drugstock*Cv1 - 
           tau2*Tv1 - 
           mu_h*Tv1               
         
         # recover with hypnozoites
         dLv1 <- p_hyp*tau2*Tv1  - 
           relr*Lv1 - 
           d_hyp*Lv1 - 
           mu_h*Lv1 +
           p_hyp2*delta1*Av1 
         
         # recover without hypnozoites
         dRv1 <- (1-p_hyp)*tau2*Tv1 - 
           rho_v*Rv1 - 
           lambdav21*Rv1*alpha1 - 
           mu_h*Rv1 +
           (1-p_hyp2)*delta1*Av1 
         
         # re-exposed or relapse
         dE2v1 <- lambdav21*Rv1*alpha1 + 
           relr*Lv1 - 
           gamma_hv2*E2v1 -
           mu_h*E2v1 
         
         #Counters
         dCincf1= lambdaf1*alpha2*Sf1 + lambdaf21*alpha2*Rf1
         dCincv1= lambdav1*alpha1*Sv1 + lambdav21*alpha1*Rv1 + relr*Lv1 
         dCtrtf1= tau_f*avail_ba*testsens_ms*drugstock*Cf1
         dCtrtv1= tau_vba*avail_ba*testsens_ms*drugstock*Cv1
         
         # PATCH 2 -Sucre 
         
         Pv2 <- (Sv2+Ev2+Av2+Cv2+Tv2+Lv2+Rv2+E2v2)
         Inf_v2 <- Cv2 + zeta_a*Av2+zeta_t*Tv2
         
         lambdav2 <- (1-itn)*relv*(a^2*b*c*m2*Inf_v2/Pv2)/(a*c*Inf_v2/Pv2+mu_m)*(gamma_m/(gamma_m+mu_m)) 
         lambdav22 <- lambdav2*pimm1 # partial immunity once recovered
         
         #P. vivax equations
         
         # susceptible
         dSv2 <- mu_h*Pv2 - 
           lambdav2*Sv2 +
           rho_v*Rv2 + 
           d_hyp*Lv2 - 
           mu_h*Sv2 - 
           mig1*Sv2 + 
           mig2*(Sf1+Sv1)
         
         # exposed
         dEv2 <- lambdav2*Sv2 - 
           gamma_hv*Ev2 - 
           mu_h*Ev2
         
         # asymptomatic
         dAv2 <- pa1*gamma_hv*Ev2 + 
           pa1*gamma_hv2*E2v2 + 
           omega_v*Cv2  - 
           delta1*Av2 - 
           mu_h*Av2 
         
         # clinical
         dCv2 <- (1-pa1)*gamma_hv*Ev2 + 
           (1-pa1)*gamma_hv2*E2v2 -
           omega_v*Cv2 - 
           tau_vs*avail_s*testsens_ms*drugstock*Cv2 - 
           mu_h*Cv2  
         
         # treated
         dTv2 <- tau_vs*avail_s*testsens_ms*drugstock*Cv2 - 
           tau2*Tv2 - 
           mu_h*Tv2               
         
         # recover with hypnozoites
         dLv2 <- p_hyp*tau2*Tv2  - 
           relr*Lv2 +
           p_hyp2*delta1*Av2 - 
           d_hyp*Lv2 - 
           mu_h*Lv2 
         
         # recover without hypnozoites
         dRv2 <- (1-p_hyp)*tau2*Tv2 + 
           (1-p_hyp2)*delta1*Av2 - 
           rho_v*Rv2 - 
           lambdav22*Rv2 - 
           mu_h*Rv2 
         
         # re-exposed or relapse
         dE2v2 <- lambdav22*Rv2+ 
           relr*Lv2 - 
           gamma_hv2*E2v2 -
           mu_h*E2v2 
         
         #Counters
         dCincv2<-lambdav2*Sv2+lambdav22*Rv2 + relr*Lv2 
         dCtrtv2<- tau_vs*avail_s*testsens_ms*drugstock*Cv2
         dCumInc <- dCincf1+dCincv1+dCincv2
         dtests <- tau_vs*avail_s*Cv2 + tau_f*avail_ba*Cf1 + tau_vba*avail_ba*Cv1
         
         # return the rate of change
         list(c(dSf1, dEf1, dAf1, dCf1, dTf1, dRf1, dE2f1, dSv1, dEv1, dAv1, dCv1,dTv1,dLv1,dRv1,dE2v1, dCincf1, dCincv1,dCtrtf1, dCtrtv1,dSv2, dEv2, dAv2, dCv2,dTv2,dLv2,dRv2,dE2v2, dCincv2, dCtrtv2, dCumInc, dtests))
       }
  ) 
}


# Run the entangled model ####
hm<-ode(times=times, y=start, func=venmodel,parms=parms)

# Store data ####
df1<-as_tibble(as.data.frame(hm)) %>% 
  mutate(# populations
         Pf1 = (Sf1+Ef1+Af1+Cf1+Tf1+Rf1+E2f1),
         Pv1 = (Sv1+Ev1+Av1+Cv1+Tv1+Lv1+Rv1+E2v1),
         Pv2 = (Sv2+Ev2+Av2+Cv2+Tv2+Lv2+Rv2+E2v2),
         
         # incidences and treatments given per timestep
         Incf1 = c(0, diff(Cincf1)), 
         Incv1 = c(0, diff(Cincv1)),
         Trtf1 = c(0, diff(Ctrtf1)),
         Trtv1 = c(0, diff(Ctrtv1)),
         Incv2 = c(0, diff(Cincv2)),
         Trtv2 = c(0, diff(Ctrtv2)),
         
         # total tests and treatments given 
         TotTests = c(0, diff(tests)),
         TotPQCQ = Trtv1 + Trtv2,
         
         # total cases at each timestep
         totcases2 = Incf1+Incv1+Incv2,
         
         # total population
         P = Pv1+Pf1+Pv2, 
         
         # total vivax population 
         totviv = Incv1+Incv2,
         
         # ratio of vivax cases to all cases
         ratios = totviv/(totviv+Incf1)) %>%  
  pivot_longer(names_to = "variable", cols = !1) %>%
  mutate(model="Venezuela")


# Plot initial graphs for all years ####

# Each plot from the model (not population or annual sums) are shown from timestep=365 
# onward to allow for the model to settle 


#Population check
df1 %>% 
  filter(variable %in% c("Pf1", "Pv1", "Pv2", "P")) %>% 
  ggplot()+
  geom_line(aes(x = time, y=value))+
  theme_minimal() +
  labs(title = "Populations", y =("population")) +
  facet_wrap(~variable)


# total cases through time predicted by model 
df1 %>% 
  filter(variable %in% c("totcases2"), time>365) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Model predicted combined incidence", y =("Population"), colour="Legend")+
  facet_wrap(~model)


# incidence through time predicted by model for each species and patch
df1 %>% 
  filter(variable %in% c("Incf1", "Incv1", "Incv2"), time>365) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Model predicted incidence for each species and location", y =("Population"), colour="Legend")+
  facet_wrap(~model)

# number of treated infections predicted by model for each species and patch 

df1 %>% 
  filter(variable %in% c("Trtf1", "Trtv1", "Trtv2"), time>365) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Model predicted total treatments given for \neach species and location", y =("population"), colour="Legend")+
  facet_wrap(~model)


# number of vivax cases compared to number of falciparum cases


df1 %>% 
  filter(variable %in% c("Incf1", "totviv"), time>365) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Model predicted incidence for \neach species", y =("population"), colour="Legend")+
  facet_wrap(~model)


# ratio of P. vivax to P. falciparum 
df1 %>% 
  filter(variable %in% c("ratios"), time>365) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Incidence ratios of P. vivax to P. falciparum", y =("population"), colour="Legend")+
  facet_wrap(~model)

prop_of_vivax <- mean((df1 %>% filter(variable %in% c("ratios"), time > 365))$value)
#0.8091381
# The actual proportion is 76% so on average, the model predicts a value
# reasonably close to it


#t <- c(0:4014)
#t1 <- seq(0,4015+1, 1)
t1<-seq(0,365*11,1) # length 4016
data_t<-data$timestep
cases<-data$cases
vencases<-approx(data_t,cases,times)$y

cases <- cases[-1]
nets <- data$nets[-1]

#Human Compartments
`%nin%` = Negate(`%in%`)

ts <- c(2008:2018)

df2 <- as.data.frame(df1 %>% filter(variable %in% c("totcases2")))$value
yearlysums <- unname(tapply(df2, (seq_along(df2)-1) %/% 365, sum))

# yearlycases predicted and actual yearly cases

ggplot() +
  geom_line(aes(x=ts, y=cases, color="Actual incidence")) + 
  geom_line(aes(x=ts, y=yearlysums, color="Model predicted incidence")) + 
  labs(title="Actual and predicted annual \nsuspected incidence (base scenario)",
       x ="Time", y = "Total incidence", color = "Legend")+ 
  scale_x_continuous(breaks= pretty_breaks())

# Sensitivity analysis on relv ####

# here I explore how much relv affects transmission and incidence

import_vector<-seq(0.8,1.05,by=0.05)
result_import<-matrix(0,nrow = 1,ncol = 3)
colnames(result_import)<-c("time","yearlysums","relv")
result_import<-as.data.frame(result_import)
yearly_sums_sa<-c()
for (i in 1:length(import_vector)){
  parms["relv"]<-import_vector[i]
  out_aux <- ode(times=times, y=start, func=venmodel,parms=parms)
  
  df_aux<-as_tibble(as.data.frame(out_aux)) %>% 
    mutate(Pf1 = (Sf1+Ef1+Af1+Cf1+Tf1+Rf1+E2f1),
           Pv1 = (Sv1+Ev1+Av1+Cv1+Tv1+Lv1+Rv1+E2v1),
           Pv2 = (Sv2+Ev2+Av2+Cv2+Tv2+Lv2+Rv2+E2v2),
           Incf1 = c(0, diff(Cincf1)),
           Incv1 = c(0, diff(Cincv1)),
           Trtf1 = c(0, diff(Ctrtf1)),
           Trtv1 = c(0, diff(Ctrtv1)),
           Incv2 = c(0, diff(Cincv2)),
           Trtv2 = c(0, diff(Ctrtv2)),
           TotTests = c(0, diff(tests)),
           TotPQCQ = Trtv1 + Trtv2,
           totcases2 = Incf1+Incv1+Incv2,
           P = Pv1+Pf1+Pv2, 
           totviv = Incv1+Incv2,
           ratios = totviv/(totviv+Incf1)) %>%  # t=1 value not t=0 for estimates for the end of that year
    pivot_longer(names_to = "variable", cols = !1) %>%
    mutate(model="Venezuela")
  
  df_sens <- as.data.frame(df_aux %>% filter(variable %in% c("totcases2")))$value
  
  #below we join together the results
  aux_mat<-cbind(as.data.frame(df_aux %>% filter(variable %in% c("time", "totcases2")))[,c(1,3)],
                 rep(import_vector[i],length(times)))
  
  yearly_sums_sa <- c(yearly_sums_sa, unname(tapply(df_sens, (seq_along(df_sens)-1) %/% 365, sum)))
  
  colnames(aux_mat)<-c("time","yearlysums","relv")
  result_import<-rbind(result_import,aux_mat)
}


sa_08<-yearly_sums_sa[1:11]
sa_085<-yearly_sums_sa[12:22]
sa_09<-yearly_sums_sa[23:33]
sa_095<-yearly_sums_sa[34:44]
sa_1<-yearly_sums_sa[45:55]
sa_105<-yearly_sums_sa[56:66]


colors <- c("Actual incidence" = "blue", "Model predicted incidence" = "red",
            "Microscopy"="blue", "RDT (same availability)"="red", 
            "RDT (100% availability)"="green",
            "0.8" = "red",
            "0.85" ="purple" ,
            "0.9" = "indigo",
            "0.95" = "light blue",
            "1" = "teal",
            "1.05" = "amber")



# sensitivity analysis of changing relv 
ggplot() + 
  geom_line(aes(x=ts, y=sa_08, color="0.8"))  + 
  geom_line(aes(x=ts, y=sa_085, color="0.85"))  + 
  geom_line(aes(x=ts, y=sa_09, color="0.9"))  + 
  geom_line(aes(x=ts, y=sa_095, color="0.95"))  + 
  geom_line(aes(x=ts, y=sa_1, color="1"))  + 
  geom_line(aes(x=ts, y=sa_105, color="1.05"))  + 
  geom_line(aes(x=ts, y=cases, color="Actual incidence")) + 
  labs(color = "relv", title="Annual incidence predicted for each value of relv",
       x="Time", y="Total incidence") + 
  scale_x_continuous(breaks= pretty_breaks())


# relv=0.9 gives best fit and so this is what was used in the analysis 



# Economic investigation into ears 2008 to 2016 ####

ts1 <- c(2008:2016)
cases <- cases[-11]
cases <- cases[-10]

# yearly cases predicted from model
df2 <- as.data.frame(df1 %>% filter(variable %in% c("totcases2")))$value
yearlysums <- unname(tapply(df2, (seq_along(df2)-1) %/% 365, sum))
yearlysums <- yearlysums[-11]
yearlysums <- yearlysums[-10]

# total yearly PQCQ treatments given and yearly costs
df3 <- as.data.frame(df1 %>% filter(variable %in% c("TotPQCQ")))$value
yearlyPQCQ <- unname(tapply(df3, (seq_along(df3)-1) %/% 365, sum))
yearlyPQCQcost <- yearlyPQCQ * (0.29+0.078) # average 14 day treatment of PQ and 3 days of CQ cost
yearlyPQCQcost <- yearlyPQCQcost[-11]
yearlyPQCQcost <- yearlyPQCQcost[-10]

# total yearly ACT treatments given and yearly costs
df4 <- as.data.frame(df1 %>% filter(variable %in% c("Trtf1")))$value
yearlyACT <- unname(tapply(df4, (seq_along(df4)-1) %/% 365, sum))
yearlyACTcost <- yearlyACT * (0.67+0.29) # 3 day treatment of AL + PQ 
yearlyACTcost <- yearlyACTcost[-11]
yearlyACTcost <- yearlyACTcost[-10]

# total yearly costs for LLINs 
yearlyNETcost <- nets * (5+2) # cost of LLIN and distribution
yearlyNETcost <- yearlyNETcost[-11]
yearlyNETcost <- yearlyNETcost[-10]

# total yearly costs of microscopy tests 
df5 <- as.data.frame(df1 %>% filter(variable %in% c("TotTests")))$value
yearlyMicroS <- unname(tapply(df5, (seq_along(df5)-1) %/% 365, sum))
yearlyMicroScost <- (yearlyMicroS * 7.79) # costs of microscopy test 
yearlyMicroScost<- yearlyMicroScost[-11]
yearlyMicroScost <- yearlyMicroScost[-10]

# total yearly costs of rdt tests 
#df5 <- as.data.frame(df1 %>% filter(variable %in% c("TotTests")))$value
#yearlyRDT <- unname(tapply(df5, (seq_along(df5)-1) %/% 365, sum))
#yearlyRDTcost <- (yearlyRDT * 2.29) # costs of rdt test 
#yearlyRDTcost <- yearlyRDTcost[-11]
#yearlyRDTcost <- yearlyRDTcost[-10]

totalyearlycosts <- yearlyPQCQcost + yearlyACTcost + yearlyNETcost + yearlyMicroScost

# Costs and for when RDT diagnosis is the routine method used ####

yearlysums2 <- c(469529.7,  488644.2,  497128.4,  511449.7,  531569.1,  562918.8,
                 624438.3,  762886.1, 1069086.5)

yearlyRDTcost <- c(328027.6,  376109.2,  409466.8,  450641.3,
                   497854.7,  552890.0,  629423.9,  766258.3,
                   1043718.1)

yearlyACTcost2 <- c(12061.25,  18618.27,  29002.33,  42345.42,  57570.54,
                    73876.02,  92463.57, 118675.28, 164589.64)

yearlyPQCQcost2 <- c(32698.12,  34998.18,  34821.77,  34493.77,  34169.66,
                     34321.48,  36026.05,  41601.30, 55239.17)

totalyearlycosts2 <- yearlyPQCQcost2 + yearlyACTcost2 + yearlyNETcost + yearlyRDTcost



# Change availability of RDT tests to 100% and compare again
# the total costs required / extra investment needed for tests and treatments

yearlysums3 <- c(414351.7, 367382.4, 330042.8, 308418.9, 299427.6, 308422.1,
                 355579.3, 491027.1, 813259.6)

yearlyRDTcost2 <- c(341743.2, 317189.1, 288276.2, 272248.5, 267298.9,
                    279681.4, 330715.4, 473280.4, 813971.0)

yearlyACTcost3 <- c(8017.832,  5508.069,  4409.120,  3903.815,  4347.682,
                    6915.218, 14765.248, 35467.648, 84158.785)

yearlyPQCQcost3 <- c(35268.75, 32721.02, 29735.81, 28008.83, 27203.47,
                     27584.89, 30369.90, 38630.43, 58298.45)


totalyearlycosts3 <- yearlyPQCQcost3 + yearlyACTcost3 + yearlyNETcost + yearlyRDTcost2

# difference in costs from having switched from microscopy in 2008-2016
# to RDT diagnosis with 100% availability

funds_diff <- totalyearlycosts3- totalyearlycosts

# Plots for years 2008 to 2016 ####

# yearlycases predicted and actual yearly cases
ggplot() +
  geom_line(aes(x=ts1, y=cases, color="Actual incidence")) + 
  geom_line(aes(x=ts1, y=yearlysums, color="Model predicted incidence")) + 
  labs(title="Actual and predicted annual \nsuspected incidence (base scenario)",
       x ="Time", y = "Total incidence", color = "Legend")+ 
  scale_x_continuous(breaks= pretty_breaks())

# yearly costs when microscopy and rdts used
ggplot() + 
  geom_line(aes(x=ts1, y=totalyearlycosts, color="Microscopy"))  + 
  geom_line(aes(x=ts1, y=totalyearlycosts2,color="RDT (same availability)")) +
  geom_line(aes(x=ts1, y=totalyearlycosts3,color="RDT (100% availability)")) + 
  labs(color = "Legend", title="Annual total cost for each diagnostic method",
       x="Time", y="Total cost") + 
  scale_x_continuous(breaks= pretty_breaks())

# difference in costs from switching microscopy to RDT diagnosis
percentdiffavg <- sum((totalyearlycosts2-totalyearlycosts)/totalyearlycosts * 100) / length((totalyearlycosts2-totalyearlycosts)/totalyearlycosts * 100)

# yearly cases for microscopy and rdts used 
ggplot() +
  geom_line(aes(x=ts1, y=yearlysums, color="Microscopy")) + 
  geom_line(aes(x=ts1, y=yearlysums2, color="RDT (same availability)")) + 
  geom_line(aes(x=ts1, y=yearlysums3, color="RDT (100% availability)")) + 
  labs(color = "Legend", title="Predicted annual incidence with each diagnostic method",
       x="Time", y="Total cases")+ 
  scale_x_continuous(breaks= pretty_breaks())

# difference in annual incidence when microscopy is switched to RDT diagnosis
percentdiffavg1 <- sum((yearlysums2-yearlysums)/yearlysums * 100) / length((yearlysums2-yearlysums)/yearlysums * 100)

# difference in costs from switching microscopy to RDT diagnosis with 100% availability
percentdiffavg2 <- sum((totalyearlycosts3-totalyearlycosts)/totalyearlycosts * 100) / length((totalyearlycosts3-totalyearlycosts)/totalyearlycosts * 100)

# difference in annual incidence when microscopy is switched to RDT diagnosis with 100% availability
percentdiffavg3 <- sum((yearlysums3-yearlysums)/yearlysums * 100) / length((yearlysums3-yearlysums)/yearlysums * 100)

