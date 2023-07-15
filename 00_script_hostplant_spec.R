# Phylogenetic comparative analyses of body size in moths
# Evolutionary relationships between phenology, body size and host-plant preferences
# Stenio Foerster
# stenio.foerster@ut.ee
# 09.07.2023

# Libraries ------------------------------------------------------------------------------

library(tidyverse)
library(ape)
library(phytools)
library(OUwie)
library(nlme)
library(car)
library(marginaleffects)
library(emmeans)
library(ggsci)

# Data -----------------------------------------------------------------------------------

# Windows
setwd("/Users/Stenio Foerster/Nextcloud/Geometridae/geom_spec/")
# MAC
setwd("~/Nextcloud/Geometridae/geom_spec/")

# data
d <- read_csv("01_geom_traits_complete.csv")
# tree
t <- force.ultrametric(read.nexus("02_geom_tree_complete.nex"))
# all good
all.equal(d$species_tree, t$tip.label)

# get species strictly associated with woody plants or herbs
d %>% filter(foodplant %in% c("W", "H")) -> dt
dt$foodplant <- ifelse(dt$foodplant == "H", "herbs", "woody")
# recode specialism into specialist or generalist
dt$specialism <- ifelse(dt$specialism == "M", "specialist", "generalist")
table(dt$foodplant, dt$specialism)
# initial plot (exploratory)
#ggplot(data = dt, mapping = aes(x = foodplant, y = log_dbm_mean, color = specialism)) + geom_boxplot()

# prune the phylogeny
tt <- drop.tip(phy = t, tip = setdiff(t$tip.label, dt$species_tree))
dt <- dt[match(tt$tip.label, dt$species_tree), ]
# all good
all.equal(dt$species_tree, tt$tip.label)


# PGLS -----------------------------------------------------------------------------------

# variance-covariance matrices
mBM <- corBrownian(value = 1, phy = tt, form = ~species_tree)
mPA <- corPagel(value = 1, phy = tt, form = ~species_tree, fixed = F)

# models
m1 <- gls(log_dbm_mean ~ foodplant * specialism, data = dt, correlation = mBM, method = "REML")
m2 <- gls(log_dbm_mean ~ foodplant * specialism, data = dt, correlation = mPA, method = "REML")

# comparison
AICcmodavg::AICc(m1)
AICcmodavg::AICc(m2)

# model 1 (BM)
summary(m1)
plot_predictions(m1, condition = c("foodplant", "specialism"), type = "response") + labs(title = "BM, AICc = 389")
emmeans(object = m1, specs = ~ foodplant|specialism, mode = "df.error") %>% pairs()
emmeans(object = m1, specs = ~ specialism|foodplant, mode = "df.error") %>% pairs()

# model 2 (Pagel)
summary(m2)
plot_predictions(m2, condition = c("foodplant", "specialism"), type = "response") + labs(title = "Pagel, lambda = 0.87, AICc = 353")
emmeans(object = m2, specs = ~ foodplant|specialism, mode = "df.error") %>% pairs()
emmeans(object = m2, specs = ~ specialism|foodplant, mode = "df.error") %>% pairs()

# model parameters (useful for plotting)
predictions(model = m1, type = "response", by = c("foodplant", "specialism"))
predictions(model = m2, type = "response", by = c("foodplant", "specialism"))

# create a common ggplot theme
th <- theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(colour = "black", size = 12, face = "bold")
        ,axis.text = element_text(colour = "black", size = 11)
        ,legend.title = element_text(face = "bold", colour = "black", size = 12)
        ,legend.text = element_text(colour = "black", size = 11))

# plot with standard errors
pp <- predictions(model = m2, type = "response", by = c("foodplant", "specialism"))
pp %>% mutate(ymin = estimate-std.error, ymax = estimate+std.error) %>% as_tibble() -> pp
ggplot(data = pp, mapping = aes(x = foodplant, y = estimate, ymin = ymin, ymax = ymax, color = specialism)) + 
  geom_pointrange(position = position_dodge2(width = 0.3), size = 0.4, linewidth = 0.3) +
  th +
  scale_y_continuous(breaks = seq(2.08, 3, 0.1)) +
  scale_color_jama()

# to save the graph as an R object (rda), first store the ggplot object into an object
# let's say "p1". Then: 
# save(p1, "plot_1.rda")
# load("plot_1.rda")

# diet specialization by host-plant growth form - compare with Davis et al. (2012)
# woody-plant feeders
m3 <- gls(log_dbm_mean ~ specialism, data = dt[dt$foodplant == "woody", ], correlation = mBM, method = "REML")
m4 <- gls(log_dbm_mean ~ specialism, data = dt[dt$foodplant == "woody", ], correlation = mPA, method = "REML")
# comparison
AICcmodavg::AICc(m3)
AICcmodavg::AICc(m4)
# summary
summary(m4)
plot_predictions(m4, condition = c("specialism"), type = "response") + labs(title = "woody-plant feeders")

# plot with standard errors
d1 <- predictions(model = m4, type = "response", by = c("specialism"))
d1 %>% mutate(ymin = estimate-std.error, ymax = estimate+std.error) %>% as_tibble() -> d1
ggplot(data = d1, mapping = aes(x = specialism, y = estimate, ymin = ymin, ymax = ymax)) + 
  geom_pointrange(position = position_dodge2(width = 0.3), size = 0.4, linewidth = 0.3) +
  th +
  scale_y_continuous(breaks = seq(2.08, 3, 0.1)) +
  labs(title = "woody-plant feeders")

# herb feeders
m5 <- gls(log_dbm_mean ~ specialism, data = dt[dt$foodplant == "herbs", ], correlation = mBM, method = "REML")
m6 <- gls(log_dbm_mean ~ specialism, data = dt[dt$foodplant == "herbs", ], correlation = mPA, method = "REML")
# comparison
AICcmodavg::AICc(m5)
AICcmodavg::AICc(m6)
# summary
summary(m6)
plot_predictions(m6, condition = c("specialism"), type = "response") + labs(title = "herb feeders")

# plot with standard errors
d2 <- predictions(model = m6, type = "response", by = c("specialism"))
d2 %>% mutate(ymin = estimate-std.error, ymax = estimate+std.error) %>% as_tibble() -> d2
ggplot(data = d2, mapping = aes(x = specialism, y = estimate, ymin = ymin, ymax = ymax)) + 
  geom_pointrange(position = position_dodge2(width = 0.3), size = 0.4, linewidth = 0.3) +
  th +
  labs(title = "herb feeders")


# Functions -----------------------------------------------------------------------------------

# Function to run all the hOUwie models
run_hOUwie_all <- function(phy, data, mserr, rate.cat.cid, discr.model, nSim) {
  
  m <- c("BM1", "BMV", "OU1", "OUA", "OUV", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA")
  #m <- c("BMV", "OUM") # test, should result in 4 models
  res <- list()
  mCD <- list()
  mCID <- list()
  
  for (i in 1:length(m)) {
    
    # CD model
    mCD[[i]] <- hOUwie(phy = phy, data = data, rate.cat = 1, discrete_model = discr.model, continuous_model = m[i], null.model = F, nSim = nSim, mserr = mserr, adaptive_sampling = F)
    # CID+ model
    mCID[[i]] <- hOUwie(phy = phy, data = data, rate.cat = rate.cat.cid, discrete_model = discr.model, continuous_model = m[i], null.model = T, nSim = nSim, mserr = mserr, adaptive_sampling = T)
  }
  
  res <- c(mCD, mCID)
  names(res) <- c(m, paste(m, "CID", sep = "_"))
  return(res)
}


# hOUwie models --------------------------------------------------------------------------

## General ----

# Woody-plant vs. herb feeders (without split by specialization)
dat_hos <- as.data.frame(subset(dt, select = c(species_tree, foodplant, log_dbm_mean, log_dbm_se)))
all.equal(dat_hos$species_tree, tt$tip.label)
# run all models
#mod_hos <- run_hOUwie_all(phy = tt, data = dat_hos, mserr = "known", rate.cat.cid = 2, discr.model = "ARD", nSim = 100)
load("output/mod_hos.rda")

# average model parameters
mod_hos_avg <- getModelAvgParams(model.list = mod_hos, type = "AICc")
# mean values
mod_hos_avg %>% group_by(tip_state) %>% summarise_all(mean)

## Specialists ----
dat_hos_spe <- as.data.frame(subset(dt[dt$specialism == "specialist", ], select = c(species_tree, foodplant, log_dbm_mean, log_dbm_se)))
phy_spe <- drop.tip(phy = tt, tip = setdiff(x = tt$tip.label, dat_hos_spe$species_tree))
all.equal(phy_spe$tip.label, dat_hos_spe$species_tree)
# run all models
#mod_hos_spe <- run_hOUwie_all(phy = phy_spe, data = dat_hos_spe, mserr = "known", rate.cat.cid = 2, discr.model = "ARD", nSim = 100)
load("output/mod_hos_spe.rda")

## Generalists ----
dat_hos_gen <- as.data.frame(subset(dt[dt$specialism == "generalist", ], select = c(species_tree, foodplant, log_dbm_mean, log_dbm_se)))
phy_gen <- drop.tip(phy = tt, tip = setdiff(x = tt$tip.label, dat_hos_gen$species_tree))
all.equal(phy_gen$tip.label, dat_hos_gen$species_tree)
# run all models
#mod_hos_gen <- run_hOUwie_all(phy = phy_gen, data = dat_hos_gen, mserr = "known", rate.cat.cid = 2, discr.model = "ARD", nSim = 100)
load("output/mod_hos_gen.rda")

## Four-states model ----
four_dat <- as.data.frame(subset(dt, select = c(species_tree, log_dbm_mean, log_dbm_se)))
four_dat <- add_column(.data = four_dat, regime = paste(dt$foodplant, dt$specialism, sep = "_"), .after = "species_tree")
all.equal(four_dat$species_tree, tt$tip.label)
# run all models
#mod_four_states <- run_hOUwie_all(phy = tt, data = four_dat, mserr = "known", rate.cat.cid = 2, discr.model = "ARD", nSim = 100)
#save(mod_four_states, file = "output/mod_four_states.rda")
load("output/mod_four_states.rda")

# summary

# four-states model (without A models)
mod_four_states[!grepl("A", names(mod_four_states), ignore.case = F)] %>% 
  getModelAvgParams(type = "AICc") %>% 
  group_by(tip_state) %>% 
  summarise_all(mean)

# generalists model (without A models)
mod_hos_gen[!grepl("A", names(mod_hos_gen), ignore.case = F)] %>% 
  getModelAvgParams(type = "AICc") %>% 
  group_by(tip_state) %>% 
  summarise_all(mean)

# specialists model (without A models)
mod_hos_spe[!grepl("A", names(mod_hos_spe), ignore.case = F)] %>% 
  getModelAvgParams(type = "AICc") %>% 
  group_by(tip_state) %>% 
  summarise_all(mean)


# hOUwie models Approach II --------------------------------------------------------------

# generalists vs. specialists within each host plant

# data
dat_spe <- as.data.frame(subset(dt, select = c(species_tree, specialism, log_dbm_mean, log_dbm_se)))
dat_woo <- as.data.frame(subset(dt[dt$foodplant == "woody", ], select = c(species_tree, specialism, log_dbm_mean, log_dbm_se)))
dat_her <- as.data.frame(subset(dt[dt$foodplant == "herbs", ], select = c(species_tree, specialism, log_dbm_mean, log_dbm_se)))

# trees
all.equal(tt$tip.label, dat_spe$species_tree)
phy_woo <- drop.tip(phy = tt, tip = setdiff(tt$tip.label, y = dat_woo$species_tree))
phy_her <- drop.tip(phy = tt, tip = setdiff(tt$tip.label, y = dat_her$species_tree))

dat_woo <- dat_woo[match(phy_woo$tip.label, dat_woo$species_tree), ]
dat_her <- dat_her[match(phy_her$tip.label, dat_her$species_tree), ]
all.equal(phy_woo$tip.label, dat_woo$species_tree)
all.equal(phy_her$tip.label, dat_her$species_tree)

# models
#mod_spe <- run_hOUwie_all(phy = tt, data = dat_spe, mserr = "known", rate.cat.cid = 2, discr.model = "ARD", nSim = 100)
#mod_woo <- run_hOUwie_all(phy = phy_woo, data = dat_woo, mserr = "known", rate.cat.cid = 2, discr.model = "ARD", nSim = 100)
#mod_her <- run_hOUwie_all(phy = phy_her, data = dat_her, mserr = "known", rate.cat.cid = 2, discr.model = "ARD", nSim = 100)

#save(mod_spe, file = "output/mod_spe.rda")
#save(mod_woo, file = "output/mod_woo.rda")
#save(mod_her, file = "output/mod_her.rda")

load("output/mod_spe.rda")
load("output/mod_woo.rda")
load("output/mod_her.rda")

mod_spe %>% getModelAvgParams() %>% group_by(tip_state) %>% summarise_all(mean)
mod_woo %>% getModelAvgParams() %>% group_by(tip_state) %>% summarise_all(mean)
mod_her %>% getModelAvgParams() %>% group_by(tip_state) %>% summarise_all(mean)


## without measurement error ----

# data
dat_spe_noME <- dat_spe[, -4]
dat_woo_noME <- dat_woo[, -4]
dat_her_noME <- dat_her[, -4]

# models
mod_spe_noME <- run_hOUwie_all(phy = tt, data = dat_spe, mserr = "none", rate.cat.cid = 2, discr.model = "ARD", nSim = 100)
mod_woo_noME <- run_hOUwie_all(phy = phy_woo, data = dat_woo, mserr = "none", rate.cat.cid = 2, discr.model = "ARD", nSim = 100)
mod_her_noME <- run_hOUwie_all(phy = phy_her, data = dat_her, mserr = "none", rate.cat.cid = 2, discr.model = "ARD", nSim = 100)

mod_spe %>% getModelAvgParams() %>% group_by(tip_state) %>% summarise_all(mean)
mod_woo %>% getModelAvgParams() %>% group_by(tip_state) %>% summarise_all(mean)
mod_her %>% getModelAvgParams() %>% group_by(tip_state) %>% summarise_all(mean)






