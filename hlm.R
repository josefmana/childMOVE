# clean the environment
rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(dplyr)
library(purrr)
library(lme4)
library(brms)
library(sjPlot)
library(XML)

# extract current directory
dir <- here()

# set ggplots theme
theme_set(theme_classic(base_size = 18))


# DATA PREPARATION ----

# read the data
dat <-
  read.table( here("_raw","data.txt"), header = T) %>%
  mutate(pohlavie = as.factor(pohlavie), intervence = as.factor(intervence))

# extract variables of interest
vars_orig <- names(dat[4:48])
vars <- unique(substr(vars_orig, 1, nchar(vars_orig)-2)) %>% sort()

# extract predictors of interest
preds <- names(dat[1:3])

# extract temporary names
tmp_names <- split(sort(vars_orig), ceiling(seq_along(vars_orig)/5))

# preapare the data set
anal_dat <-
  lapply(1:length(tmp_names),
         function(i)
           dat[ , which(colnames(dat) %in% c(preds, tmp_names[[i]]))] %>%
           pivot_longer(all_of(tmp_names[[i]]),
                        names_to = "cas",
                        names_prefix = paste0(vars[[i]], "_"),
                        values_to = vars[[i]])
         ) %>%
  # pull them together
  reduce( left_join, by = c("kod","pohlavie","intervence","cas") ) %>%
  # prepare variables
  mutate(kod = as.factor(kod),
         pohlavie = as.factor(pohlavie),
         intervence = as.factor(intervence),
         cas_mono = as.ordered(cas),
         cas = (as.numeric(cas)-1)*2)


# FIT HLMs ----

# fit REML-based LMMs
fit <-
  lapply( setNames( vars[-which(vars == "chronologicky_vek")], vars[-which(vars == "chronologicky_vek")] ),
          function(i)
            lmer(
              formula = as.formula(paste0(i, " ~ cas * intervence + chronologicky_vek + pohlavie + (1 + cas | kod)")),
              data = anal_dat,
              REML = T
            )
          )

# prepare folders for Bayesian models 
sapply( paste0( "_bayes", c("","/mono","/mono/mods") ), function(i) if( !dir.exists(i) ) dir.create(i) )

# use all cores
options(mc.cores = parallel::detectCores())

# fit monotonic predictions models
fit_mon <-
  lapply( setNames( names(fit), names(fit) ),
          function(i)
            brm(
              as.formula(paste0(i, " ~ mo(cas) * intervence + chronologicky_vek + pohlavie + (1 + cas | kod)")),
              data = anal_dat, init = 0, control = list(adapt_delta = .95), seed = 87542,
              file = paste0(dir, "/_bayes/mono/mods/", i, ".rds")
            )
          )



# summaries
sum <- list()
for (i in names(fit)) {
  sum[[i]] <- tab_model(fit[[i]])
}
for (i in names(sum)) {
  saveXML(sum[[i]], paste0(dir, "/tabs/", i, ".html"))
}
# diagnostics
dia <- list()
for (i in names(fit)) {
  dia[[i]] <- plot_model(fit[[i]], type = "diag")
}
for (i in names(dia)) {
  for (j in c(1:4)) {
    print(dia[[i]][[j]])
    ggsave(paste0(dir, "/dia/", i, "_", j, ".jpg"), device = "jpeg")
  }
}
# plots
for (i in names(fit)) {
  plot_model(fit[[i]], type = "int")
  ggsave(paste0(dir, "/plots/", i, ".jpg"), device = "jpeg")
}

# 20210414 big plot
vars <- c("hruba_motorika", "koordinacia_hornych_koncatin", "bilateralna_koordinacia", "balanc",
          "rychlost_ohybnost", "sila")
nams <- c("Celková hrubá motorika", "Koordinácia horních končatín", "Bilaterálna koordinácia",
          "Rovnováha", "Rychlosť a pohyblivosť", "Sila")
plt <- list()
for (i in 1:length(vars)) {
  plt[[vars[i]]] <- plot_model(fit[[vars[i]]], type = "int")
  plt[[vars[i]]] <- plt[[vars[i]]] +
    labs(x = "Čas (mesiace)", y = nams[[i]], color = "Skupina", title = NULL) +
    scale_fill_manual(labels = c("kontrolná sk.", "intervenčná sk."), values = c("blue", "red")) +
    scale_color_manual(labels = c("kontrolná sk.", "intervenčná sk."), values = c("blue", "red")) +
    theme_classic(base_size = 24)
}
ggpubr::ggarrange(plt[[1]], plt[[2]], plt[[3]],
                  plt[[5]], plt[[4]], plt[[6]],
                  nrow = 3, ncol = 2,
                  common.legend = T, legend = "bottom")
ggsave(paste0(dir, "/plots/the_plot.jpg"), device = "jpeg", height = 4*4.28, width = 2*9.49)

# monotonic models to control


fit_mon <- list()
for (i in vars[-which(vars == "chronologicky_vek")]) {
  fit_mon[[i]] <- brm(
    as.formula(paste0(i, " ~ mo(cas) * intervence + chronologicky_vek + pohlavie + (1 + cas | kod)")),
    data = anal_dat, inits = 0, control = list(adapt_delta = .95), seed = 87542,
    file = paste0(dir, "/_mono/", i, ".rds")
    )
}
for (i in names(fit_mon)) {
  print(conditional_effects(fit_mon[[i]], effects = "cas:intervence"))
  ggsave(paste0(dir, "/_mono/cond/", i, ".jpg"), device = "jpeg")
}
