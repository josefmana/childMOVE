setwd(dirname(rstudioapi::getSourceEditorContext()$path))
dir <- getwd()
date <- format(Sys.Date(), "%Y%m%d")
# packages
library(tidyverse)
library(dplyr)
library(lme4)
library(sjPlot)
# set ggplots theme
theme_set(theme_classic(base_size = 18))
# data get ready
dat <- read.table("_raw/data.txt", header = T)
dat <- dat %>% mutate(pohlavie = as.factor(pohlavie), intervence = as.factor(intervence))
vars_orig <- names(dat[4:48])
vars <- unique(substr(vars_orig, 1, nchar(vars_orig)-2))
vars <- sort(vars)
preds <- names(dat[1:3])
dum_names <- split(sort(vars_orig), ceiling(seq_along(vars_orig)/5))
dum_dat <- lapply(1:length(dum_names), function(i)
  dat[ , which(colnames(dat) %in% c(preds, dum_names[[i]]))] %>%
    pivot_longer(all_of(dum_names[[i]]),
                 names_to = "cas", names_prefix = paste0(vars[[i]], "_"),
                 values_to = vars[[i]])
)
names(dum_dat) <- vars
anal_dat <- dum_dat[[1]]
for (i in vars) {
  anal_dat[[i]] <- dum_dat[[i]][[i]]
}
anal_dat <- anal_dat %>% mutate(kod = as.factor(kod),
                                pohlavie = as.factor(pohlavie),
                                intervence = as.factor(intervence),
                                cas_mono = as.ordered(cas),
                                cas = (as.numeric(cas)-1)*2)
# fit hlms
fit <- list()
for (i in vars[-which(vars == "chronologicky_vek")]) {
  fit[[i]] <-
    lmer(as.formula(paste0(i, " ~ cas * intervence + chronologicky_vek + pohlavie + (1 + cas | kod)")),
         data = anal_dat)
}
# summaries
sum <- list()
for (i in names(fit)) {
  sum[[i]] <- tab_model(fit[[i]])
}
for (i in names(sum)) {
  XML::saveXML(sum[[i]], paste0(dir, "/tabs/", i, ".html"))
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
nams <- c("Celkov? hrub? motorika", "Koordin?cia horn?ch kon?at?n", "Bilater?lna koordin?cia",
          "Rovnov?ha", "R?chlos? a pohyblivos?", "Sila")
plt <- list()
for (i in 1:length(vars)) {
  plt[[vars[i]]] <- plot_model(fit[[vars[i]]], type = "int")
  plt[[vars[i]]] <- plt[[vars[i]]] +
    labs(x = "?as (mesiace)", y = nams[[i]], color = "Skupina", title = NULL) +
    scale_fill_manual(labels = c("sk. bez intervence", "interve?n? sk."), values = c("blue", "red")) +
    scale_color_manual(labels = c("sk. bez intervence", "interve?n? sk."), values = c("blue", "red")) +
    theme_classic(base_size = 24)
}
ggpubr::ggarrange(plt[[1]], plt[[2]], plt[[3]],
                  plt[[5]], plt[[4]], plt[[6]],
                  nrow = 3, ncol = 2,
                  common.legend = T, legend = "bottom")
ggsave(paste0(dir, "/plots/the_plot.jpg"), device = "jpeg", height = 4*4.28, width = 2*9.49)

# monotonic models to control
library(brms)
options(mc.cores = parallel::detectCores())
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
  ggsave(paste0(dir, "/mono/cond/", i, ".jpg"), device = "jpeg")
}
