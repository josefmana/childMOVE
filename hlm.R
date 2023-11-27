# clean the environment
rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(dplyr)
library(purrr)
library(simr) # power analysis
library(lme4) # HLMs
library(sjPlot)
library(ggpubr)

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

# prepare the data set
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
         cas = (as.numeric(cas)-1)*2,
         cas2 = cas^2 )

# add time since pre-test estimates based on age
anal_dat <- 
  
  left_join(
    anal_dat,
    anal_dat %>%
      select( kod, cas, chronologicky_vek ) %>%
      pivot_wider( values_from = chronologicky_vek, names_from = cas ) %>%
      mutate( `2` = `2`-`0`, `4` = `4`-`0`, `6` = `6`-`0`, `8` = `8`-`0`, `0` = `0`-`0` ) %>%
      pivot_longer( -kod, names_to = "cas", values_to = "odstup" ) %>%
      mutate( cas = as.numeric(cas), odstup = odstup * 12 ),
    by = c("kod","cas")

  )


# FIT HLMs ----

# fit REML-based LMMs
fit1 <-
  lapply( setNames( vars[-which(vars == "chronologicky_vek")], vars[-which(vars == "chronologicky_vek")] ),
          function(i)
            lmer(
              as.formula(paste0(i, " ~ cas * intervence + chronologicky_vek + pohlavie + (1 + cas | kod)")),
              data = anal_dat,
              REML = T
            )
          )

# go for squared model next
fit2 <- 
  lapply( setNames( vars[-which(vars == "chronologicky_vek")], vars[-which(vars == "chronologicky_vek")] ),
          function(i)
            lmer(
              as.formula(paste0(i, " ~ cas * intervence + cas2 * intervence + chronologicky_vek + pohlavie + (1 | kod)")),
              data = anal_dat,
              REML = T
            )
  )


# COMPARE HLMs ----

# compute AICs and BICs (lower = better)
comp <-
  sapply( names(fit1),
          function(i)
            sapply(
              c("AIC","BIC"),
              function(j)
                c( do.call( j, list( fit1[[i]] ) ),
                   square = do.call( j, list( fit2[[i]] ) )
                )
            )
          ) %>% t() %>% `colnames<-`( c("AIC_linear","AIC_square","BIC_linear","BIC_square") )

# save it
write.table( as.data.frame( round(comp, 2 ) ) %>% rownames_to_column("outcome"),
             "tabs/porovnani_modelu.csv",
             sep = ",", row.names = F, quote = F
             )


# POSTPROCESSING ----

# extract summary tables
sum <- 
  lapply( setNames( names(fit1), names(fit1) ),
          function(i)
            tab_model( fit1[[i]], show.stat = T, show.std = T, file = paste0(dir, "/tabs/", i, ".html") )
          )

# save all of them
for (i in names(sum) ) print( sum[[i]] )

# extract diagnostics
dia <- 
  lapply( setNames( names(fit1), names(fit1) ),
          function(i)
            plot_model(fit1[[i]], type = "diag")
          )

# plot diagnostics
for (i in names(dia)) {
  for (j in c(1:4)) {
    print(dia[[i]][[j]])
    ggsave(paste0(dir, "/dia/", i, "_", j, ".jpg") )
  }
}

# plot model-based interactions
v <-
  data.frame(
    vars = c("hruba_motorika", "koordinacia_hornych_koncatin", "bilateralna_koordinacia", "balanc",
              "rychlost_ohybnost", "sila"), # variable names
    nams = c("Celková hrubá motorika", "Koordinácia horních končatín", "Bilaterálna koordinácia",
              "Rovnováha", "Rychlosť a pohyblivosť", "Sila") # variable labels
  )

# plot it
with(
  v,
  plt <<-
    lapply( setNames(vars,vars),
            function(i)
              plot_model( fit1[[i]], type = "int" ) +
              labs(x = "Čas (mesiace)", y = nams[vars == i], color = "Skupina", title = NULL) +
              scale_fill_manual(labels = c("kontrolná sk.", "intervenčná sk."), values = c("blue", "red")) +
              scale_color_manual(labels = c("kontrolná sk.", "intervenčná sk."), values = c("blue", "red")) +
              theme_classic(base_size = 24)
    )
)

# arrange it
ggarrange( plt[[1]], plt[[2]], plt[[3]], plt[[5]], plt[[4]], plt[[6]],
           nrow = 3, ncol = 2, common.legend = T, legend = "bottom" )

# save it
ggsave( paste0(dir, "/plots/hlm_plot.jpg"), height = 4*4.28, width = 2*9.49 )

# violin plots
anal_dat %>%
  
  # prepare data
  pivot_longer( cols = names(plt), names_to = "scale", values_to = "score" ) %>%
  mutate(
    scale = sapply( 1:nrow(.), function(i) with( v, nams[ vars == scale[i] ] ) ) %>%
      factor( levels = v$nams, ordered = T ),
    Skupina = ifelse( intervence == 0, "kontrolná sk.", "intervenčná sk." ) %>%
      factor( levels = c("kontrolná sk.", "intervenčná sk."), ordered = T )
  ) %>%
  
  # plotting proper
  ggviolin( x = "cas", y = "score", color = "Skupina", add = "boxplot", trim = T ) +
  facet_wrap( ~ scale, scales = "free", nrow = 3 ) +
  labs( x = "Čas (mesiace)", y = NULL ) +
  theme( legend.position = "bottom" ) +
  scale_color_manual( values = c("blue","red") )

# save it
ggsave( paste0(dir, "/plots/violins.jpg"), height = 13.3, width = 12.4 )
