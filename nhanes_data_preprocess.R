library("sandwich")
library(R2jags)
library(parallel)
library(dplyr)
set.seed(4)

## systolic blood pressure (unit: mm Hg)
bpx <- read.csv("BPX_J.csv")[, c("SEQN", "BPXSY1", "BPXSY2")] %>% 
  mutate(SEQN = as.character(SEQN),
         BPXSY1 = as.numeric(BPXSY1),
         BPXSY2 = as.numeric(BPXSY2)) %>%
  mutate(BPXSY = (BPXSY1 + BPXSY2) / 2) %>%
  filter(., complete.cases(.))

## transform variable type
demo <- read.csv("DEMO_J.csv")[, c("SEQN", "RIDAGEYR", "RIAGENDR")] %>%
  transmute(SEQN = as.character(SEQN),
            RIDAGEYR = as.numeric(RIDAGEYR),
            MALE = as.numeric(RIAGENDR == 1)) %>%
  filter(., complete.cases(.)) 

  



## merge data by subject id
dat <- inner_join(x = bpx, y = demo, by = "SEQN")
dat <- dat[sample.int(nrow(dat), size = 200),]

saveRDS(dat, "nhanes_subset.rds")
