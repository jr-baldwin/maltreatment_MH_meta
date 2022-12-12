Childhood maltreatment and mental health problems: a meta-analysis of
quasi-experimental studies
================
Jessie Baldwin

-   [Load packages](#load-packages)
-   [Load and examine data](#load-and-examine-data)
-   [Convert effect sizes to Cohen’s
    d](#convert-effect-sizes-to-cohens-d)
    -   [Convert log odds/odds ratio to Cohen’s
        d](#convert-log-oddsodds-ratio-to-cohens-d)
    -   [Convert unstandardised betas to Cohen’s
        d](#convert-unstandardised-betas-to-cohens-d)
    -   [Convert relative risks to Cohen’s
        d](#convert-relative-risks-to-cohens-d)
    -   [Convert hazard ratios to Cohen’s
        d](#convert-hazard-ratios-to-cohens-d)
    -   [Convert correlations/standardised betas to Cohen’s
        d](#convert-correlationsstandardised-betas-to-cohens-d)
-   [Derive additional variables needed for
    meta-analysis](#derive-additional-variables-needed-for-meta-analysis)
-   [Results](#results)
    -   [Search results - descriptives](#search-results---descriptives)
    -   [Quasi-experimental evidence on the association between child
        maltreatment and mental
        health](#quasi-experimental-evidence-on-the-association-between-child-maltreatment-and-mental-health)
        -   [Quasi-experimental adjusted meta-analysis - multi-level
            random-effects meta-analysis
            model](#quasi-experimental-adjusted-meta-analysis---multi-level-random-effects-meta-analysis-model)
        -   [Unadjusted meta-analysis - random-effects multi-level
            meta-analysis
            model](#unadjusted-meta-analysis---random-effects-multi-level-meta-analysis-model)
        -   [Aggregate effect sizes across studies to derive a
            study-average effect size to include in the forest
            plot](#aggregate-effect-sizes-across-studies-to-derive-a-study-average-effect-size-to-include-in-the-forest-plot)
        -   [Forest plot for unadjusted effect
            sizes](#forest-plot-for-unadjusted-effect-sizes)
    -   [Sensitivity analyses](#sensitivity-analyses)
        -   [Publication bias analyses](#publication-bias-analyses)
        -   [Undue influence of individual cohorts/studies/effect
            sizes](#undue-influence-of-individual-cohortsstudieseffect-sizes)
        -   [Definition of maltreatment sensitivity
            analysis](#definition-of-maltreatment-sensitivity-analysis)
    -   [Moderators of the association between child maltreatment and
        mental
        health](#moderators-of-the-association-between-child-maltreatment-and-mental-health)
        -   [Type of quasi-experimental
            method](#type-of-quasi-experimental-method)
        -   [Mental health outcome](#mental-health-outcome)
        -   [Type of child maltreatment](#type-of-child-maltreatment)
        -   [Prospective vs retrospective measures of child
            maltreatment](#prospective-vs-retrospective-measures-of-child-maltreatment)
        -   [Shared rater](#shared-rater)
        -   [Longitudinal vs cross-sectional
            design](#longitudinal-vs-cross-sectional-design)
        -   [Sample characteristics](#sample-characteristics)
        -   [Study quality](#study-quality)
    -   [Supplementary tables](#supplementary-tables)
        -   [Table with descriptives on studies (eTable
            5)](#table-with-descriptives-on-studies-etable-5)
        -   [Table with study quality information (eTable
            6)](#table-with-study-quality-information-etable-6)
        -   [Table with effect sizes (eTable
            7)](#table-with-effect-sizes-etable-7)

# Load packages

``` r
library(readxl) # read in Excel data 
library(dplyr) # for data manipulation 
library(psych) # for summarising data
library(metafor) # for meta-analysis
library(effectsize) # for converting between effect sizes
library(compute.es) # for converting between effect sizes
library(MAd) # for aggregating effect sizes
library(meta) # for meta-analysis object for p-curve analysis
library(dmetar) # for p-curve analysis
library(multcomp) # for pairwise comparisons
library(Hmisc) # for capitalising text
library(mgsub) # for editing text
library(stringr) # for editing text
library(kableExtra) # for making tables
```

# Load and examine data

``` r
## Load data
setwd(Data)
data <- data.frame(read_excel("maltreatment_mh_QEdata_20220111.xlsx"))

## Examine data 
str(data)
table(data$es_type_adj) # Examine range of effect sizes for QE-adjusted estimates
table(data$es_type_unadj) # Examine range of effect sizes for unadjusted estimates
```

# Convert effect sizes to Cohen’s d

In the below sections we convert all effect sizes to Cohen’s d, in the
following order:

-   Log odds and odds ratio (adjusted: n=26 and 43, respectively;
    unadjusted: n=34 \[odds ratio\])
-   Unstandardised beta (adjusted: n=256 unadjusted: n=15)
-   Relative risk (adjusted: n=18; unadjusted: n=18)
-   Hazard ratio (adjusted: n=1; unadjusted: n=1)
-   Correlation coefficient (adjusted: n=8; unadjusted n=4) and
    standardised beta (adjusted: n=14; unadjusted: n=14)

``` r
## Derive new Cohen's d effect size and SE variables for conversions 
# Adjusted ES
data$d_adj <- NA
data$d_adj[data$es_type_adj=="cohen_d"] <- data$es_adj[data$es_type_adj=="cohen_d"]
describe(data$d_adj)
# Adjusted SE
data$d_se_adj <- NA
data$d_se_adj[data$es_type_adj=="cohen_d"] <- data$se_adj[data$es_type_adj=="cohen_d"]
describe(data$d_se_adj)

# Unadjusted ES
data$d_unadj <- NA
data$d_unadj[data$es_type_unadj=="cohen_d"] <- data$es_unadj[data$es_type_unadj=="cohen_d"]
describe(data$d_unadj)
# Unadjusted SE
data$d_se_unadj <- NA
data$d_se_unadj[data$es_type_unadj=="cohen_d"] <- data$se_unadj[data$es_type_unadj=="cohen_d"]
describe(data$d_se_unadj)
```

## Convert log odds/odds ratio to Cohen’s d

To convert log odds coefficients and odds ratios to Cohen’s d, we will:

-   convert odds ratios to log odds coefficients
-   convert log odds coefficients to Cohen’s d using the formula: d =
    log_odds \* (sqrt(3)) / 3.14159<sup>1</sup>
-   convert log odds variance to Cohen’s d variance using the formula:
    variance_d = variance_log_odds \* (3/(3.14159^2)) <sup>1</sup>
-   compare effect sizes and p-values before and after converting to
    Cohen’s d

<sup>1</sup> conversion formulae are taken from [Borenstein et
al. (2009)](https://www.meta-analysis.com/downloads/Meta-analysis%20Converting%20among%20effect%20sizes.pdf)

``` r
# View studies and number of effect sizes with log odds or odds ratios 
data %>% dplyr::filter (es_type_adj == "log_odds" | data$es_type_adj == "odds_ratio") %>% 
  group_by(author) %>%  dplyr::summarise(n_effectSizes = n())

## Convert odds ratios to log odds coefficients (SEs for odds ratios are already in log format)
# Adjusted ES
data$es_adj[data$es_type_adj=="odds_ratio"] <- log(data$es_adj[data$es_type_adj=="odds_ratio"])
# Unadjusted ES
data$es_unadj[data$es_type_unadj=="odds_ratio"] <- log(data$es_unadj[data$es_type_unadj=="odds_ratio"])

# Specify that odds ratios have been converted to log odds
data$es_type_adj[data$es_type_adj=="odds_ratio"] <- "log_odds"
data$es_type_unadj[data$es_type_unadj=="odds_ratio"] <- "log_odds"
  
## Convert log odds ratios to Cohen's d 
# Cohen's d = log_odds * (sqrt(3)) / 3.14159; from: https://www.meta-analysis.com/downloads/Meta-analysis%20Converting%20among%20effect%20sizes.pdf
# Adjusted ES
data$d_adj[data$es_type_adj == "log_odds"] <- data$es_adj[data$es_type_adj == "log_odds"]*(sqrt(3)/3.14159)
# Unadjusted ES
data$d_unadj[data$es_type_unadj == "log_odds"] <- data$es_unadj[data$es_type_unadj == "log_odds"]*(sqrt(3)/3.14159)

# Define function for conversion from log odds SE to Cohen's d SE
SE_logodds_D <- function(se_log_odds) {
  variance_log_odds <- se_log_odds^2 # calculate the variance of the log odds SE
  variance_D <- variance_log_odds * (3/(3.14159^2)) # convert the log odds variance to Cohen's d variance
  se_D <- sqrt(variance_D) # convert variance of Cohen's d to SE of Cohen's d
  return(se_D)
}

## Convert log odds SEs to Cohen's d SEs
# Adjusted SE
data$d_se_adj[data$es_type_adj == "log_odds"] <- SE_logodds_D(data$se_adj[data$es_type_adj == "log_odds"])
# Unadjusted SE
data$d_se_unadj[data$es_type_unadj == "log_odds"] <- SE_logodds_D(data$se_unadj[data$es_type_unadj == "log_odds"])

## Check effect sizes before and after converting to Cohen's d
data %>% dplyr::filter (es_type_adj == "log_odds" & es_type_unadj == "log_odds") %>% 
  dplyr::select(author, mh_outcome, es_unadj, d_unadj, es_adj, d_adj)

## Check p-values before and after converting to Cohen's d
# Derive function to calculate p-value (see box in https://www.bmj.com/content/343/bmj.d2304)
p_function <- function (es, se) {
  z <- es/se
  p <- exp(-0.717*z - 0.416*z^2)
  return(p)
} 

# P-values for adjusted effect size
p_function(data$es_adj[data$es_type_adj=="log_odds"], data$se_adj[data$es_type_adj=="log_odds"]) # log_odds
p_function(data$d_adj[data$es_type_adj=="log_odds"], data$d_se_adj[data$es_type_adj=="log_odds"]) # cohen's d
# P-values for unadjusted effect sizes  
p_function(data$es_unadj[data$es_type_unadj=="log_odds"], data$se_unadj[data$es_type_unadj=="log_odds"]) # log odds
p_function(data$d_unadj[data$es_type_unadj=="log_odds"], data$d_se_unadj[data$es_type_unadj=="log_odds"]) # cohen's d

## Specify that effect sizes for log odds have now been converted to Cohen's d 
data$es_type_adj[data$es_type_adj == "log_odds"] <- "cohen_d"
data$es_type_unadj[data$es_type_unadj == "log_odds"] <- "cohen_d"

## Check effect sizes left to convert
table(data$es_type_adj)
table(data$es_type_unadj)
```

## Convert unstandardised betas to Cohen’s d

To convert unstandardised betas to Cohen’s d, we will:

-   derive variables indexing the numbers of exposed and unexposed
    participants in each study<sup>1</sup>
-   convert unstandardised betas to Cohen’s d using an adapted formula
    from the
    [“esc_B”](https://www.rdocumentation.org/packages/esc/versions/0.5.1/topics/esc_B)
    function of the “esc” R package<sup>2</sup>
-   compare effect sizes and p-values before and after converting to
    Cohen’s d

<sup>1</sup> note: if the maltreatment variable was continuous, assume
equal sample size in N1 and N2

<sup>2</sup> note: this is adapted to first calculate the standard
deviation of the outcome from the standard error of the unstandardised
effect size reported in the studies

``` r
# View studies and number of effect sizes with unstandardised betas
data %>% dplyr::filter (es_type_adj == "unstandardised_beta" | es_type_unadj == "unstandardised_beta") %>% 
  group_by(author) %>%  summarise(n_effectSizes = n())

# View studies with unstandardised betas and continuous maltreatment variables (Alemany, Isumi, Kullberg, Lecei, Schwartz)
data %>% filter (es_type_adj == "unstandardised_beta" & maltreatment_type_var=="continuous" |
                   es_type_unadj=="unstandardised_beta" & maltreatment_type_var=="continuous") %>%
  group_by(author)  %>% count(maltreatment_type_var) %>% slice(which.max(n))  %>% arrange(desc(maltreatment_type_var))

# View studies with unstandardised betas and categorical/binary maltreatment variables (Kugler, Lynch, Ma, Voith, Thornberry)
data %>% filter (es_type_adj == "unstandardised_beta" & maltreatment_type_var!="continuous" |
                   es_type_unadj=="unstandardised_beta" & maltreatment_type_var!="continuous") %>%
  group_by(author)  %>% count(maltreatment_type_var) %>% slice(which.max(n))  %>% arrange(desc(maltreatment_type_var))

################### Derive variables reflecting the N in the exposed and unexposed groups ###################
# If the variable was continuous, assume equal group sizes (50% in each group)
# If the variable was categorical/binary, obtain the exposed/unexposed Ns from the paper
# If the design is fixed-effects and the same individuals were included at 2 time points (e.g., Voith), N_exposed and N_unexposed should both equal the total sample size
# Note: in Ma (2018), fixed-effects design was used and observations not reported. But results were extracted for coefficient "spanking more than 20 times in past year" reported by 9% of sample (222), so use 222 for N_exposed and N_unexposed

#### Derive N_exposed (maltreated) for adjusted analyses
data$n_exp_adj <- NA
# Code as half of N if maltreatment was continuous 
data$n_exp_adj[data$es_type_adj=="unstandardised_beta" & data$maltreatment_type_var=="continuous"] <- 
  data$n_adj[data$es_type_adj=="unstandardised_beta" & data$maltreatment_type_var=="continuous"]/2 
# Code as total N in Isumi, fixed-effects (all participants were examined at each time point)
data$n_exp_adj[data$es_type_adj=="unstandardised_beta" & data$author=="Isumi"] <- 2920 # fixed-effects design
# Code specific Ns for studies with categorical or binary coding of maltreatment
data$n_exp_adj[data$es_type_adj=="unstandardised_beta" & data$author=="Lynch"] <- 286 # harsh physical punishment
data$n_exp_adj[data$es_type_adj=="unstandardised_beta" & data$author=="Kugler"] <- 179 # investigated + substantiated
data$n_exp_adj[data$es_type_adj=="unstandardised_beta" & data$author=="Ma"] <- 222 # 9% of total N
data$n_exp_adj[data$es_type_adj=="unstandardised_beta" & data$author=="Voith"] <- 1978/2 # half observations included in model
data$n_exp_adj[data$es_type_adj=="unstandardised_beta" & data$author=="Thornberry" 
                    & data$maltreatment_period=="0-11"] <- 104 # childhood-limited maltreatment
data$n_exp_adj[data$es_type_adj=="unstandardised_beta" & data$author=="Thornberry" 
                    & data$maltreatment_period=="12-17"] <- 72 # any adolescent maltreatment

#### Derive N_unexposed (i.e., non-maltreated) for adjusted analyses
data$n_unexp_adj <- NA
# Code as half of N if maltreatment was continuous
data$n_unexp_adj[data$es_type_adj=="unstandardised_beta" & data$maltreatment_type_var=="continuous"] <- 
  data$n_adj[data$es_type_adj=="unstandardised_beta" & data$maltreatment_type_var=="continuous"]/2 
# Code as total N in Isumi, fixed-effects (all participants were examined at each time point)
data$n_unexp_adj[data$es_type_adj=="unstandardised_beta" & data$author=="Isumi"] <- 2920 # fixed-effects design
# Code specific Ns for studies with categorical or binary coding of maltreatment
data$n_unexp_adj[data$es_type_adj=="unstandardised_beta" &data$author=="Lynch"] <- 1640+576 # non-physical/physical mild punish
data$n_unexp_adj[data$es_type_adj=="unstandardised_beta" & data$author=="Kugler"] <- 188 # no cps investigation
data$n_unexp_adj[data$es_type_adj=="unstandardised_beta" & data$author=="Ma"] <- 222 # 9% of total N
data$n_unexp_adj[data$es_type_adj=="unstandardised_beta" & data$author=="Voith"] <- 1978/2 # half observations included in model
data$n_unexp_adj[data$es_type_adj=="unstandardised_beta" & data$author=="Thornberry" 
                      & data$maltreatment_period=="0-11"] <- 907-104 # total N minus child-limited maltreatment
data$n_unexp_adj[data$es_type_adj=="unstandardised_beta" & data$author=="Thornberry" 
                      & data$maltreatment_period=="12-17"] <- 907-72 # total N minus adolescent maltreatment


data %>% dplyr::filter (es_type_adj == "unstandardised_beta") %>% 
  dplyr::select(author, n_exp_adj, n_unexp_adj, n_adj) # Note: n_adj is not always equivalent to n_exp_adj + n_unexp_adj if a fixed-effects design was used

# View maltreatment variable type (continuous or categorical/binary) to see if need Ns for exposed and unexposed
# Note: all studies reporting unadjusted effect sizes as unstandardised betas used continuous maltreatment measures, so divide by 2
data %>% filter (es_type_unadj=="unstandardised_beta" & !is.na(es_unadj)) %>%
  group_by(author)  %>% count(maltreatment_type_var) %>% slice(which.max(n))  %>% arrange(desc(maltreatment_type_var))

#### Derive N_exposed (maltreated) for unadjusted analyses 
data$n_exp_unadj <- NA
# Code as half of N if maltreatment was continuous
data$n_exp_unadj[data$es_type_unadj=="unstandardised_beta" & data$maltreatment_type_var=="continuous"] <- 
  data$n_unadj[data$es_type_unadj=="unstandardised_beta" & data$maltreatment_type_var=="continuous"]/2 

#### Derive N_exposed (maltreated) for unadjusted analyses 
data$n_unexp_unadj <- NA
# Code as half of N if maltreatment was continuous
data$n_unexp_unadj[data$es_type_unadj=="unstandardised_beta" & data$maltreatment_type_var=="continuous"] <- 
  data$n_unadj[data$es_type_unadj=="unstandardised_beta" & data$maltreatment_type_var=="continuous"]/2 

data %>% dplyr::filter (es_type_unadj == "unstandardised_beta" & !is.na(es_adj)) %>% 
  dplyr::select(author, n_exp_unadj, n_unexp_unadj, n_unadj)

################### Convert effect sizes from unstandardised beta to Cohen's d ###################
# To convert unstandardised betas to Cohen's D, adapt the function "esc_B" from the "esc" package
# The original function requires the unstandardised beta, the SD of the outcome, the N exposed and N unexposed
# Because many studies do not report the SD of the outcome, derive from the SE of the difference
# using se / sqrt((1/grp1n) + (1/grp2n)) (from 6.5.2.3 (4) in https://training.cochrane.org/handbook/current/chapter-06#_Ref190897628)

#### Derive function to convert unstandardised beta to Cohen's d
B_to_d <- function (b, se, grp1n, grp2n) 
{
  totaln <- grp1n + grp2n
  sdy <- se / sqrt((1/grp1n) + (1/grp2n)) # convert SE of difference to average SD across groups
  sdpooled <- sqrt(abs(((sdy^2 * (totaln - 1)) - (b^2 * ((grp1n * 
                                                            grp2n)/(grp1n + grp2n))))/(totaln - 2))) # derive pooled SD
  d <- b/sdpooled # derive cohen's d
  v <- (grp1n + grp2n) / (grp1n * grp2n) + (d * d) / (2 * (grp1n + grp2n)) # derive variance of d
  se <- sqrt(v) # SE of d
  return(list(d, se))
}

### Convert B to D
# Adjusted effect size & SE
b_d_adjusted <- do.call( function(es_adj, se_adj, n_exp_adj, n_unexp_adj,...) 
  B_to_d(es_adj, se_adj,  n_exp_adj, n_unexp_adj), data[data$es_type_adj=="unstandardised_beta",])
data$d_adj[data$es_type_adj == "unstandardised_beta"] <- unlist(b_d_adjusted[1]) # Effect size
data$d_se_adj[data$es_type_adj == "unstandardised_beta"] <-  unlist(b_d_adjusted[2]) # SE

# Unadjusted effect size & SE
b_d_unadjusted <- do.call(function(es_unadj, se_unadj, n_exp_unadj, n_unexp_unadj,...) 
  B_to_d(es_unadj, se_unadj, n_exp_unadj, n_unexp_unadj), data[data$es_type_unadj=="unstandardised_beta",])
data$d_unadj[data$es_type_unadj=="unstandardised_beta"] <- unlist(b_d_unadjusted[1]) # Effect size
data$d_se_unadj[data$es_type_unadj=="unstandardised_beta"] <- unlist(b_d_unadjusted[2]) # SE

## Check effect sizes before and after converting to Cohen's d
data %>% dplyr::filter (es_type_adj == "unstandardised_beta" | es_type_unadj == "unstandardised_beta") %>% 
  dplyr::select(author, mh_outcome, es_unadj, d_unadj, es_adj, d_adj)

## Check p-values before and after converting to Cohen's d
# P-values for unadjusted effect sizes  
round(p_function(data$es_unadj[data$es_type_unadj=="unstandardised_beta"], data$se_unadj[data$es_type_unadj=="unstandardised_beta"]),5) # unstandardised beta
round(p_function(data$d_unadj[data$es_type_unadj=="unstandardised_beta"], data$d_se_unadj[data$es_type_unadj=="unstandardised_beta"]),5) # cohen's d

# P-values for adjusted effect sizes  
round(p_function(data$es_adj[data$es_type_adj=="unstandardised_beta"], data$se_adj[data$es_type_adj=="unstandardised_beta"]),5) # unstandardised beta
round(p_function(data$d_adj[data$es_type_adj=="unstandardised_beta"], data$d_se_adj[data$es_type_adj=="unstandardised_beta"]),5) # cohen's d

## Compare p-values - adjusted
data %>% dplyr::filter (es_type_adj == "unstandardised_beta" | es_type_unadj == "unstandardised_beta") %>% 
   dplyr::select(author, maltreatment_type_var, es_adj, d_adj, se_adj, d_se_adj)  %>% 
  mutate(p_original = round(p_function(data$es_adj[data$es_type_adj=="unstandardised_beta"], data$se_adj[data$es_type_adj=="unstandardised_beta"]),5), 
         p_new = round(p_function(data$d_adj[data$es_type_adj=="unstandardised_beta"], data$d_se_adj[data$es_type_adj=="unstandardised_beta"]),5),
         p_diff = p_original-p_new) %>% 
    dplyr::summarise(max = max(p_diff)) # max difference in p-value

## Compare p-values - unadjusted
data %>% dplyr::filter (es_type_unadj == "unstandardised_beta") %>% 
   dplyr::select(author, maltreatment_type_var, es_unadj, d_unadj, se_unadj, d_se_unadj)  %>% 
  mutate(p_original = round(p_function(data$es_unadj[data$es_type_unadj=="unstandardised_beta"], data$se_unadj[data$es_type_unadj=="unstandardised_beta"]),5), 
         p_new = round(p_function(data$d_unadj[data$es_type_unadj=="unstandardised_beta"], data$d_se_unadj[data$es_type_unadj=="unstandardised_beta"]),5),
         p_diff = p_original-p_new) %>% 
  dplyr::summarise(max = max(p_diff)) # max difference in p-value

## Relabel to say that effect sizes for unstandardised betas have now been converted to Cohen's d 
data$es_type_adj[data$es_type_adj == "unstandardised_beta"] <- "cohen_d"
data$es_type_unadj[data$es_type_unadj == "unstandardised_beta"] <- "cohen_d"

table(data$es_type_adj)
table(data$es_type_unadj)
```

## Convert relative risks to Cohen’s d

To convert relative risks to Cohen’s d, we will:

-   derive a variable indexing the prevalence of the mental health
    outcome (p) in non-exposed individuals (needed for the conversion)
-   convert relative risks (RR) to log odds ratios using the formula:
    logOR = log(((1 - p) \* RR) / (1 - RR \* p))
-   check that the conversion from relative risks to log odds ratios is
    consistent with
    [“riskratio_to_oddsratio”](https://easystats.github.io/effectsize/reference/oddsratio_to_riskratio.html)
    function from “effectsize” package
-   convert log odds ratios to Cohen’s d
-   compare effect sizes and p-values before and after converting to
    Cohen’s d

``` r
# View relative risk effect sizes (all from a single study by Obikane, 2018)
data %>% dplyr::filter (es_type_adj == "relative_risk" | es_type_unadj == "relative_risk") %>% 
  dplyr::select(author, year, mh_outcome, es_adj, es_unadj)

#### Derive variable indexing prevalence of the mental health outcome in non-exposed individuals
# In Obikane, can calculate prevalence in non-exposed individuals from Table 1
# as it gives the prevalence of suicidal ideation, plan, and suicide attempt in overall males and females
# and the no. and prevalence of maltreatment among cases
# So you can subtract the number of maltreated participants from the number of cases with
# suicidal ideation, plan, and suicide attempt, then divide by the number of unexposed participants

## Derive variable reflecting the prevalence of the MH outcome in the non-exposed individuals
data$prev_unexp <- NA
# Males - physical abuse, suicidal ideation, 
data$prev_unexp[data$author == "Obikane" & data$perc_female==0 & data$maltreatment_type=="physical_abuse" 
                  & data$mh_outcome=="suicidal_ideation"] <- (290-38) / (1776-99) 
# Males - child neglect, suicidal ideation 
data$prev_unexp[data$author == "Obikane" & data$perc_female==0 & data$maltreatment_type=="physical_neglect" 
                & data$mh_outcome=="suicidal_ideation"] <- (290-15) / (1776-31) 
# Males - any abuse (labelled maltreatment in dataset), suicidal ideation ######
data$prev_unexp[data$author == "Obikane" & data$perc_female==0 & data$maltreatment_type=="maltreatment" 
                & data$mh_outcome=="suicidal_ideation"] <- (290-46) / (1776-122)
# Males - physical abuse, suicidal plan 
data$prev_unexp[data$author == "Obikane" & data$perc_female==0 & data$maltreatment_type=="physical_abuse" 
                & data$mh_outcome=="suicidal_plan"] <- (96-18) / (1776-99) 
# Males - neglect, suicidal plan 
data$prev_unexp[data$author == "Obikane" & data$perc_female==0 &  data$maltreatment_type=="physical_neglect" 
                & data$mh_outcome=="suicidal_plan"] <- (96-8) / (1776-31) 
# Males - any abuse, suicidal plan 
data$prev_unexp[data$author == "Obikane" & data$perc_female==0 & data$maltreatment_type=="maltreatment" 
                & data$mh_outcome=="suicidal_plan"] <- (96-22) / (1776-122)
# Males - physical abuse, suicide attempt 
data$prev_unexp[data$author == "Obikane" & data$perc_female==0 & data$maltreatment_type=="physical_abuse" 
                & data$mh_outcome=="suicide_attempt"] <- (61-9) / (1776-99) 
# Males - neglect, suicide attempt 
data$prev_unexp[data$author == "Obikane" & data$perc_female==0 & data$maltreatment_type=="physical_neglect" 
                & data$mh_outcome=="suicide_attempt"] <- (61-3) / (1776-31) 
# Males - abuse or neglect, suicide attempt 
data$prev_unexp[data$author == "Obikane" & data$perc_female==0 & data$maltreatment_type=="maltreatment" 
                & data$mh_outcome=="suicide_attempt"] <- (61-9) / (1776-122)
# Females - child physical abuse, suicidal ideation 
data$prev_unexp[data$author == "Obikane" & data$perc_female==100 & data$maltreatment_type=="physical_abuse" 
                & data$mh_outcome=="suicidal_ideation"] <- (386-54) / (2016-116) 
# Females - child neglect, suicidal ideation
data$prev_unexp[data$author == "Obikane" & data$perc_female==100 & data$maltreatment_type=="physical_neglect" 
                & data$mh_outcome=="suicidal_ideation"] <- (386-22) / (2016-43) 
# Females - child abuse or neglect, suicidal ideation 
data$prev_unexp[data$author == "Obikane" & data$perc_female==100 & data$maltreatment_type=="maltreatment" 
                & data$mh_outcome=="suicidal_ideation"] <- (386-65) / (2016-146) 
# Females - child physical abuse, suicidal plan 
data$prev_unexp[data$author == "Obikane" & data$perc_female==100 & data$maltreatment_type=="physical_abuse" 
                & data$mh_outcome=="suicidal_plan"] <- (94-17) / (2016-116) 
# Females - child neglect, suicidal plan 
data$prev_unexp[data$author == "Obikane" & data$perc_female==100 & data$maltreatment_type=="physical_neglect" 
                & data$mh_outcome=="suicidal_plan"] <- (94-6) / (2016-43) 
# Females - abuse or neglect, suicidal plan 
data$prev_unexp[data$author == "Obikane" & data$perc_female==100 & data$maltreatment_type=="maltreatment" 
                & data$mh_outcome=="suicidal_plan"] <- (94-20) / (2016-146)  
# Females - child physical abuse, suicide attempt 
data$prev_unexp[data$author == "Obikane" & data$perc_female==100 & data$maltreatment_type=="physical_abuse" 
                & data$mh_outcome=="suicide_attempt"] <- (111-24) / (2016-116)
# Females - child neglect, suicide attempt 
data$prev_unexp[data$author == "Obikane" & data$perc_female==100 & data$maltreatment_type=="physical_neglect" 
                & data$mh_outcome=="suicide_attempt"] <- (111-13) / (2016-43) 
# Females - abuse or neglect, suicide attempt 
data$prev_unexp[data$author == "Obikane" & data$perc_female==100 & data$maltreatment_type=="maltreatment" 
                & data$mh_outcome=="suicide_attempt"] <- (111-29) / (2016-146) 

data %>% dplyr::filter (es_type_adj == "relative_risk") %>% 
  dplyr::select(author, year, perc_female, maltreatment_type, mh_outcome, es_adj, es_unadj, prev_unexp)

################### Derive function to convert RR to OR ###################
# RR = OR/ (1-p + (p * OR)) 
# OR = ((1 - p) * RR) / (1 - RR * p).
# where RR is the relative risk, OR is the odds ratio, and p is the control event rate
# e.g. the prevalence rate in the non-exposed individuals
# formula recommended by Grant et al., BMJ 2014: https://www.bmj.com/content/348/bmj.f7450
# and referenced in https://www.researchgate.net/post/Can-we-Convert-Risk-Ratio-to-Odd-Ratios-what-is-the-baseline-risk-p-of-a-selection-of-populations-to-either-suffer-die-from-coronary-heart-disease
# Note: when the incidence rate is very low (<10%), the odds ratio is very similar to the relative risk

RR_to_OR <- function(RR, RR_lowCI, RR_upCI, p) {
  OR <- ((1 - p) * RR) / (1 - RR * p)
  log_OR <- log(OR)
  OR_lowCI <- ((1 - p) * RR_lowCI) / (1 - RR_lowCI * p)
  OR_upperCI <- ((1 - p) * RR_upCI) / (1 - RR_upCI * p)
  log_OR_SE <- (log(OR_upperCI) - log(OR_lowCI))/3.92
  return(list(log_OR, log_OR_SE))
}

####### Convert RR to log OR and then log OR to Cohen's d 
## RR to log OR 
## Adjusted results
rr_logor_adjusted <- do.call( function(es_adj, lowCI_adj, upCI_adj, prev_unexp,...) 
  RR_to_OR(es_adj, lowCI_adj, upCI_adj, prev_unexp), data[data$es_type_adj=="relative_risk",])
data$d_adj[data$es_type_adj=="relative_risk"] <- unlist(rr_logor_adjusted[1]) # Effect size
data$d_se_adj[data$es_type_adj=="relative_risk"] <- unlist(rr_logor_adjusted[2]) # Standard error

## Unadjusted results
rr_logor_unadjusted <- do.call( function(es_unadj, lowCI_unadj, upCI_unadj, prev_unexp,...) 
  RR_to_OR(es_unadj, lowCI_unadj, upCI_unadj, prev_unexp), data[data$es_type_unadj=="relative_risk",])
data$d_unadj[data$es_type_unadj=="relative_risk"] <- unlist(rr_logor_unadjusted[1]) # Effect size
data$d_se_unadj[data$es_type_unadj=="relative_risk"] <- unlist(rr_logor_unadjusted[2]) # SE

## Check that ORs are the same as those obtained from using "riskratio_to_oddsratio" function from "effectsize" package
# Adjusted
log(do.call( function(es_adj, prev_unexp,...)riskratio_to_oddsratio(es_adj, prev_unexp), 
                  data[data$es_type_adj=="relative_risk",]))
data$d_adj[data$es_type_adj=="relative_risk"]
# Unadjusted
log(do.call( function(es_unadj, prev_unexp,...)riskratio_to_oddsratio(es_unadj, prev_unexp), 
                  data[data$es_type_unadj=="relative_risk",]))
data$d_unadj[data$es_type_unadj=="relative_risk"]

## Log OR to Cohen's d
# Cohen's d = log_odds * (sqrt(3)) / 3.14159; from: https://www.meta-analysis.com/downloads/Meta-analysis%20Converting%20among%20effect%20sizes.pdf
# Adjusted ES
data$d_adj[data$es_type_adj=="relative_risk"] <- data$d_adj[data$es_type_adj=="relative_risk"]*(sqrt(3)/3.14159)
# Unadjusted ES
data$d_unadj[data$es_type_unadj=="relative_risk"] <- data$d_unadj[data$es_type_unadj=="relative_risk"]*(sqrt(3)/3.14159)

## Convert log odds SE to Cohen's d SE
# Adjusted SE
data$d_se_adj[data$es_type_adj=="relative_risk"] <- SE_logodds_D(data$d_se_adj[data$es_type_adj=="relative_risk"])
# Unadjusted SE
data$d_se_unadj[data$es_type_unadj=="relative_risk"] <- SE_logodds_D(data$d_se_unadj[data$es_type_unadj=="relative_risk"])

## Check effect sizes before and after converting to Cohen's d
data %>% dplyr::filter (es_type_adj == "relative_risk") %>% 
  dplyr::select(author, perc_female, maltreatment_type, mh_outcome, es_unadj, d_unadj, es_adj, d_adj)

## Check p-values before and after converting to Cohen's d
round(p_function(log(data$es_adj)[data$es_type_adj=="relative_risk"], data$se_adj[data$es_type_adj=="relative_risk"]),5) # relative risk
round(p_function(data$d_adj[data$es_type_adj=="relative_risk"], data$d_se_adj[data$es_type_adj=="relative_risk"]),5) # cohen's d

## Relabel to say that effect sizes for relative risk have now been converted to Cohen's d 
data$es_type_adj[data$es_type_adj=="relative_risk"] <- "cohen_d"
data$es_type_unadj[data$es_type_unadj=="relative_risk"] <- "cohen_d"

table(data$es_type_adj)
table(data$es_type_unadj)
```

## Convert hazard ratios to Cohen’s d

To convert relative risks to Cohen’s d, we will:

-   derive a variable indexing the prevalence of the mental health
    outcome (p) in non-exposed individuals (needed for the conversion)
-   convert hazard ratios to relative risks using the formula: RR = (1 -
    expHR \* ln(1 - p) )/ p)
-   convert relative risks to log odds ratios using the formula: logOR =
    log(((1 - p) \* RR) / (1 - RR \* p))
-   check that the conversion from relative risks to log odds ratios is
    consistent with
    [“riskratio_to_oddsratio”](https://easystats.github.io/effectsize/reference/oddsratio_to_riskratio.html)
    function from “effectsize” package
-   convert log odds ratios to Cohen’s d
-   compare effect sizes and p-values before and after converting to
    Cohen’s d

``` r
# View relative risk effect sizes (all from a single study by Capusan, 2021)
data %>% dplyr::filter (es_type_adj == "hazard_ratio" | es_type_unadj == "hazard_ratio") %>% 
  dplyr::select(author, year, mh_outcome, es_adj, es_unadj)

#### Derive variable indexing prevalence of the mental health outcome in non-exposed individuals
data$prev_unexp[data$author=="Capusan" & data$year==2021] <- (77+531) / (1383+1979) #SUD in matched healthy controls +clinical controls

### Derive function to convert hazard ratios to relative risks 
# RR = 1 - exp (HR * ln(1 - p) )/ p)
# where HR is the hazard ratio and p is the control event rate
# e.g. the prevalence rate in the non-exposed individuals
# formula recommended by  Shor et al. (2017) doi: 10.1016/j.socscimed.2017.05.049 

HR_to_RR <- function(HR, HR_lowCI, HR_upCI, p) {
  RR <- ( 1 - exp(HR * log(1-p))) / p
  RR_lowCI <- ( 1 - exp(HR_lowCI * log(1-p))) / p
  RR_upCI <- ( 1 - exp(HR_upCI * log(1-p))) / p
  return(c(RR, RR_lowCI, RR_upCI))
}

# Convert unadjusted HR to RR
rr_Capusan_unadj <- HR_to_RR(HR = data$es_unadj[data$author=="Capusan" & data$year==2021],
                       HR_lowCI = data$lowCI_unadj[data$author=="Capusan" & data$year==2021],
                       HR_upCI = data$upCI_unadj[data$author=="Capusan" & data$year==2021],
                       p = data$prev_unexp[data$author=="Capusan" & data$year==2021])

# Convert adjusted HR to RR
rr_Capusan_adj <- HR_to_RR(HR = data$es_adj[data$author=="Capusan" & data$year==2021],
                       HR_lowCI = data$lowCI_adj[data$author=="Capusan" & data$year==2021],
                       HR_upCI = data$upCI_adj[data$author=="Capusan" & data$year==2021],
                       p = data$prev_unexp[data$author=="Capusan" & data$year==2021])

####### Convert RR to log OR and then log OR to Cohen's d 
## RR to log OR 
# Unadjusted results
or_Capusan_unadj <- RR_to_OR(RR = rr_Capusan_unadj[1],
                       RR_lowCI = rr_Capusan_unadj[2],
                       RR_upCI = rr_Capusan_unadj[3],
                       p = data$prev_unexp[data$author=="Capusan" & data$year==2021])

# Adjusted results
or_Capusan_adj <- RR_to_OR(RR = rr_Capusan_adj[1],
                       RR_lowCI = rr_Capusan_adj[2],
                       RR_upCI = rr_Capusan_adj[3],
                       p = data$prev_unexp[data$author=="Capusan" & data$year==2021])

## Check that ORs are the same as those obtained from using "riskratio_to_oddsratio" function from "effectsize" package
# Unadjusted
riskratio_to_oddsratio(rr_Capusan_unadj[1], data$prev_unexp[data$author=="Capusan" & data$year==2021])
exp(as.numeric(or_Capusan_unadj)[1])
# Adjusted
riskratio_to_oddsratio(rr_Capusan_adj[1], data$prev_unexp[data$author=="Capusan" & data$year==2021])
exp(as.numeric(or_Capusan_adj)[1])

## Log OR to Cohen's d
# Cohen's d = log_odds * (sqrt(3)) / 3.14159; from: https://www.meta-analysis.com/downloads/Meta-analysis%20Converting%20among%20effect%20sizes.pdf
# Unadjusted ES
data$d_unadj[data$author=="Capusan" & data$year==2021] <- as.numeric(or_Capusan_unadj[1])*(sqrt(3)/3.14159)
# Adjusted ES
data$d_adj[data$author=="Capusan" & data$year==2021] <- as.numeric(or_Capusan_adj[1])*(sqrt(3)/3.14159)

## Convert log odds SE to Cohen's d SE
# Unadjusted SE
data$d_se_unadj[data$author=="Capusan" & data$year==2021] <- SE_logodds_D(as.numeric(or_Capusan_unadj[2]))
# Adjusted SE
data$d_se_adj[data$author=="Capusan" & data$year==2021] <- SE_logodds_D(as.numeric(or_Capusan_adj[2]))

## Relabel to say that effect sizes for relative risk have now been converted to Cohen's d 
data$es_type_adj[data$es_type_adj=="hazard_ratio"] <- "cohen_d"
data$es_type_unadj[data$es_type_unadj=="hazard_ratio"] <- "cohen_d"

table(data$es_type_adj)
table(data$es_type_unadj)
```

## Convert correlations/standardised betas to Cohen’s d

To convert correlation coefficients and standardised betas to Cohen’s d,
we will:

-   convert effect sizes using the formula: d = 2\*r / sqrt(1-r^2)
-   convert variances using the formula: variance_d = 4\*var_r /
    ((1-r^2) ^3)
-   check that the conversion is consistent with the
    [“res”](https://www.rdocumentation.org/packages/compute.es/versions/0.2-5/topics/res)
    function from the “compute.es” package
-   compare effect sizes before and after converting to Cohen’s d

``` r
# View studies and number of effect sizes with correlations or standardised betas
data %>% dplyr::filter (es_type_adj == "correlation" | es_type_unadj == "correlation" |
                          es_type_adj == "standardised_beta" | es_type_unadj == "standardised_beta") %>% 
  group_by(author) %>%  dplyr::summarise(n_effectSizes = n())

###### Derive function to convert correlations/standardised betas to Cohen's d
# formula for variance taken from: https://www.meta-analysis.com/downloads/Meta-analysis%20Converting%20among%20effect%20sizes.pdf
cor_to_D <- function(r, r_SE) {
  d <- 2*r / sqrt(1-r^2) # convert r to Cohen's d 
  var_r <- r_SE^2 # calculate variance of r
  var_d <- 4*var_r / ((1-r^2)^3) # calculate variance of Cohen's d
  se_d <- sqrt(var_d) # calculate SE of Cohen's d
  return(list(d, se_d))
}

###### Convert correlations/standardised betas to Cohen's d 
## Adjusted effect size & SE
r_d_adjusted <- do.call( function(es_adj, se_adj,...) cor_to_D(es_adj, se_adj),
                                data [data$es_type_adj == "correlation" | data$es_type_adj == "standardised_beta",])
data$d_adj[data$es_type_adj=="correlation" | data$es_type_adj=="standardised_beta"] <- unlist(r_d_adjusted[1]) # Effect size
data$d_se_adj[data$es_type_adj=="correlation" | data$es_type_adj=="standardised_beta"] <-  unlist(r_d_adjusted[2]) # SE

## Unadjusted effect size & SE
r_d_unadjusted <- do.call( function(es_unadj, se_unadj,...) cor_to_D(es_unadj, se_unadj),
                          data [data$es_type_unadj == "correlation" | data$es_type_unadj == "standardised_beta",])
data$d_unadj[data$es_type_un == "correlation" | data$es_type_unadj=="standardised_beta"] <- unlist(r_d_unadjusted[1]) # Effect size
data$d_se_unadj[data$es_type_un == "correlation" | data$es_type_unadj=="standardised_beta"] <- unlist(r_d_unadjusted[2]) # SE

#### Check results against those obtained from the "res" function from the "compute.es" package
## Adjusted results
# Effect sizes
r_d_check_adj <- do.call( function(es_adj, se_adj,  n_adj,...) 
  res(es_adj, se_adj^2,  n_adj, dig=5), data[data$es_type_adj=="correlation" | data$es_type_adj == "standardised_beta",])
round(data$d_adj[data$es_type_adj=="correlation" | data$es_type_adj=="standardised_beta"], 5)
r_d_check_adj$d
# standard errors
data$d_se_adj[data$es_type_adj=="correlation" | data$es_type_adj=="standardised_beta"]
sqrt(r_d_check_adj$var.d)

## Unadjusted results
# Effect sizes
r_d_check_unadj <- do.call( function(es_unadj, se_unadj,  n_unadj,...) 
  res(es_unadj, se_unadj^2,  n_unadj, dig=5), data[data$es_type_unadj=="correlation" | data$es_type_unadj == "standardised_beta",])
round(data$d_unadj[data$es_type_unadj=="correlation" | data$es_type_unadj=="standardised_beta"],4)
r_d_check_unadj$d
# standard errors
data$d_se_unadj[data$es_type_unadj=="correlation" | data$es_type_unadj=="standardised_beta"] 
sqrt(r_d_check_unadj$var.d)

## Check effect sizes before and after converting to Cohen's d
data %>% dplyr::filter (es_type_adj == "correlation" | data$es_type_adj=="standardised_beta") %>% 
  dplyr::select(author, mh_outcome, es_unadj, d_unadj, es_adj, d_adj)

## Check p-values before and after converting to Cohen's d
# before converting
p_function(data$es_adj[data$es_type_adj == "correlation" | data$es_type_adj=="standardised_beta"], 
           data$se_adj[data$es_type_adj == "correlation" | data$es_type_adj=="standardised_beta"])
# after converting
p_function(data$d_adj[data$es_type_adj == "correlation" | data$es_type_adj=="standardised_beta"], 
           data$d_se_adj[data$es_type_adj == "correlation" | data$es_type_adj=="standardised_beta"])

## Relabel to say that the effect sizes are now Cohen's d 
data$es_type_adj[data$es_type_adj == "correlation" | data$es_type_adj=="standardised_beta"] <- "cohen_d"
data$es_type_unadj[data$es_type_unadj == "correlation" | data$es_type_unadj=="standardised_beta"] <- "cohen_d"

table(data$es_type_adj)
table(data$es_type_unadj)
```

# Derive additional variables needed for meta-analysis

Before moving onto analysing the results, we will derive the following
variables:

-   study reference (combining the author name and year of publication)
-   effect size ID (needed for multi-level meta-analysis models)
-   variance of the effect sizes (needed for meta-analysis)
-   mental health category variable (needed to test moderation by type
    of mental health problem)

We will also exclude effect sizes based on negative binomial regression
coefficients as they could not be converted to Cohen’s d.

``` r
# Derive study reference
data$ref <- paste0(data$author, "_", data$year)

# Derive effect size ID
data$es_id <- 1:nrow(data)

# Derive variance for all studies
data$var_adj <- data$d_se_adj^2
data$var_unadj <- data$d_se_unadj^2

# Exclude effect sizes based on negative binomial regression coefficients
data <- subset(data, es_type_adj!="negative_binomial_regression_coefficient")

# Derive specific categories for different mental health outcomes
unique(data[c("mh_outcome")]) # Show all mental health outcomes
data$mh_category <- NA
data$mh_category[grepl("depress", data$mh_outcome, ignore.case = TRUE)] <- "depression"
data$mh_category[grepl("anxiety|panic|phobia|GAD", data$mh_outcome, ignore.case = TRUE)] <- "anxiety"
data$mh_category[grepl("ADHD|inattention", data$mh_outcome, ignore.case = TRUE)]<- "ADHD"
data$mh_category[grepl("drug|cannabis|cocaine|opioid|sedative|stimulant",  data$mh_outcome, ignore.case = TRUE)] <- "drug_abuse"
data$mh_category[grepl("alcohol", data$mh_outcome, ignore.case = TRUE)] <- "alcohol_abuse"
data$mh_category[grepl("substance_use_disorder", data$mh_outcome, ignore.case = TRUE)] <- "substance_use_disorder"
data$mh_category[grepl("personality", data$mh_outcome, ignore.case = TRUE)] <- "personality_disorder"
data$mh_category[grepl("psychosis|thought_disorder|psychotic", data$mh_outcome, ignore.case = TRUE)] <- "psychosis"
data$mh_category[grepl("conduct|antisocial_behaviour|aggress|opposition|general_offending|crime|arrest", data$mh_outcome, ignore.case = TRUE)] <- "conduct_problems"
data$mh_category[grepl("ASD|autism", data$mh_outcome, ignore.case = TRUE) ] <- "autism"
data$mh_category[grepl("suicidal_ideation|suicidal_thoughts|suicidal_plan", data$mh_outcome, ignore.case = TRUE)] <- "suicidal_ideation"
data$mh_category[grepl("suicide_attempt", data$mh_outcome, ignore.case = TRUE)] <- "suicide_attempt"
data$mh_category[grepl("self-harm|self-injury|self_injury|self_harm", data$mh_outcome, ignore.case = TRUE)] <- "self_harm"
data$mh_category[grepl("emotional|internalising|trauma|psychosocial_distress", data$mh_outcome, ignore.case = TRUE)] <- "internalising_problems"
data$mh_category[grepl("externalising", data$mh_outcome, ignore.case = TRUE)] <- "externalising_problems"
data$mh_category[grepl("p-factor|psychopathology|behavioural_difficulties", data$mh_outcome, ignore.case = TRUE)] <- "psychopathology_broad"
data$mh_category[grepl("bulimia", data$mh_outcome, ignore.case = TRUE)] <- "bulimia"
data$mh_category <- as.factor(data$mh_category)
table(data$mh_category, useNA="always")

# Check new mental health category variable against original labels
data %>%  dplyr::select(ref, mh_outcome, mh_category)
table(data$mh_category)

# Derive broad categories for different mental health outcomes
data$mh_category_broad <- NA
data$mh_category_broad[data$mh_category=="ADHD" | data$mh_category=="alcohol_abuse" | 
                         data$mh_category=="conduct_problems" |  data$mh_category=="drug_abuse" |
                         data$mh_category=="substance_use_disorder" |
                         data$mh_category=="externalising_problems"] <- "externalising"
data$mh_category_broad[data$mh_category=="anxiety" | data$mh_category=="depression" | 
                         data$mh_category=="internalising_problems" |  data$mh_category=="self_harm" |
                         data$mh_category=="suicidal_ideation" |
                         data$mh_category=="suicide_attempt" |
                         data$mh_outcome=="bulimia_diagnosis"] <- "internalising"
data$mh_category_broad[data$mh_category=="psychosis"] <- "psychosis"
data$mh_category_broad[data$mh_category=="personality_disorder"] <- "personality_disorder"
data$mh_category_broad <- as.factor(data$mh_category_broad)
table(data$mh_category_broad, useNA="always")
# Show outcomes not included in broad categories
data %>% dplyr::filter (is.na(mh_category_broad)) %>% 
  dplyr::select(ref, mh_outcome)

# Add short distinct cohort identifier (e.g., distinguishing between male and female only samples in Obikane, Dinwiddie and Alvanzo)
data$cohort_short_distinct <- data$cohort_short
data$cohort_short_distinct[data$cohort_short=="ATR" & data$perc_female==100] <- "ATR (females)"
data$cohort_short_distinct[data$cohort_short=="ATR" & data$perc_female==0] <- "ATR (males)"
data$cohort_short_distinct[data$cohort_short=="J-SHINE" & data$perc_female==100] <- "J-SHINE (females)"
data$cohort_short_distinct[data$cohort_short=="J-SHINE" & data$perc_female==0] <- "J-SHINE (males)"
data$cohort_short_distinct[data$cohort_short=="NESARC" & data$perc_female==100] <- "NESARC (females)"
data$cohort_short_distinct[data$cohort_short=="NESARC" & data$perc_female==0] <- "NESARC (males)"
table(data$cohort_short_distinct)
data %>%  dplyr::select(ref, cohort_short, perc_female, cohort_short_distinct)

# View studies to check the information is all there
data %>% dplyr::filter(es_type_adj=="cohen_d" | es_type_unadj=="cohen_d")  %>%
  dplyr::select(ref, d_unadj, var_unadj, d_adj, var_adj)
```

# Results

## Search results - descriptives

We will calculate descriptives on the number of studies, number of
cohorts, number of effect sizes, sample size, and sample characteristics
(age at assessment and percentage female). The results text reporting
these descriptives is shown below the code.

``` r
# Number of studies
k_studies <- data %>%   
  group_by(ref) %>% # group by study reference
  dplyr::summarise(m = max(ref)) %>%  # select one row per study 
  nrow() # count number of rows

# Number of cohorts
k_cohort <- data %>%   
  group_by(cohort) %>% # group by study cohort
  dplyr::summarise(m = max(cohort)) %>%  # select one row per study 
  nrow() # count number of rows

# Number of distinct samples (k_cohort +3 as Obikane, Dinwiddie & Alvanzo split into male and female samples)
k_distinct_cohort <- data %>%   
  group_by(cohort_short_distinct) %>% # group by distinct study cohort
  dplyr::summarise(m = max(cohort_short_distinct)) %>%  # select one row per distinct cohort 
  nrow() # count number of rows

# Obtain total N for QE-adjusted estimate
N_adjusted <- data %>%   
  group_by(cohort_short_distinct) %>%   
  dplyr::summarise(m = max(n_adj)) %>%  # make new column indexing largest sample size per distinct cohort
  dplyr::summarise(sum = sum(m)) # sum sample sizes across distinct cohorts

# Obtain total N for unadjusted estimate
N_unadjusted <- data %>%   
  dplyr::filter(!is.na(d_unadj)) %>%  # select only rows where data is available for the unadjusted ES
  group_by(cohort_short_distinct) %>%  
  dplyr::summarise(m = max(n_unadj, na.rm=TRUE)) %>%  # make new column indexing largest sample size per distinct cohort
  dplyr::summarise(sum = sum(m)) # sum sample sizes across distinct cohorts

# Obtain k studies for unadjusted estimate
k_unadjusted <- data %>%   
  dplyr::filter(!is.na(d_unadj)) %>%  # select only rows where data is available for the unadjusted ES
  group_by(ref) %>%  # group by study reference
  dplyr::summarise(m = max(ref)) %>%  # select one row per study
  nrow() # count number of studies

# Average proportion female
perc_female <- data %>%   
  group_by(cohort_short_distinct) %>% # group by study cohort
  dplyr::summarise(mean_perc_female = mean(perc_female)) %>%  # calculate the average % female per cohort
  dplyr::summarise(overall=mean(mean_perc_female, na.rm=TRUE)) # caclulate the average % female across cohorts

# Age at mental health assessment
age <- data %>%   
  group_by(cohort_short_distinct) %>%   # group by study cohort
  dplyr::summarise(mean_age_mh = mean(mh_age_assess_r)) %>% # calculate the average age at assessment per cohort
  dplyr::summarise(overall=mean(mean_age_mh, na.rm=TRUE)) # calculate the average age at assessment across cohorts

# Number of effect sizes
n_es_adj <- length(data$es_adj[!is.na(data$es_adj)]) # no. of QE-adjusted effect sizes
n_es_unadj <- length(data$es_unadj[!is.na(data$es_unadj)]) # no. of unadjusted effect sizes

# Number of different quasi-experimental studies
n_qe_method <- data %>%
count(ref, qe_method) %>% 
group_by(qe_method) %>% count(qe_method)
```

The study selection procedure is summarized in Figure S1 in the online
supplement. We identified 34 quasi-experimental studies on the
association between child maltreatment and mental health (see Table S7
in the online supplement for study details). These studies were based on
29 distinct cohorts, comprising 54,646 participants in adjusted analyses
(56.72% female, with a mean age of 28.15 years at mental health
assessment). From these studies, we obtained 156 effect sizes for the
association between child maltreatment and mental health based on
adjusted analyses and 103 effect sizes based on unadjusted analyses.

## Quasi-experimental evidence on the association between child maltreatment and mental health

In the next section we will:

-   conduct a multi-level random-effects meta-analysis for the
    association between child maltreatment and mental health in
    quasi-experimental *adjusted* studies
-   conduct a multi-level random-effects meta-analysis for the
    association between child maltreatment and mental health in
    *unadjusted* quasi-experimental studies
-   aggregate effect sizes across each study to calculate a single
    study-average effect size (for QE-adjusted effects)
-   create a forest plot of study-average effect sizes (Figure 2)
-   create a forest plot with unadjusted effect sizes

The results text is shown below the code.

### Quasi-experimental adjusted meta-analysis - multi-level random-effects meta-analysis model

``` r
# Run multi-level random-effects meta-analysis for quasi-experimental adjusted studies
res_adjusted <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), 
                       data=data)
cluster_ci_adjusted <- metafor::robust(res_adjusted, cluster=data$ref)
(format(round(c(cluster_ci_adjusted$b, cluster_ci_adjusted$ci.lb, cluster_ci_adjusted$ci.ub), 2),nsmall=2)) # estimate & robust CIs
   
# Calculate I2 statistic
W <- diag(1/data$var_adj[!is.na(data$d_adj) & !is.na(data$var_adj)])
X <- model.matrix(res_adjusted)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2 <- 100 * sum(res_adjusted$sigma2) / (sum(res_adjusted$sigma2) + 
                                    (res_adjusted$k-res_adjusted$p)/sum(diag(P)))

# Run multi-level random-effects meta-analysis for quasi-experimental adjusted studies that also reported adjusted effects
res_adjusted_s <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), 
                       data=data, subset = !is.na(d_unadj) & !is.na(var_unadj))
cluster_ci_adjusted_s <- metafor::robust(res_adjusted_s, cluster=data$ref)
(format(round(c(cluster_ci_adjusted_s$b, cluster_ci_adjusted_s$ci.lb, cluster_ci_adjusted_s$ci.ub), 2),nsmall=2)) # estimate & robust CIs
```

### Unadjusted meta-analysis - random-effects multi-level meta-analysis model

``` r
# Run multi-level random-effects meta-analysis for unadjusted quasi-experimental studies
res_unadjusted <- rma.mv(d_unadj, var_unadj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), 
                         data=data, slab=ref)
cluster_ci_unadjusted <- robust(res_unadjusted, cluster=data$ref)
(format(round(c(cluster_ci_unadjusted$b, cluster_ci_unadjusted$ci.lb, cluster_ci_unadjusted$ci.ub), 2),nsmall=2)) # estimate & robust CIs
   
# Calculate I2 statistic
#mlm.variance.distribution(x = res_adjusted) #see: https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/fitting-a-three-level-model.html
W <- diag(1/data$var_unadj[!is.na(data$d_unadj) & !is.na(data$var_unadj)])
X <- model.matrix(res_unadjusted)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_unadj <- 100 * sum(res_unadjusted$sigma2) / (sum(res_unadjusted$sigma2) + 
                                    (res_unadjusted$k-res_unadjusted$p)/sum(diag(P)))
I2_unadj

# Derive dataframe to compare effect sizes in adjusted and unadjusted results
dat.comp <- data.frame(estimate = c(coef(res_adjusted), coef(res_unadjusted)), 
                       stderror = c(cluster_ci_adjusted$se, cluster_ci_unadjusted$se),
                       meta = c("adjusted","unadjusted"), 
                       sigma2 = round(c(sum(res_adjusted$sigma2), sum(res_unadjusted$sigma2)),3))
dat.comp

# Compare estimates
diff_adj_unadj <- rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=dat.comp, digits=3)

# Estimate decrease in effect sizes between QE-adjusted and unadjusted effect sizes
round( (coef(res_unadjusted) - coef(res_adjusted) ) / coef(res_unadjusted)*100,2)
```

A multilevel random-effects meta-analysis model showed a small
association between childhood maltreatment and mental health problems in
quasi-experimental studies (Cohen’s d=0.31, 95% CI=0.24-0.37,
![I^{2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;I%5E%7B2%7D "I^{2}")=76.27;
Figure 1). This meta-analytic association between child maltreatment and
mental health from quasi-experimental studies was 45.04% smaller than
that obtained in unadjusted analyses (k=20; Cohen’s d=0.56, 95%
CI=0.41-0.71,
![I^{2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;I%5E%7B2%7D "I^{2}")=97.29
\[see Figure S2 in the online supplement\]; p-value for
difference=0.001). This effect size difference was consistent when the
quasi-experimental meta-analysis was restricted to studies reporting
both unadjusted and adjusted effect sizes (k=20; Cohen’s d=0.26, 95%
CI=0.17-0.35)

### Aggregate effect sizes across studies to derive a study-average effect size to include in the forest plot

``` r
######## Aggregated quasi-experimental adjusted meta-analysis ######## 
data$nesting_var <- ave(data$d_adj, data$ref, FUN = seq_along) # derive a nesting variable numbering each ES in each study
data %>%  group_by(ref) %>% dplyr::select(nesting_var) %>% print(n=400) # show nesting variable

data$DupCheck <- paste0(data$ref) 
vecUniqueEst <- levels(as.factor(data$DupCheck))
listAllSample=list() # make a list to store the results in

# Run loop to aggregate multiple within-study effect size estimates into a single average
for ( i in 1:length(vecUniqueEst) ) {
  df_selected = subset(data, DupCheck==as.character(vecUniqueEst[i])) 
  
  # If there are multiple effect sizes included in the study, run the aggregate model to obtain a study-average effect
  if  ( ( length(df_selected$nesting_var)>1 )==TRUE) {
    print("run aggregate model")
    agg <- agg(id=ref, es=d_adj, var=var_adj, cor=0.6, data=df_selected)
    df_out <- data.frame(ref=df_selected$ref, cohort=df_selected$cohort, N=df_selected$n_adj,
                         qe_method=df_selected$qe_method, ES=agg$es, var=agg$var, se=sqrt(agg$var))
  }
  
  # If there is only 1 effect size included in the study, keep the original effect size estimate
  else {
    print("keep individual effect size estimate")
    df_out <- data.frame(ref=df_selected$ref, cohort=df_selected$cohort,
                         N=df_selected$n_adj,
                         qe_method=df_selected$qe_method,ES=df_selected$d_adj, var=df_selected$var_adj, se=df_selected$d_se_adj)
  }
  listAllSample[[i]]=df_out
}

datMeta_inc = do.call(rbind,listAllSample)

## Derive dataset with each study reference represented only once
agg_data <- as.data.frame(datMeta_inc %>% group_by(ref) %>% filter(row_number()==1) )
agg_data$es_id <- 1:nrow(agg_data) # Derive effect size ID variable
agg_data$year <- sapply(strsplit(agg_data$ref, "_"), `[`, 2) # Derive author variable
agg_data$author <- sapply(strsplit(agg_data$ref, "_"), `[`, 1) # Derive year of publication variable

# Order dataset according to year of publication
agg_data <- agg_data[order(agg_data$year, agg_data$author),]

# Order dataset according to effect size
agg_data <- agg_data[order(agg_data$ES),]
 
# Run meta-analysis on aggregated sample (only used to inform the forest plot for Figure 2)
agg_meta <- rma(ES, var, data=agg_data, slab=c(ref))
#agg_meta <- rma(ES, var, data=agg_data, slab=paste0(ref, "_", ref_no))

# Add reference number for forest plot
agg_data$ref_no <- NA
agg_data$ref_no[agg_data$ref=="Young-Wolff_2011"] <- 32
agg_data$ref_no[agg_data$ref=="Alvanzo_2020"] <- 48
agg_data$ref_no[agg_data$ref=="Berenz_2013"] <- 33
agg_data$ref_no[agg_data$ref=="Bornovalova_2013"] <- 30
agg_data$ref_no[agg_data$ref=="Li_2021"] <- 54
agg_data$ref_no[agg_data$ref=="Voith_2014"] <- 43
agg_data$ref_no[agg_data$ref=="Schwartz_2019"] <- 38
agg_data$ref_no[agg_data$ref=="Kullberg_2020"] <- 40
agg_data$ref_no[agg_data$ref=="Capusan_2016"] <- 31
agg_data$ref_no[agg_data$ref=="Baldwin_2019"] <- 26
agg_data$ref_no[agg_data$ref=="Magnusson_2012"] <- 27
agg_data$ref_no[agg_data$ref=="Stern_2018"] <- 24
agg_data$ref_no[agg_data$ref=="Lynch_2006"] <- 41
agg_data$ref_no[agg_data$ref=="Zvara_2017"] <- 47
agg_data$ref_no[agg_data$ref=="Kugler_2019"] <- 50
agg_data$ref_no[agg_data$ref=="Lecei_2019"] <- 52
agg_data$ref_no[agg_data$ref=="Kendler_2000"] <- 36
agg_data$ref_no[agg_data$ref=="Thornberry_2010"] <- 17
agg_data$ref_no[agg_data$ref=="Dinwiddie_2000"] <- 35
agg_data$ref_no[agg_data$ref=="Isumi_2021"] <- 55
agg_data$ref_no[agg_data$ref=="Dinkler_2017"] <- 25
agg_data$ref_no[agg_data$ref=="Alemany_2013"] <- 29
agg_data$ref_no[agg_data$ref=="Nelson_2002"] <- 37
agg_data$ref_no[agg_data$ref=="Golm_2020"] <- 45
agg_data$ref_no[agg_data$ref=="Obikane_2018"] <- 49
agg_data$ref_no[agg_data$ref=="Riggins-Caspers_2003"] <- 42
agg_data$ref_no[agg_data$ref=="Nelson_2006"] <- 34
agg_data$ref_no[agg_data$ref=="Sonuga-Barke_2017"] <- 14
agg_data$ref_no[agg_data$ref=="Gerin_2019"] <- 46
agg_data$ref_no[agg_data$ref=="Schaefer_2017"] <- 28
agg_data$ref_no[agg_data$ref=="Ma_2018"] <- 51
agg_data$ref_no[agg_data$ref=="Capusan_2021"] <- 53
agg_data$ref_no[agg_data$ref=="Barrigon_2015"] <- 39
agg_data$ref_no[agg_data$ref=="Beckett_2002"] <- 44
```

#### Forest plot depicting the study-average effects of child maltreatment on mental health from quasi-experimental studies.

``` r
# Forest plot of average estimates
par(mar=c(2, 4.1, 2, 2.1))
forest(x=agg_data$ES, sei=agg_data$se, 
       slab=paste0(gsub(".{5}$", " et al. ", agg_data$ref), "(", agg_data$ref_no, ")"),
       xlab = "Cohen's d", 
       ylim=c(-2,35.5),
       alim=c(-0.5,2.4),
       at=c(-0.5, 0, 0.5, 1, 1.5, 2),
       #order = order(agg_data$year, agg_data$author), # order by date
       header=c("Reference", "Cohen's d (95% CI)"),# Add headers on left and right side
       addfit=FALSE, #remove polygon
       top=1.5)
# add polygon for overall meta-analytic effect
abline(h=0)
addpoly(res_adjusted, row=-1, mlab="Pooled effect size (MREM)", cex=0.8)
```

![](es_conversion_meta_20221206_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

### Forest plot for unadjusted effect sizes

``` r
# Make new MH outcome short label, which is capitalized, does not have symptoms or diagnosis at the end, and no underscores
data$mh_outcome_short <- data$mh_outcome
data$mh_outcome_short <- capitalize(trimws(mgsub(data$mh_outcome_short, c("_symptoms", "_diagnosis", "_", "emotional"), 
                                      c(" ", " ", " ", "emotional problems"))))

# Add more details where outcome reads "internalising" (Schaefer & Gerin)
data$mh_outcome_short[data$author=="Schaefer" & data$mh_outcome_short=="internalising"] <- "internalising factor"
data$mh_outcome_short[data$author=="Gerin_2017" & data$mh_outcome_short=="internalising"] <- "internalising symptoms"

# Capitalize and remove underscores from child maltreatment label
data$maltreatment_type_short <- capitalize(gsub("_", " ", data$maltreatment_type))

# Shorten cohort label "Chinese longitudinal study" to be shorter
data$cohort_short_distinct[data$cohort_short_distinct=="Chinese longitudinal study"] <- "Chinese long. study"

# Make forest plot (save 900 x 1200)
par(mar=c(4,2,1,2))
forest(data$d_unadj[!is.na(data$d_unadj)],
       vi=data$var_unadj[!is.na(data$d_unadj)],
       header=c("Reference", "Cohen's d (95% CI)"),# Add headers on left and right side
       cex=0.5,
       xlab = "Cohen's d", cex.axis=0.5,
       slab=NA, # Do not show references but add this manually below
       ilab=cbind(gsub("_", " et al., ", data$ref[!is.na(data$d_unadj)]), 
                  data$cohort_short_distinct[!is.na(data$d_unadj)], 
                  data$maltreatment_type_short[!is.na(data$d_unadj)], 
                  data$mh_outcome_short[!is.na(data$d_unadj)]), # Add info on reference, cohort, mt type and mh type
       ilab.xpos=c(-2.55, -1.9, -1.4, -0.95), # x-axis position of information on ref, cohort, mt and mh type
       ilab.pos = c(4, 4, 4, 4),
       xlim=c(-2.56,2.56),
       alim=c(-0.4,1.53), at=c(-0.4,  0, 0.4, 0.8, 1.2, 1.5))
par(font=2)
text(c(-1.75, -1.15, -0.65), 105, c("Cohort name", 
                                "Maltreatment type",
                                "Mental health outcome"), cex=0.5)
abline(h=0)
addpoly(res_unadjusted, row=-1, mlab="Pooled effect size (MREM)", cex=0.5)
```

![](es_conversion_meta_20221206_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

## Sensitivity analyses

Below we will run sensitivity analyses to test for publication bias,
undue influence of individual studies, and sensitivity to the definition
of childhood maltreatment. The results text is shown below the code.

### Publication bias analyses

To test for publication bias, we use the following analyses:

-   Egger’s test
-   Leave-one-out analysis to test undue influence of individual studies
    on *publication bias*
-   P-curve analysis

#### Egger’s test

``` r
test.egger <- rma.mv(d_adj, var_adj,  mod = var_adj, 
                     random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=data)
```

#### Leave-one-out analyses to test publication bias when each cohort/study/effect size is removed

``` r
#--------------------------------------------------------------------------------------------------------#
#--------------------- Exclude each COHORT in turn to test publication bias -----------------------------#
#--------------------------------------------------------------------------------------------------------#
groupids <- unique(data$cohort) # generate groupid variable indexing each unique cohort name

leave1out_cohort_bias <- matrix(rep(NA, length(groupids)), ncol=5, nrow=length(groupids)) # generate matrix to hold results

for(i in 1:length(groupids)) {  
  dataexcl <- subset(data, cohort!=groupids[i]) # remove each cohort in turn 
  N_rows <- psych::describe(dataexcl$n_adj)$n # extract the number of rows when each cohort is removed
  mean_N <- psych::describe(dataexcl$n_adj)$mean # extract the mean sample size when each cohort is removed
  meta <- rma.mv(d_adj, var_adj,  mod = var_adj, 
                 random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=dataexcl) # run Egger's test 
  # Put results into matrix
  leave1out_cohort_bias[i,1] <- i
  leave1out_cohort_bias[i,2] <- N_rows 
  leave1out_cohort_bias[i,3] <- mean_N 
  leave1out_cohort_bias[i,4] <- meta$QMp
  leave1out_cohort_bias[i,5] <- meta$QM
  }

leave1out_cohort_bias
leave1out_cohort_bias[,1] <- as.vector(groupids) # Name column 1 as the cohort that was removed
leave1out_cohort_bias_df <- as.data.frame(leave1out_cohort_bias) # Convert to dataframe
names(leave1out_cohort_bias_df) <-c('ref', 'N_rows', 'Mean_N', 'p_value', 'Q_mod') # Rename columns
# Convert columns from character to numeric
leave1out_cohort_bias_df <- leave1out_cohort_bias_df %>% mutate_at(c('N_rows', 'Mean_N', 'p_value', 'Q_mod'), as.numeric)
# Check for non-significant p-value (indicating removing the cohort removed apparent publication bias)
leave1out_cohort_bias_df[leave1out_cohort_bias_df$p_value>0.05,]
# Save dataset
setwd(Tables)
saveRDS(leave1out_cohort_bias_df, file="leave1out_cohort_bias.Rda")

#--------------------------------------------------------------------------------------------------------#
#--------------------- Exclude each PAPER in turn to test publication bias ------------------------------#
#--------------------------------------------------------------------------------------------------------#

groupids <- unique(data$ref)  # generate groupid variable indexing each unique paper name
leave1out_paper_bias <- matrix(rep(NA, length(groupids)), ncol=5, nrow=length(groupids))  # generate matrix to hold results

for(i in 1:length(groupids)) {  
  dataexcl <- subset(data, ref!=groupids[i])  # remove each paper in turn 
  N_rows <- psych::describe(dataexcl$n_adj)$n # extract the number of rows when each paper is removed
  mean_N <- psych::describe(dataexcl$n_adj)$mean # extract the mean sample size when each paper is removed
  meta <- rma.mv(d_adj, var_adj,  mod = var_adj, 
                 random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=dataexcl) # run Egger's test 
  # Put results into matrix
  leave1out_paper_bias[i,1] <- i
  leave1out_paper_bias[i,2] <- N_rows 
  leave1out_paper_bias[i,3] <- mean_N 
  leave1out_paper_bias[i,4] <- meta$QMp
  leave1out_paper_bias[i,5] <- meta$QM
}

leave1out_paper_bias
leave1out_paper_bias[,1] <- as.vector(groupids) # Name column 1 as the paper that was removed
leave1out_paper_bias_df <- as.data.frame(leave1out_paper_bias) # Convert to dataframe
names(leave1out_paper_bias_df) <-c('ref', 'N_rows', 'Mean_N', 'p_value', 'Q_mod') # Rename columns
# Convert columns from character to numeric
leave1out_paper_bias_df <- leave1out_paper_bias_df %>% mutate_at(c('N_rows', 'Mean_N', 'p_value', 'Q_mod'), as.numeric)
# Check for non-significant p-value (indicating removing the cohort removed apparent publication bias)
leave1out_paper_bias_df[leave1out_paper_bias_df$p_value>0.05,]
# Save dataset
setwd(Tables)
saveRDS(leave1out_paper_bias_df, file="leave1out_paper_bias.Rda")

#--------------------------------------------------------------------------------------------------------#
#--------------------- Exclude each EFFECT SIZE in turn to test publication bias ------------------------#
#--------------------------------------------------------------------------------------------------------#

leave1out_es_bias <- matrix(NA, nrow(data), ncol=5) # generate matrix to hold results

for ( i in 1:nrow(data) ) {
  
  dataexcl <- data[-i,] # remove each effect size in turn 
  N_rows <- psych::describe(dataexcl$n_adj)$n # extract the number of rows when each effect size is removed
  mean_N <- psych::describe(dataexcl$n_adj)$mean # extract the mean sample size when each effect size is removed
  meta <- rma.mv(d_adj, var_adj,  mod = var_adj, 
                 random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=dataexcl) # run Egger's test
 # Put results into matrix
  leave1out_es_bias[i,1] <- i
  leave1out_es_bias[i,2] <- N_rows 
  leave1out_es_bias[i,3] <- mean_N 
  leave1out_es_bias[i,4] <- meta$QMp
  leave1out_es_bias[i,5] <- meta$QM
}

leave1out_es_bias
leave1out_es_bias[,1] <- paste0(data$ref, "_", data$nesting_var)  # Name study and ES excluded
leave1out_es_bias_df <- as.data.frame(leave1out_es_bias) # Convert to dataframe
names(leave1out_es_bias_df) <-c('ref', 'N_rows', 'Mean_N', 'p_value', 'Q_mod') # Rename columns
# Convert columns from character to numeric
leave1out_es_bias_df <- leave1out_es_bias_df %>% mutate_at(c('N_rows', 'Mean_N', 'p_value', 'Q_mod'), as.numeric)
# Check for non-significant p-value (indicating removing the cohort removed apparent publication bias)
leave1out_es_bias_df[leave1out_es_bias_df$p_value>0.05,]
# Save dataset
setwd(Tables)
saveRDS(leave1out_es_bias_df, file="leave1out_es_bias.Rda")
```

#### P-curve analysis (eFigure 1)

``` r
agg_res_mg <- metagen(TE=ES, seTE=se, data=agg_data,
                  studlab = paste(ref),comb.fixed = FALSE,
                  comb.random = TRUE, method.tau = "REML",
                  hakn = FALSE, prediction = TRUE, sm = "SMD") # run meta-analysis on aggregated data to generate "meta" object needed for pcurve
pcurve(agg_res_mg) # run pcurve analysis, save 8x5.5
```

![](es_conversion_meta_20221206_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

    ## P-curve analysis 
    ##  ----------------------- 
    ## - Total number of provided studies: k = 34 
    ## - Total number of p<0.05 studies included into the analysis: k = 27 (79.41%) 
    ## - Total number of studies with p<0.025: k = 23 (67.65%) 
    ##    
    ## Results 
    ##  ----------------------- 
    ##                     pBinomial   zFull pFull   zHalf pHalf
    ## Right-skewness test     0.000 -10.759     0 -10.489     0
    ## Flatness test           0.971   6.835     1  10.449     1
    ## Note: p-values of 0 or 1 correspond to p<0.001 and p>0.999, respectively.   
    ## Power Estimate: 93% (86.2%-97%)
    ##    
    ## Evidential value 
    ##  ----------------------- 
    ## - Evidential value present: yes 
    ## - Evidential value absent/inadequate: no

### Undue influence of individual cohorts/studies/effect sizes

To test for the undue influence of individual estimates on the
meta-analytic estimate, we will run leave-one-out analyses excluding in
turn:

-   each cohort
-   each study (i.e., paper)
-   each effect size

``` r
#--------------------------------------------------------------------------------------------------------#
#------------------ Exclude each COHORT in turn to test overall meta-analytic effect --------------------#
#--------------------------------------------------------------------------------------------------------#
groupids <- unique(data$cohort)  # generate groupid variable indexing each unique cohort name
leave1out_cohort <- matrix(rep(NA,length(groupids)), ncol=6, nrow=length(groupids))  # generate matrix to hold results
for(i in 1:length(groupids)) {  
  dataexcl <- subset(data, cohort!=groupids[i])  # remove each cohort in turn 
  N_rows <- psych::describe(dataexcl$n_adj)$n # extract the number of rows when each cohort is removed
  mean_N <- psych::describe(dataexcl$n_adj)$mean # extract the mean sample size when each cohort is removed
  meta <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=dataexcl) # run meta
  cluster_ci_leave1out_cohort <- metafor::robust(meta, cluster=dataexcl$ref) # robust standard errors
  # Put results into matrix
  leave1out_cohort[i,1] <- i # cohort number excluded
  leave1out_cohort[i,2] <- N_rows # number of rows
  leave1out_cohort[i,3] <- mean_N # mean sample size
  leave1out_cohort[i,4] <- meta$beta[1] # meta-analytic estimate
  leave1out_cohort[i,5] <- cluster_ci_leave1out_cohort$ci.lb[1] # meta-analytic lower CI
  leave1out_cohort[i,6] <- cluster_ci_leave1out_cohort$ci.ub[1] # meta-analytic upper CI
}

leave1out_cohort
leave1out_cohort[,1] <- as.vector(groupids) # Name studies
leave1out_cohort_df <- as.data.frame(leave1out_cohort) # Convert to dataframe
names(leave1out_cohort_df) <-c('cohort', 'N_rows', 'Mean_N', 'es', 'ci.lb', 'ci.ub') # Rename columns
# Convert columns from character to numeric
leave1out_cohort_df <- leave1out_cohort_df %>% mutate_at(c('N_rows', 'Mean_N', 'es', 'ci.lb', 'ci.ub'), as.numeric)
# Shorten long cohort names for forest plot
leave1out_cohort_df$cohort[leave1out_cohort_df$cohort=="Cross-sectional study of patients with psychosis and unaffected siblings from Granada and Jaen"] <- data$cohort_short[data$cohort=="Cross-sectional study of patients with psychosis and unaffected siblings from Granada and Jaen"]
# Save dataset
setwd(Tables)
saveRDS(leave1out_cohort_df, file="leave1out_cohort.Rda")

#--------------------------------------------------------------------------------------------------------#
#------------------ Exclude each PAPER in turn to test overall meta-analytic effect ---------------------#
#--------------------------------------------------------------------------------------------------------#
groupids <- unique(data$ref)  # generate groupid variable indexing each unique paper name
leave1out_study <- matrix(rep(NA, length(groupids)), ncol=6, nrow=length(groupids)) # generate matrix to hold results 

for(i in 1:length(groupids)) {  
  dataexcl <- subset(data, ref!=groupids[i])  # remove each paper in turn 
  N_rows <- psych::describe(dataexcl$n_adj)$n  # extract the number of rows when each paper is removed
  mean_N <- psych::describe(dataexcl$n_adj)$mean # extract the mean sample size when each paper is removed
  meta <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=dataexcl) # run meta
  cluster_ci_leave1out_study <- metafor::robust(meta, cluster=dataexcl$ref) # robust standard errors
  # Put results into matrix
  leave1out_study[i,1] <- i
  leave1out_study[i,2] <- N_rows # number of rows
  leave1out_study[i,3] <- mean_N # mean sample size
  leave1out_study[i,4] <- meta$beta[1] # meta-analytic estimate
  leave1out_study[i,5] <- cluster_ci_leave1out_study$ci.lb[1] # meta-analytic lower CI
  leave1out_study[i,6] <- cluster_ci_leave1out_study$ci.ub[1] # meta-analytic upper CI
   
}

leave1out_study
leave1out_study[,1] <- as.vector(groupids) # Name studies
leave1out_study_df <- as.data.frame(leave1out_study) # Convert to dataframe
names(leave1out_study_df) <-c('ref', 'N_rows', 'Mean_N', 'es', 'ci.lb', 'ci.ub') # Rename columns
# Convert columns from character to numeric
leave1out_study_df <- leave1out_study_df %>% mutate_at(c('N_rows', 'Mean_N', 'es', 'ci.lb', 'ci.ub'), as.numeric)
# Save dataset
setwd(Tables)
saveRDS(leave1out_study_df, file="leave1out_study.Rda")

#--------------------------------------------------------------------------------------------------------#
#---------------- Exclude each EFFECT SIZE in turn to test overall meta-analytic effect -----------------#
#--------------------------------------------------------------------------------------------------------#

leave1out_es <- matrix(NA, nrow(data), ncol=6) # generate matrix to hold results 

for ( i in 1:nrow(data) ) {
  
  dataexcl <- data[-i,] # remove each ES in turn 
  N_rows <- psych::describe(dataexcl$n_adj)$n # extract the number of rows when each ES is removed
  mean_N <- psych::describe(dataexcl$n_adj)$mean # extract the mean sample size when each ES is removed
  meta <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=dataexcl) # run meta
  cluster_ci_leave1out_es <- metafor::robust(meta, cluster=dataexcl$ref) # robust standard errors
  # Put results into matrix
  leave1out_es[i,1] <- i
  leave1out_es[i,2] <- N_rows # number of rows
  leave1out_es[i,3] <- mean_N # mean sample size
  leave1out_es[i,4] <- meta$beta[1] # meta-analytic estimate
  leave1out_es[i,5] <- cluster_ci_leave1out_es$ci.lb[1] # meta-analytic lower CI
  leave1out_es[i,6] <- cluster_ci_leave1out_es$ci.ub[1] # meta-analytic upper CI
}

leave1out_es
leave1out_es[,1] <- paste0(data$ref, "_", data$nesting_var)  # Name study and ES excluded
leave1out_es_df <- as.data.frame(leave1out_es) # Convert to dataframe
names(leave1out_es_df) <-c('ref', 'N_rows', 'Mean_N', 'es', 'ci.lb', 'ci.ub') # Rename columns
# Convert columns from character to numeric
leave1out_es_df <- leave1out_es_df %>% mutate_at(c('N_rows', 'Mean_N', 'es', 'ci.lb', 'ci.ub'), as.numeric)
# Save dataset
setwd(Tables)
saveRDS(leave1out_es_df, file="leave1out_es.Rda")

#--------------------------------------------------------------------------------------------------------#
#----- Show smallest/largest meta-analytic effect after excluding each cohort, paper and effect size-----#
#--------------------------------------------------------------------------------------------------------#
# Cohort removal
format(round(max(leave1out_cohort_df$es), 2),nsmall=2)
format(round(min(leave1out_cohort_df$es), 2),nsmall=2)
# Paper removal
format(round(max(leave1out_study_df$es), 2),nsmall=2)
format(round(min(leave1out_study_df$es), 2),nsmall=2)
# Effect size removal
format(round(max(leave1out_es_df$es), 2),nsmall=2)
format(round(min(leave1out_es_df$es), 2),nsmall=2)
```

``` r
# Load leave-one-out results (influence on main meta-analysis)
setwd(Tables)
leave1out_es_df <- readRDS(file="leave1out_es.Rda")
leave1out_study_df <- readRDS(file="leave1out_study.Rda")
leave1out_cohort_df <-readRDS(file="leave1out_cohort.Rda")
# Load leave-one-out results (cause of publication bias)
leave1out_cohort_bias_df <- readRDS(file="leave1out_cohort_bias.Rda")
leave1out_paper_bias_df <- readRDS(file="leave1out_paper_bias.Rda")
leave1out_es_bias_df <- readRDS(file="leave1out_es_bias.Rda")
ERA_pub_bias <- leave1out_cohort_bias_df[leave1out_cohort_bias_df$p_value>0.05,]
```

#### Forest plot showing meta-analytic effect size after omitting each cohort in turn

``` r
# Leave one out cohort (eFigure 3) save 8x6
par(mar=c(4,4,1,2)) 
forest(leave1out_cohort_df$es, 
       xlim=c(-0.3,0.56),
       cex=0.75, 
       refline=NA,
       order="obs",
       ci.lb=leave1out_cohort_df$ci.lb, 
       ci.ub=leave1out_cohort_df$ci.ub, 
       slab=gsub("\\([^\\]]*\\)","",leave1out_cohort_df$cohort, perl = TRUE),
       xlab = "Cohen's D", 
       header=c("Cohort name/sample description", "Cohen's D (95% CI)"))
```

![](es_conversion_meta_20221206_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

#### Forest plot showing meta-analytic effect size after omitting each study in turn

``` r
### Leave one out study forest plot (eFigure 4)
par(mar=c(4,4,1,2)) # save 7x6
forest(leave1out_study_df$es, 
       ci.lb=leave1out_study_df$ci.lb, 
       ci.ub=leave1out_study_df$ci.ub, 
       order="obs",
       slab=gsub("_", " et al., ",leave1out_study_df$ref),
       xlab = "Cohen's d", 
       header=c("Reference", "Cohen's d (95% CI)"))
```

![](es_conversion_meta_20221206_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

### Definition of maltreatment sensitivity analysis

``` r
# Count number of studies assessing victimisation (to be omitted in sensitivity analysis)
k_victim <- data %>%   
  group_by(ref) %>%
  filter(maltreatment_type=="victimisation") %>%
  dplyr::summarise(m = max(cohort)) %>%  nrow() 

# Count number of studies assessing ACEs (to be omitted in sensitivity analysis)
k_ACEs <- data %>%   
  group_by(ref) %>%
  filter(maltreatment_type=="ACEs") %>%
  dplyr::summarise(m = max(cohort)) %>%  nrow() 

# Count number of studies assessing institutional deprivation (to be omitted in sensitivity analysis)
k_inst_neg <- data %>%   
  group_by(ref) %>%
  filter(maltreatment_type=="institutional_neglect") %>%
  dplyr::summarise(m = max(cohort)) %>%  nrow() 

# Run meta-analysis excluding studies assessing victimisation, ACEs, and institutional deprivation
adjusted_mal_pure <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), 
                            data=data, subset=maltreatment_type!="ACEs" & maltreatment_type!="victimisation" &
                              maltreatment_type!="institutional_deprivation")
cluster_ci_mal_pure <- metafor::robust(adjusted_mal_pure, cluster=data$ref)
```

An Egger’s test suggested evidence of small-study bias
(Q_moderation=8.52, p-value=0.004; see Figure S3A in the online
supplement for funnel plot). We investigated the cause of this potential
publication bias by performing leave-one-out analyses and found that it
was due to the inclusion of data from one cohort, the ERA Study, which
comprised a comparatively small sample (N=90-148) with large effect
sizes linked to severe institutional neglect (Q_moderation after
excluding ERA=2.92, p-value=0.09; see Figure S3B in the online
supplement for funnel plot). We therefore conducted later moderation
analyses both with and without effect sizes from the ERA Study, to
ensure that results were not biased. The p-curve analysis (focusing on
one averaged effect size per study) provided evidential value that a
true effect was present (see Figure S4 in the online supplement).

Leave-one-out analyses showed that the overall meta-analytic estimate
was not unduly influenced by individual studies, cohorts, or effect
sizes. The meta-analytic effect size ranged between Cohen’s d=0.29 to
0.32 (with overlapping confidence intervals) after omitting in turn each
of the 29 cohorts, 34 studies, and 156 effect sizes (see Figures S5 and
S6 in the online supplement).

The meta-analytic effect size was also similar to the original estimate
after excluding studies which assessed maltreatment as part of broader
measures of victimization (k=3) or ACEs (k=5), or focused on
institutional neglect (k=3) (Cohen’s d=0.34, 95% CI=0.26-0.42).

## Moderators of the association between child maltreatment and mental health

We next examine whether the association between child maltreatment and
mental health in quasi-experimental studies is moderated by a-priori
selected factors, including:

-   type of quasi-experimental method
-   type of mental health outcome
-   type of child maltreatment
-   prospective vs retrospective measures of child maltreatment
-   shared rater vs different raters
-   longitudinal vs cross-sectional design
-   sample characteristics
-   study quality

When testing each moderator, we remove categories that are only included
in a single study, and where relevant, test whether the results are
sensitive to removing the ERA cohort. We do not apply robust standard
errors given evidence that this can substantially reduce power in
moderation analyses, as reported
[here](https://onlinelibrary.wiley.com/doi/pdf/10.1002/jrsm.1245?casa_token=w4rZsOaWvL0AAAAA:kjb9NbybxeYp1KnWk7VUUMeP9WqqbE8UO9A-3527lnVSsy-Mg-rg0pcDAqCY-gTfv4EGVSu4PHJD).

Text reporting the results is shown below the relevant code.

### Type of quasi-experimental method

``` r
# Run moderator analysis by QE method (exclude CoT, adoption & random-intercept cross-lagged because k=1)
qe_mod_spec <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct),
       data=data, mods= ~ factor(qe_method), subset=qe_method!="CoT" & qe_method!="adoption" & qe_method!="random_intercept_cross_lagged") 

# Order variable
data$qe_method <- factor(data$qe_method, ordered=TRUE, levels=c("MZ_twin", "twin", "sibling", 
                                                                "inverse_probability_weighting", "fixed_effects",
                                                                "propensity_score", "natural_experiment",
                                                                "CoT", "adoption", "random_intercept_cross_lagged"))

# Re-run moderator analysis removing the intercept to inform forest plot for moderation effect
# note: each estimate is the estimated average effect size for all QE methods
# the significance levels show whether the estimates are different from zero rather than comparing to each-other
qetype_sep <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct),
       data=data, mods= ~ qe_method-1, subset=qe_method!="CoT" & qe_method!="adoption" & qe_method!="random_intercept_cross_lagged")    
# Add k studies to forest plot by getting a vector linking the study type with no. studies
k_qe_method <- data %>%
  filter(qe_method!="adoption" & qe_method!="CoT" & qe_method!="random_intercept_cross_lagged") %>%
  count(ref, qe_method) %>%
  group_by(qe_method) %>% count(qe_method)

# Add N effect sizes to forest plot by getting a vector linking the study type with no. effect sizes
es_qe_method <- data %>%
  filter(qe_method!="adoption" & qe_method!="CoT" & qe_method!="random_intercept_cross_lagged") %>%
  count(qe_method) 

# Add total sample size to forest plot by getting a vector linking the study type with max sample size in each study
n_qe_method <- data %>%
  filter(qe_method!="adoption" & qe_method!="CoT" & qe_method!="random_intercept_cross_lagged") %>%
  group_by(cohort_short_distinct) %>%  
  summarise(qe_method = max(qe_method), m=max(n_adj))%>% 
  dplyr::group_by(qe_method) %>%
  dplyr::summarise(n = sum(m))
```

As shown in Figure 3, the association between child maltreatment and
mental health was present across different quasi-experimental designs
(twin differences, sibling differences, fixed-effects, and natural
experiment) and analytic approaches (propensity score matching and
inverse probability weighting) and was not significantly moderated by
the type of quasi-experimental method used (Q_moderation=3.43;
p-value=0.75).

``` r
# Forest plot (save 900x550)
forest(coef(qetype_sep), sei=qetype_sep$se, 
       slab=paste0(capitalize(mgsub(names(coef(qetype_sep)), 
                                    c("qe_method","_", "sibling", "propensity_score", "MZ_twin", "twin"), 
                                    c(""," ", "Sibling difference", "Propensity score matching", 
                                      "MZ twin difference", "Twin difference"))), 
                   " (k=", k_qe_method$n, "; n=", trimws(format(n_qe_method$n, big.mark=",")), "; ES=", es_qe_method$n, ")"),
       xlab="Cohen's d", 
       header=c("Quasi-experimental method", "Cohen's d (95% CI)"),
       #order = "obs",
       xlim=c(-1, 1.2),
       alim=c(0, 0.82),
       at=c(0, 0.2, 0.4, 0.6, 0.8))
```

![](es_conversion_meta_20221206_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

### Mental health outcome

``` r
# Run moderation analysis by specific MH outcomes (including ERA)
mod_MH_spec <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), 
       data=data, mods = ~factor(mh_category), subset=mh_category!="bulimia" & mh_category!="substance_use_disorder")  

# Run moderation analysis by specific MH outcomes (excluding ERA)
mod_MH_spec_exclERA <- rma.mv(d_adj, var_adj, random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), 
                              data=data, mods = ~factor(mh_category), subset=cohort_short_distinct!="ERA" &
                                mh_category!="bulimia" & mh_category!="substance_use_disorder")  

# Run moderation analysis for broad internalising vs externalising problems
#mod_MH_broad <- rma.mv(d_adj, var_adj,  random=list( ~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), 
       #data=data, mods = ~factor(mh_category_broad), 
       #subset=mh_category_broad!="personality_disorder" & mh_category_broad!="psychosis" & ref!="Capusan_2021")  

# Prepare to make forest plot
# Order mental health outcomes in variable
data$mh_category <- factor(data$mh_category, ordered=TRUE, 
                            levels=c("depression", "anxiety", "internalising problems",
                                     "suicidal_ideation", "self_harm", "suicide_attempt",
                                     "internalising_problems", 
                                     "ADHD", "conduct_problems", "alcohol_abuse",
                                     "drug_abuse", "externalising_problems",
                                     "psychosis", "autism", 
                                     "personality_disorder", "psychopathology_broad"))

# Forest plot for moderation effect
# remove intercept so each estimate is the estimated average log risk ratio for the three levels
# the significance levels show whether the estimates are different from zero rather than comparing to eachother
mod_mh_forest <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct),
              data=data, mods= ~ mh_category-1)    

# Add k studies to forest plot by getting a vector linking the study type with no. studies
k_mh_type <- data %>%
  dplyr::filter(!is.na(mh_category)) %>%
  dplyr::count(cohort, mh_category) %>%
  dplyr::group_by(mh_category) %>% dplyr::count(mh_category)

# Add N effect sizes to forest plot by getting a vector linking the study type with no. effect sizes
es_mh_type <- data %>%
  dplyr::filter(!is.na(mh_category)) %>%
  dplyr::count(mh_category) 

# Add total sample size to forest plot by getting a vector linking the study type with no. sample size
n_mh_type <- data %>%  
  dplyr::filter(!is.na(mh_category)) %>%
  group_by(mh_category, cohort_short_distinct) %>% 
  dplyr::count(n_adj) %>%  
  dplyr::summarise(m = max(n_adj)) %>%  print(n=100)%>% 
  dplyr::summarise(n = sum(m)) 
```

The association between child maltreatment and mental health was
generally similar across different mental health problems (Figure 3).
Though stronger associations were found for autism symptoms, 4 of 5 of
these effect sizes were from the ERA Study, and there was no overall
moderation effect when this cohort was removed in a sensitivity analysis
(Q_moderation = 22.61; p-value=0.07).

#### Figure 4. Meta-analytic effect sizes by specific mental health outcome

``` r
#9.2x8 save
par(mar=c(4,2,1,2))
## Forest plot
forest(coef(mod_mh_forest), sei=mod_mh_forest$se, 
       slab=paste0(capitalize(mgsub(names(coef(mod_mh_forest)), c("mh_category","_"), 
                         c(""," "))), 
         " (k=", k_mh_type$n, "; n=", trimws(format(n_mh_type$n, big.mark=",")), "; ES=", es_mh_type$n, ")"),
       xlab="Cohen's d", header=c("Mental health outcome", "Cohen's d (95% CI)"), 
       rows=c(16:11, 8:4, 1:-2),
       ylim=c(-3,20),
       xlim=c(-1.1,1.5), 
       alim=c(-0.2, 1.02), 
       at=c(-2, 0, 0.2, 0.4, 0.6, 0.8, 1))
# Add labels distinguishing between subgroup results and the overall result
text(-1.1, c(2,9,17), pos=4, c("Other", "Externalising",
                                "Internalising"),  font=4)
```

![](es_conversion_meta_20221206_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

### Type of child maltreatment

``` r
# Code multiple types of abuse or neglect measures as "maltreatment", because they refer to any single types of maltreatment
# and the number of studies in each category is low (e.g., k=1 for neglect, k=2 for abuse)
data$maltreatment_type_r <- as.character(data$maltreatment_type)
data$maltreatment_type_r[data$maltreatment_type=="abuse"| data$maltreatment_type=="neglect"] <- "maltreatment" 
data$maltreatment_type_r <- factor(data$maltreatment_type_r)

# Check k studies per maltreatment type
data %>% dplyr::count(ref, maltreatment_type_r) %>%
  dplyr::group_by(maltreatment_type_r)%>% dplyr::count(maltreatment_type_r)

# Run analysis testing moderation by type of maltreatment
mod_mal_type <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct),
       data=data, mods= ~ factor(maltreatment_type_r)) 

# Run pairwise comparison to see which exposures are significantly different
mod_mal_type_spec <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct),
       data=data, mods= ~ maltreatment_type_r-1) 
pairwise_mt_types <- summary(glht(mod_mal_type_spec, linfct=cbind(contrMat(rep(1,9), type="Tukey"))), test=adjusted("none"))
pvals <- pairwise_mt_types$test$pvalues
pvals[pvals<=0.05] # Show significant pairwise differences
levels(data$maltreatment_type_r) # Show levels of maltreatment variable to work out which effects differ
# 1=ACEs; 2 = emotional abuse; 3=emotional neglect; 4=institutional deprivation; # 5=maltreatment; 6=physical abuse; 7=physical neglect; 8=sexual abuse
# emotional abuse and ACEs (2-1) differed
# institutional neglect and ACEs (4-1) differed
# emotional abuse and emotional neglect (3-2) differed
# emotional abuse and maltreatment (5-2) differed
# emotional abuse and physical abuse (6-2) differed
# emotional abuse and sexual abuse (8-2) differed

# Add k studies to forest plot by getting a vector linking the study type with no. studies
k_mal_type <- data %>%
  dplyr::count(ref, maltreatment_type_r) %>%
  dplyr::group_by(maltreatment_type_r)%>% count(maltreatment_type_r)

# Add N effect sizes to forest plot by getting a vector linking the study type with no. effect sizes
es_mal_type <-data %>% count(maltreatment_type_r) 

# Add total sample size to forest plot by getting a vector linking the study type with no. sample size
n_mal_type <- data %>%  
  dplyr::filter(!is.na(maltreatment_type_r)) %>%
  group_by(maltreatment_type_r, cohort_short_distinct) %>% 
  dplyr::count(n_adj) %>% 
  dplyr::summarise(max_N = max(n_adj)) %>% print(n=100) %>% 
  dplyr::summarise(n = sum(max_N)) 
```

The association between child maltreatment and mental health was
moderated by the type of childhood maltreatment (Q_moderation=22.43;
p-value=0.0042). As shown in Figure 4, emotional abuse and institutional
neglect were more strongly associated with mental health problems than
various subtypes of maltreatment and/or broader composite measures of
maltreatment and ACEs. Specifically, pairwise comparisons showed that
emotional abuse was more strongly associated with mental health than
physical abuse, sexual abuse, emotional neglect, physical neglect, and
broader measures of maltreatment and ACEs, while institutional neglect
was more strongly associated with mental health than ACEs. However, the
estimates for emotional abuse were only based only three studies and
seven effect sizes, and should be interpreted with caution.

#### Figure 5. Effect sizes for the pooled associations between different types of child maltreatment and mental health problems.

``` r
# save 755x600
par(mar=c(4, 4, 2, 2))
forest(coef(mod_mal_type_spec), sei=mod_mal_type_spec$se, 
       slab=paste0(capitalize(mgsub(names(coef(mod_mal_type_spec)), c("maltreatment_type_r","_"), 
                         c(""," "))), 
                   " (k=", k_mal_type$n, "; n=", trimws(format(n_mal_type$n, big.mark=",")), "; ES=", es_mal_type$n, ")"),
                   #" (k=", k_type_mal$n, "; ES=", es_type_mal$n, ")"),
       order = "obs",
       rows=c(1:9),
       ylim=c(0,11),
       xlab="Cohen's d", header=c("Exposure", "Cohen's d [95% CI]"), subset=order(-coef(mod_mal_type_spec)), 
       top=2,
       xlim=c(-1.25, 1.6)) 
```

![](es_conversion_meta_20221206_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

### Prospective vs retrospective measures of child maltreatment

``` r
# Run moderation analysis by prospective vs retrospective measures
mod_pro_retro <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct),
                                       data=data, mods= ~ factor(pro_retro_assess)) 
# Get results in prospective and retrospective studies separately by removing the intercept
mod_pro_retro_spec <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct),
                                       data=data, mods= ~ factor(pro_retro_assess)-1) 
## Excluding ERA
# Run moderation analysis by prospective vs retrospective measures
mod_pro_retro_exclERA <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct),
                                       data=data, mods= ~ factor(pro_retro_assess), subset=cohort_short_distinct!="ERA")
# Get results in prospective and retrospective studies separately by removing the intercept
mod_pro_retro_spec_exclERA <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct),
                                       data=data, mods= ~ factor(pro_retro_assess)-1, subset=cohort_short_distinct!="ERA")

# k studies with retrospective measures
k_retro <- data %>% group_by(ref) %>%
  filter(pro_retro_assess=="retrospective") %>%
  dplyr::summarise(m = max(ref)) %>%  nrow() 

# k studies with prospective measures including ERA
k_pro <- data %>% group_by(ref) %>%
  filter(pro_retro_assess=="prospective") %>%
  dplyr::summarise(m = max(ref)) %>%  nrow() 

# k studies with prospective measures excluding ERA
k_pro_exclERA <- data %>% group_by(ref) %>%
  filter(pro_retro_assess=="prospective" & cohort_short_distinct!="ERA") %>%
  dplyr::summarise(m = max(ref)) %>%  nrow() 

# Total sample size prospective and retrospective
n_pro_retro <- data %>%
  group_by(cohort_short_distinct) %>%  
  summarise(pro_retro_assess = max(pro_retro_assess), m=max(n_adj))%>% 
  dplyr::group_by(pro_retro_assess) %>%
  dplyr::summarise(n = sum(m))

# Sample sizes excluding ERA
n_pro_retro_exclERA <- data %>%
  group_by(cohort_short_distinct) %>%  
  summarise(pro_retro_assess = max(pro_retro_assess), m=max(n_adj)) %>% 
  dplyr::filter(cohort_short_distinct!="ERA") %>% 
  dplyr::group_by(pro_retro_assess) %>%
  dplyr::summarise(n = sum(m))
```

The association between child maltreatment and mental health was similar
in studies assessing maltreatment prospectively (Cohen’s d=0.29, 95%
CI=0.19-0.39; k=9, ES=47) versus retrospectively Cohen’s d=0.32, 95%
CI=0.24-0.39; k=9, ES=109), Q_moderation=0.253; p-value =0.61). This was
also the case after we conducted a sensitivity analysis excluding the
(prospective) ERA Study on institutional neglect (Cohen’s d for
prospective studies \[excluding ERA\]=0.25, 95% CI=0.15-0.36; k=6,
ES=20; Q_moderation=1.01, p-value=0.31).

### Shared rater

``` r
# Run moderation analysis by shared rater
mod_rater <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=data,
                    mod=~factor(data$q_informant_same))
# Get results in shared rater and non-shared rater studies separately by removing the intercept
mod_rater_spec <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=data,
                    mod=~factor(data$q_informant_same)-1)
### Excluding ERA
# Run moderation analysis by shared rater
mod_rater_exclERA <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=data,
                    mod=~factor(data$q_informant_same),  subset=cohort_short_distinct!="ERA")
# Get results in shared rater and non-shared rater studies separately by removing the intercept
mod_rater_spec_exclERA <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=data,
                    mod=~factor(data$q_informant_same)-1, subset=cohort_short_distinct!="ERA")

# k studies using the same rater
k_sameRater <- data %>% group_by(ref) %>%
  filter(q_informant_same=="0") %>%
  dplyr::summarise(m = max(ref)) %>%  nrow() 

# k studies using different raters including ERA
k_diffRater <- data %>% group_by(ref) %>%
  filter(q_informant_same=="1") %>%
  dplyr::summarise(m = max(ref)) %>%  nrow() 

# k studies using different raters excluding ERA
k_diffRater_exclERA <- data %>% group_by(ref) %>%
  filter(q_informant_same=="1" & cohort_short_distinct!="ERA") %>%
  dplyr::summarise(m = max(ref)) %>%  nrow() 

# n studies using the same v different rater 
n_rater <- data %>%
  group_by(cohort_short_distinct) %>%  
  summarise(q_informant_same = max(q_informant_same), m=max(n_adj))%>% 
  dplyr::group_by(q_informant_same) %>%
  dplyr::summarise(n = sum(m))
```

The meta-analytic association was present when child maltreatment and
mental health outcomes were assessed through self-reports (Cohen’s
d=0.32, 95% CI=0.25-0.39; k=26, ES=115) or through different sources
(Cohen’s d=0.28, 95% CI=0.17-0.39; k=9, ES=41; Q_moderation = 0.473,
p-value=0.49), although the latter effect size was attenuated after
excluding the ERA Study, which included different informants (Cohen’s
d=0.23, 95% CI=0.1-0.35; k=6, ES=14; Q_moderation=1.775, p-value=0.18).

### Longitudinal vs cross-sectional design

``` r
# Run moderation analysis by whether the study was longitudinal or cross-sectional
mod_longit <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=data,
                    mod=~factor(data$q_longitud))
# Get results in longitudinal and cross-sectional studies separately by removing the intercept
mod_longit_spec <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=data,
                    mod=~factor(data$q_longitud)-1)

## Excluding ERA
# Run moderation analysis by whether the study was longitudinal or cross-sectional
mod_longit_exclERA <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=data,
                    mod=~factor(data$q_longitud),  subset=cohort_short_distinct!="ERA")
# Get results in longitudinal and cross-sectional studies separately by removing the intercept
mod_longit_spec_exclERA <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=data,
                    mod=~factor(data$q_longitud)-1, subset=cohort_short_distinct!="ERA")

# Longitudinal studies that control for pre-existing mental health outcomes
long_control_mh <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=data,
                    subset=q_longitud==1 & q_mh_control==1)

# k longitudinal studies 
k_long <- data %>% group_by(ref) %>%
  filter(q_longitud=="1") %>%
  dplyr::summarise(m = max(ref)) %>%  nrow() 

# k cross-sectional
k_cross <- data %>% group_by(ref) %>%
  filter(q_longitud=="0") %>%
  dplyr::summarise(m = max(ref)) %>%  nrow() 

# k longitudinal studies excluding ERA
k_long_exclERA <- data %>% group_by(ref) %>%
  filter(q_longitud=="1" & cohort_short_distinct!="ERA") %>%
  dplyr::summarise(m = max(ref)) %>%  nrow() 

# k longitudinal studies with prior MH controlled
k_long_mh <- data %>% group_by(ref) %>%
  filter(q_longitud=="1" & q_mh_control=="1") %>%
  dplyr::summarise(m = max(ref)) %>%  nrow() 
```

The meta-analytic association between child maltreatment and mental
health was consistent when mental health outcomes were assessed
longitudinally (Cohen’s d=0.36, 95% CI=0.25-0.46; k=13, ES=53) or
cross-sectionally (Cohen’s d=0.28, 95% CI=0.19-0.36; k=21, ES=103;
Q_moderation=1.41, p-value=0.23. This was also the case after excluding
the longitudinal ERA Study (Cohen’s d for longitudinal studies=0.33, 95%
CI=0.21-0.44; k=10, ES=26; Q_moderation=0.45, p-value=0.50) and when
restricting to longitudinal studies controlling for pre-existing
psychopathology (Cohen’s d=0.34, 95% CI=0.16-0.52; k=8, ES=37).

``` r
### Make table with moderation results by assessment characteristics
Moderator <- c("Retrospective measures", 
               "Prospective measures", 
               "Prospective measures (excl. ERA)",
               "Same informant", 
               "Different informants", 
               "Different informants (excl. ERA)", 
               "Cross-sectional", 
               "Longitudinal",
               "Longitudinal (excl. ERA)", 
               "Longitudinal (controlling for earlier psychopathology)")
k <- c(k_retro,
       k_pro,
       k_pro_exclERA,
       k_sameRater,
       k_diffRater,
       k_diffRater_exclERA,
       k_cross,
       k_long,
       k_long_exclERA,
       k_long_mh)

ES <- c(table(data$pro_retro_assess)["retrospective"], 
        table(data$pro_retro_assess)["prospective"],
        table(data$pro_retro_assess[data$cohort_short_distinct!="ERA"])[1],
        table(data$q_informant_same)["0"], 
        table(data$q_informant_same)["1"], 
        table(data$q_informant_same[data$cohort_short_distinct!="ERA"])[2],
        table(data$q_longitud)["0"],
        table(data$q_longitud)["1"],
        table(data$q_longitud[data$cohort_short_distinct!="ERA"])[2],
        long_control_mh$k)

n <- c(n_pro_retro$n[n_pro_retro$pro_retro_assess=="retrospective"], # retro
       n_pro_retro$n[n_pro_retro$pro_retro_assess=="prospective"], # pro
       n_pro_retro_exclERA$n[n_pro_retro_exclERA$pro_retro_assess=="prospective"], # pro excluding ERA

# Define function to combine effect estimates and CIs into a single cell
paste_res <- function(est, lowCI, upCI) {
  combined <- paste0(format(round(est, 2), nsmall=2), " (", format(round(lowCI, 2), nsmall=2, trim=TRUE),"-", format(round(upCI, 2), nsmall=2), ")" )
  return(combined)
}

estimate_CI <- c(# Prospective versus retrospective
  paste_res(mod_pro_retro_spec$beta[2],mod_pro_retro_spec$ci.lb[2], mod_pro_retro_spec$ci.ub[2]), # retro
  paste_res(mod_pro_retro_spec$beta[1],mod_pro_retro_spec$ci.lb[1], mod_pro_retro_spec$ci.ub[1]), # pro
  paste_res(mod_pro_retro_spec_exclERA$beta[1],mod_pro_retro_spec_exclERA$ci.lb[1], mod_pro_retro_spec_exclERA$ci.ub[1]), #pro excluding ERA,
  # Shared rater vs different sources
  paste_res(mod_rater_spec$beta[1],mod_rater_spec$ci.lb[1], mod_rater_spec$ci.ub[1]), # self-reports
  paste_res(mod_rater_spec$beta[2],mod_rater_spec$ci.lb[2], mod_rater_spec$ci.ub[2]), # different sources
  paste_res(mod_rater_spec_exclERA$beta[2],mod_rater_spec_exclERA$ci.lb[2], mod_rater_spec_exclERA$ci.ub[2]),#diff exclERA
  # Longitudinal versus cross-sectional
  paste_res(mod_longit_spec$beta[1],mod_longit_spec$ci.lb[1], mod_longit_spec$ci.ub[1]), # Cross-sectional
  paste_res(mod_longit_spec$beta[2],mod_longit_spec$ci.lb[2], mod_longit_spec$ci.ub[2]), # Longitudinal
  paste_res(mod_longit_spec_exclERA$beta[2], mod_longit_spec_exclERA$ci.lb[2], mod_longit_spec_exclERA$ci.ub[2]), # Longitudinal excl ERA
  paste_res(long_control_mh$beta[1], long_control_mh$ci.lb[1], long_control_mh$ci.ub[1])) # Longitudinal early MH   

mod_table <- data.frame(Moderator, k, ES, estimate_CI)
kable(mod_table, col.names = c("Moderator", "k", "ES", "Cohen's d (95% CI)")) %>% 
  kable_styling(font_size = 11) %>%
  pack_rows("Prospective vs. retrospective measures of maltreatment", 1, 3, label_row_css = "color:black") %>%
  pack_rows("Shared rater vs. different sources", 4, 6, label_row_css = "color:black") %>%
  pack_rows("Longitudinal vs. cross-sectional assessment", 7, 10, label_row_css = "color:black")
```

### Sample characteristics

``` r
# Moderation by sex
mod_sex <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort), data=data,
                       mod=~data$perc_female)

# Moderation by age at outcome
mod_age <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort), data=data,
                  mod=~data$mh_age_assess_r)
```

The meta-analytic association between child maltreatment and mental
health was not moderated by sex (Q_moderation=2.917; p-value=0.09) or
age at mental health assessment (Q_moderation=0.02; p-value=0.87). We
did not examine moderation by race or ethnicity as the majority of
studies (24 of 35) did not report the race or ethnicity of the sample.

### Study quality

``` r
# Test moderation by study quality (total score)
res_q <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=data,
                  mod=~data$q_total)
# Exclude ERA study and repeat moderation analysis by quality
res_q_era <- rma.mv(d_adj, var_adj,  random=list(~ 1 | es_id, ~ 1 | ref, ~1 | cohort_short_distinct), data=data,
                mod=~data$q_total, subset=cohort_short_distinct!="ERA")
```

The meta-analytic association between child maltreatment and mental
health was not moderated by an overall measure of study quality
(Q_moderation=0.15; p-value=0.701; see eTable 6 in the Supplement for
coding).

## Supplementary tables

### Table with descriptives on studies (eTable 5)

Here we create a table with descriptive information about the included
studies. This is done by condensing the data on each study (often
included in multiple rows, if there are multiple effect sizes) into a
single row.

``` r
#--------------------------------------------------------------------------------------------------------#
#------------------- Prepare data to create table with descriptive information --------------------------#
#--------------------------------------------------------------------------------------------------------#

# Subset to data needed for the table
cohort_details <- subset(data, select=c(ref, cohort, country, perc_female, n_adj, n_unadj, qe_method, 
                                         maltreatment_type, maltreatment_age_assess, pro_retro_assess, 
                                         maltreatment_informant, maltreatment_measure, 
                                         mh_outcome, mh_measure, mh_informant, mh_age_assess))

# Format reference column to read "author et al. (year)"
cohort_details$ref <- paste0(sapply(strsplit(cohort_details$ref, "_"), `[`, 1), " et al. (", sapply(strsplit(cohort_details$ref, "_"), `[`, 2), ")")

# Combine cohort and country information into one column (the "cohort" column)
cohort_details$cohort <- paste0(cohort_details$cohort, ", ", cohort_details$country)

# Edit QE method labels
cohort_details$qe_method <- gsub("_", " ", cohort_details$qe_method)
cohort_details$qe_method <- capitalize(mgsub(cohort_details$qe_method, c("twin", "sibling", "CoT", "propensity score", "adoption"),
                                  c("twin differences", "sibling differences", 
                                    "children of twins", "propensity score matching", "adoption design")))

# Make mini table with maltreatment type information in a single cell for each study
mal_type <- cohort_details %>% group_by(ref, cohort, qe_method) %>% 
  summarise(maltreatment_type = paste(maltreatment_type, collapse = ";")) 
mal_type$maltreatment_type <- capitalize(gsub("_", " ", sapply(strsplit(mal_type$maltreatment_type, ";"), function(x) paste(unique(x), collapse = "; ")))) # edit text for maltreatment variable

# Make mini table with mental health information in a single cell for each study
mh_type <- cohort_details %>% group_by(ref, cohort, qe_method) %>%
  summarise(mh_outcome = paste(mh_outcome, collapse = ";")) 
mh_type$mh_outcome <- trimws(gsub(" ; ", ";", capitalize(mgsub(sapply(strsplit(mh_type$mh_outcome, ";"), 
                                                                function(x) paste(unique(x), collapse = "; ")),
                                 c("_", "symptoms", "diagnosis", ";"), 
                                       c(" ", "", "", ";"))))) # edit text for MH variable

# Make mini table with maltreatment assessment in a single cell for each study
mal_measure <- cohort_details %>% group_by(ref, cohort, qe_method) %>%
  summarise(maltreatment_measure = paste(c(maltreatment_measure, maltreatment_informant), collapse = ";"))
mal_measure$maltreatment_measure <- sapply(strsplit(mal_measure$maltreatment_measure, ";"), 
                                           function(x) paste(unique(x), collapse = " ")) # report only once in a cell the type of measure and informant used
mal_measure$maltreatment_measure <- gsub("Cps", "CPS", capitalize(ifelse(grepl("record", mal_measure$maltreatment_measure), 
                                           paste0(word(mal_measure$maltreatment_measure, 2), " ",
                                                  word(mal_measure$maltreatment_measure, 1)),
                                           paste0(word(mal_measure$maltreatment_measure, 1), " (",
                                                  word(mal_measure$maltreatment_measure, 2), ")")))) # Edit label to put brackets around informant, capitalise first letter and make "Cps" read "CPS"

# Make mini table with maltreatment age of assessment in a single cell for each study
mal_age <- cohort_details %>% group_by(ref, cohort, qe_method) %>%
  summarise(mal_age = paste(c(maltreatment_age_assess, pro_retro_assess), collapse = ";"))
mal_age$mal_age <- sapply(strsplit(mal_age$mal_age, ";"), function(x) paste(unique(x), collapse = "; "))
# round up ages to 1 dp
mal_age$mal_age <- mgsub(mal_age$mal_age, 
c("33.799999999999997", "33.76", "32.65", "49.68", "18.22", "10.34, 12.16", "16.08, 23, 30", "29.79"), 
c("33.8",               "33.8", "32.7",   "49.7", "18.2",   "10.3, 12.2",   "16.1, 23, 30", "29.8"))
# Dinwiddie study: re-label because 2 ages were given (for males and females - take mean)
mal_age$mal_age[mal_age$mal_age=="44.8; 42.7; retrospective"] <- "44.1; retrospective"
# Obikane_2018 study: re-label because 2 ages were given (for males and females - take mean)
mal_age$mal_age[mal_age$mal_age=="36.6; 36.4; retrospective"] <- "36.5; retrospective"
# Put prospective or retrospective in brackets
mal_age$mal_age <- paste0(gsub("\\;*\\s*\\w*$", "", mal_age$mal_age), " (", word(mal_age$mal_age, -1), ")")
# Edit Stern study   
mal_age$mal_age[mal_age$mal_age=="5, 7, 10, 12; 18; prospective (retrospective)"] <- "5, 7, 10, 12 (prospective); 18 (retrospective)"

# Make mini table with mental health age of assessment in a single cell for each study
mh_age <- cohort_details %>%  group_by(ref, cohort, qe_method) %>%
  summarise(mh_age = paste(mh_age_assess, collapse = ";"))
mh_age$mh_age <- sapply(strsplit(mh_age$mh_age, ";"), function(x) paste(unique(x), collapse = "; "))
# round up ages to 1 dp
mh_age$mh_age <- mgsub(mh_age$mh_age, 
c("33.799999999999997", "24.92", "33.76", "49.68", "18.22", "10.34, 12.16", "25.79"), 
c("33.8",               "24.9",   "33.8",   "49.7", "18.2",   "10.3, 12.2",   "25.8"))
# Check where multiple ages reported; semi-colon means different time points were examined seperately; comma means were merged into one measure
# Dinwiddie = 2 ages are reported because of males and females; report mean age for overall sample
mh_age$mh_age[mh_age$mh_age=="44.8; 42.7"] <- "44.1"
# Kendler: 2 ages are reported: 37.76; 30, 35. 37.76 was used for all outcomes except bulimia; so for table report 37.76
mh_age$mh_age[mh_age$mh_age=="37.76; 30, 35"] <- "37.6"
# Obikane 2018: re-label because 2 ages were given (for males and females - take mean)
mh_age$mh_age[mh_age$mh_age=="36.6; 36.4"] <- "36.5"

# Make mini table with mental health assessment in a single cell for each study
mh_measure <- cohort_details %>% group_by(ref, cohort, qe_method) %>%
  summarise(mh_measure = paste(c(mh_measure, mh_informant), collapse = ";"))
mh_measure$mh_measure <- sapply(strsplit(mh_measure$mh_measure, ";"), function(x) paste(unique(x), collapse = " "))
mh_measure$mh_measure <- capitalize(paste0(word(mh_measure$mh_measure, 1), " (", gsub("_", " ", word(mh_measure$mh_measure, 2)), ")"))
# Edit for Thornberry as records were used to assess crime outcomes but interview was used for other outcomes
mh_measure$mh_measure[mh_measure$mh_measure=="Record (interview)" & mh_measure$ref=="Thornberry et al. (2010)"] <- "Interview (self) or crime record (for arrest)"
# Edit for Stern as parents/teacher used to assess ADHD up to 12; self used to assess ADHD at age 18
mh_measure$mh_measure[mh_measure$mh_measure=="Interview (mixed)" & mh_measure$ref=="Stern et al. (2018)"] <- "Interview (parent, teacher); interview (self)"

# Make mini table with sample size information in a single cell for each study
sample_size <- cohort_details %>% group_by(ref, cohort, qe_method) %>%
  summarise(mean_n_adj = paste(c(round(mean(n_adj),0), round(mean(n_unadj),0)), collapse = "; "))

## Combine mini tables into one overall table with details of the studies (1 row per study)
cohort_details_table <- list(sample_size, mal_type, mal_age, mal_measure, mh_type, mh_age, mh_measure) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("ref", "cohort", "qe_method")), .)

# Order studies according to 1) quasi-experimental method, and 2) prospective vs retrospective measures
cohort_details_table$qe_method <- factor(cohort_details_table$qe_method, ordered=TRUE, 
                                            levels = c("MZ twin differences", "Twin differences", "Sibling differences",
                                                       "Children of twins", "Adoption design", "Fixed effects", 
                                                       "Natural experiment", "Propensity score matching", 
                                                       "Inverse probability weighting"))
cohort_details_table$orderProRetro <- ifelse(grepl("prospective", cohort_details_table$mal_age), 0, 1) # make pro vs retro variable
cohort_details_table <- cohort_details_table[order(cohort_details_table$qe_method, cohort_details_table$orderProRetro),] # order studies

# Remove prospective and retrospective ordering variable
cohort_details_table <- subset(cohort_details_table, select=-orderProRetro)
```

``` r
#--------------------------------------------------------------------------------------------------------#
#------------------- Create table with descriptive information ------------------------------------------#
#--------------------------------------------------------------------------------------------------------#

options(kableExtra.auto_format = FALSE)
kable(cohort_details_table, col.names = c("Reference", "Cohort name/description, country", 
                                          "Quasi-experimental Method", 
                                          "N (QE-adjusted; unadjusted)", 
                                          "Maltreatment type",
                                          "Age at maltreatment assessment",
                                          "Maltreatment measure",
                                          "Mental health outcome(s)", 
                                          "Age at mental health assessment",
                                          "Mental health measure")) %>%  kable_styling(font_size = 13)  %>%
  column_spec(1:7,  color="black") %>%
  column_spec(9:10,  color="black") %>%
  column_spec(8,  width_min= "18em", color="black") %>%
  pack_rows("MZ twin differences", 1, 10, label_row_css = "color:black") %>%
  pack_rows("Twin differences", 11, 16, label_row_css = "color:black") %>%
  pack_rows("Sibling differences", 17, 19, label_row_css = "color:black") %>%
  pack_rows("Children of twins", 20, 20, label_row_css = "color:black") %>%
  pack_rows("Adoption design", 21, 21, label_row_css = "color:black") %>%
  pack_rows("Within-individual fixed-effects", 22, 23, label_row_css = "color:black") %>%
  pack_rows("Natural experiment", 24, 26, label_row_css = "color:black") %>%
  pack_rows("Propensity score matching", 27, 29, label_row_css = "color:black") %>%
  pack_rows("Inverse probability weighting", 30, 32, label_row_css = "color:black") 
```

<table class="table" style="font-size: 13px; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Reference
</th>
<th style="text-align:left;">
Cohort name/description, country
</th>
<th style="text-align:left;">
Quasi-experimental Method
</th>
<th style="text-align:left;">
N (QE-adjusted; unadjusted)
</th>
<th style="text-align:left;">
Maltreatment type
</th>
<th style="text-align:left;">
Age at maltreatment assessment
</th>
<th style="text-align:left;">
Maltreatment measure
</th>
<th style="text-align:left;">
Mental health outcome(s)
</th>
<th style="text-align:left;">
Age at mental health assessment
</th>
<th style="text-align:left;">
Mental health measure
</th>
</tr>
</thead>
<tbody>
<tr grouplength="10">
<td colspan="10" style="color:black">
<strong>MZ twin differences</strong>
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Dinkler et al. (2017)
</td>
<td style="text-align:left;color: black !important;">
Child and Adolescent Twin Study in Sweden (CATSS), Sweden
</td>
<td style="text-align:left;color: black !important;">
MZ twin differences
</td>
<td style="text-align:left;color: black !important;">
3568; 8166
</td>
<td style="text-align:left;color: black !important;">
Maltreatment
</td>
<td style="text-align:left;color: black !important;">
9 (prospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (parent)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
ADHD;autism
</td>
<td style="text-align:left;color: black !important;">
9
</td>
<td style="text-align:left;color: black !important;">
Interview (parent)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Stern et al. (2018)
</td>
<td style="text-align:left;color: black !important;">
E-Risk Longitudinal Twin Study (E-Risk), UK
</td>
<td style="text-align:left;color: black !important;">
MZ twin differences
</td>
<td style="text-align:left;color: black !important;">
1100; NA
</td>
<td style="text-align:left;color: black !important;">
Victimisation
</td>
<td style="text-align:left;color: black !important;">
5, 7, 10, 12 (prospective); 18 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (parent)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
ADHD
</td>
<td style="text-align:left;color: black !important;">
5, 7, 10, 12; 18
</td>
<td style="text-align:left;color: black !important;">
Interview (parent, teacher); interview (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Alemany et al. (2013)
</td>
<td style="text-align:left;color: black !important;">
Cross-sectional study of adult twins from Catalonia, Spain
</td>
<td style="text-align:left;color: black !important;">
MZ twin differences
</td>
<td style="text-align:left;color: black !important;">
170; 226
</td>
<td style="text-align:left;color: black !important;">
ACEs
</td>
<td style="text-align:left;color: black !important;">
33.8 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Positive psychotic experiences;negative psychotic experiences
</td>
<td style="text-align:left;color: black !important;">
33.8
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Baldwin et al. (2019)
</td>
<td style="text-align:left;color: black !important;">
E-Risk Longitudinal Twin Study (E-Risk), UK
</td>
<td style="text-align:left;color: black !important;">
MZ twin differences
</td>
<td style="text-align:left;color: black !important;">
1100; 2055
</td>
<td style="text-align:left;color: black !important;">
Victimisation
</td>
<td style="text-align:left;color: black !important;">
18 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Suicidal ideation; self harm; suicide attempt
</td>
<td style="text-align:left;color: black !important;">
18
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Bornovalova et al. (2013)
</td>
<td style="text-align:left;color: black !important;">
Minnesota Twin Family Study (MTFS), USA
</td>
<td style="text-align:left;color: black !important;">
MZ twin differences
</td>
<td style="text-align:left;color: black !important;">
1792; 2764
</td>
<td style="text-align:left;color: black !important;">
Abuse; emotional abuse; physical abuse; sexual abuse
</td>
<td style="text-align:left;color: black !important;">
20, 24, 29 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Borderline personality disorder
</td>
<td style="text-align:left;color: black !important;">
24.9
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Capusan et al. (2016)
</td>
<td style="text-align:left;color: black !important;">
The Study of Twin Adults: Genes and Environment (STAGE), Sweden
</td>
<td style="text-align:left;color: black !important;">
MZ twin differences
</td>
<td style="text-align:left;color: black !important;">
940; 17711
</td>
<td style="text-align:left;color: black !important;">
Maltreatment; emotional neglect; physical neglect; physical abuse;
sexual abuse; abuse; neglect
</td>
<td style="text-align:left;color: black !important;">
33.8 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
ADHD
</td>
<td style="text-align:left;color: black !important;">
33.8
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Lecei et al. (2019)
</td>
<td style="text-align:left;color: black !important;">
twinssCan Study, Belgium
</td>
<td style="text-align:left;color: black !important;">
MZ twin differences
</td>
<td style="text-align:left;color: black !important;">
266; 266
</td>
<td style="text-align:left;color: black !important;">
Maltreatment
</td>
<td style="text-align:left;color: black !important;">
18.2 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Total psychopathology;psychosis;anxiety;depression
</td>
<td style="text-align:left;color: black !important;">
18.2
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Magnusson et al. (2012)
</td>
<td style="text-align:left;color: black !important;">
The Study of Twin Adults: Genes and Environment (STAGE), Sweden
</td>
<td style="text-align:left;color: black !important;">
MZ twin differences
</td>
<td style="text-align:left;color: black !important;">
44; 13595
</td>
<td style="text-align:left;color: black !important;">
Emotional neglect; physical abuse; sexual abuse
</td>
<td style="text-align:left;color: black !important;">
33.5 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Alcohol dependence
</td>
<td style="text-align:left;color: black !important;">
33.5
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Schaefer et al. (2017)
</td>
<td style="text-align:left;color: black !important;">
E-Risk Longitudinal Twin Study (E-Risk), UK
</td>
<td style="text-align:left;color: black !important;">
MZ twin differences
</td>
<td style="text-align:left;color: black !important;">
1158; 2062
</td>
<td style="text-align:left;color: black !important;">
Victimisation
</td>
<td style="text-align:left;color: black !important;">
18 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
P-factor;internalising;externalising;thought disorder
</td>
<td style="text-align:left;color: black !important;">
18
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Young-Wolff et al. (2011)
</td>
<td style="text-align:left;color: black !important;">
Virginia Adult Twin Study of Psychiatric and Substance Use Disorders
(VATSPSUD), USA
</td>
<td style="text-align:left;color: black !important;">
MZ twin differences
</td>
<td style="text-align:left;color: black !important;">
174; 3527
</td>
<td style="text-align:left;color: black !important;">
Maltreatment
</td>
<td style="text-align:left;color: black !important;">
35 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Alcohol abuse/dependence
</td>
<td style="text-align:left;color: black !important;">
35
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
</tr>
<tr grouplength="6">
<td colspan="10" style="color:black">
<strong>Twin differences</strong>
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Berenz et al. (2013)
</td>
<td style="text-align:left;color: black !important;">
Norwegian Twin Registry (NTR), Norway
</td>
<td style="text-align:left;color: black !important;">
Twin differences
</td>
<td style="text-align:left;color: black !important;">
616; 2780
</td>
<td style="text-align:left;color: black !important;">
ACEs
</td>
<td style="text-align:left;color: black !important;">
28.2 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Paranoid personality disorder;schizoid personality disorder;schizotypal
personality disorder;histrionic personality disorder;narcissistic
personality disorder;borderline personality disorder;antisocial
personality disorder;avoidant personality disorder;obsessive compulsive
personality disorder;dependent personality disorder
</td>
<td style="text-align:left;color: black !important;">
28.2
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;color: black !important;">
Australian National Health and Medical Research Council (NH&MRC) Twin
Register, Australia
</td>
<td style="text-align:left;color: black !important;">
Twin differences
</td>
<td style="text-align:left;color: black !important;">
75; 3180
</td>
<td style="text-align:left;color: black !important;">
Sexual abuse
</td>
<td style="text-align:left;color: black !important;">
44.1 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Depression;suicidal ideation; suicide attempt; panic disorder;social
phobia;alcohol dependence;conduct disorder;psychopathology any
</td>
<td style="text-align:left;color: black !important;">
44.1
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Kendler et al. (2000)
</td>
<td style="text-align:left;color: black !important;">
Virginia Twin Registry (VTR), USA
</td>
<td style="text-align:left;color: black !important;">
Twin differences
</td>
<td style="text-align:left;color: black !important;">
133; 1403
</td>
<td style="text-align:left;color: black !important;">
Sexual abuse
</td>
<td style="text-align:left;color: black !important;">
32.7 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Depression;GAD;alcohol dependence;drug dependence;bulimia
</td>
<td style="text-align:left;color: black !important;">
37.6
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Nelson et al. (2002)
</td>
<td style="text-align:left;color: black !important;">
Australian Twin Register young adult cohort, Australia
</td>
<td style="text-align:left;color: black !important;">
Twin differences
</td>
<td style="text-align:left;color: black !important;">
64; NA
</td>
<td style="text-align:left;color: black !important;">
Sexual abuse
</td>
<td style="text-align:left;color: black !important;">
29.9 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Depression;suicide attempt; conduct disorder;alcohol dependence;social
anxiety
</td>
<td style="text-align:left;color: black !important;">
29.9
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Nelson et al. (2006)
</td>
<td style="text-align:left;color: black !important;">
Australian Twin Register young adult cohort, Australia
</td>
<td style="text-align:left;color: black !important;">
Twin differences
</td>
<td style="text-align:left;color: black !important;">
280; NA
</td>
<td style="text-align:left;color: black !important;">
Sexual abuse
</td>
<td style="text-align:left;color: black !important;">
29.9 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Cannabis dependence/abuse;opioids dependence/abuse;sedatives
dependence/abuse;stimulants dependence/abuse;cocaine
dependence/abuse;any illicit drug dependence/abuse;non-cannabis illicit
drug dependence/abuse
</td>
<td style="text-align:left;color: black !important;">
29.9
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Schwartz et al. (2019)
</td>
<td style="text-align:left;color: black !important;">
Midlife in the United States (MIDUS), USA
</td>
<td style="text-align:left;color: black !important;">
Twin differences
</td>
<td style="text-align:left;color: black !important;">
862; 862
</td>
<td style="text-align:left;color: black !important;">
ACEs
</td>
<td style="text-align:left;color: black !important;">
46, 50 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Depression;antisocial behaviour
</td>
<td style="text-align:left;color: black !important;">
50
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
</tr>
<tr grouplength="3">
<td colspan="10" style="color:black">
<strong>Sibling differences</strong>
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Capusan et al. (2021)
</td>
<td style="text-align:left;color: black !important;">
Official record study of participants from Östergötland, Sweden
</td>
<td style="text-align:left;color: black !important;">
Sibling differences
</td>
<td style="text-align:left;color: black !important;">
865; 3887
</td>
<td style="text-align:left;color: black !important;">
Maltreatment
</td>
<td style="text-align:left;color: black !important;">
0-18 (prospective)
</td>
<td style="text-align:left;color: black !important;">
Medical record
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Substance use disorder
</td>
<td style="text-align:left;color: black !important;">
29.5
</td>
<td style="text-align:left;color: black !important;">
Record (medical)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Barrigon et al. (2015)
</td>
<td style="text-align:left;color: black !important;">
Cross-sectional study of patients with psychosis and unaffected siblings
from Granada and Jaen, Spain
</td>
<td style="text-align:left;color: black !important;">
Sibling differences
</td>
<td style="text-align:left;color: black !important;">
98; NA
</td>
<td style="text-align:left;color: black !important;">
Maltreatment
</td>
<td style="text-align:left;color: black !important;">
31.7 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Psychosis
</td>
<td style="text-align:left;color: black !important;">
31.7
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Kullberg et al. (2020)
</td>
<td style="text-align:left;color: black !important;">
Netherlands Study of Depression and Anxiety (NESDA), Netherlands
</td>
<td style="text-align:left;color: black !important;">
Sibling differences
</td>
<td style="text-align:left;color: black !important;">
636; 636
</td>
<td style="text-align:left;color: black !important;">
Emotional abuse; physical abuse; sexual abuse
</td>
<td style="text-align:left;color: black !important;">
49.7 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Depression;anxiety
</td>
<td style="text-align:left;color: black !important;">
49.7
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
</tr>
<tr grouplength="1">
<td colspan="10" style="color:black">
<strong>Children of twins</strong>
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Schwartz et al. (2019)
</td>
<td style="text-align:left;color: black !important;">
Add Health, USA
</td>
<td style="text-align:left;color: black !important;">
Sibling differences
</td>
<td style="text-align:left;color: black !important;">
3112; 3112
</td>
<td style="text-align:left;color: black !important;">
ACEs
</td>
<td style="text-align:left;color: black !important;">
16.1, 23, 30 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Depression
</td>
<td style="text-align:left;color: black !important;">
30
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
</tr>
<tr grouplength="1">
<td colspan="10" style="color:black">
<strong>Adoption design</strong>
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Lynch et al. (2006)
</td>
<td style="text-align:left;color: black !important;">
Australian Twin Register children of twins, Australia
</td>
<td style="text-align:left;color: black !important;">
Children of twins
</td>
<td style="text-align:left;color: black !important;">
2502; 1926
</td>
<td style="text-align:left;color: black !important;">
Physical abuse
</td>
<td style="text-align:left;color: black !important;">
25.1 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Externalising behaviour; drug and alcohol use; internalising behaviour
</td>
<td style="text-align:left;color: black !important;">
25.1
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
</tr>
<tr grouplength="2">
<td colspan="10" style="color:black">
<strong>Within-individual fixed-effects</strong>
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Riggins-Caspers et al. (2003)
</td>
<td style="text-align:left;color: black !important;">
Cross-sectional study of adult adoptees from Iowa, USA
</td>
<td style="text-align:left;color: black !important;">
Adoption design
</td>
<td style="text-align:left;color: black !important;">
150; NA
</td>
<td style="text-align:left;color: black !important;">
Physical abuse
</td>
<td style="text-align:left;color: black !important;">
31.5 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Conduct disorder;oppositional behaviour
</td>
<td style="text-align:left;color: black !important;">
31.5
</td>
<td style="text-align:left;color: black !important;">
Interview (adoptive parent)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Isumi et al. (2021)
</td>
<td style="text-align:left;color: black !important;">
Adachi Child Health Impact of Living Difficulty (A-CHILD) , Japan
</td>
<td style="text-align:left;color: black !important;">
Fixed effects
</td>
<td style="text-align:left;color: black !important;">
2920; NA
</td>
<td style="text-align:left;color: black !important;">
Maltreatment
</td>
<td style="text-align:left;color: black !important;">
6.5, 7.5, 9.5 (prospective)
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (parent)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Behavioural difficulties
</td>
<td style="text-align:left;color: black !important;">
6.5, 7.5, 9.5
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (parent)
</td>
</tr>
<tr grouplength="3">
<td colspan="10" style="color:black">
<strong>Natural experiment</strong>
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Ma et al. (2018)
</td>
<td style="text-align:left;color: black !important;">
Fragile Families and Child Wellbeing Study (FFCWS), USA
</td>
<td style="text-align:left;color: black !important;">
Fixed effects
</td>
<td style="text-align:left;color: black !important;">
2472; NA
</td>
<td style="text-align:left;color: black !important;">
Physical abuse
</td>
<td style="text-align:left;color: black !important;">
3, 5 (prospective)
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (parent)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Aggressive behaviour
</td>
<td style="text-align:left;color: black !important;">
3,5
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (parent)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Voith et al. (2014)
</td>
<td style="text-align:left;color: black !important;">
National Survey of Child and Adolescent Well-Being (NSCAW-I), USA
</td>
<td style="text-align:left;color: black !important;">
Fixed effects
</td>
<td style="text-align:left;color: black !important;">
1022; NA
</td>
<td style="text-align:left;color: black !important;">
ACEs
</td>
<td style="text-align:left;color: black !important;">
10.3, 12.2 (prospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (mixed)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Trauma;depression
</td>
<td style="text-align:left;color: black !important;">
10.3, 12.2
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Beckett et al. (2002)
</td>
<td style="text-align:left;color: black !important;">
English and Romanian Adoptees Study (ERA), UK/Romania
</td>
<td style="text-align:left;color: black !important;">
Natural experiment
</td>
<td style="text-align:left;color: black !important;">
90; NA
</td>
<td style="text-align:left;color: black !important;">
Institutional neglect
</td>
<td style="text-align:left;color: black !important;">
0-3.6 (prospective)
</td>
<td style="text-align:left;color: black !important;">
Government record
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Self injury
</td>
<td style="text-align:left;color: black !important;">
6
</td>
<td style="text-align:left;color: black !important;">
Interview (parent)
</td>
</tr>
<tr grouplength="3">
<td colspan="10" style="color:black">
<strong>Propensity score matching</strong>
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Golm et al. (2020)
</td>
<td style="text-align:left;color: black !important;">
English and Romanian Adoptees Study (ERA), UK/Romania
</td>
<td style="text-align:left;color: black !important;">
Natural experiment
</td>
<td style="text-align:left;color: black !important;">
98; NA
</td>
<td style="text-align:left;color: black !important;">
Institutional neglect
</td>
<td style="text-align:left;color: black !important;">
0-3.6 (prospective)
</td>
<td style="text-align:left;color: black !important;">
Government record
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Depression;generalised anxiety
</td>
<td style="text-align:left;color: black !important;">
23.9
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (parent)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;color: black !important;">
English and Romanian Adoptees Study (ERA), UK/Romania
</td>
<td style="text-align:left;color: black !important;">
Natural experiment
</td>
<td style="text-align:left;color: black !important;">
148; NA
</td>
<td style="text-align:left;color: black !important;">
Institutional neglect
</td>
<td style="text-align:left;color: black !important;">
0-3.6 (prospective)
</td>
<td style="text-align:left;color: black !important;">
Government record
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
ASD;inattention overactivity;emotional;conduct problem
</td>
<td style="text-align:left;color: black !important;">
6; 11; 15; 24.1
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (parent)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Thornberry et al. (2010)
</td>
<td style="text-align:left;color: black !important;">
Rochester Youth Development Study (RYDS), USA
</td>
<td style="text-align:left;color: black !important;">
Propensity score matching
</td>
<td style="text-align:left;color: black !important;">
907; NA
</td>
<td style="text-align:left;color: black !important;">
Maltreatment
</td>
<td style="text-align:left;color: black !important;">
14, 15, 16, 17, 18 (prospective)
</td>
<td style="text-align:left;color: black !important;">
CPS record
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Arrest or incarceration; suicidal thoughts; depression
</td>
<td style="text-align:left;color: black !important;">
22.7
</td>
<td style="text-align:left;color: black !important;">
Interview (self) or crime record (for arrest)
</td>
</tr>
<tr grouplength="3">
<td colspan="10" style="color:black">
<strong>Inverse probability weighting</strong>
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Gerin et al. (2019)
</td>
<td style="text-align:left;color: black !important;">
Duke Neurogenetics Study (DNS), USA
</td>
<td style="text-align:left;color: black !important;">
Propensity score matching
</td>
<td style="text-align:left;color: black !important;">
196; NA
</td>
<td style="text-align:left;color: black !important;">
Maltreatment
</td>
<td style="text-align:left;color: black !important;">
19.5 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Internalising
</td>
<td style="text-align:left;color: black !important;">
20.5
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Zvara et al. (2017)
</td>
<td style="text-align:left;color: black !important;">
Family Life Project (FLP), USA
</td>
<td style="text-align:left;color: black !important;">
Propensity score matching
</td>
<td style="text-align:left;color: black !important;">
204; NA
</td>
<td style="text-align:left;color: black !important;">
Sexual abuse
</td>
<td style="text-align:left;color: black !important;">
29.8 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Postnatal depression
</td>
<td style="text-align:left;color: black !important;">
25.8
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;padding-left: 2em;" indentlevel="1">
Kugler et al. (2019)
</td>
<td style="text-align:left;color: black !important;">
Female Adolescent Development Study, USA
</td>
<td style="text-align:left;color: black !important;">
Inverse probability weighting
</td>
<td style="text-align:left;color: black !important;">
367; NA
</td>
<td style="text-align:left;color: black !important;">
Maltreatment
</td>
<td style="text-align:left;color: black !important;">
14, 15, 16, 17, 18 (prospective)
</td>
<td style="text-align:left;color: black !important;">
CPS record
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Drug use; depression
</td>
<td style="text-align:left;color: black !important;">
19
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;">
Alvanzo et al. (2020)
</td>
<td style="text-align:left;color: black !important;">
National Epidemiologic Survey on Alcohol and Related Conditions
(NESARC), USA
</td>
<td style="text-align:left;color: black !important;">
Inverse probability weighting
</td>
<td style="text-align:left;color: black !important;">
10396; 10396
</td>
<td style="text-align:left;color: black !important;">
ACEs
</td>
<td style="text-align:left;color: black !important;">
45.9; 46.5 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Severe alcohol problems; moderate alcohol problems
</td>
<td style="text-align:left;color: black !important;">
43.9; 44.5
</td>
<td style="text-align:left;color: black !important;">
Interview (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;">
Obikane et al. (2018)
</td>
<td style="text-align:left;color: black !important;">
The Japanese Study on Stratification, Health, Income, and Neighborhood
(J-SHINE), Japan
</td>
<td style="text-align:left;color: black !important;">
Inverse probability weighting
</td>
<td style="text-align:left;color: black !important;">
1896; 1896
</td>
<td style="text-align:left;color: black !important;">
Physical abuse; physical neglect; maltreatment
</td>
<td style="text-align:left;color: black !important;">
36.5 (retrospective)
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Suicidal ideation; suicidal plan; suicide attempt
</td>
<td style="text-align:left;color: black !important;">
36.5
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;">
Li et al. (2021)
</td>
<td style="text-align:left;color: black !important;">
Longitudinal study of students from schools in Guangdong, China, China
</td>
<td style="text-align:left;color: black !important;">
NA
</td>
<td style="text-align:left;color: black !important;">
3742; 3742
</td>
<td style="text-align:left;color: black !important;">
Emotional abuse
</td>
<td style="text-align:left;color: black !important;">
9.9; 10.4; 10.9; 11.4 (prospective)
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
<td style="text-align:left;min-width: 18em; color: black !important;">
Depression
</td>
<td style="text-align:left;color: black !important;">
10.4; 10.9; 11.4; 11.9
</td>
<td style="text-align:left;color: black !important;">
Questionnaire (self)
</td>
</tr>
</tbody>
</table>

### Table with study quality information (eTable 6)

Here we create a table reporting the results of the quality assessment
for each study. Note that occasionally the quality scores varied within
a single study (e.g., if some outcomes were assessed longitudinally and
others cross-sectionally). Here we report the maximum score for each
study, in a single row.

``` r
# Subset data to columns on reference and all quality indicators
quality_table <- subset(data, select=c(ref, q_exposed_rep,  q_control_selec,    q_exposure_assess,  q_mh_control,   q_env_confound, q_gen_confound, q_mh_assess,    q_informant_same,   q_longitud, q_attrition,    q_total))

# Double check if q_total matches the sum of the other columns (yes)
quality_table$q_total2 <- rowSums(quality_table[,c("q_exposed_rep", "q_control_selec",  "q_exposure_assess",    
                                                   "q_mh_control",  "q_env_confound",   "q_gen_confound",   
                                                   "q_mh_assess",   "q_informant_same", "q_longitud",   "q_attrition")])
identical(quality_table$q_total, quality_table$q_total2)
quality_table <- subset(quality_table, select=-c(q_total2))

# Edit reference variable
quality_table$ref <- paste0(sapply(strsplit(quality_table$ref, "_"), `[`, 1), " et al. (", sapply(strsplit(quality_table$ref, "_"), `[`, 2), ")")

# Condense so there is only one row per study representing the maximum score for that study, 
# and drop trailing zeros (e.g., 1.0) by converting variables to to character
quality_table <- quality_table %>% group_by(ref) %>% 
  summarise_at(c("q_exposed_rep","q_control_selec", "q_exposure_assess",    "q_mh_control", "q_env_confound",   "q_gen_confound",
                 "q_mh_assess", "q_informant_same", "q_longitud",   "q_attrition", "q_total"), max) %>%
  mutate_at(c("q_exposed_rep","q_control_selec",    "q_exposure_assess",    "q_mh_control",
                                               "q_env_confound",    "q_gen_confound", "q_mh_assess",    
                                               "q_informant_same",  "q_longitud",   "q_attrition", "q_total"), as.character)
```

``` r
# Make table                
kable(quality_table, col.names = c("Reference", "Represent. (exposed)", "Exposed & unexposed from same cohort",
                                   "Validated MT assessment",   "Control for pre-existing MH",  "Control for environ. confounding",
                                   "Control for genetic confounding",   "Validated MH assessment",  "Different informants for MT & MH",
                                   "Longitudinal assessment",   "Retention >70%",   "Total quality score"),
      row.names=FALSE) %>% 
  kable_styling(font_size = 9)  
```

<table class="table" style="font-size: 9px; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Reference
</th>
<th style="text-align:left;">
Represent. (exposed)
</th>
<th style="text-align:left;">
Exposed & unexposed from same cohort
</th>
<th style="text-align:left;">
Validated MT assessment
</th>
<th style="text-align:left;">
Control for pre-existing MH
</th>
<th style="text-align:left;">
Control for environ. confounding
</th>
<th style="text-align:left;">
Control for genetic confounding
</th>
<th style="text-align:left;">
Validated MH assessment
</th>
<th style="text-align:left;">
Different informants for MT & MH
</th>
<th style="text-align:left;">
Longitudinal assessment
</th>
<th style="text-align:left;">
Retention \>70%
</th>
<th style="text-align:left;">
Total quality score
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Alemany et al. (2013)
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
Alvanzo et al. (2020)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
4.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Baldwin et al. (2019)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
Barrigon et al. (2015)
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Beckett et al. (2002)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
Berenz et al. (2013)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
Bornovalova et al. (2013)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
6.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Capusan et al. (2016)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
6.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Capusan et al. (2021)
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
7.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinkler et al. (2017)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
6.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
4.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Gerin et al. (2019)
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
4.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Golm et al. (2020)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
Isumi et al. (2021)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
7.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Kendler et al. (2000)
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
4.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Kugler et al. (2019)
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
7.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Kullberg et al. (2020)
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Lecei et al. (2019)
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
Li et al. (2021)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
8.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Lynch et al. (2006)
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
4.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Ma et al. (2018)
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
8.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Magnusson et al. (2012)
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
Nelson et al. (2002)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
4.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Nelson et al. (2006)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
4.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Riggins-Caspers et al. (2003)
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
4.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Schaefer et al. (2017)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
7.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Schwartz et al. (2019)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
10
</td>
</tr>
<tr>
<td style="text-align:left;">
Stern et al. (2018)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
7.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Thornberry et al. (2010)
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
7.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Voith et al. (2014)
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
7.5
</td>
</tr>
<tr>
<td style="text-align:left;">
Young-Wolff et al. (2011)
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Zvara et al. (2017)
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
5
</td>
</tr>
</tbody>
</table>

### Table with effect sizes (eTable 7)

Here we create a table of effect sizes (in Cohen’s d format) and
confidence intervals for each study.

``` r
#--------------------------------------------------------------------------------------------------------#
#-------------------- Prepare data to create table with effect sizes ------------------------------------#
#--------------------------------------------------------------------------------------------------------#

# Subset data to columns on reference, author, year, cohort, maltreatment type, MH type, effect size and SE
es_table <- subset(data, select=c(ref, year, author, cohort_short_distinct, maltreatment_type_short, mh_outcome_short, d_adj, d_se_adj))

# Format reference column to read "author et al. (year)"
es_table$ref <- paste0(sapply(strsplit(es_table$ref, "_"), `[`, 1), " et al. (", sapply(strsplit(es_table$ref, "_"), `[`, 2), ")")

# Derive 95% CIs
es_table$lowci <- es_table$d_adj - 1.96*es_table$d_se_adj
es_table$upci <- es_table$d_adj + 1.96*es_table$d_se_adj
es_table$ci <- paste0(format(round(es_table$lowci, 3), nsmall=3, trim=TRUE), "-", format(round(es_table$upci, 3), nsmall=3, trim=TRUE))

# Order table and subset to remove unnecessary columns 
es_table <- es_table[order(es_table$year, es_table$author),]
es_table <- subset(es_table, select=-c(year, author, lowci, upci, d_se_adj))

#--------------------------------------------------------------------------------------------------------#
#-------------------------- Create table with effect sizes ----------------------------------------------#
#--------------------------------------------------------------------------------------------------------#

kable(es_table, col.names = c("Reference", "Cohort name/description", 
                                          "Maltreatment type",
                                          "Mental health outcome", 
                                          "Cohen's D",
                                          "CI"),
      row.names=FALSE, digits=3) %>% 
  kable_styling(font_size = 13)  
```

<table class="table" style="font-size: 13px; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Reference
</th>
<th style="text-align:left;">
Cohort name/description
</th>
<th style="text-align:left;">
Maltreatment type
</th>
<th style="text-align:left;">
Mental health outcome
</th>
<th style="text-align:right;">
Cohen’s D
</th>
<th style="text-align:left;">
CI
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;">
ATR (females)
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.197
</td>
<td style="text-align:left;">
-0.180-0.575
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;">
ATR (females)
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Suicidal ideation
</td>
<td style="text-align:right;">
0.252
</td>
<td style="text-align:left;">
-0.065-0.569
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;">
ATR (females)
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Suicide attempt
</td>
<td style="text-align:right;">
0.466
</td>
<td style="text-align:left;">
-0.281-1.213
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;">
ATR (females)
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Panic disorder
</td>
<td style="text-align:right;">
0.382
</td>
<td style="text-align:left;">
-0.211-0.975
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;">
ATR (females)
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Social phobia
</td>
<td style="text-align:right;">
0.224
</td>
<td style="text-align:left;">
-0.476-0.923
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;">
ATR (females)
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Alcohol dependence
</td>
<td style="text-align:right;">
0.505
</td>
<td style="text-align:left;">
-0.017-1.027
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;">
ATR (females)
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Conduct disorder
</td>
<td style="text-align:right;">
0.123
</td>
<td style="text-align:left;">
-0.599-0.845
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;">
ATR (females)
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Psychopathology any
</td>
<td style="text-align:right;">
0.242
</td>
<td style="text-align:left;">
-0.069-0.553
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;">
ATR (males)
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.224
</td>
<td style="text-align:left;">
-0.476-0.923
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;">
ATR (males)
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Suicidal ideation
</td>
<td style="text-align:right;">
0.940
</td>
<td style="text-align:left;">
0.110-1.770
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;">
ATR (males)
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Alcohol dependence
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:left;">
-0.885-0.885
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;">
ATR (males)
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Conduct disorder
</td>
<td style="text-align:right;">
0.382
</td>
<td style="text-align:left;">
-0.382-1.146
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinwiddie et al. (2000)
</td>
<td style="text-align:left;">
ATR (males)
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Psychopathology any
</td>
<td style="text-align:right;">
0.382
</td>
<td style="text-align:left;">
-0.281-1.045
</td>
</tr>
<tr>
<td style="text-align:left;">
Kendler et al. (2000)
</td>
<td style="text-align:left;">
VTR
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.186
</td>
<td style="text-align:left;">
-0.095-0.466
</td>
</tr>
<tr>
<td style="text-align:left;">
Kendler et al. (2000)
</td>
<td style="text-align:left;">
VTR
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
GAD
</td>
<td style="text-align:right;">
0.212
</td>
<td style="text-align:left;">
-0.144-0.569
</td>
</tr>
<tr>
<td style="text-align:left;">
Kendler et al. (2000)
</td>
<td style="text-align:left;">
VTR
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Alcohol dependence
</td>
<td style="text-align:right;">
0.574
</td>
<td style="text-align:left;">
0.066-1.081
</td>
</tr>
<tr>
<td style="text-align:left;">
Kendler et al. (2000)
</td>
<td style="text-align:left;">
VTR
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Drug dependence
</td>
<td style="text-align:right;">
0.382
</td>
<td style="text-align:left;">
-0.212-0.976
</td>
</tr>
<tr>
<td style="text-align:left;">
Kendler et al. (2000)
</td>
<td style="text-align:left;">
VTR
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Bulimia
</td>
<td style="text-align:right;">
0.157
</td>
<td style="text-align:left;">
-0.664-0.978
</td>
</tr>
<tr>
<td style="text-align:left;">
Beckett et al. (2002)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Self injury
</td>
<td style="text-align:right;">
1.219
</td>
<td style="text-align:left;">
0.084-2.353
</td>
</tr>
<tr>
<td style="text-align:left;">
Nelson et al. (2002)
</td>
<td style="text-align:left;">
ATR_YoungAdults
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.245
</td>
<td style="text-align:left;">
0.033-0.458
</td>
</tr>
<tr>
<td style="text-align:left;">
Nelson et al. (2002)
</td>
<td style="text-align:left;">
ATR_YoungAdults
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Suicide attempt
</td>
<td style="text-align:right;">
0.554
</td>
<td style="text-align:left;">
0.174-0.934
</td>
</tr>
<tr>
<td style="text-align:left;">
Nelson et al. (2002)
</td>
<td style="text-align:left;">
ATR_YoungAdults
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Conduct disorder
</td>
<td style="text-align:right;">
0.606
</td>
<td style="text-align:left;">
0.165-1.046
</td>
</tr>
<tr>
<td style="text-align:left;">
Nelson et al. (2002)
</td>
<td style="text-align:left;">
ATR_YoungAdults
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Alcohol dependence
</td>
<td style="text-align:right;">
0.245
</td>
<td style="text-align:left;">
0.007-0.484
</td>
</tr>
<tr>
<td style="text-align:left;">
Nelson et al. (2002)
</td>
<td style="text-align:left;">
ATR_YoungAdults
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Social anxiety
</td>
<td style="text-align:right;">
0.466
</td>
<td style="text-align:left;">
0.132-0.801
</td>
</tr>
<tr>
<td style="text-align:left;">
Riggins-Caspers et al. (2003)
</td>
<td style="text-align:left;">
Iowa adoption study
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Conduct disorder
</td>
<td style="text-align:right;">
0.387
</td>
<td style="text-align:left;">
0.059-0.715
</td>
</tr>
<tr>
<td style="text-align:left;">
Riggins-Caspers et al. (2003)
</td>
<td style="text-align:left;">
Iowa adoption study
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Oppositional behaviour
</td>
<td style="text-align:right;">
0.516
</td>
<td style="text-align:left;">
0.184-0.849
</td>
</tr>
<tr>
<td style="text-align:left;">
Lynch et al. (2006)
</td>
<td style="text-align:left;">
ATR_CoT
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Externalising behaviour
</td>
<td style="text-align:right;">
0.303
</td>
<td style="text-align:left;">
0.179-0.426
</td>
</tr>
<tr>
<td style="text-align:left;">
Lynch et al. (2006)
</td>
<td style="text-align:left;">
ATR_CoT
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Drug and alcohol use
</td>
<td style="text-align:right;">
0.343
</td>
<td style="text-align:left;">
0.220-0.467
</td>
</tr>
<tr>
<td style="text-align:left;">
Lynch et al. (2006)
</td>
<td style="text-align:left;">
ATR_CoT
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Internalising behaviour
</td>
<td style="text-align:right;">
0.165
</td>
<td style="text-align:left;">
0.042-0.289
</td>
</tr>
<tr>
<td style="text-align:left;">
Nelson et al. (2006)
</td>
<td style="text-align:left;">
ATR_YoungAdults
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Cannabis dependence/abuse
</td>
<td style="text-align:right;">
0.136
</td>
<td style="text-align:left;">
-0.152-0.424
</td>
</tr>
<tr>
<td style="text-align:left;">
Nelson et al. (2006)
</td>
<td style="text-align:left;">
ATR_YoungAdults
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Opioids dependence/abuse
</td>
<td style="text-align:right;">
1.032
</td>
<td style="text-align:left;">
0.212-1.852
</td>
</tr>
<tr>
<td style="text-align:left;">
Nelson et al. (2006)
</td>
<td style="text-align:left;">
ATR_YoungAdults
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Sedatives dependence/abuse
</td>
<td style="text-align:right;">
0.829
</td>
<td style="text-align:left;">
-0.016-1.674
</td>
</tr>
<tr>
<td style="text-align:left;">
Nelson et al. (2006)
</td>
<td style="text-align:left;">
ATR_YoungAdults
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Stimulants dependence/abuse
</td>
<td style="text-align:right;">
0.302
</td>
<td style="text-align:left;">
-0.108-0.712
</td>
</tr>
<tr>
<td style="text-align:left;">
Nelson et al. (2006)
</td>
<td style="text-align:left;">
ATR_YoungAdults
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Cocaine dependence/abuse
</td>
<td style="text-align:right;">
0.382
</td>
<td style="text-align:left;">
-0.382-1.146
</td>
</tr>
<tr>
<td style="text-align:left;">
Nelson et al. (2006)
</td>
<td style="text-align:left;">
ATR_YoungAdults
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Any illicit drug dependence/abuse
</td>
<td style="text-align:right;">
0.318
</td>
<td style="text-align:left;">
0.036-0.599
</td>
</tr>
<tr>
<td style="text-align:left;">
Nelson et al. (2006)
</td>
<td style="text-align:left;">
ATR_YoungAdults
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Non-cannabis illicit drug dependence/abuse
</td>
<td style="text-align:right;">
0.430
</td>
<td style="text-align:left;">
0.037-0.823
</td>
</tr>
<tr>
<td style="text-align:left;">
Thornberry et al. (2010)
</td>
<td style="text-align:left;">
RYDS
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Arrest or incarceration
</td>
<td style="text-align:right;">
0.243
</td>
<td style="text-align:left;">
-0.006-0.491
</td>
</tr>
<tr>
<td style="text-align:left;">
Thornberry et al. (2010)
</td>
<td style="text-align:left;">
RYDS
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Suicidal thoughts
</td>
<td style="text-align:right;">
0.369
</td>
<td style="text-align:left;">
0.067-0.672
</td>
</tr>
<tr>
<td style="text-align:left;">
Thornberry et al. (2010)
</td>
<td style="text-align:left;">
RYDS
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.209
</td>
<td style="text-align:left;">
0.004-0.413
</td>
</tr>
<tr>
<td style="text-align:left;">
Thornberry et al. (2010)
</td>
<td style="text-align:left;">
RYDS
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Arrest or incarceration
</td>
<td style="text-align:right;">
0.408
</td>
<td style="text-align:left;">
0.105-0.711
</td>
</tr>
<tr>
<td style="text-align:left;">
Thornberry et al. (2010)
</td>
<td style="text-align:left;">
RYDS
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Suicidal thoughts
</td>
<td style="text-align:right;">
0.496
</td>
<td style="text-align:left;">
0.150-0.842
</td>
</tr>
<tr>
<td style="text-align:left;">
Thornberry et al. (2010)
</td>
<td style="text-align:left;">
RYDS
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.197
</td>
<td style="text-align:left;">
-0.044-0.438
</td>
</tr>
<tr>
<td style="text-align:left;">
Young-Wolff et al. (2011)
</td>
<td style="text-align:left;">
VATSPSUD
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Alcohol abuse/dependence
</td>
<td style="text-align:right;">
-0.046
</td>
<td style="text-align:left;">
-0.491-0.399
</td>
</tr>
<tr>
<td style="text-align:left;">
Magnusson et al. (2012)
</td>
<td style="text-align:left;">
STAGE
</td>
<td style="text-align:left;">
Emotional neglect
</td>
<td style="text-align:left;">
Alcohol dependence
</td>
<td style="text-align:right;">
0.032
</td>
<td style="text-align:left;">
-0.331-0.395
</td>
</tr>
<tr>
<td style="text-align:left;">
Magnusson et al. (2012)
</td>
<td style="text-align:left;">
STAGE
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Alcohol dependence
</td>
<td style="text-align:right;">
0.201
</td>
<td style="text-align:left;">
-0.266-0.669
</td>
</tr>
<tr>
<td style="text-align:left;">
Magnusson et al. (2012)
</td>
<td style="text-align:left;">
STAGE
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Alcohol dependence
</td>
<td style="text-align:right;">
0.466
</td>
<td style="text-align:left;">
-0.060-0.993
</td>
</tr>
<tr>
<td style="text-align:left;">
Alemany et al. (2013)
</td>
<td style="text-align:left;">
Catalonia twin study
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Positive psychotic experiences
</td>
<td style="text-align:right;">
0.381
</td>
<td style="text-align:left;">
0.078-0.685
</td>
</tr>
<tr>
<td style="text-align:left;">
Alemany et al. (2013)
</td>
<td style="text-align:left;">
Catalonia twin study
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Negative psychotic experiences
</td>
<td style="text-align:right;">
0.390
</td>
<td style="text-align:left;">
0.086-0.693
</td>
</tr>
<tr>
<td style="text-align:left;">
Berenz et al. (2013)
</td>
<td style="text-align:left;">
NTR
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Paranoid personality disorder
</td>
<td style="text-align:right;">
0.063
</td>
<td style="text-align:left;">
-0.095-0.221
</td>
</tr>
<tr>
<td style="text-align:left;">
Berenz et al. (2013)
</td>
<td style="text-align:left;">
NTR
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Schizoid personality disorder
</td>
<td style="text-align:right;">
0.090
</td>
<td style="text-align:left;">
-0.069-0.248
</td>
</tr>
<tr>
<td style="text-align:left;">
Berenz et al. (2013)
</td>
<td style="text-align:left;">
NTR
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Schizotypal personality disorder
</td>
<td style="text-align:right;">
0.063
</td>
<td style="text-align:left;">
-0.095-0.221
</td>
</tr>
<tr>
<td style="text-align:left;">
Berenz et al. (2013)
</td>
<td style="text-align:left;">
NTR
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Histrionic personality disorder
</td>
<td style="text-align:right;">
0.063
</td>
<td style="text-align:left;">
-0.095-0.221
</td>
</tr>
<tr>
<td style="text-align:left;">
Berenz et al. (2013)
</td>
<td style="text-align:left;">
NTR
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Narcissistic personality disorder
</td>
<td style="text-align:right;">
0.127
</td>
<td style="text-align:left;">
-0.031-0.285
</td>
</tr>
<tr>
<td style="text-align:left;">
Berenz et al. (2013)
</td>
<td style="text-align:left;">
NTR
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Borderline personality disorder
</td>
<td style="text-align:right;">
0.155
</td>
<td style="text-align:left;">
-0.003-0.314
</td>
</tr>
<tr>
<td style="text-align:left;">
Berenz et al. (2013)
</td>
<td style="text-align:left;">
NTR
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Antisocial personality disorder
</td>
<td style="text-align:right;">
0.155
</td>
<td style="text-align:left;">
-0.003-0.314
</td>
</tr>
<tr>
<td style="text-align:left;">
Berenz et al. (2013)
</td>
<td style="text-align:left;">
NTR
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Avoidant personality disorder
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:left;">
-0.158-0.158
</td>
</tr>
<tr>
<td style="text-align:left;">
Berenz et al. (2013)
</td>
<td style="text-align:left;">
NTR
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Obsessive compulsive personality disorder
</td>
<td style="text-align:right;">
0.127
</td>
<td style="text-align:left;">
-0.031-0.285
</td>
</tr>
<tr>
<td style="text-align:left;">
Berenz et al. (2013)
</td>
<td style="text-align:left;">
NTR
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Dependent personality disorder
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:left;">
-0.158-0.158
</td>
</tr>
<tr>
<td style="text-align:left;">
Bornovalova et al. (2013)
</td>
<td style="text-align:left;">
MTFS
</td>
<td style="text-align:left;">
Abuse
</td>
<td style="text-align:left;">
Borderline personality disorder
</td>
<td style="text-align:right;">
0.090
</td>
<td style="text-align:left;">
-0.086-0.266
</td>
</tr>
<tr>
<td style="text-align:left;">
Bornovalova et al. (2013)
</td>
<td style="text-align:left;">
MTFS
</td>
<td style="text-align:left;">
Emotional abuse
</td>
<td style="text-align:left;">
Borderline personality disorder
</td>
<td style="text-align:right;">
0.190
</td>
<td style="text-align:left;">
-0.045-0.425
</td>
</tr>
<tr>
<td style="text-align:left;">
Bornovalova et al. (2013)
</td>
<td style="text-align:left;">
MTFS
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Borderline personality disorder
</td>
<td style="text-align:right;">
0.130
</td>
<td style="text-align:left;">
-0.066-0.326
</td>
</tr>
<tr>
<td style="text-align:left;">
Bornovalova et al. (2013)
</td>
<td style="text-align:left;">
MTFS
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Borderline personality disorder
</td>
<td style="text-align:right;">
-0.050
</td>
<td style="text-align:left;">
-0.344-0.244
</td>
</tr>
<tr>
<td style="text-align:left;">
Voith et al. (2014)
</td>
<td style="text-align:left;">
NSCAW-I
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Trauma
</td>
<td style="text-align:right;">
0.150
</td>
<td style="text-align:left;">
0.062-0.238
</td>
</tr>
<tr>
<td style="text-align:left;">
Voith et al. (2014)
</td>
<td style="text-align:left;">
NSCAW-I
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.118
</td>
<td style="text-align:left;">
0.030-0.206
</td>
</tr>
<tr>
<td style="text-align:left;">
Barrigon et al. (2015)
</td>
<td style="text-align:left;">
Spanish CS study
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Psychosis
</td>
<td style="text-align:right;">
1.096
</td>
<td style="text-align:left;">
0.033-2.159
</td>
</tr>
<tr>
<td style="text-align:left;">
Capusan et al. (2016)
</td>
<td style="text-align:left;">
STAGE
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
ADHD
</td>
<td style="text-align:right;">
0.180
</td>
<td style="text-align:left;">
0.105-0.255
</td>
</tr>
<tr>
<td style="text-align:left;">
Capusan et al. (2016)
</td>
<td style="text-align:left;">
STAGE
</td>
<td style="text-align:left;">
Emotional neglect
</td>
<td style="text-align:left;">
ADHD
</td>
<td style="text-align:right;">
0.190
</td>
<td style="text-align:left;">
0.115-0.265
</td>
</tr>
<tr>
<td style="text-align:left;">
Capusan et al. (2016)
</td>
<td style="text-align:left;">
STAGE
</td>
<td style="text-align:left;">
Physical neglect
</td>
<td style="text-align:left;">
ADHD
</td>
<td style="text-align:right;">
0.250
</td>
<td style="text-align:left;">
-0.040-0.540
</td>
</tr>
<tr>
<td style="text-align:left;">
Capusan et al. (2016)
</td>
<td style="text-align:left;">
STAGE
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
ADHD
</td>
<td style="text-align:right;">
0.080
</td>
<td style="text-align:left;">
-0.065-0.225
</td>
</tr>
<tr>
<td style="text-align:left;">
Capusan et al. (2016)
</td>
<td style="text-align:left;">
STAGE
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
ADHD
</td>
<td style="text-align:right;">
0.200
</td>
<td style="text-align:left;">
0.020-0.380
</td>
</tr>
<tr>
<td style="text-align:left;">
Capusan et al. (2016)
</td>
<td style="text-align:left;">
STAGE
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
ADHD
</td>
<td style="text-align:right;">
0.190
</td>
<td style="text-align:left;">
0.065-0.315
</td>
</tr>
<tr>
<td style="text-align:left;">
Capusan et al. (2016)
</td>
<td style="text-align:left;">
STAGE
</td>
<td style="text-align:left;">
Abuse
</td>
<td style="text-align:left;">
ADHD
</td>
<td style="text-align:right;">
0.240
</td>
<td style="text-align:left;">
0.015-0.465
</td>
</tr>
<tr>
<td style="text-align:left;">
Capusan et al. (2016)
</td>
<td style="text-align:left;">
STAGE
</td>
<td style="text-align:left;">
Neglect
</td>
<td style="text-align:left;">
ADHD
</td>
<td style="text-align:right;">
0.150
</td>
<td style="text-align:left;">
-0.005-0.305
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinkler et al. (2017)
</td>
<td style="text-align:left;">
CATSS
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
ADHD
</td>
<td style="text-align:right;">
0.260
</td>
<td style="text-align:left;">
-0.065-0.585
</td>
</tr>
<tr>
<td style="text-align:left;">
Dinkler et al. (2017)
</td>
<td style="text-align:left;">
CATSS
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Autism
</td>
<td style="text-align:right;">
0.500
</td>
<td style="text-align:left;">
0.200-0.800
</td>
</tr>
<tr>
<td style="text-align:left;">
Schaefer et al. (2017)
</td>
<td style="text-align:left;">
E-Risk
</td>
<td style="text-align:left;">
Victimisation
</td>
<td style="text-align:left;">
P-factor
</td>
<td style="text-align:right;">
0.644
</td>
<td style="text-align:left;">
0.385-0.904
</td>
</tr>
<tr>
<td style="text-align:left;">
Schaefer et al. (2017)
</td>
<td style="text-align:left;">
E-Risk
</td>
<td style="text-align:left;">
Victimisation
</td>
<td style="text-align:left;">
Internalising
</td>
<td style="text-align:right;">
0.655
</td>
<td style="text-align:left;">
0.449-0.862
</td>
</tr>
<tr>
<td style="text-align:left;">
Schaefer et al. (2017)
</td>
<td style="text-align:left;">
E-Risk
</td>
<td style="text-align:left;">
Victimisation
</td>
<td style="text-align:left;">
Externalising
</td>
<td style="text-align:right;">
0.676
</td>
<td style="text-align:left;">
0.488-0.863
</td>
</tr>
<tr>
<td style="text-align:left;">
Schaefer et al. (2017)
</td>
<td style="text-align:left;">
E-Risk
</td>
<td style="text-align:left;">
Victimisation
</td>
<td style="text-align:left;">
Thought disorder
</td>
<td style="text-align:right;">
0.698
</td>
<td style="text-align:left;">
0.475-0.920
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
ASD
</td>
<td style="text-align:right;">
0.772
</td>
<td style="text-align:left;">
0.340-1.204
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
ASD
</td>
<td style="text-align:right;">
1.378
</td>
<td style="text-align:left;">
0.838-1.919
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
ASD
</td>
<td style="text-align:right;">
0.937
</td>
<td style="text-align:left;">
0.289-1.586
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
ASD
</td>
<td style="text-align:right;">
0.937
</td>
<td style="text-align:left;">
0.289-1.586
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Inattention overactivity
</td>
<td style="text-align:right;">
0.827
</td>
<td style="text-align:left;">
0.287-1.367
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Inattention overactivity
</td>
<td style="text-align:right;">
0.551
</td>
<td style="text-align:left;">
0.011-1.092
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Inattention overactivity
</td>
<td style="text-align:right;">
0.882
</td>
<td style="text-align:left;">
0.342-1.422
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Inattention overactivity
</td>
<td style="text-align:right;">
1.048
</td>
<td style="text-align:left;">
0.399-1.696
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Emotional problems
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:left;">
-0.648-0.648
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Emotional problems
</td>
<td style="text-align:right;">
0.165
</td>
<td style="text-align:left;">
-0.375-0.706
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Emotional problems
</td>
<td style="text-align:right;">
0.551
</td>
<td style="text-align:left;">
-0.205-1.308
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Emotional problems
</td>
<td style="text-align:right;">
1.048
</td>
<td style="text-align:left;">
0.399-1.696
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Emotional problems
</td>
<td style="text-align:right;">
0.055
</td>
<td style="text-align:left;">
-0.377-0.487
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Emotional problems
</td>
<td style="text-align:right;">
0.221
</td>
<td style="text-align:left;">
-0.212-0.653
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Emotional problems
</td>
<td style="text-align:right;">
0.717
</td>
<td style="text-align:left;">
0.176-1.257
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Conduct problem
</td>
<td style="text-align:right;">
0.276
</td>
<td style="text-align:left;">
-0.373-0.924
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Conduct problem
</td>
<td style="text-align:right;">
0.717
</td>
<td style="text-align:left;">
0.068-1.365
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Conduct problem
</td>
<td style="text-align:right;">
0.276
</td>
<td style="text-align:left;">
-0.373-0.924
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Conduct problem
</td>
<td style="text-align:right;">
1.158
</td>
<td style="text-align:left;">
0.401-1.914
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Conduct problem
</td>
<td style="text-align:right;">
0.441
</td>
<td style="text-align:left;">
-0.099-0.981
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Conduct problem
</td>
<td style="text-align:right;">
0.055
</td>
<td style="text-align:left;">
-0.485-0.595
</td>
</tr>
<tr>
<td style="text-align:left;">
Sonuga-Barke et al. (2017)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Conduct problem
</td>
<td style="text-align:right;">
-0.221
</td>
<td style="text-align:left;">
-0.761-0.320
</td>
</tr>
<tr>
<td style="text-align:left;">
Zvara et al. (2017)
</td>
<td style="text-align:left;">
FLP
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Postnatal depression
</td>
<td style="text-align:right;">
0.290
</td>
<td style="text-align:left;">
0.016-0.564
</td>
</tr>
<tr>
<td style="text-align:left;">
Ma et al. (2018)
</td>
<td style="text-align:left;">
FFCWS
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Aggressive behaviour
</td>
<td style="text-align:right;">
0.704
</td>
<td style="text-align:left;">
0.512-0.895
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (males)
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Suicidal ideation
</td>
<td style="text-align:right;">
0.532
</td>
<td style="text-align:left;">
0.298-0.767
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (males)
</td>
<td style="text-align:left;">
Physical neglect
</td>
<td style="text-align:left;">
Suicidal ideation
</td>
<td style="text-align:right;">
0.258
</td>
<td style="text-align:left;">
-0.315-0.832
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (males)
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Suicidal ideation
</td>
<td style="text-align:right;">
0.439
</td>
<td style="text-align:left;">
0.203-0.675
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (females)
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Suicidal ideation
</td>
<td style="text-align:right;">
0.576
</td>
<td style="text-align:left;">
0.363-0.788
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (females)
</td>
<td style="text-align:left;">
Physical neglect
</td>
<td style="text-align:left;">
Suicidal ideation
</td>
<td style="text-align:right;">
0.567
</td>
<td style="text-align:left;">
0.156-0.978
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (females)
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Suicidal ideation
</td>
<td style="text-align:right;">
0.539
</td>
<td style="text-align:left;">
0.339-0.740
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (males)
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Suicidal plan
</td>
<td style="text-align:right;">
0.462
</td>
<td style="text-align:left;">
0.099-0.826
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (males)
</td>
<td style="text-align:left;">
Physical neglect
</td>
<td style="text-align:left;">
Suicidal plan
</td>
<td style="text-align:right;">
0.167
</td>
<td style="text-align:left;">
-0.437-0.771
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (males)
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Suicidal plan
</td>
<td style="text-align:right;">
0.396
</td>
<td style="text-align:left;">
0.062-0.731
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (females)
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Suicidal plan
</td>
<td style="text-align:right;">
0.447
</td>
<td style="text-align:left;">
0.111-0.782
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (females)
</td>
<td style="text-align:left;">
Physical neglect
</td>
<td style="text-align:left;">
Suicidal plan
</td>
<td style="text-align:right;">
0.314
</td>
<td style="text-align:left;">
-0.170-0.799
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (females)
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Suicidal plan
</td>
<td style="text-align:right;">
0.463
</td>
<td style="text-align:left;">
0.136-0.790
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (males)
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Suicide attempt
</td>
<td style="text-align:right;">
0.510
</td>
<td style="text-align:left;">
0.060-0.961
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (males)
</td>
<td style="text-align:left;">
Physical neglect
</td>
<td style="text-align:left;">
Suicide attempt
</td>
<td style="text-align:right;">
-0.348
</td>
<td style="text-align:left;">
-1.323-0.626
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (males)
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Suicide attempt
</td>
<td style="text-align:right;">
0.358
</td>
<td style="text-align:left;">
-0.127-0.842
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (females)
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Suicide attempt
</td>
<td style="text-align:right;">
0.651
</td>
<td style="text-align:left;">
0.378-0.924
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (females)
</td>
<td style="text-align:left;">
Physical neglect
</td>
<td style="text-align:left;">
Suicide attempt
</td>
<td style="text-align:right;">
0.961
</td>
<td style="text-align:left;">
0.400-1.521
</td>
</tr>
<tr>
<td style="text-align:left;">
Obikane et al. (2018)
</td>
<td style="text-align:left;">
J-SHINE (females)
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Suicide attempt
</td>
<td style="text-align:right;">
0.659
</td>
<td style="text-align:left;">
0.394-0.924
</td>
</tr>
<tr>
<td style="text-align:left;">
Stern et al. (2018)
</td>
<td style="text-align:left;">
E-Risk
</td>
<td style="text-align:left;">
Victimisation
</td>
<td style="text-align:left;">
ADHD
</td>
<td style="text-align:right;">
0.140
</td>
<td style="text-align:left;">
0.022-0.259
</td>
</tr>
<tr>
<td style="text-align:left;">
Stern et al. (2018)
</td>
<td style="text-align:left;">
E-Risk
</td>
<td style="text-align:left;">
Victimisation
</td>
<td style="text-align:left;">
ADHD
</td>
<td style="text-align:right;">
0.345
</td>
<td style="text-align:left;">
0.225-0.465
</td>
</tr>
<tr>
<td style="text-align:left;">
Baldwin et al. (2019)
</td>
<td style="text-align:left;">
E-Risk
</td>
<td style="text-align:left;">
Victimisation
</td>
<td style="text-align:left;">
Suicidal ideation
</td>
<td style="text-align:right;">
0.205
</td>
<td style="text-align:left;">
0.053-0.357
</td>
</tr>
<tr>
<td style="text-align:left;">
Baldwin et al. (2019)
</td>
<td style="text-align:left;">
E-Risk
</td>
<td style="text-align:left;">
Victimisation
</td>
<td style="text-align:left;">
Self harm
</td>
<td style="text-align:right;">
0.224
</td>
<td style="text-align:left;">
0.091-0.356
</td>
</tr>
<tr>
<td style="text-align:left;">
Baldwin et al. (2019)
</td>
<td style="text-align:left;">
E-Risk
</td>
<td style="text-align:left;">
Victimisation
</td>
<td style="text-align:left;">
Suicide attempt
</td>
<td style="text-align:right;">
0.136
</td>
<td style="text-align:left;">
-0.104-0.376
</td>
</tr>
<tr>
<td style="text-align:left;">
Gerin et al. (2019)
</td>
<td style="text-align:left;">
DNS
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Internalising
</td>
<td style="text-align:right;">
0.620
</td>
<td style="text-align:left;">
0.326-0.914
</td>
</tr>
<tr>
<td style="text-align:left;">
Kugler et al. (2019)
</td>
<td style="text-align:left;">
FADS
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Drug use
</td>
<td style="text-align:right;">
0.362
</td>
<td style="text-align:left;">
0.156-0.568
</td>
</tr>
<tr>
<td style="text-align:left;">
Kugler et al. (2019)
</td>
<td style="text-align:left;">
FADS
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.235
</td>
<td style="text-align:left;">
0.030-0.441
</td>
</tr>
<tr>
<td style="text-align:left;">
Lecei et al. (2019)
</td>
<td style="text-align:left;">
TwinssCan
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Total psychopathology
</td>
<td style="text-align:right;">
0.295
</td>
<td style="text-align:left;">
0.054-0.537
</td>
</tr>
<tr>
<td style="text-align:left;">
Lecei et al. (2019)
</td>
<td style="text-align:left;">
TwinssCan
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Psychosis
</td>
<td style="text-align:right;">
0.270
</td>
<td style="text-align:left;">
0.029-0.512
</td>
</tr>
<tr>
<td style="text-align:left;">
Lecei et al. (2019)
</td>
<td style="text-align:left;">
TwinssCan
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Anxiety
</td>
<td style="text-align:right;">
0.386
</td>
<td style="text-align:left;">
0.143-0.628
</td>
</tr>
<tr>
<td style="text-align:left;">
Lecei et al. (2019)
</td>
<td style="text-align:left;">
TwinssCan
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.254
</td>
<td style="text-align:left;">
0.012-0.495
</td>
</tr>
<tr>
<td style="text-align:left;">
Schwartz et al. (2019)
</td>
<td style="text-align:left;">
MIDUS
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.116
</td>
<td style="text-align:left;">
-0.018-0.249
</td>
</tr>
<tr>
<td style="text-align:left;">
Schwartz et al. (2019)
</td>
<td style="text-align:left;">
MIDUS
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Antisocial behaviour
</td>
<td style="text-align:right;">
0.226
</td>
<td style="text-align:left;">
0.093-0.360
</td>
</tr>
<tr>
<td style="text-align:left;">
Schwartz et al. (2019)
</td>
<td style="text-align:left;">
Add Health
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.107
</td>
<td style="text-align:left;">
0.036-0.177
</td>
</tr>
<tr>
<td style="text-align:left;">
Alvanzo et al. (2020)
</td>
<td style="text-align:left;">
NESARC (males)
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Severe alcohol problems
</td>
<td style="text-align:right;">
-0.090
</td>
<td style="text-align:left;">
-0.329-0.150
</td>
</tr>
<tr>
<td style="text-align:left;">
Alvanzo et al. (2020)
</td>
<td style="text-align:left;">
NESARC (females)
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Severe alcohol problems
</td>
<td style="text-align:right;">
0.005
</td>
<td style="text-align:left;">
-0.222-0.233
</td>
</tr>
<tr>
<td style="text-align:left;">
Alvanzo et al. (2020)
</td>
<td style="text-align:left;">
NESARC (males)
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Moderate alcohol problems
</td>
<td style="text-align:right;">
0.067
</td>
<td style="text-align:left;">
-0.106-0.241
</td>
</tr>
<tr>
<td style="text-align:left;">
Alvanzo et al. (2020)
</td>
<td style="text-align:left;">
NESARC (females)
</td>
<td style="text-align:left;">
ACEs
</td>
<td style="text-align:left;">
Moderate alcohol problems
</td>
<td style="text-align:right;">
0.136
</td>
<td style="text-align:left;">
-0.028-0.300
</td>
</tr>
<tr>
<td style="text-align:left;">
Golm et al. (2020)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.480
</td>
<td style="text-align:left;">
0.088-0.872
</td>
</tr>
<tr>
<td style="text-align:left;">
Golm et al. (2020)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Generalised anxiety
</td>
<td style="text-align:right;">
0.490
</td>
<td style="text-align:left;">
0.098-0.882
</td>
</tr>
<tr>
<td style="text-align:left;">
Golm et al. (2020)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.410
</td>
<td style="text-align:left;">
-0.002-0.822
</td>
</tr>
<tr>
<td style="text-align:left;">
Golm et al. (2020)
</td>
<td style="text-align:left;">
ERA
</td>
<td style="text-align:left;">
Institutional neglect
</td>
<td style="text-align:left;">
Generalised anxiety
</td>
<td style="text-align:right;">
0.380
</td>
<td style="text-align:left;">
-0.032-0.792
</td>
</tr>
<tr>
<td style="text-align:left;">
Kullberg et al. (2020)
</td>
<td style="text-align:left;">
NESDA
</td>
<td style="text-align:left;">
Emotional abuse
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.509
</td>
<td style="text-align:left;">
0.351-0.667
</td>
</tr>
<tr>
<td style="text-align:left;">
Kullberg et al. (2020)
</td>
<td style="text-align:left;">
NESDA
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
-0.002
</td>
<td style="text-align:left;">
-0.158-0.153
</td>
</tr>
<tr>
<td style="text-align:left;">
Kullberg et al. (2020)
</td>
<td style="text-align:left;">
NESDA
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.068
</td>
<td style="text-align:left;">
-0.087-0.224
</td>
</tr>
<tr>
<td style="text-align:left;">
Kullberg et al. (2020)
</td>
<td style="text-align:left;">
NESDA
</td>
<td style="text-align:left;">
Emotional abuse
</td>
<td style="text-align:left;">
Anxiety
</td>
<td style="text-align:right;">
0.292
</td>
<td style="text-align:left;">
0.135-0.448
</td>
</tr>
<tr>
<td style="text-align:left;">
Kullberg et al. (2020)
</td>
<td style="text-align:left;">
NESDA
</td>
<td style="text-align:left;">
Physical abuse
</td>
<td style="text-align:left;">
Anxiety
</td>
<td style="text-align:right;">
0.073
</td>
<td style="text-align:left;">
-0.083-0.228
</td>
</tr>
<tr>
<td style="text-align:left;">
Kullberg et al. (2020)
</td>
<td style="text-align:left;">
NESDA
</td>
<td style="text-align:left;">
Sexual abuse
</td>
<td style="text-align:left;">
Anxiety
</td>
<td style="text-align:right;">
0.104
</td>
<td style="text-align:left;">
-0.051-0.260
</td>
</tr>
<tr>
<td style="text-align:left;">
Capusan et al. (2021)
</td>
<td style="text-align:left;">
Östergötland cohort
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Substance use disorder
</td>
<td style="text-align:right;">
0.771
</td>
<td style="text-align:left;">
0.488-1.053
</td>
</tr>
<tr>
<td style="text-align:left;">
Isumi et al. (2021)
</td>
<td style="text-align:left;">
A-CHILD
</td>
<td style="text-align:left;">
Maltreatment
</td>
<td style="text-align:left;">
Behavioural difficulties
</td>
<td style="text-align:right;">
0.333
</td>
<td style="text-align:left;">
0.281-0.384
</td>
</tr>
<tr>
<td style="text-align:left;">
Li et al. (2021)
</td>
<td style="text-align:left;">
Chinese long. study
</td>
<td style="text-align:left;">
Emotional abuse
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.100
</td>
<td style="text-align:left;">
0.037-0.163
</td>
</tr>
<tr>
<td style="text-align:left;">
Li et al. (2021)
</td>
<td style="text-align:left;">
Chinese long. study
</td>
<td style="text-align:left;">
Emotional abuse
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.100
</td>
<td style="text-align:left;">
0.037-0.163
</td>
</tr>
<tr>
<td style="text-align:left;">
Li et al. (2021)
</td>
<td style="text-align:left;">
Chinese long. study
</td>
<td style="text-align:left;">
Emotional abuse
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.100
</td>
<td style="text-align:left;">
0.037-0.163
</td>
</tr>
<tr>
<td style="text-align:left;">
Li et al. (2021)
</td>
<td style="text-align:left;">
Chinese long. study
</td>
<td style="text-align:left;">
Emotional abuse
</td>
<td style="text-align:left;">
Depression
</td>
<td style="text-align:right;">
0.080
</td>
<td style="text-align:left;">
0.017-0.143
</td>
</tr>
</tbody>
</table>
