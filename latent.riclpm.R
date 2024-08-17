
load("~/Documents/Projects/Replicate-riclpm/df.Rda")

library(semTools)
library(tidyverse)
library(lavaan)
library(semTools)
library(psych)
library(MBESS)


names(df)
# names(df)
# [1] "CNOTMHon"   "r3"        "r2"        "r"         "h3"        "h2"
# [7] "h"         "age"       "gender"    "T3RegId1"  "T3RegId2"  "T3RegId3"
# [13] "T3RegId4"  "T2RegId1"  "T2RegId2"  "T2RegId3"  "T2RegId4"  "RegId1"
# [19] "RegId2"    "RegId3"    "RegId4"    "T3HonCon1" "T3HonCon2" "T3HonCon3"
# [25] "T3HonCon4" "T3HonCon5" "T2HonCon1" "T2HonCon2" "T2HonCon3" "T2HonCon4"
# [31] "T2HonCon5" "HonCon1"   "HonCon2"   "HonCon3"   "HonCon4"   "HonCon5"

psych::describe(df)
df <- df %>% ungroup() %>% mutate(age.uncentered = age)

df <- df %>% mutate(age = age - mean(age, na.rm = TRUE))

df<-df %>% mutate(gender.original = gender)
df<-df %>% mutate(gender = if_else(gender == 1, 0, 1))


# Replicate ----




RICLPM.Base <- '
#VARIANCES of latent variables

RegIdT1 ~~ 0*RegIdT1
RegIdT2 ~~ 0*RegIdT2
RegIdT3 ~~ 0*RegIdT3

HonT1 ~~ 0*HonT1
HonT2 ~~ 0*HonT2
HonT3 ~~ 0*HonT3

################
# BETWEEN PART #
################
RIRegId =~ 1*RegIdT1 + 1*RegIdT2 + 1*RegIdT3
RIHon =~ 1*HonT1 + 1*HonT2 + 1*HonT3

RIRegId ~~ psi.7_7*RIRegId
RIHon ~~ psi.8_8*RIHon
RIRegId ~~ psi.7_8*RIHon

###############
# WITHIN PART #
###############
T1WRegId =~ 1*RegIdT1
T2WRegId =~ 1*RegIdT2
T3WRegId =~ 1*RegIdT3

T1WHon =~ 1*HonT1
T2WHon =~ 1*HonT2
T3WHon =~ 1*HonT3

T1WRegId ~~ NA*T1WRegId + psi.9_9*T1WRegId
T2WRegId ~~ NA*T2WRegId + psi.11_11*T2WRegId
T3WRegId ~~ NA*T3WRegId + psi.13_13*T3WRegId

T1WHon ~~ NA*T1WHon + psi.10_10*T1WHon
T2WHon ~~ NA*T2WHon + psi.12_12*T2WHon
T3WHon ~~ NA*T3WHon + psi.14_14*T3WHon

T1WRegId ~~ NA*T1WHon + psi.9_10*T1WHon
T2WRegId ~~ NA*T2WHon + psi.11_12*T2WHon
T3WRegId ~~ NA*T3WHon + psi.13_14*T3WHon

###########
# CL + AR #
###########
T2WRegId ~ beta.11_9*T1WRegId  + beta.11_10*T1WHon
T2WHon ~ beta.12_9*T1WRegId  + beta.12_10*T1WHon

T3WRegId ~ beta.13_11*T2WRegId + beta.13_12*T2WHon
T3WHon ~ beta.14_11*T2WRegId + beta.14_12*T2WHon

##########################
# ADDITIONAL CONSTRAINTS #
##########################
RIRegId + RIHon ~~ 0*T1WRegId + 0*T1WHon

beta.11_9 == beta.13_11
beta.11_10 == beta.13_12
beta.12_9 == beta.14_11
beta.12_10 == beta.14_12
'


df.v1 <- df %>%
  rowwise() %>%
  mutate(
    RegIdT1 = mean(c(RegId1, RegId2, RegId3, RegId4), na.rm = TRUE),
    RegIdT2 = mean(c(T2RegId1, T2RegId2, T2RegId3, T2RegId4), na.rm = TRUE),
    RegIdT3 = mean(c(T3RegId1, T3RegId2, T3RegId3, T3RegId4), na.rm = TRUE),
    HonT1 = mean(c(HonCon1, HonCon2, HonCon3, HonCon4, HonCon5), na.rm = TRUE),
    HonT2 = mean(c(T2HonCon1, T2HonCon2, T2HonCon3, T2HonCon4, T2HonCon5), na.rm = TRUE),
    HonT3 = mean(c(T3HonCon1, T3HonCon2, T3HonCon3, T3HonCon4, T3HonCon5), na.rm = TRUE)
  ) %>%
  ungroup()

RICLPM.fit.v1a <- lavaan(RICLPM.Base, data = df.v1, estimator = 'ML', missing = 'ml.x', meanstructure = T, int.ov.free = T)
summary(RICLPM.fit.v1a,
        standardized = T,
        fit.measures = TRUE,
        rsquare= TRUE,
        nd=3)


# latent version----

mod.conf <- "
RegIdT1 =~ RegId1 + RegId2 + RegId3 + RegId4
RegIdT2 =~ T2RegId1 + T2RegId2 + T2RegId3 + T2RegId4
RegIdT3 =~ T3RegId1 + T3RegId2 + T3RegId3 + T3RegId4

HonT1 =~ HonCon1 + HonCon2 + HonCon3 + HonCon4 + HonCon5
HonT2 =~ T2HonCon1 + T2HonCon2 + T2HonCon3 + T2HonCon4 + T2HonCon5
HonT3 =~ T3HonCon1 + T3HonCon2 + T3HonCon3 + T3HonCon4 + T3HonCon5
"


longFacNames <- list(RegId = c("RegIdT1","RegIdT2","RegIdT3"),
                     Hon = c("HonT1","HonT2","HonT3"))


config.syntax<- measEq.syntax(configural.model = c(mod.conf),
                              data = df,
                              ID.fac = "fixed.factor",
                              longFacNames = longFacNames)


metric.syntax<- measEq.syntax(configural.model = c(mod.conf),
                              data = df,
                              ID.fac = "fixed.factor",
                              longFacNames = longFacNames,
                              long.equal  = c("loadings"))

cat(as.character(metric.syntax))
scalar.syntax<- measEq.syntax(configural.model = c(mod.conf),
                              data = df,
                              ID.fac = "fixed.factor",
                              longFacNames = longFacNames,
                              long.equal  = c("loadings","intercepts"))
cat(as.character(scalar.syntax))


# scalar.syntax<- measEq.syntax(configural.model = c(mod.conf),
#                               data = df,
#                               ID.fac = "EC",
#                               longFacNames = longFacNames,
#                               long.equal  = c("loadings","intercepts"))
# cat(as.character(scalar.syntax))


resid.syntax <- measEq.syntax(configural.model = c(mod.conf),
                              data = df,
                              ID.fac = "fixed.factor",
                              longFacNames = longFacNames,
                              long.equal  = c("loadings","intercepts","residuals"))

config.fit <- cfa(as.character(config.syntax),
                  data = df,
                  estimator = "MLR",
                  missing = "FIML")

metric.fit <- cfa(as.character(metric.syntax),
                  data = df,
                  estimator = "MLR",
                  missing = "FIML")

metric.minv.long <- semTools::compareFit(config.fit,metric.fit,nested = TRUE)
summary(metric.minv.long)

scalar.fit <- cfa(as.character(scalar.syntax),
                  data = df,
                  estimator = "MLR",
                  missing = "FIML")

scalar.minv.long <- semTools::compareFit(metric.fit,scalar.fit, nested = TRUE)
summary(scalar.minv.long)

cat(as.character(scalar.syntax))

resid.fit <- cfa(as.character(resid.syntax),
                 data = df,
                 estimator = "MLR",
                 missing = "FIML")

resid.minv.long <- semTools::compareFit(scalar.fit,resid.fit, nested = TRUE)
summary(resid.minv.long)

summary(config.fit,fit.measures = T)
summary(metric.fit,fit.measures = T)
summary(scalar.fit,fit.measures = T)
summary(resid.fit,fit.measures = T)


measure.v1 <-"
## LOADINGS:

RegIdT1 =~ NA*RegId1 + lambda.1_1*RegId1
RegIdT1 =~ NA*RegId2 + lambda.2_1*RegId2
RegIdT1 =~ NA*RegId3 + lambda.3_1*RegId3
RegIdT1 =~ NA*RegId4 + lambda.4_1*RegId4
RegIdT2 =~ NA*T2RegId1 + lambda.1_1*T2RegId1
RegIdT2 =~ NA*T2RegId2 + lambda.2_1*T2RegId2
RegIdT2 =~ NA*T2RegId3 + lambda.3_1*T2RegId3
RegIdT2 =~ NA*T2RegId4 + lambda.4_1*T2RegId4
RegIdT3 =~ NA*T3RegId1 + lambda.1_1*T3RegId1
RegIdT3 =~ NA*T3RegId2 + lambda.2_1*T3RegId2
RegIdT3 =~ NA*T3RegId3 + lambda.3_1*T3RegId3
RegIdT3 =~ NA*T3RegId4 + lambda.4_1*T3RegId4
HonT1 =~ NA*HonCon1 + lambda.13_4*HonCon1
HonT1 =~ NA*HonCon2 + lambda.14_4*HonCon2
HonT1 =~ NA*HonCon3 + lambda.15_4*HonCon3
HonT1 =~ NA*HonCon4 + lambda.16_4*HonCon4
HonT1 =~ NA*HonCon5 + lambda.17_4*HonCon5
HonT2 =~ NA*T2HonCon1 + lambda.13_4*T2HonCon1
HonT2 =~ NA*T2HonCon2 + lambda.14_4*T2HonCon2
HonT2 =~ NA*T2HonCon3 + lambda.15_4*T2HonCon3
HonT2 =~ NA*T2HonCon4 + lambda.16_4*T2HonCon4
HonT2 =~ NA*T2HonCon5 + lambda.17_4*T2HonCon5
HonT3 =~ NA*T3HonCon1 + lambda.13_4*T3HonCon1
HonT3 =~ NA*T3HonCon2 + lambda.14_4*T3HonCon2
HonT3 =~ NA*T3HonCon3 + lambda.15_4*T3HonCon3
HonT3 =~ NA*T3HonCon4 + lambda.16_4*T3HonCon4
HonT3 =~ NA*T3HonCon5 + lambda.17_4*T3HonCon5

## INTERCEPTS:

RegId1 ~ NA*1 + nu.1*1
RegId2 ~ NA*1 + nu.2*1
RegId3 ~ NA*1 + nu.3*1
RegId4 ~ NA*1 + nu.4*1
T2RegId1 ~ NA*1 + nu.1*1
T2RegId2 ~ NA*1 + nu.2*1
T2RegId3 ~ NA*1 + nu.3*1
T2RegId4 ~ NA*1 + nu.4*1
T3RegId1 ~ NA*1 + nu.1*1
T3RegId2 ~ NA*1 + nu.2*1
T3RegId3 ~ NA*1 + nu.3*1
T3RegId4 ~ NA*1 + nu.4*1
HonCon1 ~ NA*1 + nu.13*1
HonCon2 ~ NA*1 + nu.14*1
HonCon3 ~ NA*1 + nu.15*1
HonCon4 ~ NA*1 + nu.16*1
HonCon5 ~ NA*1 + nu.17*1
T2HonCon1 ~ NA*1 + nu.13*1
T2HonCon2 ~ NA*1 + nu.14*1
T2HonCon3 ~ NA*1 + nu.15*1
T2HonCon4 ~ NA*1 + nu.16*1
T2HonCon5 ~ NA*1 + nu.17*1
T3HonCon1 ~ NA*1 + nu.13*1
T3HonCon2 ~ NA*1 + nu.14*1
T3HonCon3 ~ NA*1 + nu.15*1
T3HonCon4 ~ NA*1 + nu.16*1
T3HonCon5 ~ NA*1 + nu.17*1

## UNIQUE-FACTOR VARIANCES:

RegId1 ~~ NA*RegId1 + theta.1_1*RegId1
RegId2 ~~ NA*RegId2 + theta.2_2*RegId2
RegId3 ~~ NA*RegId3 + theta.3_3*RegId3
RegId4 ~~ NA*RegId4 + theta.4_4*RegId4
T2RegId1 ~~ NA*T2RegId1 + theta.5_5*T2RegId1
T2RegId2 ~~ NA*T2RegId2 + theta.6_6*T2RegId2
T2RegId3 ~~ NA*T2RegId3 + theta.7_7*T2RegId3
T2RegId4 ~~ NA*T2RegId4 + theta.8_8*T2RegId4
T3RegId1 ~~ NA*T3RegId1 + theta.9_9*T3RegId1
T3RegId2 ~~ NA*T3RegId2 + theta.10_1c(0, 0)*T3RegId2
T3RegId3 ~~ NA*T3RegId3 + theta.11_11*T3RegId3
T3RegId4 ~~ NA*T3RegId4 + theta.12_12*T3RegId4
HonCon1 ~~ NA*HonCon1 + theta.13_13*HonCon1
HonCon2 ~~ NA*HonCon2 + theta.14_14*HonCon2
HonCon3 ~~ NA*HonCon3 + theta.15_15*HonCon3
HonCon4 ~~ NA*HonCon4 + theta.16_16*HonCon4
HonCon5 ~~ NA*HonCon5 + theta.17_17*HonCon5
T2HonCon1 ~~ NA*T2HonCon1 + theta.18_18*T2HonCon1
T2HonCon2 ~~ NA*T2HonCon2 + theta.19_19*T2HonCon2
T2HonCon3 ~~ NA*T2HonCon3 + theta.20_2c(0, 0)*T2HonCon3
T2HonCon4 ~~ NA*T2HonCon4 + theta.21_21*T2HonCon4
T2HonCon5 ~~ NA*T2HonCon5 + theta.22_22*T2HonCon5
T3HonCon1 ~~ NA*T3HonCon1 + theta.23_23*T3HonCon1
T3HonCon2 ~~ NA*T3HonCon2 + theta.24_24*T3HonCon2
T3HonCon3 ~~ NA*T3HonCon3 + theta.25_25*T3HonCon3
T3HonCon4 ~~ NA*T3HonCon4 + theta.26_26*T3HonCon4
T3HonCon5 ~~ NA*T3HonCon5 + theta.27_27*T3HonCon5

## UNIQUE-FACTOR COVARIANCES:

RegId1 ~~ NA*T2RegId1 + theta.5_1*T2RegId1
RegId1 ~~ NA*T3RegId1 + theta.9_1*T3RegId1
RegId2 ~~ NA*T2RegId2 + theta.6_2*T2RegId2
RegId2 ~~ NA*T3RegId2 + theta.10_2*T3RegId2
RegId3 ~~ NA*T2RegId3 + theta.7_3*T2RegId3
RegId3 ~~ NA*T3RegId3 + theta.11_3*T3RegId3
RegId4 ~~ NA*T2RegId4 + theta.8_4*T2RegId4
RegId4 ~~ NA*T3RegId4 + theta.12_4*T3RegId4
T2RegId1 ~~ NA*T3RegId1 + theta.9_5*T3RegId1
T2RegId2 ~~ NA*T3RegId2 + theta.10_6*T3RegId2
T2RegId3 ~~ NA*T3RegId3 + theta.11_7*T3RegId3
T2RegId4 ~~ NA*T3RegId4 + theta.12_8*T3RegId4
HonCon1 ~~ NA*T2HonCon1 + theta.18_13*T2HonCon1
HonCon1 ~~ NA*T3HonCon1 + theta.23_13*T3HonCon1
HonCon2 ~~ NA*T2HonCon2 + theta.19_14*T2HonCon2
HonCon2 ~~ NA*T3HonCon2 + theta.24_14*T3HonCon2
HonCon3 ~~ NA*T2HonCon3 + theta.20_15*T2HonCon3
HonCon3 ~~ NA*T3HonCon3 + theta.25_15*T3HonCon3
HonCon4 ~~ NA*T2HonCon4 + theta.21_16*T2HonCon4
HonCon4 ~~ NA*T3HonCon4 + theta.26_16*T3HonCon4
HonCon5 ~~ NA*T2HonCon5 + theta.22_17*T2HonCon5
HonCon5 ~~ NA*T3HonCon5 + theta.27_17*T3HonCon5
T2HonCon1 ~~ NA*T3HonCon1 + theta.23_18*T3HonCon1
T2HonCon2 ~~ NA*T3HonCon2 + theta.24_19*T3HonCon2
T2HonCon3 ~~ NA*T3HonCon3 + theta.25_2c(0, 0)*T3HonCon3
T2HonCon4 ~~ NA*T3HonCon4 + theta.26_21*T3HonCon4
T2HonCon5 ~~ NA*T3HonCon5 + theta.27_22*T3HonCon5

## LATENT MEANS/INTERCEPTS:

RegIdT1 ~ c(0, 0)*1 + alpha.1*1
RegIdT2 ~ NA*1 + alpha.2*1
RegIdT3 ~ NA*1 + alpha.3*1
HonT1 ~ c(0, 0)*1 + alpha.4*1
HonT2 ~ NA*1 + alpha.5*1
HonT3 ~ NA*1 + alpha.6*1

# ## COMMON-FACTOR VARIANCES:
#
# RegIdT1 ~~ 1*RegIdT1 + psi.1_1*RegIdT1
# RegIdT2 ~~ NA*RegIdT2 + psi.2_2*RegIdT2
# RegIdT3 ~~ NA*RegIdT3 + psi.3_3*RegIdT3
# HonT1 ~~ 1*HonT1 + psi.4_4*HonT1
# HonT2 ~~ NA*HonT2 + psi.5_5*HonT2
# HonT3 ~~ NA*HonT3 + psi.6_6*HonT3
#
# ## COMMON-FACTOR COVARIANCES:
#
# RegIdT1 ~~ NA*RegIdT2 + psi.2_1*RegIdT2
# RegIdT1 ~~ NA*RegIdT3 + psi.3_1*RegIdT3
# RegIdT1 ~~ NA*HonT1 + psi.4_1*HonT1
# RegIdT1 ~~ NA*HonT2 + psi.5_1*HonT2
# RegIdT1 ~~ NA*HonT3 + psi.6_1*HonT3
# RegIdT2 ~~ NA*RegIdT3 + psi.3_2*RegIdT3
# RegIdT2 ~~ NA*HonT1 + psi.4_2*HonT1
# RegIdT2 ~~ NA*HonT2 + psi.5_2*HonT2
# RegIdT2 ~~ NA*HonT3 + psi.6_2*HonT3
# RegIdT3 ~~ NA*HonT1 + psi.4_3*HonT1
# RegIdT3 ~~ NA*HonT2 + psi.5_3*HonT2
# RegIdT3 ~~ NA*HonT3 + psi.6_3*HonT3
# HonT1 ~~ NA*HonT2 + psi.5_4*HonT2
# HonT1 ~~ NA*HonT3 + psi.6_4*HonT3
# HonT2 ~~ NA*HonT3 + psi.6_5*HonT3
"


measure.v2 <-"
## LOADINGS:

RegIdT1 =~ NA*RegId1 + lambda.1_1*RegId1
RegIdT1 =~ NA*RegId2 + lambda.2_1*RegId2
RegIdT1 =~ NA*RegId3 + lambda.3_1*RegId3
RegIdT1 =~ NA*RegId4 + lambda.4_1*RegId4
RegIdT2 =~ NA*T2RegId1 + lambda.1_1*T2RegId1
RegIdT2 =~ NA*T2RegId2 + lambda.2_1*T2RegId2
RegIdT2 =~ NA*T2RegId3 + lambda.3_1*T2RegId3
RegIdT2 =~ NA*T2RegId4 + lambda.4_1*T2RegId4
RegIdT3 =~ NA*T3RegId1 + lambda.1_1*T3RegId1
RegIdT3 =~ NA*T3RegId2 + lambda.2_1*T3RegId2
RegIdT3 =~ NA*T3RegId3 + lambda.3_1*T3RegId3
RegIdT3 =~ NA*T3RegId4 + lambda.4_1*T3RegId4
HonT1 =~ NA*HonCon1 + lambda.13_4*HonCon1
HonT1 =~ NA*HonCon2 + lambda.14_4*HonCon2
HonT1 =~ NA*HonCon3 + lambda.15_4*HonCon3
HonT1 =~ NA*HonCon4 + lambda.16_4*HonCon4
HonT1 =~ NA*HonCon5 + lambda.17_4*HonCon5
HonT2 =~ NA*T2HonCon1 + lambda.13_4*T2HonCon1
HonT2 =~ NA*T2HonCon2 + lambda.14_4*T2HonCon2
HonT2 =~ NA*T2HonCon3 + lambda.15_4*T2HonCon3
HonT2 =~ NA*T2HonCon4 + lambda.16_4*T2HonCon4
HonT2 =~ NA*T2HonCon5 + lambda.17_4*T2HonCon5
HonT3 =~ NA*T3HonCon1 + lambda.13_4*T3HonCon1
HonT3 =~ NA*T3HonCon2 + lambda.14_4*T3HonCon2
HonT3 =~ NA*T3HonCon3 + lambda.15_4*T3HonCon3
HonT3 =~ NA*T3HonCon4 + lambda.16_4*T3HonCon4
HonT3 =~ NA*T3HonCon5 + lambda.17_4*T3HonCon5

## INTERCEPTS:

RegId1 ~ NA*1 + nu.1*1
RegId2 ~ NA*1 + nu.2*1
RegId3 ~ NA*1 + nu.3*1
RegId4 ~ NA*1 + nu.4*1
T2RegId1 ~ NA*1 + nu.1*1
T2RegId2 ~ NA*1 + nu.2*1
T2RegId3 ~ NA*1 + nu.3*1
T2RegId4 ~ NA*1 + nu.4*1
T3RegId1 ~ NA*1 + nu.1*1
T3RegId2 ~ NA*1 + nu.2*1
T3RegId3 ~ NA*1 + nu.3*1
T3RegId4 ~ NA*1 + nu.4*1
HonCon1 ~ NA*1 + nu.13*1
HonCon2 ~ NA*1 + nu.14*1
HonCon3 ~ NA*1 + nu.15*1
HonCon4 ~ NA*1 + nu.16*1
HonCon5 ~ NA*1 + nu.17*1
T2HonCon1 ~ NA*1 + nu.13*1
T2HonCon2 ~ NA*1 + nu.14*1
T2HonCon3 ~ NA*1 + nu.15*1
T2HonCon4 ~ NA*1 + nu.16*1
T2HonCon5 ~ NA*1 + nu.17*1
T3HonCon1 ~ NA*1 + nu.13*1
T3HonCon2 ~ NA*1 + nu.14*1
T3HonCon3 ~ NA*1 + nu.15*1
T3HonCon4 ~ NA*1 + nu.16*1
T3HonCon5 ~ NA*1 + nu.17*1

## UNIQUE-FACTOR VARIANCES:

RegId1 ~~ NA*RegId1 + theta.1_1*RegId1
RegId2 ~~ NA*RegId2 + theta.2_2*RegId2
RegId3 ~~ NA*RegId3 + theta.3_3*RegId3
RegId4 ~~ NA*RegId4 + theta.4_4*RegId4
T2RegId1 ~~ NA*T2RegId1 + theta.5_5*T2RegId1
T2RegId2 ~~ NA*T2RegId2 + theta.6_6*T2RegId2
T2RegId3 ~~ NA*T2RegId3 + theta.7_7*T2RegId3
T2RegId4 ~~ NA*T2RegId4 + theta.8_8*T2RegId4
T3RegId1 ~~ NA*T3RegId1 + theta.9_9*T3RegId1
T3RegId2 ~~ NA*T3RegId2 + theta.10_10*T3RegId2
T3RegId3 ~~ NA*T3RegId3 + theta.11_11*T3RegId3
T3RegId4 ~~ NA*T3RegId4 + theta.12_12*T3RegId4
HonCon1 ~~ NA*HonCon1 + theta.13_13*HonCon1
HonCon2 ~~ NA*HonCon2 + theta.14_14*HonCon2
HonCon3 ~~ NA*HonCon3 + theta.15_15*HonCon3
HonCon4 ~~ NA*HonCon4 + theta.16_16*HonCon4
HonCon5 ~~ NA*HonCon5 + theta.17_17*HonCon5
T2HonCon1 ~~ NA*T2HonCon1 + theta.18_18*T2HonCon1
T2HonCon2 ~~ NA*T2HonCon2 + theta.19_19*T2HonCon2
T2HonCon3 ~~ NA*T2HonCon3 + theta.20_20*T2HonCon3
T2HonCon4 ~~ NA*T2HonCon4 + theta.21_21*T2HonCon4
T2HonCon5 ~~ NA*T2HonCon5 + theta.22_22*T2HonCon5
T3HonCon1 ~~ NA*T3HonCon1 + theta.23_23*T3HonCon1
T3HonCon2 ~~ NA*T3HonCon2 + theta.24_24*T3HonCon2
T3HonCon3 ~~ NA*T3HonCon3 + theta.25_25*T3HonCon3
T3HonCon4 ~~ NA*T3HonCon4 + theta.26_26*T3HonCon4
T3HonCon5 ~~ NA*T3HonCon5 + theta.27_27*T3HonCon5

## UNIQUE-FACTOR COVARIANCES:

RegId1 ~~ NA*T2RegId1 + theta.5_1*T2RegId1
RegId1 ~~ NA*T3RegId1 + theta.9_1*T3RegId1
RegId2 ~~ NA*T2RegId2 + theta.6_2*T2RegId2
RegId2 ~~ NA*T3RegId2 + theta.10_2*T3RegId2
RegId3 ~~ NA*T2RegId3 + theta.7_3*T2RegId3
RegId3 ~~ NA*T3RegId3 + theta.11_3*T3RegId3
RegId4 ~~ NA*T2RegId4 + theta.8_4*T2RegId4
RegId4 ~~ NA*T3RegId4 + theta.12_4*T3RegId4
T2RegId1 ~~ NA*T3RegId1 + theta.9_5*T3RegId1
T2RegId2 ~~ NA*T3RegId2 + theta.10_6*T3RegId2
T2RegId3 ~~ NA*T3RegId3 + theta.11_7*T3RegId3
T2RegId4 ~~ NA*T3RegId4 + theta.12_8*T3RegId4
HonCon1 ~~ NA*T2HonCon1 + theta.18_13*T2HonCon1
HonCon1 ~~ NA*T3HonCon1 + theta.23_13*T3HonCon1
HonCon2 ~~ NA*T2HonCon2 + theta.19_14*T2HonCon2
HonCon2 ~~ NA*T3HonCon2 + theta.24_14*T3HonCon2
HonCon3 ~~ NA*T2HonCon3 + theta.20_15*T2HonCon3
HonCon3 ~~ NA*T3HonCon3 + theta.25_15*T3HonCon3
HonCon4 ~~ NA*T2HonCon4 + theta.21_16*T2HonCon4
HonCon4 ~~ NA*T3HonCon4 + theta.26_16*T3HonCon4
HonCon5 ~~ NA*T2HonCon5 + theta.22_17*T2HonCon5
HonCon5 ~~ NA*T3HonCon5 + theta.27_17*T3HonCon5
T2HonCon1 ~~ NA*T3HonCon1 + theta.23_18*T3HonCon1
T2HonCon2 ~~ NA*T3HonCon2 + theta.24_19*T3HonCon2
T2HonCon3 ~~ NA*T3HonCon3 + theta.25_20*T3HonCon3
T2HonCon4 ~~ NA*T3HonCon4 + theta.26_21*T3HonCon4
T2HonCon5 ~~ NA*T3HonCon5 + theta.27_22*T3HonCon5

## LATENT MEANS/INTERCEPTS:

RegIdT1 ~ NA*1 + alpha.1*1
RegIdT2 ~ NA*1 + alpha.2*1
RegIdT3 ~ NA*1 + alpha.3*1
HonT1 ~ NA*1 + alpha.4*1
HonT2 ~ NA*1 + alpha.5*1
HonT3 ~ NA*1 + alpha.6*1


## MODEL CONSTRAINTS:

lambda.1_1 == 4 - lambda.2_1 - lambda.3_1 - lambda.4_1
nu.1 == 0 - nu.2 - nu.3 - nu.4
lambda.13_4 == 5 - lambda.14_4 - lambda.15_4 - lambda.16_4 - lambda.17_4
nu.13 == 0 - nu.14 - nu.15 - nu.16 - nu.17
"




MG.RICLPM.Base <- '
#VARIANCES of latent variables

RegIdT1 ~~ c(0, 0)*RegIdT1
RegIdT2 ~~ c(0, 0)*RegIdT2
RegIdT3 ~~ c(0, 0)*RegIdT3

HonT1 ~~ c(0, 0)*HonT1
HonT2 ~~ c(0, 0)*HonT2
HonT3 ~~ c(0, 0)*HonT3

################
# BETWEEN PART #
################
RIRegId =~ c(1, 1)*RegIdT1 + c(1, 1)*RegIdT2 + c(1, 1)*RegIdT3
RIHon =~ c(1, 1)*HonT1 + c(1, 1)*HonT2 + c(1, 1)*HonT3

RIRegId ~~ c(psi.7_7.g1, psi.7_7.g2)*RIRegId
RIHon ~~ c(psi.8_8.g1, psi.8_8.g2)*RIHon
RIRegId ~~ c(psi.7_8.g1, psi.7_8.g2)*RIHon

###############
# WITHIN PART #
###############
T1WRegId =~ c(1, 1)*RegIdT1
T2WRegId =~ c(1, 1)*RegIdT2
T3WRegId =~ c(1, 1)*RegIdT3

T1WHon =~ c(1, 1)*HonT1
T2WHon =~ c(1, 1)*HonT2
T3WHon =~ c(1, 1)*HonT3

T1WRegId ~~ c(NA, NA)*T1WRegId + c(psi.9_9.g1, psi.9_9.g2)*T1WRegId
T2WRegId ~~ c(NA, NA)*T2WRegId + c(psi.11_11.g1, psi.11_11.g2)*T2WRegId
T3WRegId ~~ c(NA, NA)*T3WRegId + c(psi.13_13.g1, psi.13_13.g2)*T3WRegId

T1WHon ~~ c(NA, NA)*T1WHon + c(psi.10_10.g1, psi.10_10.g2)*T1WHon
T2WHon ~~ c(NA, NA)*T2WHon + c(psi.12_12.g1, psi.12_12.g2)*T2WHon
T3WHon ~~ c(NA, NA)*T3WHon + c(psi.14_14.g1, psi.14_14.g2)*T3WHon

T1WRegId ~~ c(NA, NA)*T1WHon + c(psi.9_10.g1, psi.9_10.g2)*T1WHon
T2WRegId ~~ c(NA, NA)*T2WHon + c(psi.11_12.g1, psi.11_12.g2)*T2WHon
T3WRegId ~~ c(NA, NA)*T3WHon + c(psi.13_14.g1, psi.13_14.g2)*T3WHon

###########
# CL + AR #
###########
T2WRegId ~ c(beta.11_9.g1, beta.11_9.g2)*T1WRegId  + c(beta.11_10.g1, beta.11_10.g2)*T1WHon
T2WHon ~ c(beta.12_9.g1, beta.12_9.g2)*T1WRegId  + c(beta.12_10.g1, beta.12_10.g2)*T1WHon

T3WRegId ~ c(beta.13_11.g1, beta.13_11.g2)*T2WRegId + c(beta.13_12.g1, beta.13_12.g2)*T2WHon
T3WHon ~ c(beta.14_11.g1, beta.14_11.g2)*T2WRegId + c(beta.14_12.g1, beta.14_12.g2)*T2WHon

##########################
# ADDITIONAL CONSTRAINTS #
##########################
RIRegId + RIHon ~~ c(0, 0)*T1WRegId + c(0, 0)*T1WHon
'

model.nons <- c(measure.v2,RICLPM.Base)

RICLPM.fit.r <- lavaan(model.nons, data = df, estimator = 'MLR', missing = 'FIML', meanstructure = T, int.ov.free = T)

summary(RICLPM.fit.r,
        standardized = T,
        fit.measures = TRUE,
        rsquare= TRUE,
        nd=3)


RICLPM.fit.o2 <- lavaan(model.nons, data = df, estimator = 'ML', missing = 'ml.x', meanstructure = T, int.ov.free = T)
summary(RICLPM.fit.o2,
        standardized = T,
        fit.measures = TRUE,
        rsquare= TRUE,
        nd=3)

df.2 <- df %>% 
  mutate(RegIdT1 = 0.25*RegId1 + 0.25*RegId2 + 0.25*RegId3 + 0.25*RegId4,
         RegIdT2 = 0.25*T2RegId1 + 0.25*T2RegId2 + 0.25*T2RegId3 + 0.25*T2RegId4,
         RegIdT3 = 0.25*T3RegId1 + 0.25*T3RegId2 + 0.25*T3RegId3 + 0.25*T3RegId4,
         HonT1 = 0.2*HonCon1 + 0.2*HonCon2 + 0.2*HonCon3 + 0.2*HonCon4 + 0.2*HonCon5,
         HonT2 = 0.2*T2HonCon1 + 0.2*T2HonCon2 + 0.2*T2HonCon3 + 0.2*T2HonCon4 + 0.2*T2HonCon5,
         HonT3 = 0.2*T3HonCon1 + 0.2*T3HonCon2 + 0.2*T3HonCon3 + 0.2*T3HonCon4 + 0.2*T3HonCon5)

write.csv(df.2, "df.2.csv")

RICLPM.fit.2r <- lavaan(RICLPM.Base, data = df.2, estimator = 'MLR', missing = 'FIML', meanstructure = T, int.ov.free = T)

summary(RICLPM.fit.2r,
        standardized = T,
        fit.measures = TRUE,
        rsquare= TRUE,
        nd=3)

RICLPM.fit.o1 <- lavaan(RICLPM.Base, data = df.2, estimator = 'ML', missing = 'ml.x', meanstructure = T, int.ov.free = T)
summary(RICLPM.fit.o1,
        standardized = T,
        fit.measures = TRUE,
        rsquare= TRUE,
        nd=3)



df.3 <- df %>%
  mutate(RegIdT1 = r,
         RegIdT2 = r2,
         RegIdT3 = r3,
         HonT1 = h,
         HonT2 = h2,
         HonT3 = h3)
RICLPM.fit.o1a <- lavaan(RICLPM.Base, data = df.3, estimator = 'ML', missing = 'ml.x', meanstructure = T, int.ov.free = T)
summary(RICLPM.fit.o1a,
        standardized = T,
        fit.measures = TRUE,
        rsquare= TRUE,
        nd=3)


