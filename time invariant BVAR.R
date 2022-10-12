library(readxl)
library(dplyr)
library(zoo)
library(forecast)
library(writexl)
library(tempdisagg)
library(rlist)
library(BVAR)

# Interpolasi data PDB Nominal dan Deflator dgn variabel referensi M2 dan CPI
setwd("")
PdbNominal <- (read_excel("r_data_bkf.xlsx", sheet = "triwulan"))$gdp
Deflator <- (read_excel("r_data_bkf.xlsx", sheet = "triwulan"))$deflator
PdbNominal <- ts(PdbNominal, start=c(2010, 1), end=c(2020, 4), frequency=4)
Deflator <- ts(Deflator, start=c(2010, 1), end=c(2020, 4), frequency=4)

Bulanan <- read_excel("r_data_bkf.xlsx", sheet = "bulanan")
m2 <- Bulanan$m2
m2 <- ts(m2, start=c(2010, 1), end=c(2020, 12), frequency=12)
cpi <- Bulanan$cpi
cpi <- ts(cpi, start=c(2010, 1), end=c(2020, 12), frequency=12)
ipi <- Bulanan$ipi
ipi <- ts(ipi, start=c(2010, 1), end=c(2020, 12), frequency=12)

pdb_interpolasi <- td(PdbNominal ~ m2,
                      conversion = "sum", method = "chow-lin-maxlog", # method="denton-cholette" or "litterman-maxlog"
                      to="monthly")
PdbBulanan <- predict(pdb_interpolasi)

pdb_interpolasi_spline = spline(PdbNominal, method = "fmm", n=nrow(Bulanan))
PdbBulananSpline = pdb_interpolasi_spline$y

deflator_interpolasi <- td(Deflator ~ cpi,
                           conversion = "mean", method = "chow-lin-maxlog", 
                           to="monthly")
DeflatorBulanan <- predict(deflator_interpolasi)

# Seasonal Adjusment
BulananNsa <- cbind(Bulanan, PdbBulanan, DeflatorBulanan)
BulananNsa <- ts(BulananNsa, start=c(2010, 1), end=c(2020, 12), frequency=12)
# BulananNsa <- window(BulananNsa, start=c(2014, 1), end=c(2020, 12), frequency=12)

SaHasil <- list()
for (col in colnames(BulananNsa)) {
  Stl <- stl(BulananNsa[,col], s.window="periodic", robust=TRUE)
  SaHasil[[col]] <- seasadj(Stl)
}
Sa <- list.cbind(SaHasil)

#============================================Estimasi Multiplier==========================================#

bkf_raw <- data.frame(Sa)
bkf_raw = select(bkf_raw, -deflator)
bkf_raw = rename(bkf_raw, gdp=PdbBulanan, deflator=DeflatorBulanan)
pdb <- bkf_raw$gdp
deflator <- bkf_raw$deflator

pct <- function(x) {
  (x-lag(x, n=1))/lag(x, n=1)*100
}

apct2 = function(x) {
  (((x/lag(x, n=1))^12)-1)*100
}

apct <- function(x) {
  (((x/lag(x, n=1))^(1/12))-1)*100
}

real <- function(x) {
  x/deflator
}

int1 <- diff(bkf_raw$int1/100, lag=1)
int2 <- diff(bkf_raw$int2/100, lag=1)
int3 <- diff(bkf_raw$int3, lag=1)

int1_raw <- bkf_raw$int1
int2_raw <- bkf_raw$int2
int3_raw <- bkf_raw$int3

bkf.real <- bkf_raw %>%
  mutate_each(funs(real), c(gdp, pegawai, barang, modal, sosial, total_spending, tax))

bkf.1 <- bkf.real %>%
  mutate_each(funs(apct), c(gdp, pegawai, barang, modal, sosial, total_spending, tax, cpi, int3))
vector_ = rep(0, length(na.approx(bkf.1$sosial))+1)
vector_[2:length(vector_)] = na.approx(bkf.1$sosial)
bkf.1$sosial = vector_

vector_ = rep(0, length(na.approx(bkf.1$total_spending))+1)
vector_[2:length(vector_)] = na.approx(bkf.1$total_spending)
bkf.1$total_spending = vector_

vector_ = rep(0, length(na.approx(bkf.1$tax))+1)
vector_[2:length(vector_)] = na.approx(bkf.1$tax)
bkf.1$tax = vector_

vector_ = rep(0, length(na.approx(bkf.1$gdp))+1)
vector_[2:length(vector_)] = na.approx(bkf.1$gdp)
bkf.1$gdp = vector_

bkf <- cbind(bkf.1$sosial, bkf.1$tax, bkf.1$gdp)
colnames(bkf) <- c("Spending", "Tax", "GDP")
bkf.ts <- ts(bkf, start=c(2010, 1), end=c(2020, 12), frequency=12)
bkf.ts <- window(bkf.ts, start=c(2010, 1), end=c(2020, 12), frequency=12)

plot.ts(bkf.ts, main= "Government Spending, Tax Revenues, and GDP 2010-2020")

# BVAR estimation
bvar_fit = bvar(bkf.ts, lags = 3, 
                n_draw = 100000, n_burn = 80000, n_thin = 100L)
# Compute + store IRF
irf_bvar_fit = irf(bvar_fit, bv_irf(horizon = 36L, fevd = TRUE))#, n_thin = 10L)
# Update the confidence bands of the IRFs
irf_bvar_fit = irf(bvar_fit, conf_bands = 0.05)

plot(irf_bvar_fit, vars_impulse = "Spending", vars_response = "GDP")

