library(bvarsv)
library(readxl)
library(dplyr)
library(zoo)
library(forecast)
library(writexl)
library(tempdisagg)
library(rlist)
library(rgl)
library(plot3D)
library(gridExtra)
library(ggplot2)

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
                           conversion = "mean", method = "litterman-maxlog", 
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

plot.ts(bkf.ts, main= "Government Spending, Tax Revenues, and GDP 2011-2020 in Percentage (%)")

# Dates in string format
tml <- paste0(floor(time(bkf.ts)),"M", (1+12*( time(bkf.ts) - floor(time(bkf.ts)))))
class(bkf.ts) # just to make sure that the dataset is time series

head(bkf.ts, n=10)

# Convenience functions for plotting
matplot2 <-function(x, y, ylim, ...){
  matplot(x =x,y =y,ylim =ylim,type ="l",xlab ="",ylab ="",
          lty =1,lwd =2,bty ="n", ...)
}
abline2 <-function(...){
  abline(...,lty =4,lwd =0.3)
}
gp <- seq(2011,2021, 6)
# Color palette, taken from http://www.cookbook-r.com/Graphs/Colors_%28ggplot2%29/
cbPalette <- c("#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
               "#D55E00","#CC79A7")
cols1 <-cbPalette[ c(2,4,2)]
cols2 <-cbPalette[ c(6,4,6)]

# Fix number of Gibbs sampler draws used in this vignette
nburn.vignette <- 80000
nrep.vignette <- 100000

# Run estimation
fit <- bvar.sv.tvp(bkf.ts, p = 3, nburn = nburn.vignette, nrep = nrep.vignette, thinfac = 100)
saveRDS(fit, file = "fit BP + total spending mom 2011 lag2 80% sosial nominal apct chow-lin.rds")
fit <- readRDS("fit BP + total spending mom 2011 lag2 80% sosial nominal apct chow-lin.rds")

# SD of Gov residual
# X axis
# lag 1: xax <- seq(2018.500,2021,1/12)
# lag 2: xax <- seq(2018.583,2021,1/12)
# lag 3: xax <- seq(2018.666,2021,1/12)
# lag 4: xax <- seq(2018.749,2021,1/12)
# lag 5: xax <- seq(2018.832,2021,1/12)

xax <- seq(2014.583,2021,1/12)
# Get posterior draws
aux <- seq(1,4*length(xax),4)
x1 <- sqrt(fit$H.draws[1, aux, ])
x1 <- t(apply(x1,1, function(z) c(mean(z), quantile(z, c(0.16,0.84)))[ c(2,1,3)]))
# Plot
matplot2(x =xax,y =x1,ylim = c(0,30),col =cols1, main ="Spending on Social")
#abline2(h = c(0.2,0.4,0.6),v =gp)

# SD of Tax residual
# Get posterior draws
x2 <- sqrt(fit$H.draws[2, aux +1, ])
x2 <- t(apply(x2,1, function(z) c(mean(z), quantile(z, c(0.16,0.84)))[ c(2,1,3)]))
# Plot
matplot2(x =xax,y =x2,ylim = c(0,30),col =cols1, main ="Spending on Social")
# abline2(h = c(0.5,1),v =gp)

# SD of GDP residual
# Get posterior draws
x3 <- sqrt(fit$H.draws[3, aux +2, ])
x3 <- t(apply(x3,1, function(z) c(mean(z), quantile(z, c(0.16,0.84)))[ c(2,1,3)]))
# Plot
matplot2(x =xax,y =x3,ylim = c(0,0.75),col =cols1, main ="GDP")
# abline2(h =1:4,v =gp)

# SD of int residual
# Get posterior draws
x4 <- sqrt(fit$H.draws[4, aux +3, ])
x4 <- t(apply(x4,1, function(z) c(mean(z), quantile(z, c(0.16,0.84)))[ c(2,1,3)]))
# Plot
matplot2(x =xax,y =x4,ylim = c(0,1.5),col =cols1,
         main ="Interest Rate")
abline2(h =1:4,v =gp)

# Impulse Response of spending
all_dts <- c("2019M1","2020M4")

tmp <- list()
for (dd in all_dts){
  t <- which(tml == dd) - 42 #lag 1=41, lag 2=42, lag 3=43, lag 4=44, lag 5=45
  t_ind <- which(all_dts == dd)
  # Compute impulse responses
  aux <- impulse.responses(fit,impulse.variable =2,
                           response.variable =2, t =t,
                           scenario =3,
                           draw.plot =FALSE)$irf
  tmp[[t_ind]] <- aux
}

# Make data for graph
gdat <- rbind(0, sapply(tmp, function(z) apply(as.matrix(z),2, median)))
# Configure and print plot
yb <- c(-1,1)

plot_title <- paste0("Impulse Response Function")
matplot2(x =0:20,y =gdat,ylim =yb,col =cols2,main =plot_title)
# abline2(v = seq(5,20,5),h = seq(yb[1], yb[2],0.05))
legend_loc <-"topright"
legend(legend_loc,legend =all_dts,col =cols2,lty =1,lwd =2,bty ="n")

# Multiplier (response of GDP)
all_dts2 <- tml[43:nrow(bkf.ts)] #lag 1=42, lag 2=43, lag 3=44, lag 4=45, lag 5=46

tmp2 <- list()
tmp2_mean <- list()
time_horizon <- seq(1:25)
sumlist <- setNames(vector(length(time_horizon), mode="list"), time_horizon)

for (dd in all_dts2) {
  t <- which(tml == dd) - 42 #lag 1=41, lag 2=42, lag 3=43, lag 4=44, lag 5=45
  t_ind <- which(all_dts2 == dd)
  
  # Compute impulse responses
  aux2 <- impulse.responses(fit,impulse.variable = 2,
                            response.variable = 4, t = t,
                            scenario =3,
                            draw.plot =FALSE, nhor = 25)$irf
  tmp2[[t_ind]] <- aux2
  for (i in 1:ncol(tmp2[[t_ind]])) {
    
    mean_try <- mean(tmp2[[t_ind]][,i])
    tmp2_mean[[i]] <- mean_try
    
    # Cumulative Impulse Response
    if (i==ncol(tmp2[[t_ind]])) {
      mean_try.df <- data.frame(tmp2_mean)
      
      sumlist[[1]][[t_ind]] <- mean_try.df[1,1]
      
      for (p in 2:25) {
        sum_p <- rowSums( mean_try.df[,1:p])
        sumlist[[p]][[t_ind]] <- sum_p
      }
      # next
    } else if (i>ncol(tmp2[[t_ind]])) {
      print("The number of loops is more than the periods estimated")
      break
    }
  }
}

# Multiplier (response of gov_spending)
all_dts2 <- tml[43:nrow(bkf.ts)] #lag 1=42, lag 2=43, lag 3=44, lag 4=45, lag 5=46

tmp3 <- list()
tmp3_mean <- list()
sumlist2 <- setNames(vector(length(time_horizon), mode="list"), time_horizon)

for (dd in all_dts2) {
  t <- which(tml == dd) - 42 #lag 1=41, lag 2=42, lag 3=43, lag 4=44, lag 5=45
  t_ind <- which(all_dts2 == dd)
  
  # Compute impulse responses
  aux3 <- impulse.responses(fit,impulse.variable = 2,
                            response.variable = 2, t = t,
                            scenario =3,
                            draw.plot =FALSE, nhor = 25)$irf
  tmp3[[t_ind]] <- aux3
  for (i in 1:ncol(tmp3[[t_ind]])) {
    
    mean_try <- mean(tmp3[[t_ind]][,i])
    tmp3_mean[[i]] <- mean_try
    
    # Cumulative Impulse Response
    if (i==ncol(tmp3[[t_ind]])) {
      mean_try.df <- data.frame(tmp3_mean)
      sumlist2[[1]][[t_ind]] <- mean_try.df[1,1]
      
      for (p in 2:25) {
        sum_p <- rowSums( mean_try.df[,1:p])
        sumlist2[[p]][[t_ind]] <- sum_p
      }
      # next
    } else if (i>ncol(tmp3[[t_ind]])) {
      print("The number of loops is more than the periods estimated")
      break
    }
  }
}


# Multiplier with conversion factor
# Pembagian response

multiplier_list <- list()
for (e in 1:25) {
  gdp_response <- as.numeric(sumlist[[e]])
  gov_response <- as.numeric(sumlist2[[e]])
  multiplier <- gdp_response/gov_response
  multiplier_list[[e]] <- multiplier
}

gdp <- bkf_raw$gdp[43:nrow(bkf.ts)] #lag 1=42, lag 2=43, lag 3=44, lag 4=45, lag 5=46
g <- bkf_raw$sosial[43:nrow(bkf.ts)]
conversion <- g/gdp
plot.ts(conversion)

timeline <- xax
multiplier_fix <- list()
for (o in 1:25) {
  multiplier <- multiplier_list[[o]]
  conversion_result <- multiplier/conversion
  conversion_result_df <- data.frame(conversion_result)
  identifier <- data.frame(rep(o, length(timeline)))
  test_df <- cbind(timeline, identifier, conversion_result_df)
  colnames(test_df) <- c("timeline", "time horizon", "multiplier")
  multiplier_fix[[o]] <- data.frame(test_df)
}

# 2D plot
rbind_multipliers <- bind_rows(multiplier_fix)
data_test2 <- rbind_multipliers
data_test2$time.horizon <- as.factor(data_test2$time.horizon)
data_test2 <- data_test2 %>%
  filter(time.horizon==3 | time.horizon==6 | time.horizon==12)

ggplot(data_test2, aes(x = timeline, y = multiplier, group=time.horizon, 
                       col=time.horizon, fill=time.horizon)) +
  geom_line(aes(linetype=time.horizon), size=1) + 
  geom_point() +
  scale_linetype_manual(values=c("solid", "solid", "dotdash")) +
  scale_color_manual(values=c('#000000','#E74C3C', '#2E86C1'))+
  theme_minimal() +
  theme(legend.position="top")
# ylim(NA, NA)

