# nonparametric book by S. Jung
# chapter 10. Part III: robust linear regression
#
#
# ----------------
library(MASS)
library(tidyverse) 
library(latex2exp) 
library(cowplot)
library(GGally)
library(WRS2)

library(showtext)
#font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()


# Use the following steps to install `WRS` package
# (Visit https://github.com/nicebread/WRS )
# # first: install dependent packages
# install.packages(c("MASS", "akima", "robustbase"))
# 
# # second: install suggested packages
# install.packages(c("akima", "cobs", "robust", "mgcv", "scatterplot3d", "quantreg", "rrcov", "lars", "pwr", "trimcluster", "mc2d", "psych", "Rfit", "DepthProc", "class", "fda", "rankFD"))
# 
# # third: install an additional package which provides some C functions
# # install.packages("devtools")
# # NOTE: This seems to be stalled and not functional any more
# # devtools::install_github("mrxiaohe/WRScpp")
# 
# # fourth: install WRS
# devtools::install_github("nicebread/WRS", subdir="pkg")
library(WRS)


# ex 1 --------------------------------------------------------------------
# mean, trimmed mean, and winsorized mean

x = c(-5,-3,-3,-3,-1,0,1,2,3,4,1000)
c(mean(x), mean(x, tr = 0.1), WRS2::winmean(x, tr = 0.1))




# figure ------------------------------------------------------------------
# Psi functions

# 함수들
x = seq(-5,5,by = 0.01)
k = 1.28
c = 4.685
dat = data.frame(x=x,LS = x, LAD = sign(x), 
                 Huber = pmax(-k, pmin(k,x)),
                 biweight = ifelse( abs(x)<=c, x*(1-(x/c)^2)^2,0) )
dat %>% pivot_longer(cols = 2:5, names_to = "Type") %>% 
  ggplot(aes(x = x, y = value, linetype = Type, color = Type)) + geom_line(size = 1)

ggsave(filename = "ch10m1.pdf", path = 'images', width = 6.5, height = 5)
ggsave(filename = "ch10m1.png", path = 'images', width = 6.5, height = 5) # two-panel

# x = seq(-1.5,1.5,by = 0.01)
# k = 1.28
# dat = data.frame(x=x,LS = x, LAD = sign(x), 
#                  Huber = pmax(-k, pmin(k,x)),
#                  biweight = ifelse( abs(x)<=1, x*(1-x^2)^2,0) )
# dat %>% pivot_longer(cols = 3:5, names_to = "Type") %>% 
#   ggplot(aes(x = x, y = value, linetype = Type)) + geom_line(size = 1.5)



# ex 2 --------------------------------------------------------------------
WRS2::mest(x,bend = 1.28)  


# inference ---------------------------------------------------------------

# Use Belgium Phone Calls 1950-1973 data 
data("phones")
as.data.frame(phones) %>% ggplot(aes(x = year, y = calls)) + geom_point()
ggsave(filename = "ch10m2.pdf", path = 'images', width = 6.5, height = 5)
ggsave(filename = "ch10m2.png", path = 'images', width = 6.5, height = 5)  

# Pretend that the "#calls" data are not a time-series, and estimate the location.
# compare mean, trimmed means, winsorized means and huber M-estimator. 
y = phones$calls

x = c(-5,-3,-3,-3,-1,0,1,2,3,4,1000)


# Inference 
# H0: true average = 15, vs H1: not H0
x = MASS::phones$calls
mean(x)
t.test(x, mu = 15)$p.value

# 윌콕슨 signed-rank 
wilcox.test(x,mu = 15, alternative = "two.sided", 
            exact = F, conf.int = T)

# 부호검정 
s = sign(x-15)
n = sum(s != 0) 
(Sobs = sum(s > 0))
(pvalue.sign.test = 2*(1- pbinom(Sobs,n,1/2))) # accept

alpha = 0.05
l = qbinom(alpha/2, size = n, prob = 1/2); u = n+ 1 - l
(ci.sign.test.x = sort(x)[c(l,u)]) # CI by sign test 

# bootstrap CI 
(mest.obs = WRS2::mest(x))

set.seed(1)
B = 1999
mest.boot = replicate(B,WRS2::mest(sample(x, length(x), replace = T)))
(se.est = sd(mest.boot))
(bias.est = mean(mest.boot) - mest.obs)
(norm.app.boot.CI = mest.obs - bias.est + se.est * qnorm(0.025) *c(1,-1))
(boot.perc.ci.mest = quantile(mest.boot, probs = c(0.025,1-0.025)))


theta.jack = vector(length = n)
for (i in 1:n){theta.jack[i] = WRS2::mest(x[-i])} 
acc = sum((mean(theta.jack) - theta.jack)^3)/
  (6 * sum(( mean(theta.jack) - theta.jack)^2)^(3/2))
acc
(b.factor = qnorm( sum(mest.boot <= mest.obs)/(B+1)))
bz = b.factor + c(qnorm(0.025), qnorm(1-0.025))
(BCa.int.probs = pnorm( (bz)/(1 - acc * (bz)) + b.factor ))
(BCa.interval = quantile(mest.boot, BCa.int.probs))


# careful with bootstrapping! Resampling might increase the number of "outliers"

# bootstrap-t CI
set.seed(1)
mest.obs = mest(x)
st.mest.boot = replicate(B,{ 
  x.boot = sample(x, length(x), replace = T)
  (WRS2::mest(x.boot)- mest.obs)/mad(x.boot)} 
) %>% t()    # m-estimator / MAD
(boot.t.ci = mest.obs - mad(x) * quantile(st.mest.boot, probs = c(1-0.025,0.025)))
# careful with studentization! Variation of MAD is often too large.
# Resulting in too large standard deviation of the bootstrapped t-distribution
# Using bootstrap-t is not recommended. 



# robust regression  ------------------------------------------------------



# ex:TS-simple  -----------------------------------------------------------

wh = read.csv("data/whappiness.csv")
x = wh$GDP.per.capita
y = wh$Score
cor.test(x,y, method = "k")$p.value

n = length(x)
pairids = gtools::combinations(n,2) 
pairbetas = data.frame(i = pairids[,1], j = pairids[,2], 
                       betaij = rep(0,nrow(pairids)))
for (i in 1:nrow(pairids)){
  betaijs = (y[pairids[i,1]] - y[pairids[i,2]]) / 
    (x[pairids[i,1]] - x[pairids[i,2]])
  pairbetas[i,3] = betaijs
}
pairbeta_rm = pairbetas %>% filter(betaij != Inf)
(betaTS = median(pairbeta_rm$betaij))
(alphaTS = median(y - betaTS * x))

wh %>% ggplot(aes(GDP.per.capita, Score)) + geom_point() + geom_smooth(method = "lm", se = F) +
  geom_abline(intercept = alphaTS, slope = betaTS)

ggsave(filename = "ch10m3.pdf", path = 'images', width = 6.5, height = 5)
ggsave(filename = "ch10m3.png", path = 'images', width = 6.5, height = 5)  
mblm::mblm(y ~ x, repeated = FALSE)


(l = (1 + SuppDists::qKendall(0.025,n)) * choose(n,2)/2)
(u = (1 + SuppDists::qKendall(1-0.025,n)) * choose(n,2)/2)
(TS.CI = sort(pairbeta_rm$betaij)[c(l,u)])



# ex:TS-multiple  ---------------------------------------------------------
h = read.csv("data/worldbank2020.csv")
dat = h %>% dplyr::select(3:9)
colnames(dat) = c("Score", "GDP", "Social", 
                  "Life", "Freedom", "Generosity", "Corruption")
 

ts.out = WRS::tsreg(as.matrix(dat[,2:7]), dat$Score,iter = 100)
ts.out$coef 

lm(Score ~ ., data = dat)$coefficients


# ex:m-regression1 --------------------------------------------------------


wh = read.csv("data/whappiness.csv")
x = wh$GDP.per.capita
y = wh$Score
MASS::rlm(y ~ x, method = "M")$coefficients
MASS::rlm(y ~ x, method = "MM")$coefficients


m1 = MASS::rlm(y ~ x, method = "M")$coefficients
m2 = MASS::rlm(y ~ x, method = "MM")$coefficients 
g1 = wh %>% ggplot(aes(GDP.per.capita, Score)) + geom_point() + geom_smooth(method = "lm", se = F) +
  geom_abline(intercept = m1[1], slope = m1[2], color = "red", linetype = "dashed", size = 1.5)+
  geom_abline(intercept = m2[1], slope = m2[2], color = "red", linetype = "dotted", size = 1.5) +
  ggtitle("x[9] = 1.537")

x[9] = 100
MASS::rlm(y ~ x, method = "M")$coefficients
MASS::rlm(y ~ x, method = "MM")$coefficients


m11 = MASS::rlm(y ~ x, method = "M")$coefficients
m21 = MASS::rlm(y ~ x, method = "MM")$coefficients
wh2 = wh
wh2$GDP.per.capita[9] = 100
g2 = wh2 %>% ggplot(aes(GDP.per.capita, Score)) + geom_point() + geom_smooth(method = "lm", se = F) +
  geom_abline(intercept = m11[1], slope = m11[2], color = "red", linetype = "dashed", size = 1.5)+
  geom_abline(intercept = m21[1], slope = m21[2], color = "red", linetype = "dotted", size = 1.5) + 
  coord_cartesian(xlim=c(1.2, 1.6)) + ggtitle("x[9] = 100")
plot_grid(g1,g2)
ggsave(filename = "ch10m4.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch10m4.png", path = 'images', width = 6.5, height = 4)  



# ex: 10.21 ex:robust_regression_inference ------------------------------------------
h = read.csv("data/worldbank2020.csv")
dat = h %>% dplyr::select(3:9)
colnames(dat) = c("Score", "GDP", "Social", 
                  "Life", "Freedom", "Generosity", "Corruption")

# 1. standard error 
MASS::rlm(Score ~ ., data = dat, method = "MM") %>% summary()

# 2. bootstrap-estimated stan.err.
n = nrow(dat)
set.seed(1)
B = 1999
MM.boot = replicate(B, {
  boot.id = sample(1:n,n,replace = TRUE)
  dat.boot = dat[boot.id,]
  MASS::rlm(Score ~ ., data = dat.boot, method = "MM")$coefficients}
) %>% t()
apply(MM.boot,2,sd)


# ex: 10.22 ---------------------------------------------------------------
# 3. bootstrap-t CI 
MM.out = summary(rlm(Score ~ ., data = dat, method = "MM"))[[4]]
beta.hat = MM.out[,1]
beta.hat.se = MM.out[,2]

B = 1999
set.seed(1)
# resampling t 
MM.boot = replicate(B, {
  boot.id = sample(1:n,n,replace = TRUE)
  dat.boot = dat[boot.id,]
  MM.out.boot = summary(rlm(Score ~ ., data = dat.boot, method = "MM"))[[4]]
  (MM.out.boot[,1] - beta.hat)/MM.out.boot[,2]}
) %>% t()

ci = matrix(NA, nrow = 7, ncol = 2); colnames(ci) = c("2.5%","97.5%")
row.names(ci) = colnames(MM.boot)
for (i in 1:7){
  ci[i,] = beta.hat[i] -  beta.hat.se[i] * quantile(MM.boot[,i],probs=c(1-0.025,0.025))
}
ci

# ex:10.23 permutation test (with WLS) ---------------------------------------------
(t.stats = beta.hat / beta.hat.se)
set.seed(1)
M = 1999
MM.perm = replicate(M,{
  perm.id = sample(1:n,n, replace = FALSE)
  dat.perm = dat
  dat.perm$Score = dat$Score[perm.id]
  MM.out.perm = summary(rlm(Score ~ ., data = dat.perm, method = "MM"))[[4]]
  (MM.out.perm[,1])/MM.out.perm[,2]}
) %>% t()

perm.pvalues = vector(length = 7)
for (i in 1:7){
  perm.pvalues[i] = mean( abs(c(t.stats[i],MM.perm[,i])) >= abs(t.stats[i]))
}
perm.pvalues 
 
 
