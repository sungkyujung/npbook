# nonparametric book by S. Jung
# chapter 9. bootstrap test
#
#
# ----------------
library(tidyverse) 
library(cowplot)
library(latex2exp)
library(extraDistr)

library(showtext)
#font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()


# Leukemia data fig 1-----------------------------------------------------------
leukemia <- read.csv("data/leukemia4.csv")
leukemia_AML <- leukemia$ge[leukemia$label == "AML"]
x = leukemia_AML 
n = length(x)
ggplot(mapping = aes(leukemia_AML)) + geom_histogram(bins = 9, color = "white") + geom_rug()
ggsave(filename = "ch9_1.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch9_1.png", path = 'images', width = 6.5, height = 4)


# test statistics ---------------------------------------------------------
leukemia <- read.csv("data/leukemia4.csv")
x <- leukemia$ge[leukemia$label == "AML"]
# hypothesis H0: mu =mu0, H1: mu > mu0
mu0 = 166

# test statistic (two test stats)
t0 = sqrt(n) * (mean(x) - mu0 ) / sd(x)
d0 = mean(x) - mu0 
m0 = median(x) / mu0
c(t0,d0,m0)



# one sample bootstrap tests ----------------------------------------------


# assume Laplace dist'n, and get MLEs
mu.hat = median(x); b.hat =  mean(abs(x - mu.hat)) # mle (general)
b.hat0 = mean(abs(x-mu0))  # mle under H0

# parametric bootstrap for (T, D, M)
set.seed(1)
boot.out = replicate(3999, {
  x.boot = rlaplace(n, mu = mu.hat, sigma = b.hat) # sample under best-fitting model
  x.boot0 = rlaplace(n, mu = mu0, sigma = b.hat0) # sample under null model 
  c(mean(x.boot), sd(x.boot), mean(x.boot0), sd(x.boot0), median(x.boot), median(x.boot0))
}) %>% t() 

t.boot0 = sqrt(n) * (boot.out[,3] - mu0) / boot.out[,4] 
t.boot = sqrt(n) *(boot.out[,1] - mu.hat) / boot.out[,2]  
d.boot0 = boot.out[,3] - mu0  
d.boot = boot.out[,1] - mu.hat  
m.boot0 = boot.out[,6] / mu0
m.boot = boot.out[,5] / mu.hat


# p-values and a figure -------------------------------------------------
p.t0 = mean(c(t.boot0,t0) >= t0)
p.tpivot = mean(c(t.boot,t0) >= t0)
p.d0 = mean(c(d.boot0,d0) >= d0)
p.dpivot = mean(c(d.boot, d0) >= d0)
p.m0 = mean(c(m.boot0,m0) >= m0)
p.mpivot = mean(c(m.boot, m0) >= m0)
g1 = ggplot(mapping = aes(x = t.boot0, y = ..density..*10)) + geom_histogram(bins = 100) + geom_vline(xintercept = t0, color = "red", linetype = "dashed") + ggtitle(paste("T, 귀무가설 재표집. p값 = ",p.t0)) 
g2 = ggplot(mapping = aes(x = t.boot, y = ..density..*10)) + geom_histogram(bins = 100) + geom_vline(xintercept = t0, color = "red", linetype = "dashed") + ggtitle(paste("T, 피벗 재표집. p값 = ",p.tpivot)) 
g3 = ggplot(mapping = aes(x = d.boot0, y = ..density..*100)) + geom_histogram(bins = 100) + geom_vline(xintercept = d0, color = "red", linetype = "dashed") + ggtitle(paste("D, 귀무가설 재표집. p값 = ",p.d0)) 
g4 = ggplot(mapping = aes(x = d.boot, y = ..density..*100)) + geom_histogram(bins = 100) + geom_vline(xintercept = d0, color = "red", linetype = "dashed") + ggtitle(paste("D, 피벗 재표집. p값 = ",p.dpivot)) 
g5 = ggplot(mapping = aes(x = m.boot0, y = ..density..)) + geom_histogram(bins = 100) + geom_vline(xintercept = m0, color = "red", linetype = "dashed") + ggtitle(paste("M, 귀무가설 재표집. p값 = ",p.m0)) 
g6 = ggplot(mapping = aes(x = m.boot, y = ..density..)) + geom_histogram(bins = 100) + geom_vline(xintercept = m0, color = "red", linetype = "dashed") + ggtitle(paste("M, 피벗 재표집. p값 = ",p.mpivot)) 
plot_grid(g1,g2,g3, g4, g5,g6 , ncol = 2)

ggsave(filename = "ch9_2.pdf", path = 'images', width = 6.5, height = 9)
ggsave(filename = "ch9_2.png", path = 'images', width = 6.5, height = 9)

# two group ---------------------------------------------------------------
source("np_supp.R")

leukemia <- read.csv("data/leukemia4.csv")
y = leukemia$ge[leukemia$label == "AML"]; m = length(x)
x = leukemia$ge[leukemia$label == "ALL"]; n = length(y)


# test statistic (two test stats)
t0 = (mean(y) - mean(x)) / pooled.se(x,y) #t.test(y,x, var.equal = T)$statistic  
m0 = median(y) - median(x)

# assume Laplace, get MLE
mu1.hat = median(x); mu2.hat = median(y)
b.hat =  mean( abs( c(x - mu1.hat, y - mu2.hat))) # mle (general)
mu0.hat = median(c(x,y))
b.hat0 = mean( abs( c(x,y) - mu0.hat))  # mle under H0

# parametric bootstrap for both (T,M)
set.seed(1)
boot.out = replicate(3999, {
  x.boot = rlaplace(m, mu = mu1.hat, sigma = b.hat) # sample under best-fitting model
  y.boot = rlaplace(n, mu = mu2.hat, sigma = b.hat)
  
  x.boot0 = rlaplace(m, mu = mu0.hat, sigma = b.hat0) # sample under null model 
  y.boot0 = rlaplace(n, mu = mu0.hat, sigma = b.hat0) 
  
  t.boot0 = (mean(y.boot0) - mean(x.boot0)) / pooled.se(x.boot0,y.boot0)
  t.boot = (mean(y.boot) - mean(x.boot) - (mu2.hat - mu1.hat) ) / pooled.se(x.boot,y.boot)
  
  m.boot0 = median(y.boot0) - median(x.boot0) 
  m.boot = median(y.boot) - median(x.boot)- (mu2.hat - mu1.hat)
  
  c(t.boot0, t.boot, m.boot0, m.boot)
}) %>% t() 

t.boot0 = boot.out[,1] 
t.boot = boot.out[,2] 
m.boot0 = boot.out[,3] 
m.boot = boot.out[,4] 


# p-values
p.t0 = mean(c(t.boot0,t0) >= t0)
p.tpivot = mean(c(t.boot,t0) >= t0) 
p.m0 = mean(c(m.boot0,m0) >= m0)
p.mpivot = mean(c(m.boot, m0) >= m0)

# p-value and a figure ----------------------------------------------------
g1 = ggplot(mapping = aes(x = t.boot0, y = ..density..)) + geom_histogram(bins = 100) + geom_vline(xintercept = t0, color = "red", linetype = "dashed") + ggtitle(paste("T, 귀무가설 재표집. p값 = ",p.t0)) 
g2 = ggplot(mapping = aes(x = t.boot, y = ..density..)) + geom_histogram(bins = 100) + geom_vline(xintercept = t0, color = "red", linetype = "dashed") + ggtitle(paste("T, 피벗 재표집. p값 = ",p.tpivot)) 
g3 = ggplot(mapping = aes(x = m.boot0, y = ..density..*100)) + geom_histogram(bins = 100) + geom_vline(xintercept = m0, color = "red", linetype = "dashed") + ggtitle(paste("M, 귀무가설 재표집. p값 = ",p.m0)) 
g4 = ggplot(mapping = aes(x = m.boot, y = ..density..*100)) + geom_histogram(bins = 100) + geom_vline(xintercept = m0, color = "red", linetype = "dashed") + ggtitle(paste("M, 피벗 재표집. p값 = ",p.mpivot)) 
plot_grid(g1,g2,g3, g4, ncol = 2)

ggsave(filename = "ch9_3.pdf", path = 'images', width = 6.5, height = 6)
ggsave(filename = "ch9_3.png", path = 'images', width = 6.5, height = 6)


# nonparametric bootstrap test 1 : ex 1 ------------------------------------------

t0 = (mean(y) - mean(x)) / pooled.se(x,y)
set.seed(1)
B = 3999
boot.out = replicate(B, {
  boot.id = sample(1:(m+n),(m+n),replace = TRUE) 
  x.boot = c(x,y)[boot.id[1:m]]
  y.boot = c(x,y)[boot.id[(m+1):(m+n)]]
  (mean(y.boot) - mean(x.boot)) / pooled.se(x.boot,y.boot)
})
boot.p = 2* min( c( mean(c(t0,boot.out) >=t0), 
                    mean(c(t0,boot.out) <=t0)))
boot.p


perm.out = replicate(B, {
  perm.id = sample(1:(m+n),(m+n),replace = FALSE) 
  x.perm = c(x,y)[perm.id[1:m]]
  y.perm = c(x,y)[perm.id[(m+1):(m+n)]]
  (mean(y.perm) - mean(x.perm)) / pooled.se(x.perm,y.perm)
})
perm.p = 2* min( c( mean(c(t0,perm.out) >=t0), 
                    mean(c(t0,perm.out) <=t0)))
perm.p



# ex 2 bootstrap two sample test by pivot resampling ---------------------------

## two group Efron and Tibs 1994 method
d0 = mean(y) - mean(x)
t0 = (mean(y) - mean(x)) / pooled.se(x,y)
set.seed(1)
B = 3999
boot.out = replicate(B, {
  x.boot = sample(x, m, replace = TRUE)
  y.boot = sample(y, n, replace = TRUE)
  t.boot = (mean(y.boot)  - mean(x.boot) - d0) / pooled.se(x,y)
})
boot.p = mean(c(t0,boot.out) >= t0)
boot.p


# ex 3 --------------------------------------------------------------------

# one sample 
x <- leukemia$ge[leukemia$label == "AML"]
n = length(x)

# hypothesis H0: mu =mu0, H1: mu > mu0
mu0 = 166

# test statistic
hat.mu = mean(x)
t0 = sqrt(n) * (mean(x) - mu0 ) / sd(x)

# nonparametric bootstrap
set.seed(1)
boot.out = replicate(3999, {
  x.boot = sample(x, n, replace = TRUE)
  sqrt(n) * (mean(x.boot) - hat.mu) / sd(x.boot)
})
boot.p = mean(c(t0,boot.out) >= t0)
boot.p


