# nonparametric book by S. Jung
# chapter 10. Part II: multiple linear regression
#
#
# ----------------
library(tidyverse) 
library(cowplot)
library(GGally)
library(latex2exp)

library(showtext)
#font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()


# ex 10.7 --------------------------------------------------------------------
h = read.csv("data/worldbank2020.csv")
dat = h %>% dplyr::select(4:9,3)
colnames(dat) = c("GDP", "Social", 
                  "Life", "Freedom", "Generosity", "Corruption","Score")
ggpairs(dat)
ggsave(filename = "ch10a_1.pdf", path = 'images', width = 6.5, height = 6.5)
ggsave(filename = "ch10a_1.png", path = 'images', width = 6.5, height = 6.5)

 

# ex 10.8 -------------------------------------------------------------
h = read.csv("data/worldbank2020.csv")
dat = h %>% dplyr::select(4:9,3)
colnames(dat) = c("GDP", "Social", 
                  "Life", "Freedom", "Generosity", "Corruption","Score")
n = nrow(dat); p = 6
lm.out = lm(Score ~ . , data =dat)

y = dat$Score
X = model.matrix(lm.out)
head(X)

s = sqrt( sum ( (lm.out$residuals)^2 ) / (n - p - 1) ) 
(se.hat = s *  ( solve(t(X) %*% X) %>% diag() %>% sqrt() ) )
beta.hat = solve(t(X) %*% X) %*% t(X) %*% dat$Score
beta.hat[6] + se.hat[6] * qt(1-0.025, n-p-1) * c(-1,1)

summary(lm.out)


# ex 10.9 -------------------------------------------------------------
## reduced model vs full model, theoretical F-test  
X0 = X[,1:5] 
X1 = X[,6:7]
q = 2
Hat = X %*% solve(t(X) %*% X) %*% t(X)
Hat0 = X0 %*% solve(t(X0) %*% X0) %*% t(X0)
SSR = t(y) %*% (Hat - Hat0) %*% y 
SSE = sum(y^2) - t(y) %*% Hat %*% y 
(F.obs = drop(SSR/SSE * (n-p-1) / q))
(p.theory = 1- pf(F.obs, q, n-p-1))
ggplot(mapping = aes(x = seq(0,10,by = 0.1), y = df(seq(0,10,by = 0.1), q, n-p-1))) +
  geom_line() + geom_vline(xintercept = F.obs, color = "red",linetype = "dashed") + 
  ggtitle(TeX("관측된 F와 그 영분포 $F_{2, 146}$")) + xlab("F") + ylab("density")

ggsave(filename = "ch10a_2.pdf", path = 'images', width = 4.8, height = 3)
ggsave(filename = "ch10a_2.png", path = 'images', width = 4.8, height = 3)




# ex 10.10 -------------------------------------------------------------
# permutation tests (three versions) 
X10 = X1 - Hat0 %*% X1
XX10 = t(X10) %*% X10
XXX10 = solve( t(X10) %*% X10) %*% t(X10)
F.pivot = function(y.perm, beta1){
  beta1.hat = XXX10 %*% y.perm
  SSE = sum(y.perm^2) - t(y.perm) %*% Hat %*% y.perm 
  SSR = t(beta1.hat - beta1) %*% XX10 %*% (beta1.hat - beta1)
  as.vector(SSR / SSE / q * (n-p-1))
}
(F.obs = F.pivot(y, c(0,0)))

## permuting observations
M = 3999
set.seed(1)
F.perm.raw = replicate(M, {
  yperm = sample(y,n, replace = FALSE)
  F.pivot(yperm, c(0,0))
})
#ggplot(mapping = aes(x = F.perm.raw)) + geom_histogram() + geom_vline(xintercept = F.obs, color = "red")  
(p.raw = mean(c(F.obs,F.perm.raw) >= F.obs))

## permuting residuals under the null 
yhat0 = Hat0 %*% y
res0 = y - yhat0
F.perm.null = replicate(M, {
  yperm = yhat0 + sample(res0,n, replace = FALSE)
  F.pivot(yperm, c(0,0))
})
#ggplot(mapping = aes(x = F.perm.null)) + geom_histogram() + geom_vline(xintercept = F.obs, color = "red")  
(p.null = mean(c(F.obs,F.perm.null) >= F.obs))

## permuting residuals under the full model 
yhat = Hat %*% y
res = y - yhat
beta.hat1 = XXX10 %*% y
F.perm.full = replicate(M, {
  yperm = yhat + sample(res,n, replace = FALSE)
  F.pivot(yperm, beta.hat1)
})
#ggplot(mapping = aes(x = F.perm.full)) + geom_histogram() + geom_vline(xintercept = F.obs, color = "red")  
(p.full = mean(c(F.obs,F.perm.full) >= F.obs))
c(p.theory,p.raw, p.null, p.full)

g1 = ggplot(mapping = aes(x = F.perm.raw)) + 
  geom_histogram(aes(y = ..density..),bins = 99) + 
  geom_vline(xintercept = F.obs, color = "red",linetype = "dashed")  + 
  geom_line(mapping = aes(x = seq(0,10,by = 0.1), y = df(seq(0,10,by = 0.1), q, n-p-1)),color = "blue") + 
  ggtitle("관측값 뒤섞기")+ scale_y_continuous(labels=NULL)

g2 = ggplot(mapping = aes(x = F.perm.null)) + 
  geom_histogram(aes(y = ..density..),bins = 99) + 
  geom_vline(xintercept = F.obs, color = "red",linetype = "dashed")  + 
  geom_line(mapping = aes(x = seq(0,10,by = 0.1), y = df(seq(0,10,by = 0.1), q, n-p-1)),color = "blue") + 
  ggtitle("귀무가설 잔차 뒤섞기")+ scale_y_continuous(labels=NULL)

g3 = ggplot(mapping = aes(x = F.perm.full)) + 
  geom_histogram(aes(y = ..density..),bins = 99) + 
  geom_vline(xintercept = F.obs, color = "red",linetype = "dashed")  + 
  geom_line(mapping = aes(x = seq(0,10,by = 0.1), y = df(seq(0,10,by = 0.1), q, n-p-1)),color = "blue") + 
  ggtitle("완전모형 잔차 뒤섞기") + scale_y_continuous(labels=NULL)

plot_grid(g1,g2,g3, nrow = 1)

ggsave(filename = "ch10a_3.pdf", path = 'images', width = 6.5, height = 3)
ggsave(filename = "ch10a_3.png", path = 'images', width = 6.5, height = 3)



# ex 10.11 -------------------------------------------------------------
## bootstrap-t interval
XX = solve(t(X) %*% X) 
Hat = X %*% XX %*% t(X)
t.pivot = function(yboot, beta){
  yhat = Hat %*% yboot
  res = yboot - yhat 
  beta.hat  = XX %*% t(X) %*% yboot
  se.beta.hat = sqrt( sum(res^2) / (n-p-1) * diag(XX))
  as.vector( (beta.hat - beta)  / se.beta.hat ) 
}

yhat = lm.out$fitted.values
adjusted.res = lm.out$residuals / sqrt(1 - hatvalues(lm.out))
beta.hat = lm.out$coefficients
set.seed(1)
B = 1999
t.boot.res = replicate(B, {
  yboot = yhat + sample(adjusted.res, n, replace = TRUE)
  t.pivot(yboot, beta.hat)
}) %>% t()

ci.beta5 = beta.hat[6] - se.hat[6] * quantile(t.boot.res[,6], probs = c(1-0.025,0.025))
ci.beta5



# ex 10.12  -------------------------------------------------------------
## bootstrap test

t5 = t.pivot(y, rep(0,p+1))[6]
(pvalue = mean(abs( c(t5, t.boot.res[,6])) >= abs(t5)))
# 귀무가설 H_0: \beta5 = 0 기각실패


ci.beta5 # \beta5의 신뢰구간이 0을 포함

alphavec = seq(0,1,length.out = B)
CImat = matrix(data = NA, nrow = B, ncol = 2)
for (ii in 1:B){
  alpha = alphavec[ii]
  CImat[ii,] = beta.hat[6] - se.hat[6] * quantile(t.boot.res[,6], probs = c(1-alpha/2,alpha/2))
}
ci.dat = data.frame(level = 1- alphavec, low = CImat[,1], upp = CImat[,2])
ci.dat %>% ggplot(aes(x = level)) + geom_ribbon(aes(ymin = low, ymax = upp))



alpha = pvalue
beta.hat[6] - se.hat[6] * quantile(t.boot.res[,6], probs = c(1-alpha/2,alpha/2))

alpha = pvalue + 0.005
beta.hat[6] - se.hat[6] * quantile(t.boot.res[,6], probs = c(1-alpha/2,alpha/2))

