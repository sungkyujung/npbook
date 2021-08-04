# nonparametric book by S. Jung
# chapter 2. one sample
#
#
# ----------------
library(tidyverse) 
library(gtools)
library(Rfit)
library(showtext)
#font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()

# example 2.1 -------------------------------------------------------------
# artificial data
iq = c(98, 121, 110, 89, 109, 108, 102, 92, 131, 114)
sort(iq)
mean(iq)
n = length(iq)
z = (mean(iq) - 100) / (10/ length(iq))
1-pnorm(z)


# example 2.2 -------------------------------------------------------------
# artificial data
x = c(1200, 2300, 5800, 1100, 1200, 1700, 1400, 1900, 3800, 1300)

data.frame(x = x) %>% 
  ggplot(aes(x = x)) +  
  geom_density() +
  geom_rug(length = unit(0.05, "npc")) + geom_vline(xintercept = 1500, color = "red", linetype = "dashed")
ggsave(filename = "ch2_sodium.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch2_sodium.png", path = 'images', width = 6.5, height = 4)


# Figure 2.2 --------------------------------------------------------------
# lost ... 


# Example 2.3 -------------------------------------------------------------
iq = c(98, 121, 110, 89, 109, 108, 102, 92, 131, 114)
ifelse(iq > 100, "+","-")
n = length(iq)
Sobs = sum(iq > 100)
plot(0:n,dbinom(0:n,n,1/2),type='l')
(pvalue = 1-pbinom(Sobs-1,n,1/2)) 
sum(dbinom(Sobs:n,n,1/2))


# example 2.4 -------------------------------------------------------------
S0 = 9
1 - pbinom(S0-1, n, 1/2) # P(S \ge S0) = 1 - P(S < S0) = 1 - P(S le S0 -1)
S0 = 8
1 - pbinom(S0-1, n, 1/2) # P(S \ge S0) 


# example 2.5 -------------------------------------------------------------
Sz = (Sobs - n/2) / sqrt(0.25*n)
1-pnorm(Sz)


# example 2.6 & 2.7-------------------------------------------------------------
sum(dbinom(8:n,n,1/2)) # type I error of sign test
1 - pnorm( ( -5 + qnorm(1-0.05)*sqrt(10) ) / sqrt(10)) # power of z-test
p1 = 1- pnorm(-5/10)  # prob. of "+"
sum(dbinom(8:n,n,p1)) # power of sign test

# H1: X ~ Laplace(105, sqrt(50))
library(extraDistr)
# power of z-test(approx)
1 - pnorm( ( -5 + qnorm(1-0.05)*sqrt(10) ) / sqrt(10))
# power of sign test 
# p(X > 100)
(p1 <- 1 - plaplace(100, mu = 105, sigma = sqrt(50)))
sum(dbinom(8:n,n,p1))



# example 2.8 -------------------------------------------------------------
# sodium data; CI.
x = sort(c(1200, 2300, 5800, 1100, 1200, 1700, 1400, 1900, 3800, 1300))
n = length(x)
alpha = 0.05
l = qbinom(alpha/2, size = n, prob = 1/2)
u = n + 1 -  l
(ci = x[c(l,u)])
(coverage = sum(dbinom(l:(u-1), n, prob = 1/2)) )
(coverage >= 1-alpha)
#sum(dbinom(l:(u-1), n, prob = 1/2)) >= 1- alpha
#sum(dbinom(0:(l-1),n,1/2))
#sum(dbinom(u:n,n,1/2))


# example 2.9 & 2.10 ------------------------------------------------------------
x = c(11,20,16,5,3,17); theta0 = 15
D = x - theta0 # 편차
sign(D)   # 부호
rank(abs(D)) # 절대편차의 순위
(SR = sum( ( sign(D)>0 ) * rank(abs(D)) ) )
n = length(x)
c = 1:(n*(n+1)/2)
plot(c,dsignrank(c,n),type = 'l') 
(pval.upper = psignrank(SR,n, lower.tail = FALSE))
(pval.lower = psignrank(SR,n, lower.tail = TRUE))
(pval.both = 2 * psignrank(SR,n,lower.tail = TRUE))
(pval.both = min( 2 * psignrank(SR,n,lower.tail = TRUE), 
                  2 * psignrank(SR,n, lower.tail = FALSE) 
))
alpha = 0.05
psignrank(18, n, lower.tail = F)
qsignrank(1-alpha,n)
psignrank( qsignrank(1-alpha,n) ,n)


# fig 2.3 -----------------------------------------------------------------
# null distribution of signrank, n = 6
# obtained by inspecting all permutations

dev = 1:6 # there is no tie; actual data values do not matter here
ranks = rank(dev)
n = length(dev)
allsignflips = gtools::permutations(2, n, v = c(1,-1), repeats.allowed = TRUE)
nperm = nrow(allsignflips) # 2^n 과 같다.
SRvec = vector(length = nperm) # 각 경우에서 순위합 계산
for (i in 1:nperm){
  SRvec[i] = sum(ranks[allsignflips[i,] > 0])
}

N = n*(n+1)/2 # 순위합이 같은 경우의 갯수가 곧 그 순위합의 빈도이다.
dsignrank.my = vector(length = N+1)
for (sr in 0:N){
  dsignrank.my[sr+1] = sum(SRvec == sr) / nperm
}
           
# resampling methods (with or without ties)
toy.data = c(-3,-2,-2,0,1,2,4,8,11,11); theta0 = 0
dev = x - 0  
reduced.dev = dev[sign(dev) != 0 ]
(ranks = rank(abs(reduced.dev))) 
n = length(ranks)
# 모든 부호의 경우를 gtools 패키지의 permutations() 함수로 생성한다. 아래 표현은 길이가 n인 벡터를 반복을 허용하면서 v = c(1,-1)의 원소들로 채우는 모든 경우
allsignflips = gtools::permutations(2, n, v = c(1,-1), repeats.allowed = TRUE)
nperm = nrow(allsignflips) # 2^n 과 같다.

SRvec = vector(length = nperm) # 각 경우에서 순위합 계산
for (i in 1:nperm){
  SRvec[i] = sum(ranks[allsignflips[i,] > 0])
}

N = n*(n+1)/2 # 순위합이 같은 경우의 갯수가 곧 그 순위합의 빈도이다.
fSR = vector(length = N+1)
for (sr in 0:N){
  fSR[sr+1] = sum(SRvec == sr) / nperm
} 
data.frame(SR = 0:N, fSR) %>% ggplot(aes(x = SR, y = fSR)) + geom_col() + ggtitle("SR_+ null distribution, n = 6, no tie")
ggsave(filename = "ch2_dsignrank.pdf", path = 'images', width = 4.8, height = 3)
ggsave(filename = "ch2_dsignrank.png", path = 'images', width = 4.8, height = 3)
 


# example 2.11 ------------------------------------------------------------

x.before = c(1.83, 0.5, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.30, 0.39, 0.94, 1.34)
x.after = c(0.88, 0.65, 0.59, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29, 0.43, 0.57, 0.92)
x = x.before - x.after 
D = x - 0    # deviation
(SR = sum( ( sign(D)>0 ) * rank(abs(D)) ))
n = length(x)
Z_SR = (SR - n*(n+1)/4)/ sqrt(n*(n+1)*(2*n+1)/24)
(p = 1-pnorm(Z_SR))


# example 2.12 -2.14 ------------------------------------------------------------
# resampling methods (with or without ties)
toy.data = c(-3,-2,-2,0,1,2,4,8,11,11); theta0 = 0
dev = toy.data - 0  
reduced.dev = dev[sign(dev) != 0 ]
(ranks = rank(abs(reduced.dev))) 
n = length(ranks)
# 모든 부호의 경우를 gtools 패키지의 permutations() 함수로 생성한다. 아래 표현은 길이가 n인 벡터를 반복을 허용하면서 v = c(1,-1)의 원소들로 채우는 모든 경우
allsignflips = gtools::permutations(2, n, v = c(1,-1), repeats.allowed = TRUE)
nperm = nrow(allsignflips) # 2^n 과 같다.

SRvec = vector(length = nperm) # 각 경우에서 순위합 계산
for (i in 1:nperm){ 
  SRvec[i] = sum(ranks[allsignflips[i,] > 0])
}

N = n*(n+1)/2 # 순위합이 같은 경우의 갯수가 곧 그 순위합의 빈도이다.
dsignrank.tie = vector(length = N+1)
for (sr in 0:N){
  dsignrank.tie[sr+1] = sum(SRvec == sr) / nperm
}
data.frame(SR = 0:N, SRdensiy = dsignrank.tie) %>% head(3)

data.frame(SR = 0:N, SRdensity = dsignrank.tie) %>% ggplot(aes(x = SR, y = SRdensity)) + geom_col()
dsignrank(0:2,n)


# example 2.15 ------------------------------------------------------------
x = x.before - x.after  
n = length(x)
w = sort(Rfit::walsh(x)) # sorted walsh averages 
median(w) # 
alpha = 0.05
u = n*(n+1)/2 + 1 - qsignrank(alpha/2,n)
l = qsignrank(alpha/2,n)
c(l,u)
(ci = c(w[l], w[u]))


# wilcox.test() -----------------------------------------------------------
wilcox.test(x, alternative = "greater", exact = F, correct =F)
wilcox.test(x, alternative = "greater", exact = T)
wilcox.test(x, conf.int = T)

