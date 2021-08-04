# nonparametric book by S. Jung
# chapter 7. bootstrap
#
#
# ----------------
library(tidyverse) 
library(cowplot)

library(showtext)
#font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()

# flights example ---------------------------------------------------------
# https://www.kaggle.com/usdot/flight-delays
load("data/flights2.RData")
set.seed(1)
small = flights2 %>% slice_sample(n = 5000) 
g1 = small %>% 
  ggplot(aes(x= SCHEDULED_DEPARTURE, y = DEPARTURE_DELAY))+
  geom_point()
g2 = small %>% ggplot(aes(x = DEPARTURE_DELAY)) + 
  geom_density() + geom_rug()
plot_grid(g1,g2)

ggsave(filename = "ch7_1.pdf", path = 'images', width = 6.5, height = 3)
ggsave(filename = "ch7_1.png", path = 'images', width = 6.5, height = 3) # two-panel


summary(flights2$DEPARTURE_DELAY)


# sampling distribution (almost exact) ------------------------------------
delay = na.omit(flights2$DEPARTURE_DELAY)
length(delay)

set.seed(2)
N = 10000
n = 100
meds = replicate(N, median(sample(delay, n, replace = TRUE)))
sd(meds)
median(meds) - median(delay)
mean(meds) - median(delay)
library(latex2exp)
data.frame(MEDs = meds) %>% 
  ggplot() + 
  geom_histogram(aes(x = MEDs, y= ..density..),binwidth  = .5, fill = "grey", color = "white") +
  geom_rug(aes(x = MEDs)) + ggtitle(TeX("$\\hat{\\theta}_n$의 표본분포"))
ggsave(filename = "ch7_2.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch7_2.png", path = 'images', width = 6.5, height = 4)

data.frame(MEDs = meds) %>% 
  ggplot() + 
  geom_density(aes(x = MEDs, y= ..density..),bandwidth  = 1, fill = "grey", color = "white") + 
  geom_rug()

sd(meds)
median(sample(delay,n))

set.seed(1)
n = 10
dat = sample(delay, n)
data.frame(delay = dat) %>% 
  ggplot() + 
  geom_histogram(aes(x = delay, y = ..density..), binwidth = 5, fill = "grey", color = "white") +
  geom_rug(aes(x = delay)) + ggtitle("관측값, n = 10")
dat
(bootstrap.dat.1 = sample(dat,n, replace = TRUE))
(bootstrap.dat.2 = sample(dat,n, replace = TRUE))
c( median(dat), median(bootstrap.dat.1), median(bootstrap.dat.2) ) 

n = 100
dat = sample(delay, n)
B = 2000
med.boot = replicate(B, median( sample(dat,n, replace = TRUE)))
data.frame(MEDs = med.boot) %>% 
  ggplot() + 
  geom_histogram(aes(x = MEDs, y = ..density..), binwidth = .5, fill = "grey", color = "white") +
  geom_rug(aes(x = MEDs)) + ggtitle(TeX("$\\hat{\\theta}_{100}$의 붓스트랩 표본분포"))
ggsave(filename = "ch7_3.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch7_3.png", path = 'images', width = 6.5, height = 4)

sd(med.boot)
(bias = mean(med.boot) - median(dat))



# windspeed example -------------------------------------------------------
set.seed(0)
v.cat = 0:9 + (0.5)
count = c(6421,14246,11953,7463,2322,659,120,11,3,0)
names(count) = v.cat
dat = DescTools::Untable( as.table(count))
v = as.numeric(as.character(dat$Var1))
n = length(v)

# formula -----------------------------------------------------------------
# estimate, and theoretical se 
pv = 1 - 2 * v + v^2 - 0.05 * v^3
theta.hat = mean(pv)
se = sd(pv)/sqrt(n)

# nonparametric bootstrap estimate of se ----------------------------------
B = 3999
set.seed(1)
theta.hat.boot = replicate(B,mean(sample(pv,n,replace = T)))
se.boot = sd(theta.hat.boot)
bias.boot = mean(theta.hat.boot) - theta.hat


# Weilbull distribution ---------------------------------------------------
library(weibullness)
mle=weibull.mle(v, threshold = 0, interval = c(0,3))
xx = seq(0,10,by = 0.1)
weibuldensity = data.frame(x = xx, density = dweibull(xx, mle$shape, mle$scale))

ggplot(data.frame(v = v.cat, Freq = count/ length(v)), aes(v,Freq)) + geom_col() + 
  geom_line(data = weibuldensity, aes(x = xx, y = density)) + ggtitle("Windspead @Seoul, March, 2001-2011.")
ggsave(filename = "ch7_4.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch7_4.png", path = 'images', width = 6.5, height = 4)
data.frame(v = v.cat, Freq = count)
mle

# parametric bootstrap ----------------------------------------------------
mle=weibull.mle(v, threshold = 0, interval = c(0,3))
k.hat = mle$shape
lambda.hat = mle$scale

# theta.hat
theta.hat = mean(1 - 2 * v + v^2 - 0.05 * v^3)

set.seed(1)
theta.boot = replicate(B,{
  v.boot = rweibull(n, k.hat, lambda.hat)
  pv.boot = 1 - 2 * v.boot + v.boot^2 - 0.05 * v.boot^3
  mean(pv.boot)
}) 
head(theta.boot)

#theta.hat.mle = 1 - 2 * lambda.hat * gamma(1 + 1 / k.hat) + 
#  lambda.hat^2 * gamma(1 + 2 / k.hat) - 
#  0.05 * lambda.hat^3 * gamma(1 + 3 / k.hat) # This is the parameter for f( ; hat/eta)

se.boot.par = sd(theta.boot)
bias.boot.par = mean(theta.boot) - theta.hat

matrix(c(se, NA, 
         se.boot,bias.boot,
         se.boot.par, bias.boot.par),nrow = 2)
 


# ex: 7.8 two-sample -------------------------------------------------------------


leukemia4 = read.csv("data/leukemia4.csv")
ALL = leukemia4$ge[leukemia4$label == "ALL"]
AML = leukemia4$ge[leukemia4$label == "AML"]
n1 = length(ALL)
n2 = length(AML)

md = mean(ALL) - mean(AML)
set.seed(1)
B = 2000
md.boot = replicate(B, 
           mean(ALL[sample(1:n1, n1, replace = TRUE)]) - 
           mean(AML[sample(1:n2, n2, replace = TRUE)])
           )
data.frame(md.boot = md.boot) %>% ggplot(aes(x = md.boot)) + geom_histogram(aes(y = ..density..), color = "white")
ggsave(filename = "ch7_5.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch7_5.png", path = 'images', width = 6.5, height = 4)
sd(md.boot)


# ex 7.9 bivariate ---------------------------------------------------------------
set.seed(1)
wh = read.csv("data/whappiness.csv")
score = wh$Score; gdp = wh$GDP.per.capita
(k.obs = cor(gdp,score, method = "k"))

# random permutation test
M = 1999
n = length(gdp)
k.perm = replicate(M, {
  id.perm = sample(1:n, n, replace = FALSE)
  cor(gdp,score[id.perm], method = "k")
})
(pvalue.perm = mean(c(k.obs,k.perm) >= k.obs))

# 붓스트랩 
B = 1999
k.boot = replicate(B, {
  id.boot = sample(1:n,n, replace = TRUE)
  cor(gdp[id.boot], score[id.boot], method = "k")
})
sd(k.boot)
data.frame( k = c(k.boot,k.perm), 
            type = c(rep("bootstrap",B), rep("permutation",M))) %>% 
  ggplot(aes(k, color = type, fill = type)) + 
  geom_density(alpha = 0.15) + 
  xlab("Kendall's tau") + 
  geom_vline(xintercept = k.obs, color = "red",linetype = "dashed")
ggsave(filename = "ch7_6.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch7_6.png", path = 'images', width = 6.5, height = 4)


