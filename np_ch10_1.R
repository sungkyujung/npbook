# nonparametric book by S. Jung
# chapter 10. Part I: simple regression
#
#
# ----------------
library(tidyverse) 
library(cowplot)
library(latex2exp)

library(showtext)
#font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()


# ex 1 --------------------------------------------------------------------
h = read.csv("data/worldbank2020.csv")
lm(Score ~ GDP.per.capita, data = h) %>% summary()

h %>% ggplot(aes(x = GDP.per.capita, y = Score)) + geom_point()
ggsave(filename = "ch10_1.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch10_1.png", path = 'images', width = 6.5, height = 4)


# ex 10.2 --------------------------------------------------------------------
h = read.csv("data/worldbank2020.csv")
y = h$Score
x = h$GDP.per.capita


(lm.out = lm(y ~ x))
adjusted.res = lm.out$residuals / sqrt(1 - hatvalues(lm.out))
#plot(lm.out)
y = h$Score
x = h$GDP.per.capita
n = length(x)
Sxx = sum( (x - mean(x))^2 ) 
Sxy = sum( (x - mean(x)) * (y - mean(y)) ) 
beta.hat = Sxy / Sxx
res = y - mean(y) - beta.hat * ( x - mean(x))
s = sqrt( sum( res^2 ) / (n - 2) ) 
se.beta.hat = s / sqrt(Sxx)
t = beta.hat / se.beta.hat
g1 = h %>% ggplot(mapping = aes(x = GDP.per.capita, y = adjusted.res)) + geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = c(-1,1)*2*s, linetype = "dashed", color = "red")
g2 = ggplot(mapping = aes(x = adjusted.res)) + geom_histogram(color = "white", bins = 11) + geom_rug()
plot_grid(g1,g2)
ggsave(filename = "ch10_2.pdf", path = 'images', width = 6.5, height = 3)
ggsave(filename = "ch10_2.png", path = 'images', width = 6.5, height = 3)


# ex 10.3 --------------------------------------------------------------------
y = h$Score
x = h$GDP.per.capita
n = length(x)
Sxx = sum( (x - mean(x))^2 ) 
Sxy = sum( (x - mean(x)) * (y - mean(y)) ) 
beta.hat = Sxy / Sxx

M = 1999
beta.hat.perm = replicate(M, {
  yperm = sample(y,n, replace = FALSE)
  sum( (x - mean(x)) * (yperm - mean(yperm))) / Sxx 
})
(p = mean( c(beta.hat, beta.hat.perm) >= beta.hat ))
ggplot(mapping = aes(x = beta.hat.perm)) + geom_histogram(color = "white", bins = 30) + 
  geom_rug() + geom_vline(xintercept = beta.hat, color = "red",linetype = "dashed") 
ggsave(filename = "ch10_3.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch10_3.png", path = 'images', width = 6.5, height = 4)


# ex 10.4  -------------------------------------------------------------
# bootstrap CI
t.pivot = function(x,y, beta0){
  n = length(x)
  Sxx = sum( (x - mean(x))^2 ) 
  Sxy = sum( (x - mean(x)) * (y - mean(y)) ) 
  beta.hat = Sxy / Sxx
  res = y - mean(y) - beta.hat * ( x - mean(x))
  se.beta.hat = sqrt( sum( res^2 ) / (n - 2) ) / sqrt(Sxx)
  t = (beta.hat - beta0)  / se.beta.hat
  t
}
set.seed(1)
B = 1999
boot.t = replicate(B, {
  boot.id = sample(1:n,n, replace = TRUE) # case resampling
  xboot = x[boot.id]
  yboot = y[boot.id]
  t.pivot(xboot, yboot, beta.hat)
})

res = y - mean(y) - beta.hat * ( x - mean(x))
s = sqrt( sum( res^2 ) / (n - 2) ) 
se.beta.hat = s / sqrt(Sxx)

boot.t.ci = beta.hat - se.beta.hat * quantile(boot.t, probs = c(1-0.025, 0.025))
boot.t.ci
(norm.appr.ci = beta.hat + se.beta.hat * qt(c(0.025, 1-0.025), df = n-2))


# 잔차 붓스트랩

(lm.out = lm(y ~ x))
adjusted.res = lm.out$residuals / sqrt(1 - hatvalues(lm.out))
set.seed(1)
boot.t.res = replicate(B, {
  boot.res = sample(adjusted.res,n, replace = TRUE) # residual resampling
  yboot = lm.out$fitted.values + boot.res
  t.pivot(x,yboot, beta.hat)
}) 
(res.boot.t.ci = beta.hat - se.beta.hat * quantile(boot.t.res, probs = c(1-0.025, 0.025)))


beta.hat + se.beta.hat * qt(c(0.025, 1-0.025), df = n-2)
#beta.hat - se.beta.hat * qt(c(1-0.025, 0.025), df = n-2)



# ex 10.5 & ex 10.6  -------------------------------------------------------------
# bootstrap tests 
# 1. 관측값  붓스트랩, 귀무가설에서  재표집 
set.seed(1)

t0 = t.pivot(x,y,0)
B = 1999
boot.t = replicate(B, {
  xboot = sample(x,n, replace = TRUE)
  yboot = sample(y,n, replace = TRUE)
  t.pivot(xboot,yboot, 0)
})
p.null.case = mean(c(t0,boot.t ) >= t0) 
boot.t.null.case = boot.t
g1 = ggplot(mapping = aes(x = boot.t.null.case)) + geom_histogram(color = "white", bins = 30) + 
  geom_rug() + ggtitle(paste("귀무가설 재표집, 관측값 붓스트랩")) 


# 2. 잔차 붓스트랩, 귀무가설에서 재표집
adjusted.res = (y - mean(y))/ sqrt( 1- 1/n)
boot.t = replicate(B, {
  boot.res = sample(adjusted.res, n, replace = TRUE)
  yboot = mean(y) + boot.res
  t.pivot(x,yboot, 0)
})
p.null.res = mean(c(t0,boot.t ) >= t0) 
c(p.null.case, p.null.res)

boot.t.null.res = boot.t
g2 = ggplot(mapping = aes(x = boot.t.null.res)) + geom_histogram(color = "white", bins = 30) + 
  geom_rug() + ggtitle(paste("귀무가설 재표집, 잔차 붓스트랩")) 


# 3. 피벗-관측값 붓스트랩 
boot.t = replicate(B, {
  boot.id = sample(1:n,n, replace = TRUE)
  xboot = x[boot.id]
  yboot = y[boot.id]
  t.pivot(xboot,yboot, beta.hat)
}) 
p.pivot.case = mean(c(t0,boot.t ) >= t0)
boot.t.pivot.case = boot.t

g3 = ggplot(mapping = aes(x = boot.t.pivot.case)) + geom_histogram(color = "white", bins = 30) + 
  geom_rug() + ggtitle(paste("피벗 재표집, 관측값 붓스트랩")) 


# 4. 피벗-잔차 붓스트랩
(lm.out = lm(y ~ x))
adjusted.res = lm.out$residuals / sqrt(1 - hatvalues(lm.out))
boot.t = replicate(B, {
  boot.res = sample(adjusted.res,n, replace = TRUE)
  yboot = lm.out$fitted.values + boot.res
  t.pivot(x,yboot, beta.hat)
}) 
p.pivot.res = mean(c(t0,boot.t ) >= t0) 
boot.t.pivot.res = boot.t

g4 = ggplot(mapping = aes(x = boot.t.pivot.res)) + geom_histogram(color = "white", bins = 30) + 
  geom_rug() + ggtitle(paste("피벗 재표집, 잔차 붓스트랩")) 


plot_grid(g1,g2,g3,g4)
ggsave(filename = "ch10_4.pdf", path = 'images', width = 6.5, height = 6)
ggsave(filename = "ch10_4.png", path = 'images', width = 6.5, height = 6)
