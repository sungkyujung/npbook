# nonparametric book by S. Jung
# chapter 8. bootstrap confidence intervals
#
#
# ----------------
library(tidyverse) 
library(cowplot)
library(latex2exp)

library(showtext)
#font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()


# Game data ---------------------------------------------------------------
# Originally from https://www.kaggle.com/tamber/steam-video-games
# Preprocessed to extract data regarding TF2 
TFtime = read.csv("data/TeamFortress2.csv")$x 
qplot(TFtime) + geom_histogram(bins = 50) +  
  xlab("time (hours)") + ggtitle("Team Fortress 2 play time at Steam")
ggsave(filename = "ch8_1.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch8_1.png", path = 'images', width = 6.5, height = 4)


# ex 2 --------------------------------------------------------------------
TFtime = read.csv("data/TeamFortress2.csv")$x 
(hat.theta = mean(TFtime))
n = length(TFtime) 
B = 1999
set.seed(1)
theta.boot = replicate(B, mean(sample(TFtime, n, replace = TRUE)))
(se.est = sd(theta.boot))
(bias.est = mean(theta.boot) - hat.theta)
(norm.app.boot.CI = hat.theta - bias.est + se.est * qnorm(0.025) *c(1,-1))


# ex 2 figure  ------------------------------------------------------------

data.frame( theta.boot = theta.boot) %>% 
  ggplot(aes(theta.boot)) + geom_histogram(aes(y = ..density..), color = "white",bins = 100) + xlab("time") + 
  geom_vline(xintercept = hat.theta, color = "red", linetype = "dashed")

ggsave(filename = "ch8_2.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch8_2.png", path = 'images', width = 6.5, height = 4)

pdf(file="images/ch8_2qq.pdf",width = 6.5, height = 4)
car::qqPlot(theta.boot)
dev.off() 


# ex 8.3 --------------------------------------------------------------------

set.seed(1) 
boot.out = replicate(B, {
  TFtime.b = sample(TFtime, n, replace = TRUE)
  c(mean(TFtime.b), sd(TFtime.b)/sqrt(n)) 
})
boot.out = t(boot.out); colnames(boot.out) = c("theta*","se*")
head(boot.out)
hat.theta = mean(TFtime)
t.pivot = sort( (boot.out[,"theta*"] - hat.theta) / boot.out[,"se*"] )


se.hat.theta = sd(TFtime) / sqrt(n)
(boot.t.ci = hat.theta - se.hat.theta * c(t.pivot[1950], t.pivot[50])) 

hat.theta - se.hat.theta * quantile(t.pivot, probs= c(1-0.025,0.025),type = 1)

boot.dist = sort(boot.out[,"theta*"])
(basic.ci = 2*hat.theta - c(boot.dist[1950], boot.dist[50]))
 



# ex 8.4 --------------------------------------------------------------------
# studentized, second level bootstrap
B = 1999
C = 25
boot.out = replicate(B, {
  TFtime.b = sample(TFtime, n, replace = TRUE)
  theta.boot = mean(TFtime.b)
  theta.boot.se = sd( replicate(C, mean (sample(TFtime.b, n, replace = TRUE)) )) 
  c(theta.boot, theta.boot.se)
}) 
boot.out = t(boot.out); colnames(boot.out) = c("theta*","se*")
head(boot.out)
hat.theta = mean(TFtime)
t.pivot = sort( (boot.out[,"theta*"] - hat.theta) / boot.out[,"se*"] )

se.hat.theta.boot = sd(boot.out[,"theta*"])
hat.theta - se.hat.theta * c(t.pivot[1950], t.pivot[50]) # bootstrap-t


# ex 8.5 --------------------------------------------------------------------

set.seed(1) 
B = 1999
boot.out = replicate(B, mean(sample(TFtime,n,replace=TRUE)))
(perc.int = quantile(boot.out, probs = c(0.025, 1-0.025)))

set.seed(1) 
sort.b = replicate(B, mean(sample(TFtime, n, replace = TRUE))) %>% sort()
(perc.int = c(sort.b[50],sort.b[1950]))


l = sort.b[50]
u = sort.b[1950]
c(l,u)
data.frame( theta.boot = sort.b) %>% 
  ggplot(aes(theta.boot)) + geom_histogram(color = "white",bins = 100) + xlab("") + 
  geom_vline(xintercept = mean(TFtime), color = "red",linetype = "dashed") + 
  geom_segment(aes( x = l, xend = u, y = 1, yend = 1), color = "blue" ,
               lineend = "round", size = 3) + 
  annotate("text", x = l, y = -3, label = TeX("$\\hat{\\theta}^*_{(50)}$")) + 
  annotate("text", x = u, y = -3, label = TeX("$\\hat{\\theta}^*_{(1950)}$"))

ggsave(filename = "ch8_3.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch8_3.png", path = 'images', width = 6.5, height = 4)


# ex 6 --------------------------------------------------------------------
set.seed(1)
hat.theta = mean(TFtime)
n = length(TFtime); B = 1999
theta.boot = replicate(B, mean(sample(TFtime, n, replace = TRUE)))
(b.factor = qnorm( sum(theta.boot <= hat.theta)/(B+1)))
(BC.int.probs = pnorm( 2*b.factor + c(qnorm(0.025), qnorm(1-0.025)) ) )
(BC.interval = quantile(theta.boot, BC.int.probs))



# ex 7 --------------------------------------------------------------------

theta.jack = vector(length = n)
for (i in 1:n){theta.jack[i] = mean(TFtime[-i])} 
acc = sum((mean(theta.jack) - theta.jack)^3)/
  (6 * sum(( mean(theta.jack) - theta.jack)^2)^(3/2))
acc
bz = b.factor + c(qnorm(0.025), qnorm(1-0.025))
(BCa.int.probs = pnorm( (bz)/(1 - acc * (bz)) + b.factor ))
(BCa.interval = quantile(theta.boot, BCa.int.probs))

mean(TFtime) + sd(TFtime) * qnorm(c(0.025,1-0.025))/sqrt(n)

rbind( mean(TFtime) + sd(TFtime) * qnorm(c(0.025,1-0.025))/sqrt(n), 
       norm.app.boot.CI,
       basic.ci,
       boot.t.ci,
       perc.int,
       BC.interval,
       BCa.interval)


 

# not included ------------------------------------------------------------
# using packages nptest and boot

library(nptest)
set.seed(1)
out = np.boot(TFtime, mean, R = 1999, level = 0.95)
out$normal
out$basic
out$percent

library(boot)
x = TFtime
mean.sd = function(x,index){
  x.boot = x[index]
  n = length(x)
  c(mean(x.boot),var(x.boot)/n)
}
set.seed(1) 
boot.out = boot::boot(x, mean.sd, 1999)
boot.ci(boot.out, conf = 0.95, type = c("norm","basic", "stud", "perc"), var.t0 = var(x)/n)
# boot.ci(boot.out, conf = 0.95, type = "bca", var.t0 = var(x)/n)
# # problem in estimating acc.