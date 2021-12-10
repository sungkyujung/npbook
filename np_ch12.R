# nonparametric book by S. Jung
# chapter 12: local polynomial regression
#
#
# ----------------
library(tidyverse) 
library(latex2exp) 
library(cowplot)
library(KernSmooth) 
source("np_supp.R")
library(showtext)
#font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()
 

# 12.1 ex:biden ----------------------------------------------------------------

poll <- read.table("data/US-Biden-Joe-Favorable-poll-responses-clean.tsv",
                    sep = "\t",header = T)
library(lubridate) 
poll %>% 
  ggplot(aes(as_date(start_date),y = Favorable)) + 
  geom_point() + 
  geom_smooth(se = FALSE,color = "red") + 
  geom_smooth(method = "lm", se= TRUE)

ggsave(filename = "ch11r-biden.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch11r-biden.png", path = 'images', width = 6.5, height = 4)


# fig:mxatx ---------------------------------------------------------------
# fix x
x = as.numeric(poll$start_date %>% as_date(), "days")
y = poll$Favorable
xgrid = round(seq(min(x),max(x), length.out = 100))
xp = median(x)-11
as_date(xp)
h = 300
kx = Kh(xp-x,h)
w = kx / sum(kx)
Kh = function(x,h){ exp( - x^2 /h^2)/ sqrt(2*pi)/h}
kernreg = function(xgrid, x, y, h){
  kx = sapply(x,function(xi){ Kh(xgrid-xi,h)})
  W <- kx / rowSums(kx)
  drop(W %*% y) # xgrid의 각 점에서 sum w * y 계산
} 
dat = data.frame(xgrid, yy = kernreg(xgrid,x,y,300))
data.frame(x,y,w) %>% ggplot() + geom_point(aes(x,y,color = w)) + 
  scale_color_gradient(low="grey", high="red") +
  geom_line(data = dat, aes(x = xgrid, y = yy),color = "blue") + 
  scale_x_continuous(labels =NULL) + geom_vline(xintercept = xp, linetype = "dashed") + ggtitle("커널회귀곡선 m(x) 추정")
ggsave(filename = "ch11r2.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch11r2.png", path = 'images', width = 6.5, height = 4)



# kernreg definition ------------------------------------------------------
## manual implementation of Nadaraya-Watson estimator
dat = data.frame(start_date = as.numeric(poll$start_date %>% as_date(), "days"),
                 poll$Favorable)

x = as.numeric(poll$start_date %>% as_date(), "days")
y = poll$Favorable
xgrid = round(seq(min(x),max(x), length.out = 100))
h = 100
Kh = function(x,h){ exp( - x^2 /h^2)/ sqrt(2*pi)/h}
kernreg = function(xgrid, x, y, h){
  kx = sapply(x,function(xi){ Kh(xgrid-xi,h)})
  W <- kx / rowSums(kx)
  drop(W %*% y) # xgrid의 각 점에서 sum w * y 계산
}  
mhat = kernreg(xgrid, x,y,100)

# see also lcreg() at np_supp.R
source("np_supp.R")
lcreg.out = lcreg(xgrid,x,y,100)
max(abs(lcreg.out$mhat - mhat))



# fig 12.3 ----------------------------------------------------------------


xgrid = round(seq(min(x),max(x), length.out = 1000))
dat = data.frame(x = as_date(xgrid),h100 = kernreg(xgrid, x,y, 100), 
                 h10 = kernreg(xgrid, x,y,10), 
                 h300 = kernreg(xgrid, x,y,300)) %>% 
  pivot_longer(cols = 2:4, values_to = "y", names_to = "bandwidth")
poll %>% ggplot(aes(as_date(start_date),y = Favorable)) + geom_point() + 
  geom_line(data = dat, aes(x=x,y=y,linetype = bandwidth,color= bandwidth),size=1.2)
ggsave(filename = "ch11r3.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch11r3.png", path = 'images', width = 6.5, height = 4)




# fig:12.4 loclinearatx --------------------------------------------------------
# local linear
source("np_supp.R")
ll.out = llreg(xgrid,x,y,300)
ll = data.frame(x = xgrid, y = ll.out$mhat)
xp = 16000
h = 300
kx = Kh(xp-x,h)
w = kx / sum(kx)
W = diag(kx)
X = model.matrix(y ~ x)
beta.xp = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
dat1 = data.frame(x,y,kx,yy = beta.xp[1] + beta.xp[2] * x) 
g1 = dat1 %>% ggplot() + geom_point(aes(x,y,color = kx)) + 
  scale_color_gradient(low="grey", high="red") + geom_vline(xintercept = xp, linetype = "dashed") +
  geom_line(aes(x = x, y = yy,color = kx)) +  annotate(
    geom = "curve", xend = xp, yend = beta.xp[1] + beta.xp[2] *xp, x = 16000, y = 50, 
    curvature = -.3, arrow = arrow(length = unit(2, "mm"))
  ) + guides(color = FALSE)+
  annotate("text", x= 16000, y = 50, label = TeX(r'($(x, \\hat{m}_1(x))$)', output="character"),hjust = "right",parse=T) +
  scale_x_continuous(labels =NULL) + ggtitle("국소선형회귀추정, x = 2013/10/22") + guides(fill=FALSE)+ 
  geom_line(data = ll, aes(x=as_date(x),y=y),color = "blue")

xp = round(16750)
h = 300
kx = Kh(xp-x,h)
w = kx / sum(kx)
W = diag(kx)
X = model.matrix(y ~ x)
beta.xp = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
dat2 = data.frame(x,y,kx,yy = beta.xp[1] + beta.xp[2] * x) 

g2 = dat2 %>% ggplot() + geom_point(aes(x,y,color = kx)) + 
  scale_color_gradient(low="grey", high="red") + geom_vline(xintercept = xp, linetype = "dashed") +
  geom_line(aes(x = x, y = yy,color = kx)) +  annotate(
    geom = "curve", xend = xp, yend = beta.xp[1] + beta.xp[2] *xp, x = 16000, y = 50, 
    curvature = -.3, arrow = arrow(length = unit(2, "mm"))
  ) + guides(color = FALSE)+ 
  annotate("text", x= 16000, y = 50, label = TeX(r'($(x, \\hat{m}_1(x))$)', output="character"),hjust = "right",parse=T) +
  scale_x_continuous(labels =NULL) + ggtitle("국소선형회귀추정, x = 2015/11/11 ") + 
  geom_line(data = ll, aes(x=as_date(x),y=y),color = "blue")
plot_grid(g1,g2)
ggsave(filename = "ch11r4.pdf", path = 'images', width = 6.7, height = 4)
ggsave(filename = "ch11r4.png", path = 'images', width = 6.7, height = 4)
 

# example:12.2 ex:loc-linear --------------------------------------------------


library(KernSmooth)
p = 1
h = 50
ll = KernSmooth::locpoly(x,y,degree = p, bandwidth = h)
head(as.data.frame(ll),3)

xgrid = ll$x
llregout = llreg(xgrid, x, y, h)
cbind(xgrid, llregout$mhat) %>% head(3)


# figure 12.5 ------------------------------------------------------------------


degree = 0 
bw = c(50,300)
ll = KernSmooth::locpoly(x,y,degree = degree, bandwidth = bw[1],gridsize =1000)
llbw = data.frame(x = ll$x, y = ll$y, bw = as.factor(bw[1]))
ll = KernSmooth::locpoly(x,y,degree = degree, bandwidth = bw[2],gridsize =1000)
l0bw = rbind(llbw, data.frame(x = ll$x, y = ll$y, bw = as.factor(bw[2])))
g1= ggplot() + geom_point(aes(x,y)) + 
  geom_line(aes(x,y,color = bw,linetype = bw), data = l0bw,size=1.5)  + 
  ggtitle("나다라야-왓슨 (local constant)")+scale_x_continuous(labels =NULL) 
degree = 1  
ll = KernSmooth::locpoly(x,y,degree = degree, bandwidth = bw[1],gridsize =1000)
llbw = data.frame(x = ll$x, y = ll$y, bw = as.factor(bw[1]))
ll = KernSmooth::locpoly(x,y,degree = degree, bandwidth = bw[2],gridsize =1000)
l1bw = rbind(llbw, data.frame(x = ll$x, y = ll$y, bw = as.factor(bw[2])))
g2= ggplot() + geom_point(aes(x,y)) + 
  geom_line(aes(x,y,color = bw,linetype = bw), data = l1bw,size=1.5)  + 
  ggtitle("국소선형 (local linear)")+scale_x_continuous(labels =NULL) 
degree = 2  
ll = KernSmooth::locpoly(x,y,degree = degree, bandwidth = bw[1],gridsize =1000)
llbw = data.frame(x = ll$x, y = ll$y, bw = as.factor(bw[1]))
ll = KernSmooth::locpoly(x,y,degree = degree, bandwidth = bw[2],gridsize =1000)
l2bw = rbind(llbw, data.frame(x = ll$x, y = ll$y, bw = as.factor(bw[2])))
g3 = ggplot() + geom_point(aes(x,y)) + 
  geom_line(aes(x,y,color = bw,linetype = bw), data = l2bw,size=1.5)  + 
  ggtitle("국소 2차 (local quadratic)") +scale_x_continuous(labels =NULL) 
plot_grid(g1,g2,g3, align = "h", nrow = 1)
#ggsave(filename = "ch11r5.pdf", path = 'images', width = 6.5, height = 3)
#ggsave(filename = "ch11r5.png", path = 'images', width = 6.5, height = 3)

lbwdat = rbind(data.frame(l0bw,degree = "나다라야-왓슨 (local constant)"),
      data.frame(l1bw,degree = "국소선형 (local linear)"),
      data.frame(l2bw,degree = "국소 2차 (local quadratic)"))
lbwdat$degree = factor(lbwdat$degree,levels = c("나다라야-왓슨 (local constant)",
                                                "국소선형 (local linear)",
                                                "국소 2차 (local quadratic)"))
head(lbwdat)
ggplot() + facet_grid(rows = vars(degree)) + 
  geom_point(aes(x,y), data = data.frame(x=x,y=y)) + 
  geom_line(aes(x,y,color = bw,linetype = bw), data = lbwdat ,size=1.5)
ggsave(filename = "ch11r5.pdf", path = 'images', width = 6.5, height = 9)
ggsave(filename = "ch11r5.png", path = 'images', width = 6.5, height = 9)



# bandwidth selection -----------------------------------------------------

# ex:biden-dpill ----------------------------------------------------------
(h.dpi = KernSmooth::dpill(x,y))
ll = KernSmooth::locpoly(x,y,degree = 1, bandwidth = h.dpi ,gridsize =1000)
ll = data.frame(ll)
ggplot() + geom_point(aes(x,y)) + 
  geom_line(aes(x,y), data = ll,size=1.5, color = "blue")  +scale_x_continuous(labels =NULL) 


# cross validation --------------------------------------------------------
h.range = seq(40, 120, by = 2)
K = length(h.range)
cv = vector(length = K)
for (k in 1:K){
  h = h.range[k]
  llreg.out = llreg(x,x,y,h)
  cv[k] = mean( ( (y - llreg.out$mhat) / 
                    (1- diag(llreg.out$weight)) )^2 )
}
(h.cv = h.range[which.min(cv)])

qplot(h.range, cv) + geom_point()
ggsave(filename = "ch11r6.pdf", path = 'images', width = 4.8, height = 3)
ggsave(filename = "ch11r6.png", path = 'images', width = 4.8, height = 3)

ll = KernSmooth::locpoly(x,y,degree = 1, bandwidth = h.dpi,gridsize =1000)
llbw = data.frame(x = ll$x, y = ll$y, bw = "dpi")
ll = KernSmooth::locpoly(x,y,degree = 1, bandwidth = h.cv,gridsize =1000)
l0bw = rbind(llbw, data.frame(x = ll$x, y = ll$y, bw = "cv"))
ggplot() + geom_point(aes(x,y)) + 
  geom_line(aes(x,y,color = bw,linetype = bw), data = l0bw,size=1.5)  + 
  ggtitle("Local linear: bandwidth by h.dpi and h.cv")+scale_x_continuous(labels =NULL) 

ggsave(filename = "ch11r7.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch11r7.png", path = 'images', width = 6.5, height = 4)


# 신뢰구간 1 ------------------------------------------------------------------

n = length(x)
llreg.x = llreg(x, x, y, h = h.cv)
residuals = y - llreg.x$mhat 
tr = sum( diag( t(diag(n) - llreg.x$weight) %*% (diag(n) - llreg.x$weight)))
(hat.sigma2 = sum(residuals^2) / tr)

xgrid = seq(min(x), max(x), length.out = 100)
llreg.out = llreg(xgrid, x, y, h = h.cv)
mhat = llreg.out$mhat
se = sqrt( rowSums(llreg.out$weight^2) * hat.sigma2 ) 
ci.theory = data.frame(x = xgrid, 
                       ci.low = mhat - qnorm(1-0.025)*se,
                       mhat = mhat,
                       ci.up = mhat + qnorm(1-0.025)*se)


# 신뢰구간 2 ------------------------------------------------------------------

xgrid = seq(min(x), max(x), length.out = 100)
mhat = llreg(xgrid, x, y, h = h.cv)$mhat
set.seed(1)
B = 1999
m.boot = replicate(B,{
  boot.id = sample(1:n,n,replace = TRUE)
  x.boot = x[boot.id]
  y.boot = y[boot.id]
  llreg(xgrid, x.boot, y.boot, h = h.cv)$mhat
}) %>% t()
se.boot = apply(m.boot,2,sd)

ci.pi = data.frame(x = xgrid, 
                   ci.low = mhat - qnorm(1-0.025)*se.boot,
                   mhat = mhat,
                   ci.up = mhat + qnorm(1-0.025)*se.boot)

q1 = apply(m.boot,2,function (x) {quantile(x,0.025)})
q2 = apply(m.boot,2,function (x) {quantile(x,1-0.025)})

ci.perc = data.frame(x = xgrid, 
                     ci.low =q1,
                     mhat = mhat,
                     ci.up = q2)

# ggplot() + geom_point(aes(x,y)) + 
#   geom_line(data=ci.theory, aes(x,mhat), color = "blue") + 
#   geom_ribbon(data=ci.theory, aes(x,mhat,ymin=ci.low,ymax=ci.up),alpha=0.3, 
#               fill = "blue", linetype = "dashed") +
#   geom_ribbon(data=ci.pi, aes(x,mhat,ymin=ci.low,ymax=ci.up),alpha=0.3, 
#               fill = "red", linetype = "dashed") +
#   geom_ribbon(data=ci.perc, aes(x,mhat,ymin=ci.low,ymax=ci.up),alpha=0.3, 
#               fill = "yellow", linetype = "dashed") 

g1=ggplot() + geom_point(aes(x,y)) + 
  geom_line(data=ci.theory, aes(x,mhat), color = "blue") + 
  geom_ribbon(data=ci.theory, aes(x,mhat,ymin=ci.low,ymax=ci.up),alpha=0.3, 
              fill = "blue", linetype = "dashed") + ggtitle("이론적인 신뢰구간")

g2=ggplot() + geom_point(aes(x,y)) + 
  geom_line(data=ci.pi, aes(x,mhat), color = "blue") + 
  geom_ribbon(data=ci.pi, aes(x,mhat,ymin=ci.low,ymax=ci.up),alpha=0.3, 
              fill = "blue", linetype = "dashed") + ggtitle("붓스트랩-플러그인 신뢰구간")


g3=ggplot() + geom_point(aes(x,y)) + 
  geom_line(data=ci.theory, aes(x,mhat), color = "blue") + 
  geom_ribbon(data=ci.perc, aes(x,mhat,ymin=ci.low,ymax=ci.up),alpha=0.3, 
              fill = "blue", linetype = "dashed") + ggtitle("붓스트랩-퍼센타일 신뢰구간")
plot_grid(g1,g2,g3, align = "h", nrow = 1)

ggsave(filename = "ch11r8.pdf", path = 'images', width = 6.5, height = 3)
ggsave(filename = "ch11r8.png", path = 'images', width = 6.5, height = 3)

##

lbwdat = rbind(data.frame(ci.theory,type = "정규근사이론"),
               data.frame(ci.pi,type = "붓스트랩-플러그인"),
               data.frame(ci.perc,type = "붓스트랩-퍼센타일"))
lbwdat$type = factor(lbwdat$type,levels = c("정규근사이론",
                                                "붓스트랩-플러그인",
                                                "붓스트랩-퍼센타일"))
head(lbwdat)
ggplot() + facet_grid(rows = vars(type)) + 
  geom_point(aes(x,y), data = data.frame(x=x,y=y)) + 
  geom_line(data=lbwdat, aes(x,mhat), color = "blue") + 
  geom_ribbon(data=lbwdat, aes(x,mhat,ymin=ci.low,ymax=ci.up),alpha=0.3, 
              fill = "blue", linetype = "dashed")
ggsave(filename = "ch11r8.pdf", path = 'images', width = 6.5, height = 9)
ggsave(filename = "ch11r8.png", path = 'images', width = 6.5, height = 9)
