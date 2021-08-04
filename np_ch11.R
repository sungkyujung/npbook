# nonparametric book by S. Jung
# chapter 11.  kernel density estimate
#
#
# ----------------
library(tidyverse) 
library(latex2exp) 
library(cowplot)
library(showtext)
#font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()


# Hidalgo stamp thickness data --------------------------------------------

# 1872 Hidalgo stamp data (이달고)
# Izenman, A. J. and Sommer, C. J. Philatelic mixtures and multimodal densities. Journal of the American Statistical association, 83(404):941-953, 1988.
thickness = read.csv("data/stamp.csv")$x  


# Fig 1 -------------------------------------------------------------------

summary(thickness)
xlim <- c(min(thickness), max(thickness))
par(mfrow=c(2,2))
h = 0.015;
binselect1 = seq(xlim[1],xlim[2]+h,h)
hist(thickness,binselect1, main="h = 0.015", xlim=xlim,col="lightblue")

h = 0.010;
binselect1 = seq(xlim[1],xlim[2]+h,h)
hist(thickness,binselect1, main="h = 0.010", xlim=xlim,col="lightblue")

h = 0.005;
binselect1 = seq(xlim[1],xlim[2]+h,h)
hist(thickness,binselect1, main="h = 0.005", xlim=xlim,col="lightblue")

h = 0.002;
binselect1 = seq(xlim[1],xlim[2]+h,h)
hist(thickness,binselect1, main="h = 0.002", xlim=xlim,col="lightblue")
dev.copy(pdf,'images/ch11_hidalgo_stamp.pdf')
dev.off()


# figure 2 ----------------------------------------------------------------


par(mfrow=c(2,2))
h = 0.012;
loc = 0
binselect1 = seq(xlim[1]-loc,xlim[2]+h,h)
hist(thickness,binselect1, xlim=xlim,col="lightblue", main = "start = 0.06")
loc = 0.002
binselect1 = seq(xlim[1]-loc,xlim[2]+h,h)
hist(thickness,binselect1, xlim=xlim,col="lightblue", main = "start = 0.058")
loc = 0.004
binselect1 = seq(xlim[1]-loc,xlim[2]+h,h)
hist(thickness,binselect1, xlim=xlim,col="lightblue", main = "start = 0.056")
loc = 0.006
binselect1 = seq(xlim[1]-loc,xlim[2]+h,h)
hist(thickness,binselect1, xlim=xlim,col="lightblue", main = "start = 0.054")
dev.copy(pdf,'images/ch11_hidalgo_stamp_loc.pdf')
dev.off()
par(mfrow=c(1,1))


# ex 2 soduim data --------------------------------------------------------

dat = c(1200, 2300, 5800, 1100, 1200, 1700, 1400, 1900, 3800, 1300)
n = length(dat)
h = 500
x = seq(500,7000, by = 50)
Kern_each = matrix(nrow = length(x), ncol = n)
for (i in 1:n){
  Kern_each[,i] = dnorm( (x -dat[i])/h )/(n*h)
}
kde = rowSums(Kern_each)
ggplot() + geom_rug(aes(x = dat)) + geom_line(aes(x = x, y = kde)) + 
  geom_line(aes(x = x, y = Kern_each[,9]), linetype = "dashed", color = "blue")+
  geom_line(aes(x = x, y = Kern_each[,2]), linetype = "dashed", color = "blue")

ggsave(filename = "ch11sod.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch11sod.png", path = 'images', width = 6.5, height = 4)





# ex 11.3: Hidalgo thickness -------------------------------------------------

a = density(thickness)
b = density(thickness, adjust = 1/4)
a$bw
b$bw
temp.dat = data.frame( x = c(a$x,b$x), kde = c(a$y,b$y), h = as.factor( c( rep(0.004, 512),rep(0.001, 512)) ) )
ggplot() + geom_histogram(aes(x = thickness, y = ..density..), bins = 30, color="white", fill = "grey") + 
  geom_line(data = temp.dat, mapping = aes(x,kde, color = h, linetype = h), size = 1)
ggsave(filename = "ch11hi-1.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch11hi-1.png", path = 'images', width = 6.5, height = 4)


# fig: kernel -------------------------------------------------------------
x = seq(-3,3,by = 0.01)
Gaussian = dnorm(x)
Epanechinikov = 3/4 * ( 1-x^2) * ifelse( abs(x) <= 1, 1, 0)
Uniform = 1/2 * ifelse( abs(x) <= 1, 1, 0)
data.frame(x, Gaussian, Epanechinikov, Uniform) %>% pivot_longer(cols = 2:4, names_to = "type", values_to = "kernel") %>% 
  ggplot(aes(x, kernel, color = type, linetype = type)) + geom_line()
ggsave(filename = "ch11kern.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch11kern.png", path = 'images', width = 6.5, height = 4)


# fig: variance-bias ------------------------------------------------------

h = seq(1,5,by = 0.01)
bias2 = h^4
variance = 300/h
mse = bias2 + variance 
data.frame(bandwidth = h,bias2,variance,mse) %>% pivot_longer(cols = 2:4, names_to = "terms", values_to = "amount") %>% 
  ggplot(aes(bandwidth, amount, color = terms, linetype = terms)) + 
  geom_line() + 
  annotate("text", 4.5, -1, label = "More smooth") +annotate("text", 1.5, -1, label = "Less smooth")
ggsave(filename = "ch11varbias.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch11varbias.png", path = 'images', width = 6.5, height = 4)


# 11.4 ex:kdebandwidth ---------------------------------------------------------

bw.SJ(thickness)
kde.nrd = density(thickness, bw = "nrd")
kde.nrd0 = density(thickness, bw = "nrd0") # default bw = "nrd0"
kde.SJ = density(thickness, bw = "SJ")
temp.dat = data.frame( x = c(kde.nrd$x,kde.nrd0$x,kde.SJ$x), kde = c(kde.nrd$y,kde.nrd0$y,kde.SJ$y), 
                       bw = as.factor( c( rep("nrd (정규분포 참조)", 512),rep("nrd0 (어림짐작)", 512),rep("SJ",512)) ) )
ggplot() + geom_histogram(aes(x = thickness, y = ..density..), bins = 30, color="white", fill = "grey") + 
  geom_line(data = temp.dat, mapping = aes(x,kde, color = bw, linetype = bw), size = 1)
ggsave(filename = "ch11hibw.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch11hibw.png", path = 'images', width = 6.5, height = 4)


# multi-dimensional kernels -----------------------------------------------
library(bivariate) 
limits = c(-5,5)
prec = 50
par(mfrow=c(2,2))
plot(nbvpdf (),ref.arrows=F, main = "(a)",
     n = prec, xlab = "x1", ylab = "x2", xlim = limits, ylim = limits)
plot(nbvpdf(sd.X = 2, sd.Y = 1) ,n = prec,  main = "(b)",
     ref.arrows=F, xlab = "x1", ylab = "x2", xlim = limits, ylim = limits)
plot(nbvpdf(sd.X = 2, sd.Y = 2) , n = prec,main = "(c)",
     ref.arrows=F, xlab = "x1", ylab = "x2", xlim = limits, ylim = limits)
plot(nbvpdf(sd.X = 1, sd.Y = 1, cor = -0.8) , n = prec, main = "(d)",
     ref.arrows=F, xlab = "x1", ylab = "x2", xlim = limits, ylim = limits)
dev.copy(pdf,'images/ch11-bivariatekernel.pdf')
dev.off()


# ex 11.5 penguins data -----------------------------------------------------------

#remotes::install_github("allisonhorst/palmerpenguins")
#library(palmerpenguins)
#head(penguins)
#data(penguins)
#temp = na.omit(penguins)
#write.csv(temp, "data/penguins.csv",row.names = FALSE)
penguins = read.csv("data/penguins.csv")
#penguins %>% ggplot(aes(bill_length_mm,flipper_length_mm, color = sex)) + geom_point()
penguins %>% ggplot(aes(bill_length_mm,flipper_length_mm, color = species, shape = species )) + geom_point() 
ggsave(filename = "ch11pen1.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch11pen1.png", path = 'images', width = 6.5, height = 4)


# now use ks package
penguins.b = penguins %>% 
  select(bill_length_mm, flipper_length_mm)

library(ks)
H <- ks::Hpi(penguins.b)
H
kde.p <- ks::kde(x = penguins.b, H = H)

par(mfrow=c(1,1))
# contour of density estimate. 
plot(kde.p, display = "filled.contour", cont = seq(5, 95, by = 10), 
     xlab = "bill length", ylab = "flipper length")
points(kde.p$x, pch = 1, cex = 0.5) 
dev.copy(pdf,'images/ch11pen2.pdf')
dev.off()

# 3D plot
plot(kde.p, display = "persp", xlab = "bill_length", ylab = "flipper length")
dev.copy(pdf,'images/ch11pen3.pdf')
dev.off()


#kde.by.group <- ks::kda(x = penguins.b, x.group = penguins$species)
#plot(kde.by.group, lwd = 2, col.pt = 1, cont = seq(5, 95, by = 10), drawpoints = TRUE)


# Mollified density -------------------------------------------------------
set.seed(1)
x = c(rlnorm(281,0,1),rnorm(19,0,4))
plot(x)
xx = seq(0,2,by = .001)
f.wigg = function(xx) {abs(sin(xx) + 5*sin(xx/100) + sin(xx*100)/10)*xx*(2-xx)}
f.k = function(y) {dnorm(y,0,0.1)}
# convolution
f.z = function(z) integrate(function(x,z) f.k(z-x)*f.wigg(x), -Inf,Inf,z)$value
f.Z = Vectorize(f.z)
data.frame(xx = xx, f = f.wigg(xx), mollified = f.Z(xx)) %>% 
  pivot_longer(cols = 2:3, names_to = "type", values_to = "density")  %>%  ggplot(aes(xx,density,color = type )) + geom_line()

ggsave(filename = "ch11mol.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch11mol.png", path = 'images', width = 6.5, height = 4)


# ex:hidalgo-pointwiseci --------------------------------------------------


## CI
h = bw.SJ(thickness)
n = length(thickness)
ks.out = ks::kde(thickness, h = h)
f.hat = ks.out$estimate
ci.lower = f.hat - qnorm(1-0.025) * sqrt( f.hat / (2*sqrt(pi)) / (n*h))
ci.upper = f.hat + qnorm(1-0.025) * sqrt( f.hat / (2*sqrt(pi)) / (n*h))
data.frame(x = ks.out$eval.points, ci.low = ci.lower, ci.up = ci.upper) %>% head()


temp.dat.norm = data.frame(x = ks.out$eval.points, kde = f.hat, ci.lower, ci.upper)
ggplot() + 
  geom_histogram(aes(x = thickness, y = ..density..), bins = 30, color="white", fill = "grey") + 
  geom_line(data = temp.dat.norm, mapping = aes(x,kde), size = 1) + 
  geom_ribbon(data=temp.dat.norm, aes(x,kde,ymin=ci.lower,ymax=ci.upper),alpha=0.3, fill = "blue")

ggsave(filename = "ch11hidci-1.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch11hidci-1.png", path = 'images', width = 6.5, height = 4)


# ex: hidalgo-bootstrap-plugin --------------------------------------------------------

x = ks.out$eval.points
kde0 = ks.out$estimate 
n = length(thickness)

B = 1999
kde.boot = replicate(B,{
  ks::kde(sample(thickness,n, replace = TRUE),# 붓스트랩 샘플
          h = h,          # 같은 너비 h 사용 
          eval.points = x # 같은 점들 X에서
  )$estimate      #  kde(x) 계산
})
se.fhat = apply(kde.boot,1,sd)
ci.lower = f.hat - qnorm(1-0.025) * se.fhat
ci.upper = f.hat + qnorm(1-0.025) * se.fhat 

temp.dat.boot1 = data.frame(x = ks.out$eval.points, kde = f.hat, ci.lower, ci.upper)
# ggplot() + 
#   geom_histogram(aes(x = thickness, y = ..density..), bins = 30, color="white", fill = "grey") + 
#   geom_line(data = temp.dat.boot1, mapping = aes(x,kde), size = 1) + 
#   geom_ribbon(data=temp.dat.boot1, aes(x,kde,ymin=ci.lower,ymax=ci.upper),alpha=0.3, fill = "blue")+ 
#   geom_ribbon(data=temp.dat.norm, aes(x,kde,ymin=ci.lower,ymax=ci.upper),alpha=0.3, fill = "red")


# ex:hidalgo-bootstrap-last -----------------------------------------------

# standard bootstrap interval
kde.quantiles = apply(kde.boot, 1, 
                     function(x){quantile(x,probs = c(1-0.025,0.025))})
ci.lower = 2*f.hat - kde.quantiles[1,] 
ci.upper = 2*f.hat - kde.quantiles[2,]
temp.dat.boot2 = data.frame(x = ks.out$eval.points, kde = f.hat, ci.lower, ci.upper)
ggplot() + 
  geom_histogram(aes(x = thickness, y = ..density..), bins = 30, color="white", fill = "grey") + 
  geom_line(data = temp.dat.boot2, mapping = aes(x,kde), size = 1) + 
  geom_ribbon(data=temp.dat.boot2, aes(x,kde,ymin=ci.lower,ymax=ci.upper),alpha=0.3, fill = "yellow")

# absolute deviation - bootstrap 
abs.dev.boot = apply(kde.boot, 2, function(x){ abs(x - f.hat) })
abs.dev.boot.quantiles = apply(abs.dev.boot, 1, 
                               function(x){quantile(x,probs = 1-0.05)})
ci.lower = f.hat - abs.dev.boot.quantiles
ci.upper = f.hat + abs.dev.boot.quantiles

temp.dat.boot3 = data.frame(x = ks.out$eval.points, kde = f.hat, ci.lower, ci.upper)

g0 = ggplot() + 
  geom_histogram(aes(x = thickness, y = ..density..), bins = 30, color="white", fill = "grey") + 
  geom_line(data = temp.dat.boot3, mapping = aes(x,kde), size = 1)
g1 = g0 + geom_ribbon(data=temp.dat.norm, aes(x,kde,ymin=ci.lower,ymax=ci.upper),alpha=0.3, fill = "blue") + 
  ggtitle("플러그인 대표본 신뢰구간")
g2 = g0 + geom_ribbon(data=temp.dat.boot1, aes(x,kde,ymin=ci.lower,ymax=ci.upper),alpha=0.3, fill = "blue")+
  ggtitle("붓스트랩-플러그인 신뢰구간")
g3 = g0 + geom_ribbon(data=temp.dat.boot2, aes(x,kde,ymin=ci.lower,ymax=ci.upper),alpha=0.3, fill = "blue")+
  ggtitle("표준 붓스트랩 신뢰구간")
g4 = g0 + geom_ribbon(data=temp.dat.boot3, aes(x,kde,ymin=ci.lower,ymax=ci.upper),alpha=0.3, fill = "blue")+
  ggtitle("절대편차 붓스트랩 신뢰구간")
plot_grid(g1,g2,g3,g4)

ggsave(filename = "ch11hidci-2.pdf", path = 'images', width = 6.5, height = 5)
ggsave(filename = "ch11hidci-2.png", path = 'images', width = 6.5, height = 5)


# ex:confband-hidalgo -----------------------------------------------------

Delta = apply(abs.dev.boot,2,max) 
q.Delta = quantile(Delta, 1-0.05)
conf.band.lower = f.hat - q.Delta
conf.band.upper = f.hat + q.Delta
temp.dat.boot4 = data.frame(x = ks.out$eval.points, kde = f.hat, conf.band.lower, conf.band.upper)
g0 + geom_ribbon(data=temp.dat.boot4, aes(x,kde,ymin=conf.band.lower,ymax=conf.band.upper),
                 alpha=0.3, fill = "blue")+
  geom_ribbon(data=temp.dat.boot3, aes(x,kde,ymin=ci.lower,ymax=ci.upper),alpha=0.3, 
              fill = "white", linetype = "dashed") + 
  ggtitle("붓스트랩 신뢰밴드") 
ggsave(filename = "ch11hidci-3.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch11hidci-3.png", path = 'images', width = 6.5, height = 4)


conf.band.lower = pmax(f.hat - q.Delta,0)
temp.dat.boot4a = data.frame(x = ks.out$eval.points, kde = f.hat, conf.band.lower, conf.band.upper)
g0 + geom_ribbon(data=temp.dat.boot4a, 
                 aes(x,kde,ymin=conf.band.lower,ymax=conf.band.upper),
                 alpha=0.3, fill = "blue")
ggsave(filename = "ch11hidci-4.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch11hidci-4.png", path = 'images', width = 6.5, height = 4)

# 
# # bootstrap-t interval 
# t.boot = apply(kde.boot, 2, function(x){
#   (x - f.hat)/sqrt( x / (2*sqrt(pi)) / (n*h)) }
# )
# plot(x,kde.boot[,3])
# t.quantiles = apply(t.boot, 1, function(x){
#   quantile(x,probs = c(1-0.025,0.025))}
# )
# ci.lower = f.hat - t.quantiles[1,] * sqrt( f.hat / (2*sqrt(pi)) / (n*h))
# ci.upper = f.hat - t.quantiles[2,] * sqrt( f.hat / (2*sqrt(pi)) / (n*h))
# temp.dat.boot2 = data.frame(x = ks.out$eval.points, kde = f.hat, ci.lower, ci.upper)
# ggplot() + 
#   geom_histogram(aes(x = thickness, y = ..density..), bins = 30, color="white", fill = "grey") + 
#   geom_line(data = temp.dat.boot2, mapping = aes(x,kde), size = 1) + 
#   geom_ribbon(data=temp.dat.boot2, aes(x,kde,ymin=ci.lower,ymax=ci.upper),alpha=0.3, fill = "yellow")+ 
#   geom_ribbon(data=temp.dat.boot1, aes(x,kde,ymin=ci.lower,ymax=ci.upper),alpha=0.3, fill = "blue")+ 
#   geom_ribbon(data=temp.dat.norm, aes(x,kde,ymin=ci.lower,ymax=ci.upper),alpha=0.3, fill = "red")
# 
# 
# 
# 
