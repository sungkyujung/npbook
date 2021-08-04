# nonparametric book by S. Jung
# chapter 6. anova
#
#
# ----------------
library(tidyverse) 
library(gtools)
library(cowplot)
library(latex2exp)
library(car)
source("np_supp.R")  
library(showtext)
#font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()


# kickboard ---------------------------------------------------------------
# artificial data
kick = read.csv("data/kickboard.csv", header = TRUE) 
ggplot(kick, aes(x = brand, y = time)) + geom_boxplot() +
  geom_dotplot(binaxis = "y",stackdir = "center", dotsize =0.5, alpha = 0.5)
ggsave(filename = "ch6_kick.pdf", path = 'images', width = 4.8, height = 3) # smaller output
ggsave(filename = "ch6_kick.png", path = 'images', width = 4.8, height = 3)

# pizza -------------------------------------------------------------------
# data from https://dasl.datadescription.com/datafile/activating-baking-yeast/
pizza <- read.delim("data/pizza.txt")
pizza %>% ggplot(aes(x = Recipe, y = Activation.Times)) + geom_boxplot() +
  geom_dotplot(binaxis = "y",stackdir = "center", dotsize =0.5, alpha = 0.5)



# aov ---------------------------------------------------------------------
# anova
kick = read.csv("data/kickboard.csv", header = TRUE) 
aov(time ~ brand, data = kick) %>% summary()


# permutation test -------------------------------------------------------------
pizza <- read.delim("data/activating-baking-yeast.txt")
pizza
factorial(16) / (factorial(4))^4 
n = nrow(pizza)
k = 4 

(Fobs = F.stat(Activation.Times, Recipe, data = pizza))
set.seed(0)
M = 9999
Fperm = replicate(M, {
  perm.x = sample(pizza$Activation.Times, replace = FALSE)
  F.stat(perm.x, pizza$Recipe)
  })
(p.val = mean(c(Fobs,Fperm) >= Fobs))

p.val * (M+1)

# figure 6.1
x = seq(0,50,by = 0.01)
densityf = data.frame( x = x, densityF = df(x,k-1, n-k))
fig1 = data.frame( Fperm = Fperm) %>% ggplot(aes(x = Fperm)) + 
  geom_histogram(aes(y = ..density..), binwidth = 1, boundary = 0, color = "white") + 
  geom_vline(xintercept = Fobs, color = "red",linetype = "dashed") + 
  geom_line(aes(x =x ,y = densityF), data = densityf, color = "blue")
fig2 = data.frame( Fperm = Fperm) %>% ggplot(aes(x = Fperm)) + 
  geom_histogram(aes(y = ..density..), binwidth = 0.25,color = "white") + 
  geom_vline(xintercept = Fobs, color = "red",linetype = "dashed") + 
  geom_line(aes(x =x ,y = densityF), data = densityf, color = "blue") + xlim(0,5) 
                 
plot_grid(fig1, fig2)
ggsave(filename = "ch6_pizza_F.pdf", path = 'images', width = 6.5, height = 3) # two-panel
ggsave(filename = "ch6_pizza_F.png", path = 'images', width = 6.5, height = 3) # two-panel
 

pdf(file="images/ch6_pizza_F_qq.pdf",width = 6.5, height = 4)
car::qqPlot(Fperm, dist = "f", df1 = k-1, df2 = n-k)
dev.off() 




# rank_transform ----------------------------------------------------------
(pizza = pizza %>% mutate(rank = rank(Activation.Times)))
# rank-based F statistic
F.stat(rank, Recipe, data = pizza)
KW.stat(rank, Recipe, data = pizza)
KW.stat(Activation.Times, Recipe, data = pizza)
kruskal.test(pizza$Activation.Times, pizza$Recipe)$statistic


# example 3 and a figure --------------------------------------------------
(KWobs = KW.stat(Activation.Times, Recipe, data = pizza))
n = nrow(pizza)
k = 4 
set.seed(0)
M = 9999
KWperm = replicate(M, {
  perm.id = sample(1:n,n, replace = FALSE)
  KW.stat(pizza$Activation.Times[perm.id], pizza$Recipe)
})
(p.val = mean(c(KWobs,KWperm) >= KWobs))
p.val * (M+1)

x = seq(0,15,by = 0.01)
densityf = data.frame( x = x, densityF = dchisq(x,k-1))
data.frame( KW = KWperm) %>% ggplot(aes(x = KW)) + 
  geom_histogram(aes(y = ..density..), binwidth = 0.2,color = "white") + 
  geom_line(aes(x =x ,y = densityF), data = densityf, color = "blue")

ggsave(filename = "ch6_pizza_KW.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch6_pizza_KW.png", path = 'images', width = 6.5, height = 4)

pdf(file="images/ch6_pizza_KW_qq.pdf",width = 6.5, height = 4)
df = k-1 
car::qqPlot(KWperm, dist = "chisq", df =df)
dev.off()
 


# 예제 4 킥보드 ----------------------------------------------------------------
(KWobs = KW.stat(time, brand, data = kick))
n = nrow(kick)
k = 3
set.seed(0)
M = 9999
KWperm = replicate(M, {
  perm.id = sample(1:n,n, replace = FALSE)
  KW.stat(kick$time[perm.id], kick$brand)
})


x = seq(0,15,by = 0.01)
densityf = data.frame( x = x, densityF = dchisq(x,k-1))
data.frame( KW = KWperm) %>% ggplot(aes(x = KW)) + 
  geom_histogram(aes(y = ..density..), binwidth = 0.2,color = "white") + 
  geom_line(aes(x =x ,y = densityF), data = densityf, color = "blue")

ggsave(filename = "ch6_ex1_KW.pdf", path = 'images', width = 3, height = 3)
ggsave(filename = "ch6_ex1_KW.png", path = 'images', width = 3, height = 3)

ggsave(filename = "ch6_ex1_KWl.pdf", path = 'images', width = 5.5, height = 3) 

pdf(file="images/ch6_ex1_KW_qq.pdf",width = 3.5, height = 4)
df = k-1 
car::qqPlot(KWperm, dist = "chisq", df =df )
dev.off()


pdf(file="images/ch6_ex1_KW_qql.pdf",width = 6, height = 4)
df = k-1 
car::qqPlot(KWperm, dist = "chisq", df =df )
dev.off()


# 예제 6.6 --------------------------------------------------------------------
(KWobs = KW.stat(time, brand, data = kick))
(pvalue.chisq = 1 - pchisq(KWobs, df = 2))
set.seed(0)
M = 9999
KWperm = replicate(M, {
  x.perm = sample(kick$time, replace = FALSE)
  KW.stat(x.perm, kick$brand)
})
(p.val = mean(c(KWobs,KWperm) >= KWobs))

# function `kruskal.test`
kruskal.test(time ~ brand, kick)





# 순서가 있는 대립가설 -------------------------------------------------------------
# For ordered alternative ---------------------------------------
# tree data; artificial data
tree = read.csv("data/tree.csv")
summary(tree)
summary(aov(increment ~ slope, data = tree))
kruskal.test(increment ~ slope, data = tree)

head(tree)
# this is needed for JT test
(tree$slope = factor(tree$slope, ordered = TRUE))

# two sided test (n <= 100 exact, n > 100 asymptotic)
library(DescTools)
JonckheereTerpstraTest(increment ~ slope, data = tree)

# one sided test (n <= 100 exact, n > 100 asymptotic)
with(tree, JonckheereTerpstraTest(increment , slope,  alternative = "increasing"))

# not needed here (because n < 100), but if desired, random permutation can be applied
set.seed(0)
with(tree, JonckheereTerpstraTest(increment , slope,  alternative = "increasing", nperm= 5000))





# CBD ---------------------------------------------------------------------

# ex 6.12
farm = read.csv("data/smartfarm1.csv")
# classical F-test
(Fobs = F.stat.CBD(x, treatment, block, data = farm))
# summary(aov(x ~ treatment + block, data = farm))
k = 3;b = 3
(pval = 1-pf(Fobs,k-1, (k-1)*(b-1)))

# exact permutation using F statistic (for b = 3 only)
Fpermtest = F.test2.exact(x,treatment,block, data = farm)
Fpermtest$pvalue


# figure (comparing the permutation distribution  with theoretical null distribution)
xx = seq(0,60,by = 0.1)
densityf = data.frame( x = xx, densityF = df(xx,k-1, (k-1)*(b-1)))
fig1 = data.frame( Fperm = Fpermtest$Fperms) %>% ggplot(aes(x = Fperm)) + 
  geom_histogram(aes(y = ..density..), binwidth = .5,color = "white") + geom_rug() + 
  geom_vline(xintercept = Fobs, color = "red", linetype = "dashed") + 
  geom_line(aes(x =xx ,y = densityF), data = densityf, color = "blue")
fig1
ggsave(filename = "ch6_smartfarm_F.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch6_smartfarm_F.png", path = 'images', width = 6.5, height = 4) 



# Friedman ----------------------------------------------------------------

# rank transformation  
farm = read.csv("data/smartfarm1.csv")
as.data.frame(farm %>% arrange(block) %>% group_by(block) %>% mutate(ranks = rank(x)) )


farmr = as.data.frame(farm %>%  group_by(block) %>% mutate(ranks = rank(x)) )
farmr
# or use rankinblock(x,block)
with(farm, rankinblock(x,block))


# FM statistic
(FMobs = FM.stat.CBD(x,treatment,block,data = farm))

# pvalue (asymptotic)
(pval1 = 1-pchisq(FMobs,k-1))
# this is what friedman.test returns
friedman.test(x ~ treatment | block,data = farm)
with(farmr, friedman.test(x,treatment,block))
with(farmr, friedman.test(ranks,treatment,block))

# pvalue (exact permutation test)
Friedman.test.exact(x,treatment,block, data = farm)$pvalue


# figure ------------------------------------------------------------------
out = Friedman.test.exact(x,treatment,block, data = farm)
FMperm = out$Fperms
FMobs = out$Fobs

pdf(file="images/ch6_tree_FM_qq.pdf",width = 6.5, height = 4)
car::qqPlot(FMperm,dist="chisq", df=2)
dev.off()





# random permutation within block -----------------------------------------
randperm.inblock = function(x, block){
  block.id = unique(block)
  b = length(block.id)
  xperm = x
  for (j in 1:b){
    xperm[block == block.id[j]] = sample( xperm[block == block.id[j]], replace = FALSE)
  }
  return(xperm)
}
randperm.inblock(x,block)

# random permutation Friedman test ----------------------------------------
(FMobs = FM.stat.CBD(x, treatment, block, data = farm))
set.seed(0)
M = 9999
FMperm = with(farm, 
              replicate(M, {
    xperm = randperm.inblock(x,block)
    FM.stat.CBD(xperm, treatment, block)}) )
pval2 = mean(c(FMobs,FMperm) >= FMobs)

set.seed(0)
Friedman.test.randperm(x, treatment, block, M = 9999, data = farm)$pvalue




# Page trend test ---------------------------------------------------------

#######################
### ordered alternative, H1: t1 < t2 < t3 < t4
x = c(128, 128, 160, 142, 157, 179,
      122, 137, 177, 177, 160, 138,
      207, 188, 181, 164, 155, 175,
      120, 208, 199, 194, 177, 195)
treatment = rep(1:4, each=6)
block = rep(1:6, times = 4)
dat = data.frame(x, treatment, block)
write.csv(dat,"data/toydat_page.csv",row.names = F)

datmatrix = matrix(x, byrow=T, ncol=6)
rownames(datmatrix)=paste("trt",1:4)
colnames(datmatrix)=paste("block",1:6)
datmatrix  

crank::page.trend.test(t(datmatrix), rank = FALSE)   # for this function, columns are treatments # not recommended.

dat = read.csv("data/toydat_page.csv")
(Lobs = Page.L.stat(x,treatment,block, data =dat))
set.seed(0)
Page.trend.randperm(x,treatment,block, data = dat)$pvalue
 



# Skillings and Mack ------------------------------------------------------
smartfarm = read.csv("data/smartfarm.csv")
ggplot(smartfarm, aes(x = block, y = x, color = treatment)) + geom_boxplot()
with(smartfarm, rankinblock(x,block))


(SMobs = SM.stat(x,treatment,block,data = smartfarm))
(pvalue.asymp = 1 - pchisq(SMobs, df = 2))

SM.perm.test(x,treatment,block,M=9999,data = smartfarm)$pvalue








