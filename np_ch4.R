# nonparametric book by S. Jung
# chapter 4. 두 모집단에 대한 추론
#
#
# ----------------
library(tidyverse) 
library(DescTools)
library(showtext)
#font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()

# example 4.1 -------------------------------------------------------------
#  From Golub et al. 1999
leukemia4= read.csv(file = "data/leukemia4.csv")
head(leukemia4,3)

leukemia4 %>% ggplot(aes(x = ge, fill = label,color = label)) +  
  geom_density(color="#e9ecef", alpha=0.6, position = 'identity') +
  geom_rug(length = unit(0.05, "npc"))
ggsave(filename = "leukemia4.pdf", path = 'images', width = 4.8, height = 3)
ggsave(filename = "leukemia4.png", path = 'images', width = 4.8, height = 3)
sum(leukemia4$label == "ALL") # x, m 
sum(leukemia4$label == "AML") # y, n 


# example 4.2 & 4.3 -------------------------------------------------------------
# example 2: cow data (artificial)
cowdat = read.csv("data/cow.csv")
ranks = rank(cowdat$Marb)
data.frame(cowdat, ranks)
(W.obs = sum( ranks[cowdat$RFI == "Low"] ))


# fig 4.1 -----------------------------------------------------------------------

allperms.id = gtools::combinations(7,4)
W.perm = apply(allperms.id,1,
               function (perm.y.id){ sum(ranks[perm.y.id]) }) 
table(W.perm)


data.frame(W.perm = W.perm) %>% ggplot(aes(x =W.perm)) + geom_histogram(binwidth = 0.5,color = "white") + 
  geom_dotplot(binwidth = .3) + geom_vline(xintercept = 14, color = "red", linetype = "dashed") + 
  ggtitle("윌콕슨 순위합 검정통계량의 영분포; m = 3, n =4") 
ggsave(filename = "ch4_wilcox.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch4_wilcox.png", path = 'images', width = 6.5, height = 4)


# example 4.4 -------------------------------------------------------------
sum(W.perm <= W.obs)
?dwilcox



# example 4.5 -------------------------------------------------------------

# back to leukemia data
# large sample approximation
# take ALL as x, AML as y 
m = sum(leukemia4$label == "ALL") # 47
n = sum(leukemia4$label == "AML") # 25
N = n+m
ew = n*(N+1)/2
vw = n*m*(N+1)/12
Wobs = sum( rank(leukemia4$ge)[leukemia4$label == "AML"] )

w = 500:1300 
df = data.frame(w,norm = dnorm(w,ew,sqrt(vw)))
# random permutation
B = 1999
set.seed(1)
Wperm = vector(length = B)
for (b in 1:B){
  Wperm[b] = sum(sample(1:N, n, replace = FALSE))
}
data.frame(W = Wperm) %>% ggplot() + 
  geom_histogram(aes(x = W, y = ..density..), binwidth = 50, fill = "grey",color = "black") + 
  geom_line(data = df, aes(x = w, y = norm), color = "blue") +
  geom_vline(xintercept = Wobs, color = "red", linetype = "dashed") 
ggsave(filename = "ch4_wilcox_clt.pdf", path = 'images', width = 4.8, height = 3) # smaller output
ggsave(filename = "ch4_wilcox_clt.png", path = 'images', width = 4.8, height = 3) # smaller output

sum(Wperm > Wobs)
(p_rp = 2*mean(c(Wobs,Wperm) >= Wobs))
z = ( Wobs - ew ) / sqrt(vw)
(p_z = 2*pnorm(z, lower.tail = FALSE))


# example 4.6 -------------------------------------------------------------
# leukemia4 data
# Dbar null distribution

ge = leukemia4$ge 
Dobs = mean(ge[leukemia4$label == "AML"])
tobs = -t.test(ge ~ leukemia4$label,var.equal = TRUE)$statistic
Dvec = vector(length = B)
Tvec = vector(length = B)
set.seed(1)
perm.out = replicate(B, {
  perm.label = sample(leukemia4$label,N, replace = FALSE)
  D.perm = mean(ge[perm.label == "AML"])
  t.perm = -t.test(ge ~ perm.label,var.eqaul = FALSE)$statistic
  c(D.perm, t.perm)  
})
D.perm = perm.out[1,]
t.perm = perm.out[2,]


dd = 120:270
meand = mean(ge)
vard = var(ge)*(N-1)/N /n /(N-1) * m 
df = data.frame(dd,norm = dnorm(dd, meand, sqrt(vard)))
data.frame(Dbar = D.perm) %>% ggplot() + 
  geom_histogram(aes(x = Dbar, y = ..density..), binwidth = 10, fill = "grey",color = "black") + 
  geom_line(data = df, aes(x = dd, y = norm), color = "blue") +
  geom_vline(xintercept = Dobs, color = "red", linetype = "dashed")
ggsave(filename = "ch4_dbar_null_s.pdf", path = 'images', width = 4.8, height = 3) # smaller output
ggsave(filename = "ch4_dbar_null_s.png", path = 'images', width = 4.8, height = 3) # smaller output

sum(Dobs < D.perm) 
(p = 2 * mean(Dobs <= c(Dobs,D.perm)))

# fig 4.4 -----------------------------------------------------------------

df2 = data.frame(dd,norm = dnorm(dd, mean(D.perm), sd(D.perm)))
data.frame(Dbar = D.perm) %>% ggplot() + 
  geom_histogram(aes(x = Dbar, y = ..density..), binwidth = 10, fill = "grey",color = "black") + 
  geom_line(data = df, aes(x = dd, y = norm), color = "blue") +
  geom_line(data =df2, aes(x = dd, y = norm), color = "black", linetype = "dotted",size =1.3) +
  geom_vline(xintercept = Dobs, color = "red", linetype = "dashed")
ggsave(filename = "ch4_dbar_null_2.pdf", path = 'images', width = 6.5, height = 4) 
ggsave(filename = "ch4_dbar_null_2.png", path = 'images', width = 6.5, height = 4) 


# fig 4.5 -----------------------------------------------------------------

tt = seq(-4,4,by=0.1)
df2 = data.frame(tt,norm = dnorm(tt, mean(t.perm), sd(t.perm)))

data.frame(t = t.perm) %>% ggplot() + 
  geom_histogram(aes(x = t, y = ..density..), binwidth = 0.5, fill = "grey",color = "black") +
  geom_line(data =df2, aes(x = tt, y = norm), color = "black", linetype = "dotted",size =1.3) + 
  geom_vline(xintercept = tobs, color = "red", linetype = "dashed") 
ggsave(filename = "ch4_t_null_2.pdf", path = 'images', width = 6.5, height = 4) 
ggsave(filename = "ch4_t_null_2.png", path = 'images', width = 6.5, height = 4) 

sum(tobs <= t.perm) 
mean( abs(c(t.perm,tobs) >= abs(tobs)) )
2*(1-pnorm((tobs - mean(t.perm)) /  sd(t.perm)))

# example 4.7 -------------------------------------------------------------
# confidence interval 

x = cowdat$Marb[cowdat$RFI == "High"]
y = cowdat$Marb[cowdat$RFI == "Low"]
m = length(x); n = length(y)
D  = sapply(y, function(y){y-x}) # Delta = E(Y) - E(X)
(D.ord = sort(D))
alpha = 0.2 
kl = qwilcox(alpha/2,m,n)
ku = m*n - kl + 1
c(kl,ku)
(ci = D.ord[c(kl,ku)])
median(D.ord)
(coverage = sum(dwilcox(kl:(ku-1),m,n)))

#sum(dwilcox(0:(kl-1),m,n))
#sum(dwilcox(ku:(m*n),m,n))
#sum(dwilcox(0:kl,m,n))
wilcox.test(y,x, exact = TRUE, conf.int =TRUE, conf.level = 1-alpha)


# example 4.8 -------------------------------------------------------------
# confidence interval large sample  
x = leukemia4$ge[leukemia4$label == "ALL"]
y = leukemia4$ge[leukemia4$label == "AML"]
D = sapply(y, function(y){y-x}) 
D.ord = sort(D)

m = length(x); n = length(y)
ew = n * (n+m + 1)/2 # W 통계량의 모평균
vw = n*m*(n+m+1)/12  # W 통계량의 모분산
eu = ew - n*(n+1)/2  # U 통계량의 모평균 
vu = vw              # U 통계량의 모분산 
(ci.limits = round( eu + 1.96 * sqrt(vu) *c(-1,1)))
(ci = D.ord[ci.limits])


wilcox.test(y,x, exact = FALSE, conf.int =TRUE, conf.level = 0.95)

# wilcox.test -------------------------------------------------------------
# 양측검정, 
wilcox.test(ge ~ label, exact = FALSE, conf.int = TRUE, data = leukemia4)
# exact = TRUE 로 하면 (동점이 없는 경우의 검정)
wilcox.test(ge ~ label, exact = TRUE, data = leukemia4)

x = leukemia4$ge[leukemia4$label == "ALL"]
y = leukemia4$ge[leukemia4$label == "AML"]
wilcox.test(y,x, exact = FALSE, correct = FALSE, alternative = "greater")




# 척도모수 --------------------------------------------------------------------


# fig 4.4 -----------------------------------------------------------------
x = seq(0, 10, by = 0.1)
f = dnorm(x, 5, sd = 1)
g = dnorm(x, 5, sd = 2)
data.frame(x, f, g) %>% pivot_longer(2:3, names_to = "Population",values_to = "density") %>% 
  ggplot(aes(x =x , y = density, color = Population, linetype = Population)) +
  geom_line()
ggsave(filename = "ch4_scale0.pdf", path = 'images', width = 4.8, height = 3) # smaller output
ggsave(filename = "ch4_scale0.png", path = 'images', width = 4.8, height = 3) # smaller output




# fig 4.5 -----------------------------------------------------------------
# 카레공장
x = c(110, 115, 120, 105, 103, 117)
y = c(109, 118, 113, 112, 119)
dat = data.frame(label = c(rep("Old",6),rep("New",5)),x = c(x,y), rank = rank(c(x,y)))
dat2 = cbind(dat %>% arrange(x),AB = c(1:6, 5:1), ST = c(1,4,5,8,9,11,10,7,6,3,2)) %>% mutate(Mood = (rank-6)^2)
dat2 %>% pivot_longer(4:6, names_to = "Type", values_to = "score") %>% 
  ggplot(aes(x = rank, y =score, color = Type, linetype = Type)) + 
  geom_line()
ggsave(filename = "ch4_scale.pdf", path = 'images', width = 4.8, height = 3) # smaller output
ggsave(filename = "ch4_scale.png", path = 'images', width = 4.8, height = 3) # smaller output


# example 4.11 ------------------------------------------------------------


# 카레공장 순열검정

N= 11
m = 6
n = 5
choose(11,5)

x = c(110, 115, 120, 105, 103, 117)
y = c(109, 118, 113, 112, 119)
curry = data.frame(label = c(rep("Old",6),rep("New",5)),x = c(x,y),rank = rank(c(x,y)))
curry2 = cbind(curry %>% arrange(x),
             AB = c(1:6, 5:1), 
             ST = c(1,4,5,8,9,11,10,7,6,3,2)) %>% 
  mutate(Mood = (rank-6)^2)
head(curry2,3)

stats.obs = curry2 %>% filter(label == "New") %>% summarise(sum(AB),sum(ST),sum(Mood)) %>% unlist()
stats.obs
N = 11; n = 5
allperms = gtools::combinations(N,n) # perm. ids for "y" group
stats.perm = apply(allperms,1,
      function (perm.y.id){ 
        curry2[perm.y.id,] %>% summarise(sum(AB),sum(ST),sum(Mood)) %>% unlist()
      } )
p.AB = mean(stats.perm[1,] >= stats.obs[1])   # the lerger the extreme
p.ST = mean(stats.perm[2,] >= stats.obs[2])   
p.Mood = mean(stats.perm[3,] <= stats.obs[3]) # the smaller the extreme
(p.vals = c(p.AB = p.AB, p.ST = p.ST, p.Mood = p.Mood))

g1 = data.frame(AB = stats.perm[1,]) %>% ggplot(aes(x = AB)) +  geom_histogram(binwidth = 0.5,color = "white") + 
  geom_dotplot(binwidth = .01) + geom_vline(xintercept = stats.obs[1], color = "red") 
g2 = data.frame(ST = stats.perm[2,]) %>% ggplot(aes(x = ST)) +  geom_histogram(binwidth = 1,color = "white") + 
  geom_dotplot(binwidth = .01) + geom_vline(xintercept = stats.obs[2], color = "red") 
g3 = data.frame(Mood = stats.perm[3,]) %>% ggplot(aes(x = Mood)) +  geom_histogram(binwidth = 4,color = "white") + 
  geom_dotplot(binwidth = .01) + geom_vline(xintercept = stats.obs[3], color = "red") 
library(cowplot)
plot_grid(g1,g2,g3, nrow = 1, align = "v")
ggsave(filename = "ch4_curry.pdf", path = 'images', width = 6.5, height = 3) # three-panel
ggsave(filename = "ch4_curry.png", path = 'images', width = 6.5, height = 3) # three-panel



# example 4.12 ------------------------------------------------------------
# equal location 
x = c(110, 115, 120, 105, 103, 117)
y = c(109, 118, 113, 112, 119)
ansari.test(y,x,alternative = "less")
DescTools::SiegelTukeyTest(y,x,alternative = "less",exact = TRUE)
mood.test(y,x,alternative = "less")

DescTools::SiegelTukeyRank(curry$x,curry$label)


# example 4.13 ------------------------------------------------------------
## unequal location
leukemia4= read.csv(file = "data/leukemia4.csv")
x = leukemia4$ge[leukemia4$label == "ALL"]
y = leukemia4$ge[leukemia4$label == "AML"]

var.test(x,y)$p.value
ansari.test(x,y)
ansari.test(x - median(x), y - median(y))
mood.test(x - median(x), y - median(y))
DescTools::SiegelTukeyTest(x,y,adjust.median = TRUE)
DescTools::SiegelTukeyTest(x - median(x), y - median(y))


# KS ----------------------------------------------------------------------

# fig 4.7 -----------------------------------------------------------------
set.seed(1)
x = rnorm(40)
pdf(file="images/ch4_ks1.pdf",width = 6.5, height = 4)
plot(ecdf(x))
lines(seq(-2.5,3,by = 0.1),pnorm(seq(-2.5,3,by = 0.1)),col = "red")
dev.off()


# fig 4.8 -----------------------------------------------------------------
pdf(file="images/ch4_ks2.pdf",width = 6.5, height = 4)
leukemia4= read.csv(file = "data/leukemia4.csv")
x = leukemia4$ge[leukemia4$label == "ALL"]
y = leukemia4$ge[leukemia4$label == "AML"]
plot(ecdf(x),main = "Empirical CDFs",xlim = c(-50,450))
lines(ecdf(y), col = "red")
legend(1, .95, legend=c("ALL", "AML"),
       col=c("black", "red"), lty=1, cex=0.8)
dev.off()



# example 4.14 ------------------------------------------------------------

ks.test(x,y)
out <- ks.test(x,y)
m = length(x)
n = length(y)
K = sqrt(m*n / (m+n)) * out$statistic

Kdistn = function(x){ 
  a = 0
  for (j in 1:200){
    a = a + (-1)^(j-1) * exp(-2 * j^2 * x^2)
  }
  1 - 2 * a
}
kseq = seq(0.01,2.5,by = 0.01)
data.frame(k = kseq,CDF = Kdistn(kseq)) %>% 
  ggplot(aes(k,CDF)) + geom_line() + geom_vline(xintercept = K, color = "red", linetype = "dashed")
ggsave(filename = "ch4_ks3.pdf", path = 'images', width = 4.8, height = 3) # smaller output
ggsave(filename = "ch4_ks3.png", path = 'images', width = 4.8, height = 3) # smaller output

