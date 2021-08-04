# nonparametric book by S. Jung
# chapter 5. association
#
#
# ----------------
library(tidyverse) 
library(gtools)
library(latex2exp)
library(showtext)
#font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()

# example 5.1 ----------------------------------------------------------------
##### Reading ability data set
##
ranks = 1:10  
score = c(10,20,25,50,80,82,85,90,95,100)
plot(ranks,score)
cor(ranks,score)
cor(x,y)
cor(x,y, method = "s")


# example 5.2 ----------------------------------------------------------------
# fig 5.1 ----------------------------------------------------------------
# data from https://databank.worldbank.org/
wh = read.csv("data/whappiness.csv")
wh %>% ggplot(aes(GDP.per.capita, Score)) + geom_point()
ggsave(filename = "ch5_happy.pdf", path = 'images', width = 4.8, height = 3) # smaller output
ggsave(filename = "ch5_happy.png", path = 'images', width = 4.8, height = 3)



# example 5.3 -------------------------------------------------------------

wh %>% ggplot(aes(GDP.per.capita, Score)) + geom_point() + 
  geom_smooth(method = "lm", se = F)
ggsave(filename = "ch5_happy2.pdf", path = 'images', width = 4.8, height = 3) # smaller output
ggsave(filename = "ch5_happy2.png", path = 'images', width = 4.8, height = 3)





# example 5.4 & fig 5.2 -------------------------------------------------------------
# reading ability plot --

ranks = 1:10  
score = c(10,20,25,50,80,82,85,90,95,100)
data.frame(rank = ranks, exam = score) %>% ggplot(aes(rank,exam)) + geom_point()
ggsave(filename = "ch5_read1.pdf", path = 'images', width = 4.8, height = 3) # smaller output
ggsave(filename = "ch5_read1.png", path = 'images', width = 4.8, height = 3)
cor(ranks,score)
cor(ranks,score, method = "spearman")
r = ranks
(s = rank(score))
cor(r,s)

# example 5.5 -------------------------------------------------------------
# happiness corr -
wh = read.csv("data/whappiness.csv")
cor(wh$Score,wh$GDP.per.capita, method = "pearson")
cor(wh$Score,wh$GDP.per.capita, method = "spearman")
wh = read.csv("data/whappiness.csv")
(r = rank(wh$Score))
(s = rank(wh$GDP.per.capita))
cor(r,s, method = "pearson")
cor(wh$Score,wh$GDP.per.capita, method = "spearman")



# example 5.6 -------------------------------------------------------------
# t.test
ranks = 1:10  
score = c(10,20,25,50,80,82,85,90,95,100)
n = 10
(r.obs = cor(ranks, score))
(t = r.obs * sqrt((n-2)/(1-r.obs^2)))
2*(1-pt(t, df = n-2))
cor.test(ranks,score)

# example 5.7 -------------------------------------------------------------
# rand perm test correlation -
n = 10
factorial(n)
r.obs = cor(ranks, score)

# rand perm.
M = 1999
set.seed(1)
r.perm = replicate(M, {
  perm.id = sample(1:n,n,replace = FALSE)
  cor(ranks, score[perm.id])})
(p.perm = 2*mean(c(r.obs,r.perm) >= r.obs))

data.frame(r.perm = r.perm) %>% ggplot(aes(x = r.perm)) +  geom_histogram(binwidth = 0.05,color = "white") + 
   geom_vline(xintercept = r.obs, color = "red", linetype = "dashed") +
  ggtitle("시험점수-평가순위. 표본상관계수의 랜덤순열분포")
ggsave(filename = "ch5_read2.pdf", path = 'images', width = 6.5, height = 4)  
ggsave(filename = "ch5_read2.png", path = 'images', width = 6.5, height = 4)



# fig 5.3 -----------------------------------------------------------------
# rand.perm by z-test  
data.frame(r.perm = r.perm) %>% ggplot(aes(x = r.perm)) +  
  geom_histogram(aes(y = stat(density)),binwidth = 0.05,color = "white")  +
  stat_function(
    fun = dnorm,
    args = list(mean = 0, sd = 1/3),
    lwd = 1, 
    col = 'red')+  ggtitle( "r의 순열분포 vs N(0,1/9)")
ggsave(filename = "ch5_read3.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch5_read3.png", path = 'images', width = 6.5, height = 4)
 



# 순위상관 --------------------------------------------------------------------


# omitted from the book ---------------------------------------------------
ranks = 1:10  
score = c(10,20,25,50,80,82,85,90,95,100)
n = 10
(rs.obs = cor(ranks, score, method = "spearman"))
rs.perm = replicate(1999, 
  cor(ranks, score[sample(1:n,n,replace = F)], method = "s")
)
(p.perm = 2*mean(c(rs.obs,rs.perm) >= rs.obs))



# example 5.8 -------------------------------------------------------------
n = nrow(wh)
r = rank(wh$Score)
s = rank(wh$GDP.per.capita)
(rs.obs = cor(r,s))
set.seed(1)
rs.perm = replicate(1999, cor(r, s[sample(1:n,n,replace = F)]))
(p.perm = 2*mean(c(rs.obs,rs.perm) >= rs.obs))
cor.test(r,s, method = "spearman")

data.frame(rs.perm = rs.perm) %>% ggplot(aes(x = rs.perm)) +  geom_histogram(binwidth = 0.05,color = "white") + 
  geom_dotplot(binwidth = .01) + geom_vline(xintercept = rs.obs, color = "red", linetype = "dashed") +
  ggtitle("GDP-행복도. 순위상관계수의 랜덤순열분포")
ggsave(filename = "ch5_wh4.pdf", path = 'images', width = 6.5, height = 4)  
ggsave(filename = "ch5_wh4.png", path = 'images', width = 6.5, height = 4)




# 켄달의 타우 ------------------------------------------------------------------

# example 5.9 -------------------------------------------------------------
n = nrow(wh)
allpairs = gtools::combinations(n,2,1:n)
head(allpairs,3)
U = apply(allpairs,1,function(pair.id){
  sign( diff(wh$Score[pair.id]) * diff(wh$GDP.per.capita[pair.id]) )
})
sum(U == 1)
sum(U == -1)
(tau = sum(U) / choose(n,2))

# due to the equal values, the following is not the same as tau
# 2* (sum(U == 1)+ sum(U ==0)*1/2) / choose(n,2) - 1 

cor(wh$Score,wh$GDP.per.capita, method = "kendall")




# fig 5.4 켄달의 타우 그림 ------------------------------------------------------------------
library(SuppDists)
n = 10
N = choose(n,2)
xgrid = 2 * (0:N) / (N) - 1
fktau = SuppDists::dKendall(xgrid,n)  # Kendall's tau dist'n

yy = dnorm(xgrid, mean = 0, sd = sqrt(4*n +10)/sqrt(9*n*(n-1))) / sum(dnorm(xgrid, mean = 0, sd = sqrt(4*n +10)/sqrt(9*n*(n-1))) )
ggplot(mapping = aes(x = perm.tau)) +  
  geom_col(aes(x = xgrid, y = fktau / sum(fktau)), color = "white", alpha = 0.2) + 
  geom_line(aes(x = xgrid, y = yy),color = "red") + xlab("켄달의 타우 통계량") + ylab("영분포 (질량/밀도함수)")
ggsave(filename = "ch5_kendalltau.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch5_kendalltau.png", path = 'images', width = 6.5, height = 4)




# example 5.10 ----------------------------------------------------------------------

n = nrow(wh)
r = rank(wh$Score)
s = rank(wh$GDP.per.capita)
(tau.obs = cor(r,s, method = "k"))

# (1) 랜덤순열
set.seed(1)
tau.perm = replicate(1999, cor(r, s[sample(1:n,n,replace = F)], method = "k"))
p.perm = 2*mean(c(tau.obs,tau.perm) >= tau.obs)

# (2) 켄달의 타우 분포
p.tau = 2*(1- SuppDists::pKendall(tau.obs,n))

# (3) 정규분포
p.z = 2*(1 - pnorm(tau,sd = sqrt(4*n +10)/sqrt(9*n*(n-1)) ))

c(p.perm, p.tau, p.z)



# last section ------------------------------------------------------------
# fig 5.5
# example 5.11

## three cases
library(mvtnorm)
library(cowplot)
n = 100
set.seed(1)
x = rmvnorm(n, sigma = matrix(c(1,0.9,0.9,1),nrow = 2))
g1= as.data.frame(x) %>% ggplot(aes(V1,V2)) + geom_point()+ geom_rug()+ ggtitle("(a)") + labs(x = NULL, y= NULL)
y = x
y[,1] = x[,1]^4 * sign(x[,1])
y[,2] = x[,2]^4 * sign(x[,2])
y[,1] = x[,1]^3
y[,2] = x[,2]^3

g2 = as.data.frame(y) %>% ggplot(aes(V1,V2)) + geom_point()+ geom_rug()+ ggtitle("(b)") + labs(x = NULL, y= NULL)
z = x

#z[,2] = (x[,2])^2  * sign(x[,2])
z[,2] = (x[,2])^2  * sign(x[,2])
z[,1] = -(x[,1]-3)^4
g3 = as.data.frame(z) %>% ggplot(aes(V1,V2)) + geom_point()+ geom_rug()+ ggtitle("(c)") + labs(x = NULL, y= NULL)

w = cbind(x[,1], sin(x[,1]*3)+rnorm(n,sd = 0.2))
g4 = as.data.frame(w) %>% ggplot(aes(V1,V2)) + geom_point()+ geom_rug() + ggtitle("(d)") + labs(x = NULL, y= NULL)
cowplot::plot_grid(g1,g2,g3,g4)
ggsave(filename = "ch5_summary.pdf", path = 'images', width = 6.5, height = 5)
ggsave(filename = "ch5_summary.png", path = 'images', width = 6.5, height = 5)


o1 = cor.test(x[,1],x[,2],method = "p")
o2 = cor.test(x[,1],x[,2],method = "s")
o3 = cor.test(x[,1],x[,2],method = "k")
call = c(o1$estimate, o1$p.value, o2$estimate, o2$p.value, o3$estimate, o3$p.value)


o1 = cor.test(y[,1],y[,2],method = "p")
o2 = cor.test(y[,1],y[,2],method = "s")
o3 = cor.test(y[,1],y[,2],method = "k")
call = rbind(call, 
             c(o1$estimate, o1$p.value, o2$estimate, o2$p.value, o3$estimate, o3$p.value))


o1 = cor.test(z[,1],z[,2],method = "p")
o2 = cor.test(z[,1],z[,2],method = "s")
o3 = cor.test(z[,1],z[,2],method = "k")
call = rbind(call, 
             c(o1$estimate, o1$p.value, o2$estimate, o2$p.value, o3$estimate, o3$p.value))

o1 = cor.test(w[,1],w[,2],method = "p")
o2 = cor.test(w[,1],w[,2],method = "s")
o3 = cor.test(w[,1],w[,2],method = "k")
call = rbind(call, 
             c(o1$estimate, o1$p.value, o2$estimate, o2$p.value, o3$estimate, o3$p.value))
rownames(call) = c("(a)","(b)","(c)","(d)")
call
format(round(call, 2))

