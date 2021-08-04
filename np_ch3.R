# nonparametric book by S. Jung
# chapter 3. permutation tests
#
#
# ----------------
library(tidyverse) 
library(gtools)
library(showtext)
#font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()


# example 3.1 -------------------------------------------------------------
N = 7
n = 4
m = 3
x = c(19, 22, 24, 26)
y = c(23, 33, 40)
dat = c(x,y)


# example 3.2 -------------------------------------------------------------
tsetse = read.csv("data/tsetse.csv")
ggplot(tsetse, aes(y = AL, x= label)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis = "y",stackdir = "center")


# Table 3.1 --------------------------------------------------------------- 
allperms <- gtools::combinations(7,4,dat)
allperms_data = t(apply(allperms, 1, function(x){c(x, setdiff(dat,x))}))
allperms_data = t(apply(allperms_data, 1, 
                        function(x){
                          c(x, mean(x[5:7]) - mean(x[1:4]))
                        }))
colnames(allperms_data) = c(rep("new",4),rep("old",3),"D")
nrow(allperms_data)
choose(7,4)
allperms_data


# Fig 3.1 -----------------------------------------------------------------
D.perm = allperms_data[,8]


data.frame(D.perm = allperms_data[,8]) %>% ggplot(aes(x = D.perm)) +  
  geom_histogram(binwidth = 4,color = "white") + 
  scale_y_continuous(breaks = c(0,2,4,)) + 
  geom_dotplot(binwidth = .5) + 
  geom_vline(xintercept = 9.25, color = "red", linetype = "dashed") 
ggsave(filename = "ch3_Dperm.pdf", path = 'images', width = 4.8, height = 3) # smaller output
ggsave(filename = "ch3_Dperm.png", path = 'images', width = 4.8, height = 3)
 

# example 3.3 & 3.4 ------------------------------------------------------------
Dobs = mean(y) - mean(x)

(p.right = mean(D.perm >= Dobs)) # H1: mu(y) > mu(x)
sum(D.perm >= Dobs)

(p.left = mean(D.perm <= Dobs)) # H1: mu(y) < mu(x)
sum(D.perm <= Dobs)

 (p.both = mean(abs(D.perm) >= abs(Dobs))) # H1: mu(y) != mu(x) 
# (p.both = min(p.right,p.left) *2)# correct if the null distribution is symmetric
p.both * choose(7,4)


# section 3.2 -------------------------------------------------------------
# computations of D = ybar - xbar, t, trimmed mean difference, median difference
# omitted....


# section 3.3 -------------------------------------------------------------
# tsetse fly example
N = nrow(tsetse)
m = sum(tsetse$label == "forest")
n = sum(tsetse$label == "river")
choose(N,m)
choose(N,n)

# a permutation
set.seed(0)
permuted.label = sample(c(rep("forest",9),rep("river",6)), 15, replace = FALSE)
cbind(tsetse, permuted.label)


# fig 3.2 -----------------------------------------------------------------
t.obs = t.test(AL ~ label,data = tsetse, var.equal = TRUE)$statistic 
allperms <- gtools::combinations(N,m)
t.perm = apply(allperms,1,
                 function (perm.x.id){ 
                    permuted.x = tsetse$AL[perm.x.id]
                    permuted.y = tsetse$AL[-perm.x.id]
                    t.test(permuted.x, permuted.y, var.equal = T)$statistic} )

data.frame(t.perm = t.perm) %>% 
  ggplot(aes(x = t.perm)) +  
  geom_histogram(binwidth = 0.5,color = "white") + 
  geom_dotplot(binwidth = .01) + 
  geom_vline(xintercept = t.obs, color = "red",linetype = "dashed") 
ggsave(filename = "ch3_tperm_tsetset.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch3_tperm_tsetset.png", path = 'images', width = 6.5, height = 4)

sum(t.perm >= t.obs)
sum(abs(t.perm) >= abs(t.obs))
mean(abs(t.perm) >= abs(t.obs))

# section 3.3.2 -----------------------------------------------------------

## random permutation
N = 7
rand.permutation = sample(1:N, N, replace = FALSE)
rand.permutation


N = 7 
label = c("New","New","New","New","Old","Old","Old")
set.seed(0)
perm.id = sample(1:N, N, replace = FALSE)
perm.id
label[perm.id]
set.seed(0)
permuted.label = sample(label, N, replace = FALSE)
permuted.label



# fig 3.3 -----------------------------------------------------------------
# alternatives
library(tidyverse)
x = seq(-3,3,by = 0.01) 
library(cowplot)
a1 <- tibble(x, F = pnorm(x,0,1),
             G = pnorm(x,1,1)) %>% pivot_longer(2:3, names_to = "Popn", values_to = "Prob") %>% 
  ggplot(aes(x = x, y = Prob,color = Popn)) + geom_line() + ggtitle("F(x) > G(x)")
a2 <- tibble(x, F = dnorm(x,0,1),
             G = dnorm(x,1,1)) %>% pivot_longer(2:3, names_to = "Popn", values_to = "Density") %>% 
  ggplot(aes(x = x, y = Density,color = Popn)) + geom_line() + ggtitle("F(x) > G(x)")
cowplot::plot_grid(a1,a2)
ggsave(filename = "ch3_twosamplediff.pdf", path = 'images', width = 6.5, height = 3)
ggsave(filename = "ch3_twosamplediff.png", path = 'images', width = 6.5, height = 3) # two-panel



# section 3.6.1 -----------------------------------------------------------
dat = read.csv("data/tsetse.csv"); N = 15; m = 9; n = 6
head(dat,3)
allperms <- gtools::combinations(N,m)
head(allperms,3)
nrow(allperms)
D.perm = vector(length = choose(N,m))
for (i in 1:choose(N,m)){
  permuted.x = dat$AL[allperms[i,]]
  permuted.y = dat$AL[-allperms[i,]]
  D.perm[i] = mean(permuted.y) - mean(permuted.x)
  }
head(D.perm,3)
D.perm = apply(allperms,1,
               function (perm.x.id){
                 permuted.x = dat$AL[perm.x.id]
                 permuted.y = dat$AL[-perm.x.id]
                 mean(permuted.y) - mean(permuted.x)
               } )
head(D.perm,3)
Dobs = with(dat,
            mean(AL[label == "river"]) - mean(AL[label == "forest"]))
(p = mean(abs(D.perm) >= abs(Dobs)))


# random permutations
M = 1999
D.perm = replicate(M, {
  x.id = sample(1:N,m, replace =FALSE)
  permuted.x = dat$AL[x.id]
  permuted.y = dat$AL[-x.id]
  mean(permuted.y) - mean(permuted.x)
})
(p = mean(c(Dobs,D.perm) <= Dobs))                        
                      
# ex 3.8 ------------------------------------------------------------------
CVD.salt.tab <- tribble(
  ~Cardio, ~Salt, ~Freq,
  #------------------
  "Yes", "H", 5,
  "Yes", "L", 25,
  "No", "H", 2,
  "No", "L", 28,
) 

library(DescTools)
data = Untable(xtabs(Freq ~ Cardio + Salt, data = CVD.salt.tab))
choose(60,7)
fisher.test(xtabs(Freq ~ Cardio + Salt, data = CVD.salt.tab))
head(data,10)


# example 3.9 -------------------------------------------------------------

nulldist = dhyper(0:7, m = 30, n = 30, k = 7)
round(nulldist,2)
data.frame(d = 0:7, nulldist) %>% ggplot(aes(x = d, y = nulldist)) + 
  geom_col() + ggtitle("Fisher's exact test, null distribution") 
ggsave(filename = "ch3_fisher.pdf", path = 'images', width = 4.8, height = 3) 
ggsave(filename = "ch3_fisher.png", path = 'images', width = 4.8, height = 3) # smaller output
# 

