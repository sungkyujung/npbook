# nonparametric book by S. Jung
# chapter 1. introduction
#
#
# ----------------
library(tidyverse)
library(lubridate)  # needed for as_date()
library(showtext)
#font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()



# Biden data --------------------------------------------------------------
# https://elections.huffingtonpost.com/pollster/joe-biden-favorable-rating
# 
biden <- read.table("data/US-Biden-Joe-Favorable-poll-responses-clean.tsv",
                sep = "\t",header = T)
qplot(as_date(start_date), Favorable, data = biden )
ggsave(filename = "ch0_biden.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch0_biden.png", path = 'images', width = 6.5, height = 4)


# ngram -------------------------------------------------------------------
library(ngramr)
data <- as.data.frame(matrix(ncol=1, nrow=100))
data$V1 <- seq(from=1919, to=2018)

names(data)[names(data)=="V1"] <- "Year"
search_terms <- c("nonparametric statistics",
                  "nonparametric regression",
                  "permutation test","bootstrapping")
# Create a loop
for(i in 1:length(search_terms)){
  
  # Get each search term and store those in objects
  term <- search_terms[i]
  
  # Search for the term in the English 2012 corpus, starting from the year 1900 to 2008
  # Then house the output in a dataframe
  temp <- ngram(term, corpus = "eng_2019", year_start = 1919, year_end = 2018, smoothing = 3, 
                case_ins = T, aggregate = T, count = T)
  
  # Merge NYT data with dataframe created step 1, matching by years
  data <- merge(data, temp[,c("Year", "Count")], by ="Year", all.x=TRUE)
  
  # Reaname column by search term
  colname <- paste(term, sep="")
  
  # Rename added column with ID
  names(data)[names(data)=="Count"] <- colname
  
  # Remove temporary dataframe
  rm(temp)
  
}
tail(data)

data %>% pivot_longer(2:5, names_to = "ngram", values_to = "count") %>% 
  ggplot(aes(Year, count, color = ngram, linetype = ngram)) + 
  geom_line(size = 1)
ggsave(filename = "ch0_ngram.pdf", path = 'images', width = 6.5, height = 4)
ggsave(filename = "ch0_ngram.png", path = 'images', width = 6.5, height = 4)



# normal and binomial density ---------------------------------------------
# lost ....



# example 1.1 & 1.2 -------------------------------------------------------------
x = c(3, 2, 0, 8, 12, 11, 5, -2)
mean(x); sd(x)
n = length(x); 
(ci = mean(x) + c(-1,1) * qt(1-0.025,n-1) * sd(x) / sqrt(n))
t.test(x, alternative = "greater")$p.value



# example 1.4 -------------------------------------------------------------
x = c(2.20, 2.20, 2.40, 2.40, 2.50, 2.70, 2.80, 2.90, 
      3.03, 3.03, 3.10, 3.37, 3.40, 3.40, 3.40, 3.50, 
      3.60, 3.70, 3.70, 3.70, 3.70, 3.77, 5.28, 28.50)
g <- data.frame(x = x) %>% ggplot(aes(x=x)) 
g + geom_histogram(aes(y = ..density..), color = "white") + 
  geom_dotplot(binwidth = 0.2)
mean(x)
sd(x)
y = x[1:23]
mean(y)
sd(y)



