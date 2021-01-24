#### File of functions & datasets ####
#### This file accumulates functions and datasets for future use

######## Datasets ########
#### CH 1-3
kidney <- read.table("https://web.stanford.edu/~hastie/CASI_files/DATA/kidney.txt", header=T)

leukemia_big <- read.csv("http://web.stanford.edu/~hastie/CASI_files/DATA/leukemia_big.csv")
leukemia_big <- data.frame(t(leukemia_big))
leukemia_big$group <- sapply(rownames(leukemia_big), function(word){substr(word, 1, 3)})

score <- read.table("https://web.stanford.edu/~hastie/CASI_files/DATA/student_score.txt",
                    header = T)

#### CH 4 ####
gfr <- read.table("https://web.stanford.edu/~hastie/CASI_files/DATA/gfr.txt")$V1

######## Functions ########
#### CH 1-3
slr_se <- function(newx, lm_fit){
  n <- length(lm_fit$residuals)
  sigma_hat <- sqrt(sum((lm_fit$residuals)^2)/n)
  x <- lm_fit$model[,2]
  xbar <- mean(x)
  ssx <- sum((x - xbar)^2)
  sigma_hat * sqrt(1/n + (newx - xbar)^2/ssx)
}

get_t <- function(d){
  x1 <- d %>% filter(group == "AML") %>% pull(1)
  x2 <- d %>% filter(group == "ALL") %>% pull(1)
  sigma_hat <- (sum((x1 - mean(x1))^2) + sum((x2 - mean(x2))^2))/70
  #print(x1)
  sd_hat <- sqrt(sigma_hat * (1/25 + 1/47))
  t <- (mean(x1) - mean(x2))/sd_hat
  return(t)
}


