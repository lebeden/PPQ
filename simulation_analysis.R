library(tidyverse)
library(reshape2)

# a is the result of simsalapar for this experiment
summarizeSim <- function(a){


res <- array2df(getArray(a))
names(res)
dat <- res %>% dcast(formula = n+gamma+theta+s+n.sim~Var1)
dat <- data.frame(apply(dat,2,as.numeric))
dat <-  dat %>% rename("B_lambda" = B_gamma, "B_gamma" = B_lambda) # names were wrong in sim - fixed

comparer <- function(x,y){
  ans <- numeric(length(x))
  for (i in 1:length(x)){
    ans[i] <- all.equal(x[i],y[i])
  }
  ans <- ifelse(ans == "TRUE",0,1) # if they're equal - the algorithm didn't move
  return(ans)
}

boundary <- function(x,y){
    ans <- numeric(length(x))
    for (i in 1:length(x)){
      ans[i] <- all.equal(x[i],y[i]*2) # if the answer is half the original parameter - converged on boundary
    }
    ans <- ifelse(ans == "TRUE",1,0) # if x = 2y : boundary solution
    return(ans)

}
f_dat <- dat %>%
  mutate(converged = comparer(gamma,B_gamma) & comparer(theta,B_theta) & comparer(10,B_lambda),
         boundary = boundary(gamma,B_gamma) | boundary(theta,B_theta)) %>%
  mutate(liron_converged = comparer(L_lambda,500)) %>%
  select(-n.sim)
return(f_dat)


}

#full_res_mac <- summarizeSim(a) %>% bind_rows(summarizeSim(a1))
#write.csv(full_res_mac,"res_macbook.csv",row.names = FALSE)
#full_res_pc <- summarizeSim(a) %>% bind_rows(summarizeSim(a1))
#write.csv(full_res_pc,"res_pc.csv",row.names = FALSE)
#full_res_laptop <- summarizeSim(a) %>% bind_rows(summarizeSim(a1))
#write.csv(full_res_laptop,"res_laptop.csv",row.names = FALSE)

r1 <- read.csv("res_macbook.csv")
r2 <- read.csv("res_pc.csv")
r3 <- read.csv("res_laptop.csv")
full_res %>% View
full_res <- rbind(r1,r2,r3)

full_res %>%
  filter(converged & !boundary & liron_converged) %>%
  mutate(biasB = B_theta - theta,
         biasL = L_theta - theta

  ) %>%
  group_by(gamma,theta,s) %>%
  summarize(biasB = mean(biasB),
            biasL = mean(biasL),
            varB = var(B_theta),
            varL = var(L_theta),
            rho_mean = mean((gamma/2 + 10)/s)) %>%
  View# the mean is redundant for rho, just for technical reasons


full_res %>%
  group_by(gamma,theta,s) %>%
  summarise(boris_success_rate = mean(converged),
            boris_boundary_rate = mean(boundary),
            rho_mean = mean((gamma/2 + 10)/s),
            n.sim = n()) %>%
  round(3) %>%
  View("results of some simulations")

