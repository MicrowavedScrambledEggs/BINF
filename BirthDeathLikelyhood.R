library(ape)
finches <- read.nexus("http://www.r-phylo.org/w/images/0/02/Geospiza.nex")

## Birth death model likelyhood function
# Equation from Nee et al 1994 The reconstructed evolutionary process (model 21)
# Takes a birth rate (brate), a death rate (drate), and an object of class phylo 
# Outputs the likelyhood of the tree with the given parameters under the birth death model
bd.likelyhood.function = function(brate, drate, phy)
{
  N = length(phy$tip.label) # Number of lineages
  # times since each lineage's birth
  # Trees start from the birth of the second lineage so we put NA for the first
  t = c(NA, branching.times(phy)) 
  
  a = drate / brate # extinction rate
  r = brate - drate # speciation rate
  
  # P(t, T) Likelyhood of a lineage alive at time t is still alive at later time T
  # r /(r*exp(-r*(T-t)))
  
  # Likelyhood that the ith extant lineage to be born was born at time t_i
  # (i - 1)*brate*P(t_i, T)
  # The likelyhood that the 1st and 2nd lineages were born is 1 as otherwise we wouldn't have a tree
  
  # Likelyhood that a lineage does not give birth to another lineage within time period t
  # (1 - u_t) where u_t is the likelyhood that a lineage DOES give birth to another lineage
  # within time t
  # u_t = (brate*(1-exp(-r*t)))/(brate - drate*exp(-r*t))
  
  # So the likelyhood of a tree with exant tips is the likelyhood that each of the lineages were born 
  # when the tree says they were born times the likelyhood that each lineage did not give birth to 
  # extra lineages times the likelyhood that the first 2 lineages are extant
  
  # Using algebra this works out to:
  # likelyhood = (N-1)!*r^(N-2)*exp(r*sum(t[3:N]))*(1-a)^N*prod(1/(exp(r*t[2:N])-a)**2)
  
  cat("b:", brate, "d:", drate, "N:", N, "a", a, "r", r, "\n")
  if(r < 0 || a > 1) # if drate greater than birth rate return a really small number
    {likelyhood = 1e-100}
  else
    {likelyhood = factorial(N-1)*r**(N-2)*exp(r*sum(t[3:N]))*(1-a)**N*prod(1/(exp(r*t[2:N])-a)**2)}
  cat("likelyhood", likelyhood, "\n")
  return(likelyhood)
}

## Birth death model likelyhood function
# Equation is the log form of model 21 from Nee et al 1994 "The reconstructed evolutionary process"
# Takes a birth rate (brate), a death rate (drate), and an object of class phylo 
# Outputs the log likelyhood of the tree with the given parameters under the birth death model
bd.log.likelyhood.function <- function(brate, drate, phy)
{
  N = length(phy$tip.label) # Number of lineages
  # times since each lineage's birth
  # Trees start from the birth of the second lineage so we put NA for the first
  t = c(NA, branching.times(phy)) 
  
  a = drate / brate # extinction rate
  r = brate - drate # speciation rate
  
  cat("b:", brate, "d:", drate, "N:", N, "a", a, "r", r, "\n")
  if(r < 0 || a > 1) # if drate greater than birth rate return a really small number
    {loglikelyhood = log(1e-100)}
  else
    {loglikelyhood = lfactorial(N-1)+log(r)*(N-2)+r*sum(t[3:N])+log(1-a)*N-2*sum(log(exp(r*t[2:N])-a))}
  cat("log likelyhood", loglikelyhood, "\n")
  return(loglikelyhood)
}

# Test the log likelihood matches the output of the function birthdeath() in the package ape 
# Birth death finds the birthrate and deathrate via maximum likelyhood so we will use optim
# on my function to find the birthrate and the death rate that produce the highest log likelihood
my_max_like_bd <- optim(c(2.5,2.25), function(p) bd.log.likelyhood.function(p[1],p[2],finches),
      lower=c(0,0), control = list(fnscale=-1))

print(my_max_like_bd)
cat("Maximum Log likelihood found by my function: ", my_max_like_bd$value, "\n")
cat("Birth Rate: ", my_max_like_bd$par[1], "\tDeath Rate: ", my_max_like_bd$par[2], "\n")

# Check ape birth death maximum likelihood and parameters
ape_max_like_bd <- birthdeath(phy = finches)
print(ape_max_like_bd)

# ape birthdeath objects only store the d / b and the b - d so to compare we will have to do some
# maths
cat("My max likelyhood birth death model:\nd / b = ", my_max_like_bd$par[2], " / " , my_max_like_bd$par[1],
    " = ", my_max_like_bd$par[2] / my_max_like_bd$par[1], "\nb - d = ", my_max_like_bd$par[1], " - ", 
    my_max_like_bd$par[2], " = ", my_max_like_bd$par[1] - my_max_like_bd$par[2], "\n")
cat("APE max likelyhood birth death model:\n")
print(ape_max_like_bd$para)




