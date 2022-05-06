library(lme4)
library(arm)

## Multi-level power calculations from Gelman's ARM book:
BHET.fake <- function (mu_b3, mu_g3, mu_g4){
  wave <- rep (seq(1,3,1), 19*50) 
  person <- rep (1:19, each=3)           # person ID’s
  village <- rep (1:50, each=19*3)
  treatment <- sample (rep (0:1, 50/2))
  treated <- treatment[village]
  post <- ifelse(wave > 1, 1, 0)
  policy <- ifelse(post==1 & treated==1, 1, 0)

  # # hyperparameters:
  mu.b0.true <- 0        # more generally, these could
  sigma.b0.true <- 0.1
  
  mu.b1.true <- 0.25       # be specified as additional
  mu.b2.true <- 0.25       # arguments to the function
  mu.b3.true <- mu_b3
  
  mu.g0.true <- 0
  sigma.g0.true <- 0.1
  mu.g1.true <- 0.25
  mu.g2.true <- 0.25
  mu.g3.true <- mu_g3
  mu.g4.true <- mu_g4
  
  sigma.m.true <- 1
  sigma.y.true <- 1

  # # village-level parameters
  b0.true <- rnorm(50, mu.b0.true, sigma.b0.true) 
  b1.true <- rnorm(50, mu.b1.true, 0)
  b2.true <- rnorm(50, mu.b2.true, 0)
  b3.true <- rnorm(50, mu.b3.true, 0)
  
  g0.true <- rnorm(50, mu.g0.true, sigma.g0.true) 
  g1.true <- rnorm(50, mu.g1.true, 0)
  g2.true <- rnorm(50, mu.g2.true, 0)
  g3.true <- rnorm(50, mu.g3.true, 0)
  g4.true <- rnorm(50, mu.g4.true, 0)

  # # data
  mm <- b0.true[village] + b1.true[village]*treated
     + b2.true[village]*post + b3.true[village]*policy
  s <- sigma.m.true
  location <- log(mm^2 / sqrt(s^2 + mm^2))
  shape <- sqrt(log(1 + (s^2 / mm^2)))

  m <- rnorm (50*19*3, b0.true[village] + b1.true[village]*treated
    + b2.true[village]*post + b3.true[village]*policy, sigma.m.true)
  
  # m <- rlnorm(50*19*3, location, shape)
  
  y <- rnorm (50*19*3, g0.true[village] + g1.true[village]*treated 
    + g2.true[village]*post + g3.true[village]*policy 
    + g4.true[village]*m, sigma.y.true)
  
  return (data.frame (y, m, wave, person, village, treated, post, policy))
}

# power calculation function
BHET.power <- function (mu_b3, mu_g3, mu_g4, n.sims=100){
  te.signif <- rep (NA, n.sims)
  cde.signif <- rep (NA, n.sims)
  for (s in 1:n.sims){
    fake <- BHET.fake (mu_b3, mu_g3, mu_g4)
    lme.power.te <- lmer (y ~ treated + post + policy +
      (1 | village), data=fake)
    te.hat <- fixef(lme.power.te)["policy"]
    te.se <- se.fixef(lme.power.te)["policy"]
    te.signif[s] <- (te.hat - 2*te.se) > 0 # returns TRUE or FALSE
    lme.power.cde <- lmer (y ~ treated + post + policy + m + 
                             (1 | village), data=fake)
    cde.hat <- fixef(lme.power.cde)["policy"]
    cde.se <- se.fixef(lme.power.cde)["policy"]
    cde.signif[s] <- (cde.hat - 2*cde.se) > 0 # returns TRUE or FALSE
  }
power.te <- mean (te.signif) 
power.cde <- mean (cde.signif) # proportion of TRUE
cat("TE power:", power.te , "CDE power:", power.cde)
#return(list(power.te = power.te, power.cde = power.cde))
	#p.j <- power1*power2
	#return (p.j)
return (power.cde)
}

BHET.power (mu_b3=1, mu_g3=0.5, mu_g4=0.5, n.sims=100)

# 3.  Loop for power calculations

g3.values <- seq( from=0 , to=1 , length.out=21 )
b3.values <- c(0.1, 0.25, 0.5)
power.values <- array (NA, c(length(g3.values),length(b3.values)))
for (i1 in 1:length(g3.values)){
  for (i2 in 1:length(b3.values)){
    cat ("computing power calculation for g3 =", g3.values[i1], ", b3 =", b3.values[i2])
    power.values[i1,i2] <- BHET.power (mu_g3=g3.values[i1], mu_b3=b3.values[i2], mu_g4=0.2, n.sims=1000)
    cat ("power =", power.values[i1,i2])
  }
}

png("power-plot.png") 
plot( NULL , xlim=range(g3.values) , ylim=c(0,1) , xlab="Controlled direct effect (SD units)" , ylab="Power", cex.lab=1.2, cex.axis=1.2)
abline(h=0.8, col="black", lty=2)
lines( g3.values , power.values[,1], lwd=1.5)
lines( g3.values , power.values[,2], col="red", lwd=1.5)
lines( g3.values , power.values[,3], col="blue", lwd=1.5)
legend(0.5, 0.5, title = "Policy mediator coefficient:", 
  legend=c("T → M = 0.1", "T → M = 0.25",   "T → M = 0.5"), 
  col=c("black","red", "blue"), lty=c(1,1,1), lwd=2, cex=1.2)
dev.off()

plot(x, y1, type="b", pch=19, col="red", xlab="x", ylab="y")
# Add a line
lines(x, y2, pch=18, col="blue", type="b", lty=2)
# Add a legend
legend(1, 95, legend=c("Line 1", "Line 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8)


g3.values <- seq( from=0 , to=1 , length.out=11 )
power.values <- array (NA, 1)
{
    cat ("computing power calculation for g3 =", 0.1)
    power.values[0.1] <- BHET.power (mu_g3=0.1, mu_b3=0.05, mu_g4=0.2, n.sims=10)
    cat ("power =", power.values[0.1])
  }


# plot all the curves

plot (c(0,max(J.values)), c(0,1), xaxs="i", yaxs="i", xlab="number of children", ylab="power", type="n")
for (i2 in 1:length(K.values)){
  lines (c(0,J.values), c(.025,power.values[,i2]))
}

# just plot the curve for K=7

postscript ("c:/books/multilevel/power.ps", height=3.6, width=4)
plot (c(0,J.values), c(.025,power.values[,3]), ylim=c(0,1), yaxs="i", xaxs="i", type="l",
  xlab="number of children", ylab="power", yaxt="n", mgp=c(1.6,.5,0))
axis (2, seq(0,1,.1),rep("",11), mgp=c(1.6,.5,0), tcl=-.2)
axis (2, seq(0,1,.5),c("0","0.5","1"), mgp=c(1.6,.5,0))
dev.off ()



villages <- 50
vid <- seq(1:50)
u_v0 <- rnorm(50,0,1)
u_v1 <- rnorm(50,0,1)
u_v2 <- rnorm(50,0,1)

treatment <- sample (rep (0:1, villages/2))
treated <- treatment[vid]

v <- data.frame(vid, u_v0, u_v1, u_v2, treated)

vl <- reshape(data=v, idvar=c("vid","treated"), 
  varying = c("u_v0","u_v1","u_v2"), v.name=c("u_v"), times=c("0","1","2"), 
  new.row.names = 1:1000, direction="long")

bysort village: generate pid = 1000*village + _n

person <- rep (1000*vid + 1:19, each=19)


m <- 130
s <- 75
location <- log(m^2 / sqrt(s^2 + m^2))
shape <- sqrt(log(1 + (s^2 / m^2)))
print(paste("location:", location))
print(paste("shape:", shape))
draws3 <- rlnorm(n=1000000, location, shape)



