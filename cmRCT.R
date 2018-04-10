###########################################
# Simulation code dependence cmRCT trials #
# Author: R.H.H. Groenwold                #
###########################################

### dependencies ###
require(ggplot2)
### ------------ ###

# ====================================
# function to round matrix, d is vector with the number of decimal places per column
rnd <- function(X,d){
	X.t <- X
	for(i in 1:length(d)){
		X.t[,i] <- round(X[,i],d[i])
		}
	return(X.t)}

# function to provide descriptive statistics
smrz <- function(x) c(mean(x),sd(x)/sqrt(length(x)))	

# ===================================
# simulation function:

sim.fun <- function(n.sim, n, block, R, delta.1, delta.2){
	RESULTS.t <- array(NA, dim=c(n.sim,8,length(R)) )
	RESULTS.m <- array(NA, dim=c(n.sim,8,length(R)) )

	for (j in 1:length(R)){
		for (i in 1:n.sim){

		if(block==TRUE){
			n10 <- n
			n11 <- n
			n20 <- n
			n21 <- n}
		else{
			n10 <- rbinom(1,2*n,.5)
			n11 <- 2*n-n10
			n20 <- rbinom(1,2*n,.5)
			n21 <- 2*n-n20}

		# option 1.	subjects can only participate in 1 trial
			N <- round(4*n/R[j])	
			ids <- 1:N 
			dat <- rnorm(N)
				sample.10 <- sample(ids,n10)
				used.id <- sample.10
				sample.11 <- sample(ids[-used.id],n11)
				used.id <- c(used.id,sample.11)
				sample.20 <- sample(ids[-used.id],n20)
				used.id <- c(used.id,sample.20)
				sample.21 <- sample(ids[-used.id],n21)

				used.id <- c(used.id,sample.21)
				dat[sample.11] <- dat[sample.11] + delta.1
				dat[sample.21] <- dat[sample.21] + delta.2
				RESULTS.t[i,1,j]  <- t.test(dat[sample.11],dat[sample.10], var.equal=TRUE)$statistic
				RESULTS.t[i,2,j]  <- t.test(dat[sample.21],dat[sample.20], var.equal=TRUE)$statistic
				RESULTS.m[i,1,j] 	<- mean(dat[sample.11]) - mean(dat[sample.10])
				RESULTS.m[i,2,j] 	<- mean(dat[sample.21]) - mean(dat[sample.20])

		# option 2. control subjects can participate in multiple trials
		# treated subjects only in one trial
			N <- round(3*n/R[j])	
			ids <- 1:N 
			dat <- rnorm(N)
				sample.11 	<- sample(ids,n11)
				sample.21 	<- sample(ids[-sample.11],n21)
				used.id 	<- c(sample.11,sample.21)
				sample.10 	<- sample(ids[-used.id],n10)
				sample.20 	<- sample(ids[-used.id],n20)
				dat[sample.11] <- dat[sample.11] + delta.1
				dat[sample.21] <- dat[sample.21] + delta.2
				RESULTS.t[i,3,j]  <- t.test(dat[sample.11],dat[sample.10], var.equal=TRUE)$statistic
				RESULTS.t[i,4,j]  <- t.test(dat[sample.21],dat[sample.20], var.equal=TRUE)$statistic
				RESULTS.m[i,3,j] 	<- mean(dat[sample.11]) - mean(dat[sample.10])
				RESULTS.m[i,4,j] 	<- mean(dat[sample.21]) - mean(dat[sample.20])

		# option 3. subjects can participate in multiple trials, 
		# but only in 1 trial can they receive treatment
			N <- round(2*n*(1+sqrt(1-R[j]))/R[j])	
			ids <- 1:N 
			dat <- rnorm(N)
				sample.11 	<- sample(ids,n11)
				sample.21 	<- sample(ids[-sample.11],n21)
				sample.10 	<- sample(ids[-sample.11],n10)
				sample.20 	<- sample(ids[-sample.21],n20)
				dat[sample.11] <- dat[sample.11] + delta.1
				dat[sample.21] <- dat[sample.21] + delta.2
				RESULTS.t[i,5,j]  <- t.test(dat[sample.11],dat[sample.10], var.equal=TRUE)$statistic
				RESULTS.t[i,6,j]  <- t.test(dat[sample.21],dat[sample.20], var.equal=TRUE)$statistic
				RESULTS.m[i,5,j] 	<- mean(dat[sample.11]) - mean(dat[sample.10])
				RESULTS.m[i,6,j] 	<- mean(dat[sample.21]) - mean(dat[sample.20])
	
		# option 4. all subjects can participate in multiple trials
			N <- round(2*n/R[j])	
			ids <- 1:N 
			dat <- rnorm(N)
				sample.10 <- sample(ids,n10)
				sample.11 <- sample(ids[-sample.10],n11)
				sample.20 <- sample(ids,n20)
				sample.21 <- sample(ids[-sample.20],n21)
				dat[sample.11] <- dat[sample.11] + delta.1
				dat[sample.21] <- dat[sample.21] + delta.2
				RESULTS.t[i,7,j]  <- t.test(dat[sample.11],dat[sample.10], var.equal=TRUE)$statistic
				RESULTS.t[i,8,j]  <- t.test(dat[sample.21],dat[sample.20], var.equal=TRUE)$statistic
				RESULTS.m[i,7,j] 	<- mean(dat[sample.11]) - mean(dat[sample.10])
				RESULTS.m[i,8,j] 	<- mean(dat[sample.21]) - mean(dat[sample.20])
			}
		}
	output <- list(t = RESULTS.t, m = RESULTS.m)	
	return(output)
	}

# ====================================================================
# === SCENARIO 1 ===
# no treatment effect in trial A (delta.1 = 0)
# no treatment effect in trial B (delta.2 = 0)

set.seed(2505)
n.sim 	<- 1e5				# no. simulations
n 		<- 50					# sample size per treatment arm
R		<- seq(.05,.95, by=0.05) 	# proportion of cohort included in at least on of two trials
delta.1 	<- 0					# treatment effect in trial 1
delta.2 	<- 0					# treatment effect in trial 2
block		<- TRUE				# block randomisation (TRUE/FALSE)

# run simulation:
sim.out <- sim.fun(n.sim, n, block, R, delta.1, delta.2)

# analyze results:
CORS <- matrix(ncol=4,nrow=length(R))
for (j in 1:length(R)){
	M <- sim.out$t[,,j]
	CORS[j,1] <- cor(M[,1], M[,2])
	CORS[j,2] <- cor(M[,3], M[,4])
	CORS[j,3] <- cor(M[,5], M[,6])
	CORS[j,4] <- cor(M[,7], M[,8])
	}

cor.lower <- CORS - qnorm(0.975)*sqrt((1-CORS^2) / (n.sim-2)) 
cor.upper <- CORS + qnorm(0.975)*sqrt((1-CORS^2) / (n.sim-2)) 

df <- data.frame(option=rep(c("1","2","3","4"),each=19), R = rep(R,4), CORS = as.vector(CORS), 
		cor.lower = as.vector(cor.lower), cor.upper = as.vector(cor.upper))

ggplot(data=df, aes(x=R,y=CORS, group=option)) +
	geom_line(aes(linetype=option)) +
	geom_ribbon(aes(ymin=cor.lower, ymax=cor.upper), alpha=0.15) + 
	scale_linetype_manual(values=c("solid","longdash","twodash", "dotted")) +
	labs(x='proportion of cohort included in at least one of two trials',y='correlation between trial results') +
	scale_y_continuous(limits=c(-.6, .6)) +
	annotate("text", x = rep(0.9,4), y = c(.05,.3,-.3,-.05), label = c("1","2","3","4")) +
	theme_bw() 

# ggsave("Figure1.tiff", device="tiff",dpi=200)

# ==============================================================
# ==============================================================
# rerun n.sim = 1e6
# for R = 0.5 
# no treatment effect in trial A (delta.1 = 0)
# no treatment effect in trial B (delta.2 = 0)
# to assess conditional type I error:

TABLE.2 		<- matrix(ncol=5,nrow=12)
rownames(TABLE.2) <- c("option 1. unconditional", "p<0.05", "p>=0.05",
				"option 2. unconditional", "p<0.05", "p>=0.05",
				"option 3. unconditional", "p<0.05", "p>=0.05",
				"option 4. unconditional", "p<0.05", "p>=0.05")
colnames(TABLE.2) <- c("type I error","bias","ESE","bias/ESE","type II error")

set.seed(2505)
n.sim 	<- 1e6
n 		<- 50
R		<- 0.5
delta.1 	<- 0
delta.2 	<- 0
block		<- TRUE

# run simulation:
sim.out <- sim.fun(n.sim, n, block, R, delta.1, delta.2)

# analyze results:
t.values 	<- sim.out$t[,,1]
df 		<- 2*n - 2
p.values 	<- 1-pt(t.values,df)
sign.p 	<- 1*(p.values < 0.05)
m.values	<- sim.out$m[,,1]

TABLE.2[1,1] <- mean( sign.p[,2])
TABLE.2[2,1] <- mean( sign.p[,2][sign.p[,1]==1])
TABLE.2[3,1] <- mean( sign.p[,2][sign.p[,1]==0])
TABLE.2[4,1] <- mean( sign.p[,4])
TABLE.2[5,1] <- mean( sign.p[,4][sign.p[,3]==1])
TABLE.2[6,1] <- mean( sign.p[,4][sign.p[,3]==0])
TABLE.2[7,1] <- mean( sign.p[,6])
TABLE.2[8,1] <- mean( sign.p[,6][sign.p[,5]==1])
TABLE.2[9,1] <- mean( sign.p[,6][sign.p[,5]==0])
TABLE.2[10,1] <- mean( sign.p[,8])
TABLE.2[11,1] <- mean( sign.p[,8][sign.p[,7]==1])
TABLE.2[12,1] <- mean( sign.p[,8][sign.p[,7]==0])

TABLE.2[1,2:3] <- smrz(m.values[,2])
TABLE.2[2,2:3] <- smrz(m.values[,2][sign.p[,1]==1])
TABLE.2[3,2:3] <- smrz(m.values[,2][sign.p[,1]==0])
TABLE.2[4,2:3] <- smrz(m.values[,4])
TABLE.2[5,2:3] <- smrz(m.values[,4][sign.p[,3]==1])
TABLE.2[6,2:3] <- smrz(m.values[,4][sign.p[,3]==0])
TABLE.2[7,2:3] <- smrz(m.values[,6])
TABLE.2[8,2:3] <- smrz(m.values[,6][sign.p[,5]==1])
TABLE.2[9,2:3] <- smrz(m.values[,6][sign.p[,5]==0])
TABLE.2[10,2:3] <- smrz(m.values[,8])
TABLE.2[11,2:3] <- smrz(m.values[,8][sign.p[,7]==1])
TABLE.2[12,2:3] <- smrz(m.values[,8][sign.p[,7]==0])

TABLE.2[,4] <- TABLE.2[,2]/TABLE.2[,3]

# ====================================================================
# === SCENARIO 2 ===
# rerun n.sim = 1e6
# for R = 0.5 
# no treatment effect in trial A (delta.1 = 0)
# treatment effect in trial B (delta.2 = 0.56)
# to assess conditional type II error:

set.seed(2505)
n.sim 	<- 1e6
n 		<- 50
R		<- 0.5
delta.1 	<- 0
delta.2 	<- sqrt(2*(qnorm(0.95)+qnorm(.8))^2/n)
block		<- TRUE

# run simulation:
sim.out <- sim.fun(n.sim, n, block, R, delta.1, delta.2)

# analyze results:
t.values 	<- sim.out$t[,,1]
df 		<- 2*n - 2
p.values 	<- 1-pt(t.values,df)
sign.p 	<- 1*(p.values < 0.05)

TABLE.2[1,5] <- mean(1- sign.p[,2])
TABLE.2[2,5] <- mean(1- sign.p[,2][sign.p[,1]==1])
TABLE.2[3,5] <- mean(1- sign.p[,2][sign.p[,1]==0])
TABLE.2[4,5] <- mean(1- sign.p[,4])
TABLE.2[5,5] <- mean(1- sign.p[,4][sign.p[,3]==1])
TABLE.2[6,5] <- mean(1- sign.p[,4][sign.p[,3]==0])
TABLE.2[7,5] <- mean(1- sign.p[,6])
TABLE.2[8,5] <- mean(1- sign.p[,6][sign.p[,5]==1])
TABLE.2[9,5] <- mean(1- sign.p[,6][sign.p[,5]==0])
TABLE.2[10,5] <- mean(1- sign.p[,8])
TABLE.2[11,5] <- mean(1- sign.p[,8][sign.p[,7]==1])
TABLE.2[12,5] <- mean(1- sign.p[,8][sign.p[,7]==0])

TABLE.2

rnd(TABLE.2,c(3,4,4,2,3))

# ========================================================================
# END
# ========================================================================
