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

sim.fun <- function(n.sim, n, R, delta.1, delta.2){
	RESULTS.t <- array(NA, dim=c(n.sim,8,length(R)) )
	RESULTS.m <- array(NA, dim=c(n.sim,8,length(R)) )

	for (j in 1:length(R)){
		for (i in 1:n.sim){
		
		# option 1.	subjects can only participate in 1 trial
			N <- round(4*n/R[j])	
			ids <- 1:N 
			dat <- rnorm(N)
				sample.10 <- sample(ids,n)
				used.id <- sample.10
				sample.11 <- sample(ids[-used.id],n)
				used.id <- c(used.id,sample.11)
				sample.20 <- sample(ids[-used.id],n)
				used.id <- c(used.id,sample.20)
				sample.21 <- sample(ids[-used.id],n)
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
				sample.11 	<- sample(ids,n)
				sample.21 	<- sample(ids[-sample.11],n)
				used.id 	<- c(sample.11,sample.21)
				sample.10 	<- sample(ids[-used.id],n)
				sample.20 	<- sample(ids[-used.id],n)
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
				sample.11 	<- sample(ids,n)
				sample.21 	<- sample(ids[-sample.11],n)
				sample.10 	<- sample(ids[-sample.11],n)
				sample.20 	<- sample(ids[-sample.21],n)
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
				sample.10 <- sample(ids,n)
				sample.11 <- sample(ids[-sample.10],n)
				sample.20 <- sample(ids,n)
				sample.21 <- sample(ids[-sample.20],n)
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
n.sim 	<- 1e4
n 		<- 50
R		<- seq(.05,.95, by=0.05) 
delta.1 	<- 0
delta.2 	<- 0

# run simulation:
sim.out <- sim.fun(n.sim, n, R, delta.1, delta.2)

# analyze results:
CORS <- matrix(ncol=4,nrow=length(R))
for (j in 1:length(R)){
	M <- sim.out$t[,,j]
	CORS[j,1] <- cor(M[,1], M[,2])
	CORS[j,2] <- cor(M[,3], M[,4])
	CORS[j,3] <- cor(M[,5], M[,6])
	CORS[j,4] <- cor(M[,7], M[,8])
	}

df <- data.frame(option=rep(c("1","2","3","4"),each=19), R = rep(R,4), CORS = as.vector(CORS))
ggplot(data=df, aes(x=R,y=CORS, group=option)) +
	geom_line(aes(linetype=option)) +
	scale_linetype_manual(values=c("solid","longdash","twodash", "dotted")) +
	labs(x='proportion of cohort included in trials',y='correlation between trial results') +
	scale_y_continuous(limits=c(-.55, .55))

# ggsave("Figure1.tiff", device="tiff",dpi=200)

# ==============================================================
# rerun n.sim = 1e6
# for R = 0.5 
# no treatment effect in trial A (delta.1 = 0)
# no treatment effect in trial B (delta.2 = 0)
# to assess conditional type I error:

set.seed(2505)
n.sim 	<- 1e6
n 		<- 50
R		<- 0.5
delta.1 	<- 0
delta.2 	<- 0

# run simulation:
sim.out <- sim.fun(n.sim, n, R, delta.1, delta.2)

# analyze results:
t.values 	<- sim.out$t[,,1]
df 		<- 2*n - 2
p.values 	<- 1-pt(t.values,df)
sign.p 	<- 1*(p.values < 0.05)
m.values	<- sim.out$m[,,1]

TABLE.2 	<- matrix(ncol=4,nrow=8)
TABLE.2[1,1] <- mean( sign.p[,2][sign.p[,1]==1])
TABLE.2[2,1] <- mean( sign.p[,2][sign.p[,1]==0])
TABLE.2[3,1] <- mean( sign.p[,4][sign.p[,3]==1])
TABLE.2[4,1] <- mean( sign.p[,4][sign.p[,3]==0])
TABLE.2[5,1] <- mean( sign.p[,6][sign.p[,5]==1])
TABLE.2[6,1] <- mean( sign.p[,6][sign.p[,5]==0])
TABLE.2[7,1] <- mean( sign.p[,8][sign.p[,7]==1])
TABLE.2[8,1] <- mean( sign.p[,8][sign.p[,7]==0])

TABLE.2[1,2:3] <- smrz(m.values[,2][sign.p[,1]==1])
TABLE.2[2,2:3] <- smrz(m.values[,2][sign.p[,1]==0])
TABLE.2[3,2:3] <- smrz(m.values[,4][sign.p[,3]==1])
TABLE.2[4,2:3] <- smrz(m.values[,4][sign.p[,3]==0])
TABLE.2[5,2:3] <- smrz(m.values[,6][sign.p[,5]==1])
TABLE.2[6,2:3] <- smrz(m.values[,6][sign.p[,5]==0])
TABLE.2[7,2:3] <- smrz(m.values[,8][sign.p[,7]==1])
TABLE.2[8,2:3] <- smrz(m.values[,8][sign.p[,7]==0])

# ====================================================================
# === SCENARIO 2 ===
# no treatment effect in trial A (delta.1 = 0)
# treatment effect in trial B (delta.2 = 0.56)
# for R = 0.5: 
# to assess conditional type II error:

set.seed(2505)
n.sim 	<- 1e6
n 		<- 50
R		<- 0.5
delta.1 	<- 0
delta.2 	<- sqrt(2*(qnorm(0.95)+qnorm(.8))^2/n)

# run simulation:
sim.out <- sim.fun(n.sim, n, R, delta.1, delta.2)

# analyze results:
t.values 	<- sim.out$t[,,1]
df 		<- 2*n - 2
p.values 	<- 1-pt(t.values,df)
sign.p 	<- 1*(p.values < 0.05)

TABLE.2[1,4] <- mean(1- sign.p[,2][sign.p[,1]==1])
TABLE.2[2,4] <- mean(1- sign.p[,2][sign.p[,1]==0])
TABLE.2[3,4] <- mean(1- sign.p[,4][sign.p[,3]==1])
TABLE.2[4,4] <- mean(1- sign.p[,4][sign.p[,3]==0])
TABLE.2[5,4] <- mean(1- sign.p[,6][sign.p[,5]==1])
TABLE.2[6,4] <- mean(1- sign.p[,6][sign.p[,5]==0])
TABLE.2[7,4] <- mean(1- sign.p[,8][sign.p[,7]==1])
TABLE.2[8,4] <- mean(1- sign.p[,8][sign.p[,7]==0])

# =========================================================================
rownames(TABLE.2) <- c("option 1. p<0.05", "p>=0.05",
				"option 2. p<0.05", "p>=0.05",
				"option 3. p<0.05", "p>=0.05",
				"option 4. p<0.05", "p>=0.05")
colnames(TABLE.2) <- c("type I error","bias","ESE","type II error")

rnd(TABLE.2,c(3,4,4,3))

# ========================================================================
# END
# ========================================================================
