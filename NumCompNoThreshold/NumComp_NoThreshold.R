
#### Functions #### 

ANS_sample<-function(n1,n2,w){
  #Inputs: 
  #   n1 and n2 are the number of dots to be compared
  #   w is the weber fraction e.g. estimated (mle) with subjects data
  #Output:
  #   Sample from the approximate number system (ANS)

  rnorm(1, (n2 - n1), (w * sqrt(n1^2 + n2^2))) 
}

jerkTrajectory <-function(pos,vel,acc,target,tf,t0){
  # This function is an optimal feedback trajectory minimizing jerk (Hoff & Arbib, 1993)
  # Inputs:
  #   pos, vel, & acc: current position, velocity, and accelaration (scalar)
  #   target: end position (scalar)
  #   tf and t0: final and current time, respectively (in ms) (scalar)
  # Output:
  #   update to the current motor state:
  #       row 1 = change in pos (i.e speed) 
  #       row 2 = change in vel (i.e acceleration) 
  #       row 3 = change in acc (i.e. jerk)   
  
  D = tf - t0
  update = matrix(c(0, 1, 0, 0, 0, 1, -60/(D^3), -36/(D^2), -9/D),
                nrow=3,ncol=3,byrow=T) %*% matrix(c(pos, vel, acc)) + 
           matrix(c(0, 0, 60/(D^3)))*target
  update 
}

my.erfc<-function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)

weberFraction<-function(x,w){
  #Inputs:
  #   x: number ratio (scalar)
  #   w: weber ratio (scalar)
  
  1 - 1/2 * my.erfc(abs(x - 1) / (sqrt(2) * w * sqrt(x^2 + 1))) 
}

my.entropy <- function(p) {
  # p: probability (scalar)
  
  H = ifelse(p != 0 & p != 1,-(p*log(p)+(1-p)*log(1-p)),
             ifelse(p == 0,-(1-p)*log(1-p),-p*log(p)))
  H
}

modelPos<-function(w, W, R, D, r){
  ### This function computes horz. trajectories per num. ratio 
  ##  based on evidence coming from the ANS (approx. number system)

  #Inputs:
  #   w = Weber fraction (scalar) 
  #   W = known likelihood noise (scalar)
  #   R = sampling rate (scalar)
  #   D = EMA (scalar)
  #   r = num. ratios (vector)
  #Output:
  # see list at the end

  positions = list()
  plan = list()
  targetcorrect = 113 #(mm)
  RT = rep(0, length(r)) #Response time
  FT = rep(0, length(r)) #Flight time
  maxn = 25 #max number of dots shown to subjects
  n = rbind(r,rep(1,length(r))) 
  for (i in 1:length(n[1,])){ #This generates random integers for the num. ratios
    ntemp=1.5
    while (sum(ntemp)%%1!=0) {
      ntemp = t(t((1+round(runif(1,0,(maxn-1))))*n[,i])) 
    }
    n[,i] = ntemp
  }
  sideCorrect = round(runif(length(n[1,]),1,2)) #1 left, 2 right
  R = max(1,round(R))
  LFrt = rep(0.28,length(r)) #Detectable lift - off point (mm) (as defined by mean displacement at button release) 
  for (i in 1:length(r)){
    PlanHorz = matrix(rep(0,3),nrow=3,ncol=1)
    tc = 1 #Current time
    tf = 0 #Current flight time
    RTRT = 0
    flighttime = sample(flightTimes[flightTimes[,'Ratio'] == r[i], 'FlightTime'], 1) #from data
    probCorrect = 0.5
    probMistake = 0.5
    certainty = 0
    n1 = n[1,i]
    n2 = n[2,i]
    checkLO = 0 #checks lift-off
    while (tf < flighttime){

      #Evidence & Inference
      if (tc %% R == 1 | R == 1) {
        
        #Evidence
        evidence = ANS_sample(n1, n2, w) #if sample is less than zero, then the model is having the wrong perception 
        if(tc == 1){
          Samples = evidence 
        }else {
          EMA = D * evidence + (1 - D) * Samples[length(Samples)]
          Samples = c(Samples, EMA)
        }
        
        #Inference 
        meanS = mean(Samples)
        nBayes = length(Samples)      
        probMistake = pnorm(0,meanS,(W)/nBayes)  
        probCorrect = 1 - probMistake
        certainty = (my.entropy(1/2) - my.entropy(probCorrect))/(my.entropy(1/2)) #(0 uncertain, 1 certain). 
      } 
      
      
      #Motor plan
      if (probCorrect >= 0.5){
        weight2 = certainty #Weight correct
        weight1 = 0 #Weight incorrect
      } else {
        weight2 = 0 
        weight1 = certainty  
      }
      update1 = jerkTrajectory( PlanHorz[1, length(PlanHorz[1,])], 
                                PlanHorz[2, length(PlanHorz[2,])],
                                PlanHorz[3, length(PlanHorz[3,])], 
                                -targetcorrect, flighttime, tf) 
      update2 = jerkTrajectory( PlanHorz[1, length(PlanHorz[1,])], 
                                PlanHorz[2, length(PlanHorz[2,])],
                                PlanHorz[3, length(PlanHorz[3,])], 
                                targetcorrect, flighttime, tf) 
      PlanHorz = cbind(PlanHorz,
                       PlanHorz[,length(PlanHorz[1,])] + (weight1*update1 + weight2*update2)) #row 1= pos; row 2= vel; row 3 = accel.
      

      #Increase in flight time
      tf = tf + 1 
      
      #Increase in time 
      tc = tc + 1
      
      #Check: Lift-off 
      if (!checkLO && abs(PlanHorz[1,length(PlanHorz[1,])]) > LFrt[i]) { 
        RTRT = tc #Response time
        checkLO = 1
      }
    }
    
    plan[[i]] = PlanHorz
    if (RTRT == 0){
      RT[i] = length(PlanHorz[1,])
    } else {
      RT[i] = RTRT
    }
    FT[i] = flighttime
  }
  
  op = list(plan = plan, #Each element of the list is a ratio (easy to hardest) 
                         #with a matrix with position, vel., and acc. (rows) at each time point (columns)
            RT = RT, #Vector with response times for each ratio (easy to hardest)
            FT = FT, #Vector with flight times for each ratio (easy to hardest)
            numbers = n, #number used. Each column is a ratio (easy to hardest)
            sideCorrect = sideCorrect, #Random number simulating where the largest number is placed (1 left, 2 right) 
            LFrt = LFrt #lift-off detection point (see above)
            )
  op
}




# Example using best fit parameters to a data set #### 
# This plots one example trial for each ratio in r
load('FT.RData') # flight times (to sample)
r = c(0.1, 0.25, 0.5, 0.75, 0.9) # ratios
nDT = 400 # non-decision times (in ms)
nPoints = 101 # normalized time points
target = 113 # horz. position of target (mm)
R = 32 # rate accumulation i.e. a sample every X ms
W = 7.716049 # likelihood noise 
D = 0.06368313 # EMA parameter (0 full memory decay/linear ballistic accumulation; 1 no memory decay/traditional accumulation)
w = 0.1746911 # Weber fraction

infoByRatio = modelPos(w, W, R, D, r) #list (see function)
par(mfrow = c(2,2))

#raw positions (note: negative positions are incorrect directions)
pos = list()
for (i in 1:length(r)) {
  pos_temp = infoByRatio$plan[[i]][1,] 
  if (max(abs(pos_temp) >= (0.99*target))){ #this accounts for early stop in the recording (see paper for details)
    temp = min( which(abs(pos_temp) >= (0.99*target)) )
  } else {
    temp = infoByRatio$FT[i]
  }
  pos[[i]] = pos_temp[infoByRatio$RT[i]:temp]
}
plot(pos[[1]], col = 'cyan4', type = 'l', lwd = 2, 
     xlim = c(0, max(length(pos[[1]]), length(pos[[2]]), length(pos[[3]]),
                     length(pos[[4]]), length(pos[[5]]))),
     ylim = c(-120,120),
     xlab = 'Raw time (ms)', ylab = 'Horz. pos (mm)') 
lines(pos[[2]], col = 'black', lwd = 2)
lines(pos[[3]], col = 'purple', lwd = 2)
lines(pos[[4]], col = 'red', lwd = 2)
lines(pos[[5]], col = 'forestgreen', lwd = 2)

#time-normalized positions (note: negative positions are incorrect directions)
r1_norm = spline(pos[[1]], n = nPoints)$y 
r2_norm = spline(pos[[2]], n = nPoints)$y 
r3_norm = spline(pos[[3]], n = nPoints)$y 
r4_norm = spline(pos[[4]], n = nPoints)$y 
r5_norm = spline(pos[[5]], n = nPoints)$y 
plot(r1_norm, col = 'cyan4', type = 'l', lwd = 2, 
     xlim = c(0, nPoints), 
     ylim = c(-120,120),
     xlab = 'Normalized time', ylab = 'Horz. pos (mm)')
lines(r2_norm, col = 'black', lwd = 2)
lines(r3_norm, col = 'purple', lwd = 2)
lines(r4_norm, col = 'red', lwd = 2)
lines(r5_norm, col = 'forestgreen', lwd = 2)

#RT
plot(r, infoByRatio$RT + nDT, cex = 1.5, 
     ylim = c(0.98*min(infoByRatio$RT + nDT), 1.02*max(infoByRatio$RT + nDT)),
     col = c('cyan4','black','purple','red','forestgreen'),
     xlab = 'Ratio', ylab = 'RT (ms)')

