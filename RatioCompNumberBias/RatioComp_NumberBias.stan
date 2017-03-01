data { 
  int<lower = 0>		ntrials;			//Number of trials (i.e. type of trials)
  real<lower = 0> 		weberD;				//Weber fraction (from data)
  
  int<lower = 0>   		r_s[ntrials]; 			//number of correct responses
  int<lower = 0>   		n_r[ntrials]; 			//total number of responses

  int<lower = 0>  		unique_winners[3];		//amount of unique winners [prob 1, prob 2, across probs] 
  int<lower = 0>  		unique_losers[3];		//amount of unique losers [prob 1, prob 2, across probs] 
  int<lower = 0>  		unique_numbers;			//amount of unique numbers (across winners and losers) 
  int<lower = 0>  		unique_probs;			//amount of unique probabilities

 
  int <lower = 0>		winners[unique_winners[3]]; 	//unique winners (across probs.)
  int <lower = 0>		losers[unique_losers[3]]; 	//unique losers (across probs.)
  matrix[unique_probs,4]	probs; 				//col:1 unique probs.; 2: Numerator; 3: denominator; 4: Info.
	

  int <lower = 0>		w1[ntrials]; 			//winner presented in trial i for prob. 1 (left) 
  int <lower = 0>		w2[ntrials]; 			//for prob 2. (right)
  int <lower = 0>		l1[ntrials]; 			//loser presented in trial i for prob. 1 (left)
  int <lower = 0>		l2[ntrials]; 			//for prob 2. (right)
  int <lower = 0>		p1[ntrials]; 			//prob1 presented in trial i(left) 
  int <lower = 0>		p2[ntrials]; 			//prob 2 (right)    		
}

parameters {
  						
  
  //Weights 
  real<lower = -5, upper = 5> wgt_w;  			//Weight on winners
  real<lower = -5, upper = 5> wgt_l;  			//Weight on losers	     
  real<lower = -5, upper = 5> wgt_p;  			//Weight on prob.

  
  //Probs.
  real<lower = 0.001, upper = 0.999> p[unique_probs];	//probabilities 

  
  //Nums.	
  real<lower = 0.001, upper = 100> numW[unique_winners[3]];	//winners (perception)
  real<lower = 0.001, upper = 100> numL[unique_losers[3]];	//losers (perception) 

}
model {  

  real bp[ntrials];	//binomial probability (correct)		
  real sm[2];		//softmax  

  
 for (i in 1:unique_winners[3]) {
	numW[i] ~ normal(winners[i], weberD * winners[i])T[0,];
 }
 
 for (i in 1:unique_losers[3]) {
	numL[i] ~ normal(losers[i], weberD * losers[i])T[0,];
 }
  
 for (i in 1:unique_probs) {
        
	p[i] ~ beta(probs[i,2] + 1, probs[i,3] - probs[i,2]  + 1);

 }


  for (i in 1:ntrials){
	
				
	//Decision system 
	sm[1] = wgt_w * numW[w1[i]] + wgt_l * numL[l1[i]] + wgt_p * p[p1[i]];
	sm[2] = wgt_w * numW[w2[i]] + wgt_l * numL[l2[i]] + wgt_p * p[p2[i]];


	bp[i] = exp(sm[2])/(exp(sm[2]) + exp(sm[1])); //softmax
	if (probs[p1[i],1] > probs[p2[i],1]) {
		bp[i] = 1 - bp[i];
	} 
	bp[i] = fmax(10^(-10),bp[i] - 10^(-10)); //the decision rule may be 0s or 1s and produce a log(0) in the log likelihood.
	
  }     
  r_s ~ binomial( n_r, bp ); 	//likelihood correct
}
generated quantities {
  real log_lik[ntrials];
  real bp[ntrials];	//binomial probability 
  real sm[2];
    

  for (i in 1:ntrials){      

	//Decision system
	sm[1] = wgt_w * numW[w1[i]] + wgt_l * numL[l1[i]] + wgt_p * p[p1[i]];
	sm[2] = wgt_w * numW[w2[i]] + wgt_l * numL[l2[i]] + wgt_p * p[p2[i]];
	

	bp[i] = exp(sm[2])/(exp(sm[2]) + exp(sm[1])); //softmax for right side
	if (probs[p1[i],1] >= probs[p2[i],1]) {
		bp[i] = 1 - bp[i]; 
	}
	bp[i] = fmax(10^(-10),bp[i] - 10^(-10)); //the decision rule may be 0s or 1s and produce a log(0) in the log likelihood.
	
	log_lik[i]  = binomial_lpmf(r_s[i] | n_r[i] , bp[i]);
  }
  
}
