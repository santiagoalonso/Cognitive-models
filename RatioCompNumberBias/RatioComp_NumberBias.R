#This cognitive model compares two probabilities by weighting the ratios, 
#number of winners, and number of losers on each side. It uses a softmax to decide.

load("/Users/sadiaz/Desktop/Mind/MyCodes/RatioCompNumBias.RData")
dataStan <- list(ntrials, weber, r_s, n_r,
                 unique_winners, unique_losers,
                 unique_numbers, unique_probs,
                 winners, losers, probs,
                 w1, w2,
                 l1, l2,
                 p1, p2)

nchains = 1
iter = 10000 
warmup = iter/2
fileStan = 'RatioComp_NumberBias.stan'
fitstan <- stan(file = fileStan, 
                data = dataStan, 
                iter = iter, 
                chains = nchains, 
                warmup = warmup,
                control = list(adapt_delta = 0.8,
                               stepsize = 5, #default is 1; default max_treedepth = 10
                               max_treedepth = 12))