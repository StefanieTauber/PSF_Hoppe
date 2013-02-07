
Hoppe_Urn = function(sigma, theta, n){

### two dimensional hoppe urn, sigma and theta are parameters from the Pitman Sampling Formula
### n is the number of ball that shall be drawn)
	
### urn rules are as follows: initially there is one black ball (the "mutator") in the urn with weight theta. Then you always sample with replacement. If a black ball is drawn, it is put back, and 2 additional balls are added. One black ball with weight sigma and one ball with a new color with weight (1-sigma). If a colored ball is drawn it is put back and another ball of the same color with unit weight is added.

	urn_updated = c(1,1,2)
	prob_updated = c(theta, sigma, 1-sigma)
	
	urn_orig = urn_updated
	
	prob_updated = as.vector(tapply(prob_updated,urn_updated,sum))
	urn_updated = unique(urn_updated)
	
	
	### until n ball (none of them black) are drawn
	while (length(urn_orig[urn_orig!=1]) < n)
	{
		
	ball = sample(urn_updated,1,replace = T, prob = prob_updated)
	
	### black ball is drawn
	if (ball == 1){
		
		urn_orig = c(urn_orig, 1, max(urn_updated) + 1)
		
		urn_updated = c(urn_updated, 1, max(urn_updated) + 1)
		prob_updated = c(prob_updated, sigma, 1-sigma)
		
		index = order(urn_updated, decreasing = F)
		prob_updated = prob_updated[index]
		urn_updated = urn_updated[index]
		
		prob_updated = as.vector(tapply(prob_updated,urn_updated,sum))
		urn_updated = unique(urn_updated)
		
		if (prob_updated[1] < 0) prob_updated[1] <- 0
	
	}
	
	### colored ball is drawn
	if (ball !=1 ){
		
		urn_orig = c(urn_orig, ball)
		
		urn_updated = c(urn_updated, ball)
		prob_updated = c(prob_updated, 1)
		
		index = order(urn_updated, decreasing = F)
		prob_updated = prob_updated[index]
		urn_updated = urn_updated[index]
		
		prob_updated = as.vector(tapply(prob_updated,urn_updated,sum))
		urn_updated = unique(urn_updated)
		
		if (prob_updated[1] < 0) prob_updated[1] <- 0
	
		}
	
}
 ### return value: 
 ### orig_urn: all balls which are not black
return(list(orig_urn = urn_orig[urn_orig>1], urn_updated = urn_updated, prob_updated = prob_updated))
	
}

set.seed(123)
# infinite universe, low innovation rate,draw 10 balls
Hoppe_Urn(0.5,0.2,10) # although we draw 10 times, we just observe 5 different balls
# infinite universe, high innovation rate, draw 10 balls
Hoppe_Urn(0.5,100,10) # 10 draws, 10 colored balls

# finite universe, low innovation rate,draw 10 balls, at maximum 4 colored balls possible
Hoppe_Urn(-0.5,2,10)
max(Hoppe_Urn(-0.5,2,10000)$urn_updated) # if we sample really a lot, we see the highest number is 5, which means that there are 4 colored balls (as 1 is the black ball) 
# infinite universe, high innovation rate, draw 10 balls, at maximum 8 colored balls possible
Hoppe_Urn(0.5,4,10)
max(Hoppe_Urn(-0.5,4,10000)$urn_updated) # if we sample really a lot, we see the highest number is 9, which means that there are 8colored balls (as 1 is the black ball) 