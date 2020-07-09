



#  Below this line: Section 1, hypothetical crime scene sample, Cars. 

#  //////////////////////////////////////////////////////////////////

#  In the following code, a hypothetical example is used:
#  There are two different genes, VW and Ford.
#  In the entire human population, the VW gene has six alleles: red, orange, yellow, green, blue, purple.
#  The Ford gene has only one allele in the population: black. Thus, all individuals in the population are always homozygous black at the Ford gene.
#  For a single crime scene sample, Cars, 33% of the crime scene DNA is from Individual 1 and 66% DNA is from Individual 2.
#  Individual 1 is Homozygous yellow at the VW gene, and Individual 2 is heterozygous (yellow, green) at the VW gene. 



#  Below this line: Section 2: gammaDist functions and arguments as defined by R 

#   /////////////////////////////////////////////////////
    
#  shape = kappa_s  #wikipedia entry on gamma distribution 
#  rate = theta_r   #wikipedia entry on gamma distribution
    
#  dgamma(x, shape, rate = 1, scale = 1/rate, log = FALSE)
#  pgamma(q, shape, rate = 1, scale = 1/rate, lower.tail = TRUE,
#         log.p = FALSE)
#  qgamma(p, shape, rate = 1, scale = 1/rate, lower.tail = TRUE,
#         log.p = FALSE)
#  rgamma(n, shape, rate = 1, scale = 1/rate)

#  Arguments

#  x, q:	  vector of quantiles.
#  p:       vector of probabilities.
#  n:     	number of observations. If length(n) > 1, the length is taken to be the number required.
#  rate:    an alternative way to specify the scale. Rate = 1/scale
#  shape, scale:	shape and scale are parameters. Must be positive, scale strictly.
#  log, log.p: 	logical; if TRUE, probabilities/densities p are returned as log(p).
#  lower.tail: logical; if TRUE (default), probabilities are P[X ≤ x], otherwise, P[X > x].
 
 
 
 #      Below this line: Section 3: definitions of parameters for application of gamma distributions to probabilistic genotyping.
 
#      /////////////////////////////////////////////////////
   
   # The following definitions connect model building in Cowell et al 2015 to BulletProof software input/output.
   
   #  Q = question = forensic sample, such as Q3 in BulletProof, or Cars in Section 1, above.
   #  i = a specific single human individual, equal to a specific known (K) or unknown (U) trace in Cowell 2015, or Individual 1 or 2 in Cars. 
   #  I = Total number of individuals, {1, ..., I}, potentially contributing to a single specific crime scene sample. I = 2 in Cars sample. 
    #  phi_i = specific fraction (phi) of total DNA from a specific individual (i) at a specific single crime scene sample, before PCR.  The sum of all phi_i = sum{phi_1, phi_2, ..., phi_I}  = 1. Phi_1 = 0.33 in Cars.
   
   #  m = a single specific marker, equal to system number such as D8S1179 in Q3 on BulletProof. m = VW or Ford in Cars sample. 
   #  M = total number of markers {1, ..., M} analyzed for a single specific sample. M = 15 for Q3 in BulletProof. M = 2 (Ford and VW) in Cars.
   #  a = specific single allele for a specific given marker. Allele 1 of D8S1179 would be m1a1. In Cars at the VW marker, a1 is yellow and a2 = green. 
    #  Each marker has a unique STR sequence at a specific chromosome locus, with different alleles for a marker having different number of repeats of that sequence.
    #  For example, Allele 1 (a1) for D8S1179 (m1) has 9 repeats. A single set of primers for a marker amplifies all the alleles for that given marker. 
    #  A = total types of alleles for a given marker {1, ..., A}. For marker B8S1179 in BulletProof, A = 4 = {Alelle 1, Allele 2, Allele 3, Allele 4}. In the Cars sample, A = 1 at the Ford marker and A = 2 at the VW marker. 
  
   #  n_a_i = number (n) of a specific single allele (a) from a specific single individual (i). For Individual 2 in Cars, who is heterozygous at VW locus, n_a_i would be = 1 for VW yellow and n_a_i = 2 for Ford black.
    #  N = copy number of all alleles {1, ..., A} of a single marker (m) for a crime scene for all individuals {1, ..., I}. This is "bold n" in Cowell 2015.
    # For sample Cars (section 1, above) N = 4 for VW, and N = 4 for Ford (N = "bold n" in Cowell 2015, pg 10)
    # B_a(N, phi) = effective number of alleles for a specific marker (m) in a single crime scene sample, weighted per DNA fraction, phi_i 
    # B_a = sum{i*phi_i*n_a_i} for each i in {1, ..., I} and each a in {1, ..., A}
    # For Cars, B_a = sum{0.33*2(VW yellow) + 0.66*1(VW yellow) + 0.66*1(VW green)} = 3*0.66 = approx 1.98?
    # Since B_a is defined as covering all alleles A for a given marker m, shouldn't it be B_A or B_m, and not B_a? 
    # B_a will always be equal to the ploidy number (i.e. 2 for diploid humans), no matter the number of individuals contributing to a given sample. 
   
   #  h_a_i = peak height (h) for a single specific allele (a) contributed by a single specific individual (i). 
    # h_a_i is a gamma-distributed random variable (pg 5 in Cowell 2015), because PCR is a random process that amplifies a different amt of DNA in each PCR reaction, even with identical templates/primers, etc.
    # h_a_i <- dgamma(rho*phi_i*n_a_i, eta),  where dgamma(shape,scale) is formula for a gamma pdf in R. See next for additional definitions of these parameters.

 
  # The following parameters are used in the teaching similation code, section 4 below. 
 
   #  rho = proportional to total amount of DNA of a single crime scene sample, before PCR amplification, usually measured in nanograms.
   #  eta = scale of the gamma distribution. Scale is the literal scale of the gamma distribution, 
    #  meaning the range of the x-axis values (peak height RFUs). Scale = 1/rate. Rate is Theta, the second parameter for the R function (see Section 1).
   
   #  h_m_a = a random variable equal to the sum of all peak heights h_a_i from all individuals I {1, ..., I}, across a single allele (a), for a single marker (m). h_m_a does not appear in Cowell 2015.
    #  h_m_a = height of peaks for Yellow VW from indiv 1 and indiv 2 in Cars sample. h_m_a = 345 RFU for m1h1 (marker D8S1179, Allele 1) in BulletProof.
    #  Cowell assumes that peak heights do not vary by the length in bp of an allele. So, assume that the VW gene with a length of 900 bp (yellow) or 1,200 bp (green) will give same peak height.
     #  Similarly, Allele 1 for marker D8S1179 has 9 repeats whereas Allele 4 of marker D8S1179 has 14 repeats. Cowell assumes Allele 1 and Allele 4 would have the same peak heights, assuming the same amount of template DNA for Allele 1 and Allele 4 in the crime scene sample..
######  Therefore, peak heights per number n of a specific allele are the same across alleles for a single marker in a sample, h_m1_a1 = h_m1_a2 = h_m1_a3 if n_a1_i1 = n_a2_i1 = n_a3_i1?
######  Therefore, the distribution of peak heights, H_a, for one single specific allele, a, is the same across all alleles for a given marker?
######  H_a = dgamma(rho*B_a(phi,N), eta)???? Why does this not appear anywhere in Cowell 2015? Why is H_a expressed as a set, instead of the sum of a set (eq 3, pg 10)? 
######  H_a = {h_m1_a1, h_m1_a2} = {h_VW_yellow, h_VW_green}?
######  H_a1 = H_a2 = H_a3 for a given marker, {h_VW_yellow, h_VW, green} = {h_Ford_black}?


    #  Cowell 2015 assumes that the fraction of DNA in a single crime scene sample is constant across markers; i.e. phi_i1_m1 = phi_i1_m2 etc.
     #  This assumption would be violated if DNA is degraded, for example such that the VW alleles are missing from an individual's contribution to the Cars sample crime scene DNA.
######  Therefore, the heights for all alleles A {1, ..., A} for each marker are the same in a sample, h_m1_A = h_m2_A = h_m3_A, also assuming the same amount of starting DNA for each allele in the sample?
######  Best guess, H_a = row sum of peak heights for a single specific marker, since Cowell 2015 states in words that H_a is a sum, although the formula provided (eq 3, pg 11) uses set notation but does NOT include a summation symbol?
     #  The H_a for each marker is unique, but they are all converging on a similar value.
     #  For Q3 in BulletProof, the row sums of peak heights are all pretty similar for the fifteen markers:
     #  Row sums of peak heights in RFU for Q3: {5907, 3583, 3766, 4501, ..., 3653, 4162, 3238
     
     #  The expectation of any random variable is the average value of that variable, after infinite individual measurements. 
     #  mu_H_a = expectation of H_a = (kappa_s*theta_r) = rho*eta, when B_a = 2.  
     #  Because of model assuming H_a_m1 = H_a_m2 etc, therefore mu1 = mu2 etc. 
    
   #  Mu = mean(mu1,mu2, ..., mu_M)/2. Mu is an estimate of the average peak height for the sample, divided by 2 to account for diploidy. 
   #  Mu = Mu in BulletProof. Mu in BulletProof is an estimator of the true mu for the real gamma distribution of a crime scene sample. 
   #  s = sample standard deviation = eta*sqrt(rho). See NIST website. https://www.itl.nist.gov/div898/handbook/eda/section3/eda366b.htm
   #  variance = s^2 = (eta^2)*rho. Agrees with Cowell 2015.
   #  Sigma = coefficient of variation = s/Mu = 1/sqrt(rho)

#####  There seems to be different correction factors to estimate s depending on which distribution. (Wikipedia article, "standard deviation")
#####  Therefore, compare to see if s = sigma, where sigma = stdev(H_a), and s = eta*sqrt(rho). 

  #  The probability P of allelic dropout given a threshold C is given using the cdf of a gamma distribution:
  #  P(Z_a = 0|mu) = dgamma(C, shape=mu/eta, scale=eta) where Z_a = H_a if H_a >= C. (Cowell 2015 pg 12 eq 6)

  
#       Below this line: Section 4: show different rates of allelic drop-out between homozyg and heterozyg individual trace

#    //////////////////////////////////////////////////////////////
  
  #       Objective: plot dgamma with horizontal line indicating quartile for cdf mass calc, pgamma.

  #   

      
rho <- 2   # Chose a random integer to represent the amount of DNA in sample. 
eta <- 28   # For Cowell 2015 case study, typical value for eta = 28
phi <- 1   # for our simulation, only one contributer per trace so phi_i = 1 for all i.
 
n_a_i1 <- 0.25  # number of alleles of type a for heterozygous individual i1; should be 1
n_a_i2 <- 0.5  # number of alleles of type a for a homozygous individual i2; should be 2
B_a <-  2    # effective number of alleles for a given marker for diploid organisms

    
H_a <- rgamma(100, shape=rho*phi*B_a, scale=eta)  # generate a vector H_a of 100 peak heights.
# Because of the model simplifications across different alleles and markers, 
# H_a could be 100 obs of one marker, or one obs of 100 markers, or 50 obs each of 2 markers, all from a single individual contributer to a specific crime scene.

Heterozyg <- dgamma(H_a, shape=rho*phi*n_a_i1, scale = eta)  # Let Heterozyg be a trace of a single individual that is heterozygous for all alleles.

Homozyg <- dgamma(H_a, shape=rho*phi*n_a_i2, scale=eta)  # Let Homoz be a trace of a single individual that is homozygous for all alleles

C <- (median(H_a))/2.7  # I picked a visually convenient spot to put a vertical line based on random sample H_a. Cowell 2015 uses C = 50 for typical threshold cut-off.

xtotal <- range(H_a)     #figure out range of x axis for plot
ytotal <- range(Heterozyg)     #figure out range of y axis for plot

plot(Heterozyg~H_a, pch=20, main="Gamma PDFs", xlim=xtotal, ylim=ytotal, xlab="Peak Heights, RFU", ylab="Probability")     # plot "Heterozyg" pdf
points(Homozyg~H_a, pch=0)   # overlay points of "Homozyg" pdf onto the plot

abline(v=C, col="red", lwd=3, lty=2) #enter a dashed red line at cutoff threshold, C

cdf_het <- round(dgamma(C, shape=rho*phi*n_a_i1, scale=eta),4)  #calculates the mass for Heterozyg at cutoff C of H_a
cdf_hom <- round(dgamma(C, shape=rho*phi*n_a_i2, scale=eta),4)  #calculates the mass for Homozyg at cutoff C of H_a

spot <- (max(xtotal)-100)  #put upper rght-hand corner of legend in the right spot for chnging x axis scales between simulations

legend(spot, 0.008, legend=c(cdf_het, cdf_hom), pch=c(20,0), title="Drop-Out Prob.")  #plot masses and mass IDs as fig legend

Mu <-  rho*eta/2
Sigma <-  1/sqrt(rho)

print(Mu)
print(Sigma)


