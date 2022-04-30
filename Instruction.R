##The version of R needs to be 4.1.1 or later.
##Rtools needs to be downloaded and installed in advance.
##You can download and install the latest Rtools at the official website (https://cran.r-project.org/bin/windows/Rtools/).
##Install the latest version of "rstan" (install.packages("rstan")).

################################################################################
################################################################################
##############################  Example  #######################################
################################################################################
################################################################################

################################################################################
##############################  example of G_Bayes_XCI  ########################
################################################################################

##example 1:
##quantitative trait with covariate
##estimate the gamma of the gene using the GBN method
##Bayesian method with the prior distribution of gamma being a truncated normal distribution specified in our paper
###the prior distributions of other unknown parameters are consistent with those in our paper

###########################  step 1: installing the R package  ######################

##setwd to the root of the folder path of "GEXCIS_0.1.0.zip" and run the following command
install.packages("GEXCIS_0.1.0.zip",repos = NULL)
##or setwd to the root of the folder path of "GEXCIS_0.1.0.tar.gz" and run the following command
install.packages("GEXCIS_0.1.0.tar.gz",type='source')


############################  step 2: loading the R package  ########################

library("GEXCIS")


##############  step 3: viewing the structure of the example data  #############

head(phenotype1)
##  pid iid fid mid sex         x1          y
##1   1   1   0   0   2 -0.1687666  1.0536866
##2   2   1   0   0   2 -0.1976122  1.2608629
##3   3   1   0   0   2 -2.3744448 -1.2320881
##4   4   1   0   0   2  1.5185388  0.4924657
##5   5   1   0   0   2  0.1481116  0.6393415
##6   6   1   0   0   2 -0.7957047 -0.6921572

##"pid", "iid", "fid", "mid" and "sex" in the first five columns are the pedigree structure information as described in the "GEXCIS-manual". "x1" in the sixth column is a covariate and "y" in the last column is the trait value.

head(genotype1)
##   pid iid fid mid sex rs1 rs2 rs3 rs4 rs5 rs6 rs7 rs8 rs9 rs10 rs11 rs12 rs13 rs14 rs15 rs16 rs17 rs18 rs19 rs20 rs21 rs22 rs23 rs24 rs25 rs26 rs27 rs28 rs29
##1   1   1   0   0   2   0   0   2   1   0   0   1   1   1    1    0    0    1    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0
##2   2   1   0   0   2   0   0   0   0   0   0   0   0   1    0    0    1    1    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0
##3   3   1   0   0   2   0   1   1   1   0   0   1   1   1    0    0    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    2    0    0
##4   4   1   0   0   2   0   0   0   0   0   0   1   0   1    1    0    0    1    0    0    0    0    1    0    0    0    0    0    0    0    1    1    0    0
##5   5   1   0   0   2   0   0   0   1   0   0   1   1   1    0    0    2    0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    0    0
##6   6   1   0   0   2   0   0   0   1   0   0   1   0   0    0    0    1    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
##   rs30 rs31 rs32 rs33 rs34 rs35 rs36 rs37 rs38 rs39 rs40 rs41 rs42 rs43 rs44 rs45 rs46 rs47 rs48 rs49 rs50
##1    0    0    1    0    0    1    0    0    1    2    0    2    0    0    0    0    0    2    2    0    0
##2    0    0    1    1    0    0    0    0    1    1    1    0    0    2    2    0    0    1    2    0    0
##3    0    1    1    0    0    0    0    0    1    1    0    1    0    0    1    0    0    1    2    0    0
##4    0    0    1    0    0    1    0    0    1    1    0    0    0    0    0    0    0    1    1    0    0
##5    0    1    1    0    0    1    0    0    0    1    0    0    0    0    0    0    0    1    2    0    0
##6    0    0    1    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    1    0    0

##"pid", "iid", "fid", "mid" and "sex" in the first five columns are the pedigree structure information as described in the "GEXCIS-manual". "rs1-rs50" are genotype codes for 50 SNPs, and each genotype is coded as 0, 1 or 2, indicating the number of the minor alleles.


################  step 4: explaining the corresponding functions  ##############

##Users can use the Bayesian methods to estimate the degree of the skewness of XCI for the genes through the function "G_Bayes_XCI". 
##When the input variable "prior"="normal", the Bayesian method is GBN; when the input variable "prior"="uniform", the Bayesian method is GBU.
##When the input variable "prior"="customize", users could specify the prior distributions of ¦Ã and other unknown parameters according to their own research background.
##Here let "prior" = "normal"


########################  step 5: running the functions  #######################

##If the function "G_Bayes_XCI" is used, the R package "rstan" need to be loaded  first
library(rstan)
##hide run log
rstan_options(javascript=FALSE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

##Special note: The example here is for demonstration only.
##In practical applications, "chains_num", "iter_num" and "warmup_num" should be appropriately increased to ensure the convergence.
##Otherwise, a series of warnings will appear, prompting to increase the number of iterations.
##To ensure the convergence, we recommend using the corresponding values listed in the "GEXCIS-manual".

##The results may be different for different runs, because of the sampling randomness of the HMC algorithm. 
##If the fixed results are wanted, the seed number should be set before running the "G_Bayes_XCI". 
##Note that different version of R may lead to different results under the same seed number. The results of the examples are obtained under the R with version 4.1.1.

##The following errors and warmings can be ignored:

##Error in open.connection(con, open = mode):
##schannel: next InitializeSecurityContext failed: SEC_E_UNTRUSTED_ROOT (0x80090325) -The certificate chain was issued by an authority that is not trusted.

##Warning message:
##  In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
##  '-E' not found
##recompiling to avoid crashing R session

set.seed(123456)
G_Bayes_XCI(phenotype=phenotype1,genotype=genotype1,trait_type="quantitative",phenotype_missing=NA,genotype_missing=NA,prior="normal",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)
##$Point_Estimate
##[1] 0.4687626

##$HPDI_Lower
##[1] 0.1133476

##$HPDI_Upper
##[1] 1.235775

######################  step 6: interpreting the results  ######################

##The point estimate of gamma obtained by the GBN method is 0.4687626, and the 95% HPDI of gamma is (0.1133476, 1.235775). 
##The XCI pattern of the gene on the trait may be XCI-R or XCI-E.
##If the estimate of gamma is significantly different from 1, the gene is statistically inferred to undergo XCI-S, For example, gamma = 0.4 represents that the XCI-S on this gene is overall skewed towards the minor alleles, where only about 20% (0.4/2) of the cells have the minor allele active overall, and the other 80% of the cells have the major allele active overall. 



######################  Examples 2-10 have the similar steps  ##################

##the structures of the remaining sample data are similar to phenotype1 and genotype1, also see the "GEXCIS-manual" for the details
head(phenotype2)
head(phenotype3)
head(phenotype4)
head(phenotype5)
head(phenotype6)
head(phenotype7)
head(genotype2)
head(genotype3)

##example 2:
##quantitative trait with covariate
##the prior distribution of gamma is a uniform distribution specified in our paper
##the prior distributions of other unknown parameters are consistent with those in our paper
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype1,genotype=genotype1,trait_type="quantitative",phenotype_missing=NA,genotype_missing=NA,prior="uniform",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 3:
##quantitative trait with covariate
##the prior distributions of gamma and other unknown parameters are "customize"
##users are required to define the prior distribution of each parameter according to their own research background, for example:
model_customize="
data {
    int<lower=0> N ;
    vector[N] y;
    vector[N] x1;
    vector[N] gg_1;
    vector[N] gg_2;
}
parameters {
    real beta_0;
    real beta_1;
    real beta_c;
    real<lower=0,upper=2> gamma;
    real<lower=0> sigma;
}
model {
    vector[N] theta;
    theta = beta_0 + beta_1*x1 + beta_c*gamma*gg_1 + beta_c*(2-gamma)*gg_2;
    target += normal_lpdf( beta_0 | 0, 100 );
    target += normal_lpdf( beta_1 | 0, 10 );
    target += normal_lpdf( beta_c | 0, 20 );
    target += normal_lpdf( gamma | 1, 2 );
    target += exponential_lpdf(sigma | 1);
    for(i in 1:N)
      target += normal_lpdf( y[i] | theta[i], sigma );
}"
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype1,genotype=genotype1,trait_type="quantitative",phenotype_missing=NA,genotype_missing=NA,
            prior="customize",model_customize=model_customize,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 4:
##quantitative trait without covariate
##the prior distribution of gamma is a truncated normal distribution specified in our paper
##the prior distributions of other unknown parameters are consistent with those in our paper
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype2,genotype=genotype1,trait_type="quantitative",phenotype_missing=NA,genotype_missing=NA,prior="normal",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 5:
##qualitative trait with covariate
##the prior distribution of gamma is a truncated normal distribution specified in our paper
##the prior distributions of other unknown parameters are consistent with those in our paper
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype5,genotype=genotype1,trait_type="qualitative",phenotype_missing=NA,genotype_missing=NA,prior="normal",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 6:
##qualitative trait without covariate
##the prior distribution of gamma is a uniform distribution specified in our paper
##the prior distributions of other unknown parameters are consistent with those in our paper
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype6,genotype=genotype1,trait_type="qualitative",phenotype_missing=NA,genotype_missing=NA,prior="uniform",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)


##example 7:
##qualitative trait without covariate
##the prior distributions of gamma or other unknown parameters are "customize"
##users are required to define the prior distribution of each parameter according to their own research background, for example:
model_customize = "
        data {
          int<lower=0> N;
          int y[N] ;
          vector[N] gg_1;
          vector[N] gg_2;
        }
        parameters {
          real beta_0;
          real beta_c;
          real<lower=0,upper=2> gamma;
        }
        model {
          vector[N] theta;
          theta = beta_0 + beta_c*gamma*gg_1 + beta_c*(2-gamma)*gg_2;
          target += normal_lpdf( beta_0 | 0, 1 );
          target += normal_lpdf( beta_c | 0, 20 );
          target += normal_lpdf( gamma | 0.5, 2 );
          target += bernoulli_logit_lpmf( y | theta );
        }"
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype6,genotype=genotype1,trait_type="qualitative",phenotype_missing=NA,genotype_missing=NA,
            prior="customize",model_customize=model_customize,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 8:
##quantitative trait with covariate and missing values
##both "phenotype" and "genotype" contain the missing values (denoted by NA)
##the prior distribution of gamma is a uniform distribution specified in our paper
##the prior distributions of other unknown parameters are consistent with those in our paper
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype3,genotype=genotype2,trait_type="quantitative",phenotype_missing=NA,genotype_missing=NA,prior="uniform",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 9:
##qualitative trait with covariate and missing values
##both "phenotype" and "genotype" contain the missing values (denoted by NA)
##the prior distribution of gamma is a uniform distribution specified in our paper
##the prior distributions of other unknown parameters are consistent with those in our paper
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype7,genotype=genotype2,trait_type="qualitative",phenotype_missing=NA,genotype_missing=NA,prior="uniform",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 10:
##quantitative trait with covariate and missing values
##both "phenotype" and "genotype" contain the missing values (denoted by 9)
##the prior distribution of gamma is a truncated normal distribution specified in our paper
##the prior distributions of other unknown parameters are consistent with those in our paper
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype4,genotype=genotype3,trait_type="quantitative",phenotype_missing=9,genotype_missing=9,prior="normal",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)



################################################################################
############################  example of G_Frequen_XCI  ########################
################################################################################

##example 11:
##quantitative trait with covariate
##obtain the point estimate and the penalized point estimate of gamma for the genes, as well as the 95% CIs derived by the Fieller's and PF methods 

###########################  step 1: installing the R package  ######################

##setwd to the root of the folder path of "GEXCIS_0.1.0.zip" and run the following command
install.packages("GEXCIS_0.1.0.zip",repos = NULL)
##or setwd to the root of the folder path of "GEXCIS_0.1.0.tar.gz" and run the following command
install.packages("GEXCIS_0.1.0.tar.gz",type='source')


############################  step 2: loading the R package  ########################

library("GEXCIS")


##############  step 3: viewing the structure of the example data  #############

head(phenotype1)
##  pid iid fid mid sex         x1          y
##1   1   1   0   0   2 -0.1687666  1.0536866
##2   2   1   0   0   2 -0.1976122  1.2608629
##3   3   1   0   0   2 -2.3744448 -1.2320881
##4   4   1   0   0   2  1.5185388  0.4924657
##5   5   1   0   0   2  0.1481116  0.6393415
##6   6   1   0   0   2 -0.7957047 -0.6921572

##"pid", "iid", "fid", "mid" and "sex" in the first five columns are the pedigree structure information as described in the "GEXCIS-manual". "x1" in the sixth column is a covariate and "y" in the last column is the trait value.

head(genotype1)
##   pid iid fid mid sex rs1 rs2 rs3 rs4 rs5 rs6 rs7 rs8 rs9 rs10 rs11 rs12 rs13 rs14 rs15 rs16 rs17 rs18 rs19 rs20 rs21 rs22 rs23 rs24 rs25 rs26 rs27 rs28 rs29
##1   1   1   0   0   2   0   0   2   1   0   0   1   1   1    1    0    0    1    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0
##2   2   1   0   0   2   0   0   0   0   0   0   0   0   1    0    0    1    1    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0
##3   3   1   0   0   2   0   1   1   1   0   0   1   1   1    0    0    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    2    0    0
##4   4   1   0   0   2   0   0   0   0   0   0   1   0   1    1    0    0    1    0    0    0    0    1    0    0    0    0    0    0    0    1    1    0    0
##5   5   1   0   0   2   0   0   0   1   0   0   1   1   1    0    0    2    0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    0    0
##6   6   1   0   0   2   0   0   0   1   0   0   1   0   0    0    0    1    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
##   rs30 rs31 rs32 rs33 rs34 rs35 rs36 rs37 rs38 rs39 rs40 rs41 rs42 rs43 rs44 rs45 rs46 rs47 rs48 rs49 rs50
##1    0    0    1    0    0    1    0    0    1    2    0    2    0    0    0    0    0    2    2    0    0
##2    0    0    1    1    0    0    0    0    1    1    1    0    0    2    2    0    0    1    2    0    0
##3    0    1    1    0    0    0    0    0    1    1    0    1    0    0    1    0    0    1    2    0    0
##4    0    0    1    0    0    1    0    0    1    1    0    0    0    0    0    0    0    1    1    0    0
##5    0    1    1    0    0    1    0    0    0    1    0    0    0    0    0    0    0    1    2    0    0
##6    0    0    1    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    1    0    0

##"pid", "iid", "fid", "mid" and "sex" in the first five columns are the pedigree structure information as described in the "GEXCIS-manual". "rs1-rs50" are genotype codes for 50 SNPs, and each genotype is coded as 0, 1 or 2, indicating the number of the minor alleles.


################  step 4: explaining the corresponding functions  ##############

##Users can obtain the point estimate and the penalized point estimate of gamma for the genes, as well as the 95% CIs derived by the Fieller's and PF methods through "G_Frequen_XCI".


########################  step 5: running the functions  #######################

G_Frequen_XCI(phenotype=phenotype1,genotype=genotype1,trait_type="quantitative",phenotype_missing=NA,genotype_missing=NA,alpha=0.05)
##$penalized_point_estimate
##[1] 0.4763586

##$PF_lower
##[1] 0.06569505

##$PF_upper
##[1] 1.23757

##$PF_length
##[1] 1.171875

##$point_estimate
##[1] 0.4786828

##$F_lower
##[1] 0.06919327

##$F_upper
##[1] 1.347783

##$F_length
##[1] 1.27859

##$F_D
##[1] 0

######################  step 6: interpreting the results  ######################

##The penalized point estimate of gamma is 0.4763586, the 95% CI obtained by the PF method is (0.06569505, 1.23757), and the length of this CI is 1.171875;
##The point estimate of gamma is 0.4786828, the 95% CI obtained by the Fieller's method is (0.06919327, 1.347783). This interval is a continuous interval and has a length of 1.27859;
##The XCI pattern of the gene on the trait may be XCI-R or XCI-E.
##If the estimate of gamma is significantly different from 1, the gene is statistically inferred to undergo XCI-S, For example, gamma = 0.4 represents that the XCI-S on this gene is overall skewed towards the minor alleles, where only about 20% (0.4/2) of the cells have the minor allele active overall, and the other 80% of the cells have the major allele active overall. 


######################  Examples 12-17 have the similar steps  #################

##the structures of the remaining sample data are similar to phenotype1 and genotype1, also see the "GEXCIS-manual" for details
head(phenotype2)
head(phenotype3)
head(phenotype4)
head(phenotype5)
head(phenotype6)
head(phenotype7)
head(genotype2)
head(genotype3)

##example 12:
##quantitative trait without covariate
G_Frequen_XCI(phenotype=phenotype2,genotype=genotype1,trait_type="quantitative",phenotype_missing=NA,genotype_missing=NA,alpha=0.05)

##example 13:
##qualitative trait with covariate
G_Frequen_XCI(phenotype=phenotype5,genotype=genotype1,trait_type="qualitative",phenotype_missing=NA,genotype_missing=NA,alpha=0.05)

##example 14:
##qualitative trait without covariate
G_Frequen_XCI(phenotype=phenotype6,genotype=genotype1,trait_type="qualitative",phenotype_missing=NA,genotype_missing=NA,alpha=0.05)

##example 15:
##quantitative trait with covariate and missing values.
##both "phenotype" and "genotype" contain the missing values (denoted by NA)
G_Frequen_XCI(phenotype=phenotype3,genotype=genotype2,trait_type="quantitative",phenotype_missing=NA,genotype_missing=NA,alpha=0.05)

##example 16:
##qualitative trait with covariate and missing values.
##both "phenotype" and "genotype" contain the missing values (denoted by NA)
G_Frequen_XCI(phenotype=phenotype7,genotype=genotype2,trait_type="qualitative",phenotype_missing=NA,genotype_missing=NA,alpha=0.05)

##example 17:
##quantitative trait with covariate and missing values.
##both "phenotype" and "genotype" contain the missing values (denoted by 9)
G_Frequen_XCI(phenotype=phenotype4,genotype=genotype3,trait_type="quantitative",phenotype_missing=9,genotype_missing=9,alpha=0.05)