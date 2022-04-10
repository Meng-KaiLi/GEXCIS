##The version of R needs to be 4.1.1 or later.
##Rtools needs to be downloaded and installed in advance.
##You can download and install the latest Rtools at the official website (https://cran.r-project.org/bin/windows/Rtools/).
##Install the latest version of "rstan" (install.packages("rstan")).

################################################################################
##############################  Example  #######################################
################################################################################
##install the package
##setwd to the root of the folder path of "GEXCIS_0.1.0.zip" and run the following command
install.packages("GEXCIS_0.1.0.zip",repos = NULL)

##or setwd to the root of the folder path of "GEXCIS_0.1.0.tar.gz" and run the following command
install.packages("GEXCIS_0.1.0.tar.gz",type='source')

##load the package
library("GEXCIS")

##view the sample data
head(phenotype1)
head(phenotype2)
head(phenotype3)
head(phenotype4)
head(phenotype5)
head(phenotype6)
head(phenotype7)
head(genotype1)
head(genotype2)
head(genotype3)
################################example of G_Bayes_XCI##########################
library(rstan)
##hide run log
rstan_options(javascript=FALSE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

##Special note: The example here is for demonstration only.
##In practical applications, chains_num, iter_num and warmup_num should be appropriately increased to ensure that the sampled samples are convergent.
##Otherwise, a series of warnings will appear, prompting to increase the number of iterations

##The following errors and warming can be ignored：

##Error in open.connection(con, open = mode) :
##  schannel: next InitializeSecurityContext failed: SEC_E_UNTRUSTED_ROOT (0x80090325) - 证书链是由不受信任的颁发机构颁发的。
##此外: Warning message:
##  In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
##  '-E' not found
##recompiling to avoid crashing R session


##example 1:
##quantitative phenotype with covariate
##the prior distribution of gamma is a normal distribution
##and the prior distribution of other unknown parameters is consistent with that in the article
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype1,genotype=genotype1,phenotype_type="quantitative",phenotype_missing=NA,allele_missing=NA,prior="normal",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 2:
##quantitative phenotype with covariate
##the prior distribution of gamma is a uniform distribution
##and the prior distribution of other unknown parameters is consistent with that in the article
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype1,genotype=genotype1,phenotype_type="quantitative",phenotype_missing=NA,allele_missing=NA,prior="uniform",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 3:
##quantitative phenotype with covariate
##the prior distribution of gamma or other unknown parameters is customize.
##Users are required to define the prior distribution of each parameter according to the research background, for example:
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
G_Bayes_XCI(phenotype=phenotype1,genotype=genotype1,phenotype_type="quantitative",phenotype_missing=NA,allele_missing=NA,
            prior="customize",model_customize=model_customize,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 4:
##quantitative phenotype without covariate
##the prior distribution of gamma is a normal distribution
##and the prior distribution of other unknown parameters is consistent with that in the article
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype2,genotype=genotype1,phenotype_type="quantitative",phenotype_missing=NA,allele_missing=NA,prior="normal",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 5:
##qualitative phenotype with covariate
##the prior distribution of gamma is a normal distribution
##and the prior distribution of other unknown parameters is consistent with that in the article
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype5,genotype=genotype1,phenotype_type="qualitative",phenotype_missing=NA,allele_missing=NA,prior="normal",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 6:
##qualitative phenotype without covariate
##the prior distribution of gamma is a uniform distribution
##and the prior distribution of other unknown parameters is consistent with that in the article
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype6,genotype=genotype1,phenotype_type="qualitative",phenotype_missing=NA,allele_missing=NA,prior="uniform",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)


##example 7:
##qualitative phenotype without covariate
##the prior distribution of gamma or other unknown parameters is customize.
##Users are required to define the prior distribution of each parameter according to the research background, for example:
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
G_Bayes_XCI(phenotype=phenotype6,genotype=genotype1,phenotype_type="qualitative",phenotype_missing=NA,allele_missing=NA,
            prior="customize",model_customize=model_customize,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 8:
##quantitative phenotype with covariate and missing values. Both phenotype and genotype contain the missing value (denoted by NA)
##the prior distribution of gamma is a uniform distribution
##and the prior distribution of other unknown parameters is consistent with that in the article
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype3,genotype=genotype2,phenotype_type="quantitative",phenotype_missing=NA,allele_missing=NA,prior="uniform",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 9:
##qualitative phenotype with covariate and missing values. Both phenotype and genotype contain the missing value (denoted by NA)
##the prior distribution of gamma is a uniform distribution
##and the prior distribution of other unknown parameters is consistent with that in the article
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype7,genotype=genotype2,phenotype_type="qualitative",phenotype_missing=NA,allele_missing=NA,prior="uniform",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)

##example 10:
##quantitative phenotype with covariate and missing values. Both phenotype and genotype contain the missing value (denoted by 9)
##the prior distribution of gamma is a normal distribution
##and the prior distribution of other unknown parameters is consistent with that in the article
set.seed(123456)
G_Bayes_XCI(phenotype=phenotype4,genotype=genotype3,phenotype_type="quantitative",phenotype_missing=9,allele_missing=9,prior="normal",model_customize=NULL,
            chains_num=2,iter_num=1000,warmup_num=500,acceptance_rate=0.99)



################################example of G_Frequen_XCI########################
##example 11:
##quantitative phenotype with covariate
G_Frequen_XCI(phenotype=phenotype1,genotype=genotype1,phenotype_type="quantitative",phenotype_missing=NA,allele_missing=NA,alpha=0.05)

##example 12:
##quantitative phenotype without covariate
G_Frequen_XCI(phenotype=phenotype2,genotype=genotype1,phenotype_type="quantitative",phenotype_missing=NA,allele_missing=NA,alpha=0.05)

##example 13:
##qualitative phenotype with covariate
G_Frequen_XCI(phenotype=phenotype5,genotype=genotype1,phenotype_type="qualitative",phenotype_missing=NA,allele_missing=NA,alpha=0.05)

##example 14:
##qualitative phenotype without covariate
G_Frequen_XCI(phenotype=phenotype6,genotype=genotype1,phenotype_type="qualitative",phenotype_missing=NA,allele_missing=NA,alpha=0.05)

##example 15:
##quantitative phenotype with covariate and missing values. Both phenotype and genotype contain the missing value (denoted by NA)
G_Frequen_XCI(phenotype=phenotype3,genotype=genotype2,phenotype_type="quantitative",phenotype_missing=NA,allele_missing=NA,alpha=0.05)

##example 16:
##qualitative phenotype with covariate and missing values. Both phenotype and genotype contain the missing value (denoted by NA)
G_Frequen_XCI(phenotype=phenotype7,genotype=genotype2,phenotype_type="qualitative",phenotype_missing=NA,allele_missing=NA,alpha=0.05)

##example 17:
##quantitative phenotype with covariate and missing values. Both phenotype and genotype contain the missing value (denoted by 9)
G_Frequen_XCI(phenotype=phenotype4,genotype=genotype3,phenotype_type="quantitative",phenotype_missing=9,allele_missing=9,alpha=0.05)
