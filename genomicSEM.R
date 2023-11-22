# Install GenomicSEM package first. Follow instructions in https://github.com/GenomicSEM/GenomicSEM

require(GenomicSEM)

# Munge GWAS summary statistics first

trait_files <-c('~/usr/ldsc/regenie/SSB/Bisexuality_male.regenie',
                '~/usr/ldsc/regenie/N_children/N_children_all.regenie',
                '~/usr/ldsc/regenie/Risk_taking/risk_taking_male.regenie')

hm3 <- "~/usr/ldsc/eur_w_ld_chr/w_hm3.snplist"

# name the traits
trait.names<-c("B","NC","RT")

N=c(409746,450748,434907)

munge(files=trait_files,hm3=hm3,trait.names=trait.names,N=N)

# munged files are saved as 'B.sumstats.gz','NC.sumstats.gz', and 'RT.sumstats.gz' in the working directory.


# enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
# prevalance of 
# Bisexual behavior in male+female populaiton: 7864/(361154+7864)
# risktaking behavior in male+female populaiton: 113968/436485
# Bisexual behavior in male populaiton: 3431/(163978+3431)
# risktaking behavior in male populaiton: 68945/200379
# Bisexual behavior in female populaiton: 4433/(197176+4433)
# risktaking behavior in female populaiton: 45023/236106

# specify appropriate prevalance, modify this for different sex groups.
sample.prev<-c(3431/(163978+3431),NA,68945/200379)

# vector of population prevalences
population.prev<-sample.prev

# the folder of LD scores
ld<-"~/usr/ldsc/eur_w_ld_chr/"

# the folder of LD weights [typically the same as folder of LD scores]
wld<-"~/usr/ldsc/eur_w_ld_chr/"

munged_files <- c(
    'B.sumstats.gz',
    'NC.sumstats.gz',
    'RT.sumstats.gz'
)

# run LDSC
LDSCoutput<-ldsc(traits=munged_files,sample.prev=sample.prev,population.prev=population.prev,
                 ld=ld,wld=wld,trait.names=trait.names)

CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")

MODEL <- 'NC ~ B + RT
          B ~ RT'

result <- usermodel(LDSCoutput, estimation = "DWLS", model = MODEL,  std.lv = TRUE, imp_cov = FALSE)

#print the results
result