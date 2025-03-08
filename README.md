# SNP freq simulator

This a simple python program that simulates the allele proportion distribution that may be observed from a sequencing experiment. The program doesnt simulate the effect of sequencing or mapping error on the proportions. 


## How to use

```
usage: snp_freq_sim_exp.py [-h] [-n_snps N_SNPS] [-p P] [-f F] [-d {poisson}] [-lam LAM] [-reps REPS] [-k_bins K_BINS]

Simulate SNP data for multiple parameter combinations and repetitions.

options:
  -h, --help      show this help message and exit
  -n_snps N_SNPS  Number of SNPs to simulate (integer >= 1), comma-separated list for multiple values
  -p P            Ploidy (integer >= 1), comma-separated list for multiple values
  -f F            Correction factor (positive integer), comma-separated list for multiple values
  -d {poisson}    Sequence depth distribution (default: poisson)
  -lam LAM        Lambda for Poisson distribution (integer >= 1), comma-separated list for multiple values
  -reps REPS      Number of repetitions per parameter combination (positive integer)
  -k_bins K_BINS  Number of bins for the histogram of observed proportions (positive integer)
```

### Parameter explanation

- -n_snps 




## Simulation: how does it work?




### Single simulation step by step

The simulator first samples the number of biallelic, triallelic.... sites. Up to the maximum multipliticity allowed by the ploidy. In a realistic scenario each increasing multiplicity is increasingly rare. The most comon multiplicity would be biallelic. The odds ratio of biallelic vs triallelic sites is controlled by the "-f" correction factor. Each multiplicity after triallelic is assumed to be about 10 times less probable than the previous. 

For each site (a total of "-n_snps" sites) a read depth is sampled from a poisson distribution with mean defined by "-lam".

Then to simulate the read counts per allele per site, a random sampling of n-alleles (the n is the allele multiplicity, biallelic, triallelic and so on...) a total of "read depth" times. So for example if i want to simulate a biallelic site with a sampled read depth of 10. I may end up sampling 3 reads form one allele and 7 from the other. The simulator also samples randomly the bases of each allele. 

The probabilities for the random sampling of n-alleles are chosen randomly from the pool of all possible unique probabilities that the sampling of n classes, "depth" times allow. This sampling may not be biologically accurate.

After that the allele proportion is calculated: the number of reads observed for an allele at a site "s" / The total number of reads observed for the site "s". 

The distribution of "-n_snps" allele proportions is normalized into a histogram of "k_bins". So that a model has an easier time handling the data. 



