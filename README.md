# Pbalsamifera_stomataEvol

This repo contains the data and scripts used to investigate the climatic drivers and genetic architecture of stomatal traits in **Populus balsamifera** 

Collaborators: Karl Fetter (UConn), Baxter Worthing (UVM) and Stephen Keller (UVM)

## Table of Contents
- [Pbalsamifera\_stomataEvol](#pbalsamifera_stomataevol)
  - [Table of Contents](#table-of-contents)
- [Quantiative Genetics](#quantiative-genetics)
  - [Kinship Matrix](#kinship-matrix)
    - [Create ind sets](#create-ind-sets)
    - [Filter SNPs](#filter-snps)
    - [Make Kinship matrix](#make-kinship-matrix)
    - [View kinship matrix](#view-kinship-matrix)

# Quantiative Genetics

## Kinship Matrix

Estimate a kinship matrix, more precisely a genomic relatedness matrix, with `AGHmatrix`. The matrix will be used structure the intercepts of the random effects, which are correlated as defined by the kinship matrix. To build the matrix, use the `Gmatrix` function from `AGHmatrix`. For details see `https://pubmed.ncbi.nlm.nih.gov/27902800/`.

The first hurdle is to choose individuals as input. Create a vector of balsamifera only `ind_codes`, removing BxT and other hybrids. Filter the 227K vcf file for missingness using the balsam-only individuals, filter for missingness, and filter for a MAF of 1/N. Ideally, the kinship matrix is made from putatively neutral sites, and I'll use LD filtering of the 227K SNPs to get an LD pruned set of independent loci.

The second set of hurdles is to convert the data from it's current format to the diploid binary (i.e. 11, 12, 21, 22) format that AGHmatrix expects. This can be done with tidyverse and/or plink.

### Create ind sets

Using the NewHybrids filial calls, create a list of `ind_codes` for the balsam-only genotypes.

```{R, eval = F}
# Filter for balsam-only
library(tidyverse)

# Read in the trait data for the test genotypes.
nh <- read.table("~/Dropbox/landscape_stomata/gitlab/data/NewHybrids_Ancestry.txt", T, '\t')

# Take a look at a contingency table of the ancestry and filial groups (sanity check)
table(nh$filial_call, nh$ancestry)
         
          BxA BxB BxD BxT
  F1       13   0  15   0
  F2        3   0   0   0
  P1.F1     0   0   0   2
  P1.F2     0   0   0  31
  P1.P1F1   0   0   0  11
  P1.P1F2   0   0   0  32
  Pure_1    0 351   0   0

# Read in file of known angustifolia & deltoides hybirds
ang_del = read.table("~/Dropbox/admap/results/global_ancestry/ang_del_inds",F,'\t')
names(ang_del) <- "ind_code"

# Filter for balsam and trichocarpa inds
ind_set <- nh %>% 
  filter(ancestry == "BxB") %>%
  filter(!ind_code %in% ang_del$ind_code) %>% # filter for ang_del
  select(ind_code, ancestry, prop_B) %>%
  distinct(ind_code, .keep_all = TRUE) %>%
  mutate(keep = paste(lapply(str_split(ind_code, "_"), '[[', 1 ),
                      lapply(str_split(ind_code, "_"), '[[', 2 ),
                      sep=" ")
         ) %>%
  select(keep)

ind_set <- nh %>% 
  filter(ancestry == "BxB") %>%
  filter(!ind_code %in% ang_del$ind_code) %>% # filter for ang_del
  select(ind_code, ancestry, prop_B) %>%
  distinct(ind_code, .keep_all = TRUE)

# Check that all the balsams have expected ancestry of 1
table(ind_set$prop_B)

# This yields 351 balsam-only genotypes
write.table(ind_set %>% select(ind_code), 
            "~/Dropbox/landscape_stomata/gitlab/data/balsam_inds",
            sep="\t", quote = F, row.names = F, 
            col.names = F)
```
### Filter SNPs

Convert SNP data from the entire genotyped data set (N = 485) and the reference datasets to a diploid genpop format. Total sample size for the file is 531. Inds in rows and genotypes in columns. Genotype is diploid and binary (1,2). I also want to filter for LD among the SNPs, which wil yield an optimal set of SNPs for estimating Fst.

To filter sites, first filter for the inds to keep and for monomorphic sites. Then output a list of sites to keep based on LD filter. Finally, input the list into an --extract command in plink.

```{unix, eval = F}
# Execute these commands in
cd ~/Dropbox/landscape_stomata/gitlab/data/genotypes

# Filter for inds, filter for complete sites. First use without --recode for a dry run
plink2 --bfile nsf_poplar_485inds_227607snps --keep ../balsam_inds_N351 --maf 0.01 --geno 0 --make-bed --out nsf_poplar_351inds_Nsnps.NoMissing

# Filter and create the new snp file
plink2 --bfile nsf_poplar_485inds_227607snps --keep ../balsam_inds_N351 --maf 0.01 --geno 0 --make-bed --export vcf --out nsf_poplar_351inds_214714snps.NoMissing

# Create list of LD filtered sites
# See how many sites are retained with 10 kb windows that move 100 sites at a time with r2 = 0, or complete absence of LD
plink --bfile pb_pt_418inds_33289snps.NoMissing --noweb --indep-pairwise 10 100 0 --out pb_pt_418inds_33289snps.NoMissing.LD0

plink2 --bfile nsf_poplar_351inds_12893snps.NoMissing --indep-pairwise 200 50 0.1 --out nsf_poplar_351inds_Nnps.NoMissing.LD10
_stop here. I think we need to call SNPs for just the balsam set_.


# Filter sites for LD and output bfiles
plink --bfile pb_pt_418inds_33289snps.NoMissing --noweb --extract pb_pt_418inds_33289snps.NoMissing.LD0.prune.in --make-bed --recode --out pb_pt_418inds_30209snps.NoMissing.LD0

# Convert the bed files to diploid genotypes with plink2.
# These will be coded as 0,1,2
plink2 --bfile pb_pt_418inds_30209snps.NoMissing.LD0 --make-pgen --export A --out pb_pt_418inds_30209snps.NoMissing.LD0.geno

# Output a vcf file too
plink2 --bfile pb_pt_418inds_30209snps.NoMissing --export vcf --out pb_pt_418inds_30209snps.NoMissing

```

The file pb_pt_418inds_30209snps.NoMissing.LD0.geno.raw contains the genotypes formatted in 0/1/2.

### Make Kinship matrix

I'll compute the genomic relatedness matrix based on VanRaden (2009) METHOD. 

```{r, eval = T,  message=F}
# Read in the data
snps_agh = read.table("~/Dropbox/admap/results/genotypes/pb_pt_pheno/pb_pt_418inds_30209snps.NoMissing.LD0.geno.raw",T,'\t')

# Convert the family id's to ind_codes & keep trailing zero
snps_agh <- snps_agh %>%
  mutate(ind_code = paste(snps_agh$FID, formatC(snps_agh$IID, width=2, flag="0"), sep="_")) %>%
  relocate(ind_code, .before = FID)

# Computing the additive relationship matrix based on VanRaden 2008
library(AGHmatrix)
G <- Gmatrix(SNPmatrix = as.matrix(snps_agh[-c(1:7)]), 
             missingValue=-9, 
             maf=2/418, 
             method="VanRaden")

# Output
# Initial data: 
# 	Number of Individuals: 418 
# 	Number of Markers: 30209 
# 
# Missing data check: 
# 	Total SNPs: 30209 
# 	 0 SNPs dropped due to missing data threshold of 1 
# 	Total of: 30209  SNPs 
# MAF check: 
# 	No SNPs with MAF below 0.004784689 
# Monomorphic check: 
# 	 9 monomorphic SNPs 
# 	Total: 30200 SNPs 
# Summary check: 
# 	Initial:  30209 SNPs 
# 	Final:  30200  SNPs ( 9  SNPs removed) 
#  
# Completed! Time = 11.752  seconds

# Set the row names to the ind_codes
row.names(G) <- snps_agh$ind_code

# Save it
write.table(round(G, 3),
            file = "~/Dropbox/admap/results/kinship/G.matrix",
            sep='\t',
            row.names = T,
            col.names = T,
            quote = F)

```

### View kinship matrix

```{r}

# View the kinship matrix & save it
# The names need fixing up, (or turn them off & use population labels.)
heatmap(G)

# View a histogram of the self-variance components (the diagonal)
h0 <- ggplot(data.frame(V1=diag(G)), aes(x=V1))  + 
  geom_histogram(bins = 25, color="black", fill="grey", lwd=0.3) +
  theme_minimal() + 
  labs(x="Self-variances of G")

# View a histogram of kinship (i.e. co-variance estimates).
h1 <- ggplot(data.frame(V1=G[upper.tri(G)]), aes(x=V1))  + 
  geom_histogram(bins = 25, color="black", fill="grey", lwd=0.3) +
  theme_minimal() + 
  labs(x="Kinship values (covariances of G)")

# Absolute value of covariance estimates
h2 <- ggplot(data.frame(V1=abs(G[upper.tri(G)])), aes(x=V1))  + 
  geom_histogram(bins = 25, color="black", fill="grey", lwd=0.3) +
  theme_minimal() + 
  labs(x="Absolute value of kinship")

# plot them together & save
pout <- cowplot::plot_grid(h0, h1, h2, ncol = 1, nrow = 3)
pout

setwd("~/Dropbox/admap/results/kinship/figs")
cowplot::save_plot(
  filename = "Kinship_histograms_v0.pdf",
  plot = pout,
  ncol = 1,
  nrow = 1,
  base_height = 7,
  base_asp = 0.75,
  base_width = NULL
  )

# Combine the kinship matrix and the histograms into a single plot.
cowplot::plot_grid(grid::grid.draw(pheat), cowplot::plot_grid(h0, h1, h2, ncol = 1), ncol = 2)

# Not working, I'll have to plot the kinship matrix in ggplot.

```
 
