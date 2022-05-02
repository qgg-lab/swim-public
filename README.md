## SWIM analysis and server codes

The codes are largely organized based on programming languages. The main workflows are always programmed using `bash` scripts, which may call additional `bash`, SLURM `sbatch` (for use on a cluster), `R`, or `perl` scripts.

### Construction and evaluation of the haplotype reference panel

#### 1. prepare sequence files

codes in (bash/prepSeq.bash), which includes downloading and processing SRA sequences to conform to formats (directory structures, file names) that are compatible with downstream steps.

#### 2. mapping, post-mapping processing, and call variants in GVCF

codes in (bash/gvcfs.bash) which calls (bash/fq2gvcf.bash).

#### 3. combine GVCFs and genotype GVCFs

codes in (bash/combineGenotype.bash) and (bash/genotype.bash), including variant quality recalibration.

#### 4. QC, filtering

codes in (bash/QC.bash).

#### 5. evaluation of imputation software combinations, reference panel composition

codes in (bash/acc_ind550_imputation.bash)

### SWIM server

codes in the directory (server/)

### Figures

codes in the directory (figureCode/) and data needed to generate the figures in (figureData/)