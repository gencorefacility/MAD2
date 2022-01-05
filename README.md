# MAD2

#### Clone the git repository and edit the following parameters in the nextflow.config file:

`params.reads`: Path to input reads in fastq format. Set as `"[!*]"` if none or when simulating reads. 

`params.bams`: Path to aligned reads in bamformat. Set as `"[!*]"` if none or when simulating reads. 

`params.ref`: Path to reference genome fasta. BWA index, fasta index, and picard reference dictionary must exist in the same dir.

`params.adapters`: Path to adapter file. Can be empty when simulating reads. 

`params.pl`: Platform. Ex: Illumina

`params.pm`: Machine: Ex: NextSeq

`params.fcid`: Unique name for this analysis (alphanumeric, no spaces)

`params.outdir`: Output path

`params.read_pair_regex`: `"_read[12]"` if simulating reads. `"_n0[12]"` if using gencore data.

`params.do_sim_reads`: `true` if simulating reads, else `false`. 

`params.mut_model_vcf`: Path to VCF file to build mutation model

`params.error_model_fq_read1`: Read 1 of paired end fastq reads for error model

`params.error_model_fq_read2`: Read 2 of paired end fastq reads for error model

`params.readsim_model_bam`: Bam file for modeling

`params.mut_rate`: Mutation rate (between 0 and 0.3)

`params.readsim_cov`: Simulation coverage

`params.m_rates`: PCR error rates. Ex: 0.0005. This is an array with 2 elements, one for each PCR replicate. Do not change the length of this array or modify the `"m1"` or `"m2"` values. (Todo: get rid of this array and make this a single floating value parameter representing one error rate that will be applied to both PCR replicates)

`params.seed`: Random seed for PCR errors. This parameter ensures 'random' AF's between PCR replicates are equivilent. 

`params.readsim_downsample_fracs`: Simulation downsampling fractions [random seed, fraction]

`params.readsim_allele_fracs`: Simulation allele frequencies. Accepts floating points values and `'random'`

`params.freebayes_configs`, `params.m2_configs`, `params.hc_configs`, `params.lofreq_configs`, `params.ivar_configs`: Variant caller parameter configurations [name, params]

`params.vs_configs`, `params.timo_configs`: Variant caller parameter configurations including samtool configuration options [name, samtools params, variant caller params]

`workDir`: Nextflow workdir. Lots of temporary/intermediary files. Can delete once final analysis is complete. 

`process`: Configuration for scheduler. Replace 'slurm' with execution scheduler. Resources required will vary greatly depending on input, however reasonable defaults have been provided by default. 

#### Run the main.nf script, providing the path to the config, and specifying the `-with-singularity` parameter and the repo name:

`nextflow run main.nf -c <path_to_config> -with-singularity gencorefacility/mad:2`
