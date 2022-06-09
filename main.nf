/*  Minor Allele Simulation + Detection Pipeline
 *  Usage: nextflow run /path/to/main.nf
 *
 *  Author: Mohammed Khalfan < mkhalfan@nyu.edu >
 *  NYU Center for Genetics and System Biology 2020
 */

// Setting some defaults here,
// can be overridden in config or via command line
params.out = "${params.outdir}/out"
params.tmpdir = "${params.outdir}/gatk_temp"

println "ref: $params.ref"
println "outdir: $params.out"

// Stage some files we will need
ref = file(params.ref)
error_model_fq_read1 = file(params.error_model_fq_read1)
error_model_fq_read2 = file(params.error_model_fq_read2)
mut_model_vcf = file(params.mut_model_vcf)
readsim_model_bam = file(params.readsim_model_bam)

// Prepare the fastq read pairs for input.
// Use the size parameter to not auto-group, and instead
// use the mapping through getBaseName() and subtract
// two regexs to get the ID.
// This enables support for CGSB sequence data file naming format
Channel
    .fromFilePairs( params.reads, size: -1)
    { file -> file.getBaseName() - ~/${params.read_pair_regex}/ - ~/.f*q/ }
    .set { read_pairs_ch }

Channel
    .fromFilePairs( params.bams, size: 1) 
    { file -> file.getBaseName() - ~/.bam/ }
    .set { bams_in_ch }

process buildIndex{
    input:
    path genome from params.ref

    output:
    path '*' into index_ch

    script:
    """
    bwa index $genome 
    samtools faidx $genome
    java -jar \$PICARD_JAR CreateSequenceDictionary R=${genome} O=${genome.baseName}.dict
    """
}

process genMutModel{
    output:
    file('MutModel.p') into mut_model_ch

    when:
    params.do_sim_reads
    

    script:
    """
    python2 /apps/neat-genreads/2.0/utilities/genMutModel.py \
        -r $ref \
        -m $mut_model_vcf \
        -o MutModel.p \
	--no-whitelist
    """
}

process seqErrorModel{
    output:
    file('SeqErrorModel.p') into seq_err_model_ch
    file('SeqErrorModel.p') into seq_err_model_ch_2

    when:
    params.do_sim_reads

    script:
    """
    python2 /apps/neat-genreads/2.0/utilities/genSeqErrorModel.py \
        -i $error_model_fq_read1 \
        -o SeqErrorModel.p \
        -i2 $error_model_fq_read2    
    """
}

process gcModel{
    output:
    file('gc_model.p') into gc_model_ch, gc_model_ch_2, gc_model_ch_3

    when:
    params.do_sim_reads

    script:
    """
    bedtools genomecov -d -ibam $readsim_model_bam > ${readsim_model_bam}.genomecov
    python2 /apps/neat-genreads/2.0/utilities/computeGC.py -r $ref -i ${readsim_model_bam}.genomecov -o gc_model.p
    """
}

process fraglenModel{
    output:
    file('fraglen.p') into fraglen_model_ch

    when:
    params.do_sim_reads

    script:
    """
    samtools view $readsim_model_bam | python2 /apps/neat-genreads/2.0/utilities/computeFraglen.py    
    """
}

process simulate_golden_snps{
    publishDir "${params.out}/vcfsim_1", mode:'copy'

    input:
    file(mut_model) from mut_model_ch
    file(seq_err_model) from seq_err_model_ch
    file(gc_model) from gc_model_ch

    output:
    set val(pair_id),
	file("${pair_id}_golden.vcf") into vcfsim_1_out_ch
    file("${pair_id}_golden.vcf") into vcfsim_1_out_ch_2
    file('*') into vcfsim_1_out

    when:
    params.do_sim_reads

    script:
    pair_id = params.fcid
    """
    python2 /apps/neat-genreads/2.0/genReads.py \
	-r $ref \
	-p 100 \
	-R 151 \
	-o $pair_id \
	-e $seq_err_model \
	--gc-model $gc_model \
	-m $mut_model \
	-M $params.mut_rate \
	--vcf \
	--no-fastq \
	-c $params.readsim_cov
    """
}

process simulate_pcr_snps{
    publishDir "${params.out}/vcfsim_2", mode:'copy'

    input:
    set val(pair_id),
	file(golden_vcf) from vcfsim_1_out_ch
    file(gc_model) from gc_model_ch_3
    each m_rate from params.m_rates

    output:
    set val(pair_id),
	file("*.vcf") into vcfsim_2_out_ch

    script:
    name = m_rate[0]
    m = m_rate[1]
    pair_id = pair_id + "_${name}"
    """
    python2 /apps/neat-genreads/2.0/genReads.py \
	-r $ref \
        -p 100 \
        -R 151 \
        -o $pair_id \
        -v $golden_vcf \
        --vcf \
	--no-fastq \
	-c $params.readsim_cov \
        -M $m
    """

}

process set_vcf_afs{
    publishDir "${params.out}/set_vcf_afs", mode:'copy'

    input:
    file (golden_vcf) from vcfsim_1_out_ch_2
    set val(pair_id), file(vcf) from vcfsim_2_out_ch
    each af from params.readsim_allele_fracs

    output:
    set val(pair_id),
        file("${pair_id}.vcf") \
	into set_vcf_afs_ch
    set val(pair_id),
	file("${pair_id}_golden.vcf") \
	into golden_vcf_ch_out, golden_vcf_comp_ch, \
	golden_vcf_bcftools_stats_ch 
    file("${pair_id}_golden.vcf") \
	into analyze_af_report_vcf

    script:
    pair_id = pair_id + "_AF_${af}"
    """
    prepare_neat_vcf.py $vcf $golden_vcf $ref $af $pair_id $params.seed
    """
}

process reorder_model_bam{
    output:
    file("reordered_model_bam.bam") into reorder_model_bam_ch

    when:
    params.do_sim_reads

    script:
    """
    java -jar \$PICARD_JAR ReorderSam \
	INPUT=$readsim_model_bam \
	OUTPUT=reordered_model_bam.bam \
	REFERENCE=$ref
    """
}

process simulate_reads{
    publishDir "${params.out}/readsim_2", mode:'copy'

    input:
    set val(pair_id), 
 	file(vcf) from set_vcf_afs_ch
    file(seq_err_model) from seq_err_model_ch_2
    file(fraglen_model) from fraglen_model_ch
    file(gc_model) from gc_model_ch_2
    file(model_bam) from reorder_model_bam_ch

    output:
    set val(pair_id),
        file("${pair_id}_read1.fq"),
	file("${pair_id}_read2.fq") \
	into readsim_out_ch
    file("*") into readsim_out

    script:
    """
    # Simulate reads inserting snps from the output
    # of the above step directly into the reads
    python2 /apps/neat-genreads/2.0/genReads.py \
	-r $ref \
	-p 100 \
	-R 151 \
	-o ${pair_id} \
	-e $seq_err_model \
	-v ${vcf} \
	--vcf \
	--pe-model $fraglen_model \
	--gc-model $gc_model \
	-c $params.readsim_cov \
	-M 0
    
    # Simulate reads inserting snps from the output
    # of the above step directly into the reads 
    # using ReSeq
    #reseq illuminaPE \
    #    -r $ref \
    #    -b $model_bam \
    #    -V ${pair_id}.vcf \
    #    -1 ${pair_id}_read1.fq \
    #    -2 ${pair_id}_read2.fq \
    #    -c $params.readsim_cov \
    # 	--noBias
    """
}

process downsample_readsim_fq{
    publishDir "${params.out}/downsampled_fastqs", mode:'copy'

    input:
    set pair_id, 
	file(read1), 
	file(read2),
	file(vcf) \
	from readsim_out_ch
	.join(golden_vcf_ch_out)
    each seed_frac_pair from params.readsim_downsample_fracs

    output:
    set val(pair_id),
        file("${pair_id}_read[12].fq") \
        into readsim_downsampled_ch
    //file("${pair_id}_golden.vcf") into downsample_bzip_tabix_vcf_ch
    //file("${pair_id}_golden.vcf") into golden_vcf_comp_ch
    //set val(pair_id), 
    //	file("${pair_id}_golden.vcf") \
    //	into golden_vcf_comp_ch
    //set val(mx_id), 
    //	file("${mx_id}_golden.vcf") \
    //	into golden_vcf_comp_ch_mx
    
    script:
    seed = seed_frac_pair[0]
    frac = seed_frac_pair[1]
    pair_id = pair_id + "_frac_" + frac
    //clean_id = pair_id.replaceFirst(/_m[12]_/, '_') 
    mx_id = pair_id.replaceFirst(/_m[12]_/, '_mx_')
    downsampled_dp = params.readsim_cov * frac


    if (frac < 1.0)
    """
    seqtk sample -s${seed} ${read1} $frac > ${pair_id}_read1.fq
    seqtk sample -s${seed} ${read2} $frac > ${pair_id}_read2.fq
    #modify_neat_dp.py $vcf $downsampled_dp > ${pair_id}_golden.vcf
    #cp ${pair_id}_golden.vcf ${mx_id}_golden.vcf
    """

    else
    """
    cp ${read1} ${pair_id}_read1.fq
    cp ${read2} ${pair_id}_read2.fq
    #modify_neat_dp.py $vcf $downsampled_dp > ${pair_id}_golden.vcf
    #cp ${pair_id}_golden.vcf ${mx_id}_golden.vcf
    """
}

process trim {
    publishDir "${params.out}/trimmed", mode:'copy'

    input:
    set pair_id,
        file(reads) from read_pairs_ch
	.mix(readsim_downsampled_ch)

    output:
    set val(pair_id),
	file("${pair_id}_trimmed_1.fq.gz"),
	file("${pair_id}_trimmed_2.fq.gz") \
	into trimmed_ch_bwa, trimmed_ch_star

    script:
    trim_adapters = ''
    if (params.adapters != '') {
	trim_adapters = "ILLUMINACLIP:${params.adapters}:2:30:10:8:true"
    } 
    """
    java -jar \$TRIMMOMATIC_JAR \
	PE \
	-phred33 \
	-threads ${task.cpus} \
	${reads[0]} \
	${reads[1]} \
	${pair_id}_trimmed_1.fq.gz \
	${pair_id}.unpair_trimmed_1.fq.gz \
	${pair_id}_trimmed_2.fq.gz \
	${pair_id}.unpair_trimmed_2.fq.gz \
	${trim_adapters} \
	LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20
    """
}

process star{

    container 'docker://gencorefacility/star:2.7.6a'

    publishDir "${params.out}/star", mode:'copy'

    input:
    set pair_id,
        file(read_1),
        file(read_2) from trimmed_ch_star

    output:
    set val(pair_id),
        file("${pair_id}_star.Aligned.out.sam") \
        into star_aligned_reads_ch

    when:
    params.aligner_star

    script:
    pair_id = pair_id + "_STAR"
    
    """
    zcat $read_1 | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ${pair_id}_trimmed_1.sorted.fq
    zcat $read_2 | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ${pair_id}_trimmed_2.sorted.fq

    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${params.star_ref} \
        --outFilterScoreMinOverLread 0.3 \
        --outFilterMatchNminOverLread 0.3 \
        --readFilesIn ${pair_id}_trimmed_1.sorted.fq ${pair_id}_trimmed_2.sorted.fq \
        --outReadsUnmapped Fastx \
        --outFileNamePrefix ${pair_id}_star.
    """
}

process addReadGroups {
    //publishDir "${params.out}/star_readgroups_added", mode:'copy'

    input:
    set val(sample_id),
        file(sam) from star_aligned_reads_ch
	.mix(bams_in_ch)

    output:
    set val(sample_id),
        file("${sample_id}_star.Aligned.out.RG.Sorted.bam") \
        into rg_added_ch

    script:
    """
    gatk AddOrReplaceReadGroups \
        -I ${sam} \
        -O ${sample_id}_star.Aligned.out.RG.Sorted.bam \
        --SORT_ORDER coordinate \
        -RGID ${sample_id} \
        -RGLB ${sample_id} \
        -RGPL ${params.pl} \
        -RGPU ${params.fcid} \
        -RGSM ${sample_id}
    """
}

process bwa {
    publishDir "${params.out}/aligned_reads", mode:'copy'
	
    input:
    file genome from ref
    file index from index_ch
    set pair_id, 
	file(read_1),
	file(read_2) from trimmed_ch_bwa
     
    output:
    set val(pair_id), 
	file("${pair_id}_aligned_reads.bam") \
	into aligned_reads_ch

    when:
    params.aligner_bwa || (!params.aligner_star and !params.aligner_bbmap)
	
    script:
    pair_id = pair_id + "_BWA"
    readGroup = "@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:${params.pl}\\tPM:${params.pm}\\tSM:${pair_id}"
    """
    bwa mem \
	-K 100000000 \
	-v 3 -t ${task.cpus} \
	-Y \
	-R \"${readGroup}\" \
	$genome \
	$read_1 \
	$read_2 \
	> ${pair_id}_aligned_reads.sam

    java -jar \$PICARD_JAR SortSam \
	I=${pair_id}_aligned_reads.sam \
	O=${pair_id}_aligned_reads.bam \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true
    """
}

process check_for_mapped_reads{

    input:
    set val(sample_id),
        file(bam) from aligned_reads_ch
        .mix(rg_added_ch)

    output:
    file ("${sample_id}.txt") optional true into no_mapped_reads_ch
    set val(sample_id),
        file("${bam.baseName}_mapped.bam") optional true into mapped_reads_ch

    script:
    """
    x=(\$(samtools view -F 4 $bam | wc -l))

    if [ \$x -eq 0 ]
    then
        echo $sample_id " FAILED mapped reads check"
        echo $sample_id > ${sample_id}.txt
    else
        mv $bam ${bam.baseName}_mapped.bam
    fi
    """
}

process markDuplicatesSpark  {
    publishDir "${params.out}/sorted", mode:'copy'

    input:
    set val(sample_id), 
	file(bam) from mapped_reads_ch

    output:
    set val(sample_id),
	file("${sample_id}_sorted_dedup.bam") \
	into sorted_dedup_bam_ch, sorted_dedup_ch_for_metrics, \
	downsample_bam_ch, pilon_ch, bcftools_ch, mutect2_ch, \
	tims_pipeline_ch, varscan_ch, ivar_ch, cliquesnv_ch
    set val(sample_id),
        file("${sample_id}_sorted_dedup.bam"),
        file("${sample_id}_sorted_dedup.bam.bai") \
	into bw_ch, freebayes_ch, lofreq_ch, qualimap_ch
    set val(sample_id),
	file("${sample_id}_dedup_metrics.txt") into dedup_qc_ch
    val(sample_id) into pair_id_ch
    set val(sample_id), file(bam) into genomecov_ch

    script:
    """
    gatk \
	MarkDuplicatesSpark \
	-I ${bam} \
	-M ${sample_id}_dedup_metrics.txt \
	-O ${sample_id}_sorted_dedup.bam
    """ 
}

process genomecov{
    input:
    set val(sample_id),
	file(bam) from genomecov_ch

    output:
    file("${sample_id}.tsv") into cov_plot_ch

    script:
    """
    bedtools genomecov -d -ibam $bam > ${sample_id}.tsv
    sed -i "s/^/${sample_id}\t/" ${sample_id}.tsv
    """
}

process cov_plot{
    publishDir "${params.out}/reports", mode:'copy'

    input:
    file(tsv) from cov_plot_ch.collect()

    output:
    file("*") into cov_plot_out_ch

    script:
    """
    cat *.tsv > cov_data.tsv
    sed -i '1i name\tsegment\tntpos\ttotalcount' ${params.fcid}_${workflow.runName}_cov_data.tsv
    cov_plots.R ${params.fcid}-${workflow.runName}
    """
}

process getMetrics {
    publishDir "${params.out}/metrics", mode:'copy'

    input:
    path index from index_ch
    set val(sample_id),
	file(sorted_dedup_reads) from sorted_dedup_ch_for_metrics

    output:
    set val(sample_id), 
            file("${sample_id}_alignment_metrics.txt"),
            file("${sample_id}_insert_metrics.txt"),
            file("${sample_id}_insert_size_histogram.pdf"),
            file("${sample_id}_depth_out.txt") \
            into metrics_output, metrics_multiqc_ch

    script:
    """
    java -jar \$PICARD_JAR \
        CollectAlignmentSummaryMetrics \
	R=${ref} \
        I=${sorted_dedup_reads} \
	O=${sample_id}_alignment_metrics.txt
    java -jar \$PICARD_JAR \
        CollectInsertSizeMetrics \
        INPUT=${sorted_dedup_reads} \
	OUTPUT=${sample_id}_insert_metrics.txt \
        HISTOGRAM_FILE=${sample_id}_insert_size_histogram.pdf 
    samtools depth -a ${sorted_dedup_reads} > ${sample_id}_depth_out.txt
    """
}

process timo{
    publishDir "${params.out}/timo", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam) from tims_pipeline_ch
    each timo_config from params.timo_configs

    output:
    file("${sample_id}_timo_${name}.vcf") \
	into tims_bzip_tabix_vcf_ch, timo_reps
    set val(sample_id),
        file("${sample_id}_timo_${name}.vcf") \
        into tims_vcf_ch, tims_bcftools_stats_ch
    set val("${sample_id}"),
	val("${sample_id}_timo_${name}") \
	into timo_rep_ids

    file("${sample_id}_timo_${name}_no-binom-check.vcf") \
        into tims_bzip_tabix_vcf_ch_2, timo_reps_2
    set val(sample_id),
        file("${sample_id}_timo_${name}_no-binom-check.vcf") \
        into tims_vcf_ch_2, tims_bcftools_stats_ch_2
    set val("${sample_id}"),
        val("${sample_id}_timo_${name}_no-binom-check") \
	into timo_rep_ids_2

    file ("${sample_id}_*") into tims_out_ch

    script:
    name = timo_config[0]
    samtools_params = timo_config[1]
    timo_params = timo_config[2]
    """
    samtools view -b $samtools_params -o filtered.bam $preprocessed_bam
    samtools index filtered.bam
    timo.v2.py $timo_params --infile filtered.bam --ref $ref

    ## parse_tims_output.py will look for all files created by 
    ## timo.v2.py in the working directory and convert
    ## them into a single VCF file
    parse_tims_output.py $ref ${sample_id}
    mv ${sample_id}.vcf ${sample_id}_timo_${name}.vcf
    parse_tims_output.py $ref ${sample_id} true
    mv ${sample_id}.vcf ${sample_id}_timo_${name}_no-binom-check.vcf
    mkdir ${sample_id}_${name}
    mv FILES/fullvarlist/*.csv ${sample_id}_${name}/.
    """

}

process test_pileup{

    //publishDir "${params.out}/failed", mode:'copy', pattern: '*.txt'

    input:
    set val(sample_id),
        file(bam) from varscan_ch

    output:
    file ("${sample_id}.txt") optional true into failed_ch
    set val(sample_id),
        file("${bam.baseName}_passed.bam") optional true into pileup_passed_ch

    script:
    """
    x=(\$(samtools mpileup  -f $ref $bam | wc -l))

    if [ \$x -eq 0 ]
    then
	echo $sample_id " FAILED pileup check"
	echo $sample_id > ${sample_id}.txt
    else
	mv $bam ${bam.baseName}_passed.bam
    fi
    """
}

process varscan {
    publishDir "${params.out}/varscan", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam) from pileup_passed_ch
    each vs_config from params.vs_configs

    output:
    file("${sample_id}_varscan_${name}.vcf") \
	into varscan_bzip_tabix_vcf_ch, vs_reps
    set val(sample_id),
        file("${sample_id}_varscan_${name}.vcf") \
        into varscan_vcf_ch, varscan_bcftools_stats_ch
    //file("${sample_id}_varscan_${name}.vcf") into vs_reps
    set val("${sample_id}"),
	val("${sample_id}_varscan_${name}") into vs_rep_ids

    script:
    name = vs_config[0]
    samtools_params = vs_config[1]
    vs_params = vs_config[2]
    """
    samtools mpileup $samtools_params -f $ref --max-depth 0 $preprocessed_bam |\
	java -jar \$VARSCAN_JAR mpileup2snp $vs_params \
	--output-vcf 1 > ${sample_id}_varscan_${name}.vcf
    """
}

process ivar{
    publishDir "${params.out}/ivar", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam) from ivar_ch
    each ivar_config from params.ivar_configs

    output:
    file("${sample_id}_ivar_${name}.vcf") \
	into ivar_bzip_tabix_vcf_ch, ivar_reps
    set val(sample_id), 
	file("${sample_id}_ivar_${name}.vcf") \
	into ivar_vcf_ch, ivar_bcftools_stats_ch
    //file("${sample_id}_ivar_${name}.vcf") into ivar_reps
    set val("${sample_id}"),
	val("${sample_id}_ivar_${name}") into ivar_rep_ids

    script:
    name = ivar_config[0]
    ivar_params = ivar_config[1]
    """
    samtools mpileup -aa -A -d 0 -B -Q 0 ${preprocessed_bam} |\
	ivar variants $ivar_params \
	-p ${sample_id}_ivar_${name} \
	-r $ref
    ivar_to_vcf.py ${sample_id}_ivar_${name}.tsv
    """
}

process lofreq{
    publishDir "${params.out}/lofreq", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam),
	file(preprocessed_bam_index) from lofreq_ch
    each lofreq_config from params.lofreq_configs

    output:
    file("${sample_id}_lofreq_${name}.vcf") \
	into lofreq_bzip_tabix_vcf_ch, lofreq_reps
    set val(sample_id),
        file("${sample_id}_lofreq_${name}.vcf") \
        into lofreq_vcf_ch, lofreq_bcftools_stats_ch
    file '*' into lofreq_out_ch
    //file("${sample_id}_lofreq_${name}.vcf") into lofreq_reps
    set val("${sample_id}"),
	val("${sample_id}_lofreq_${name}") into lofreq_rep_ids

    script:
    name = lofreq_config[0]
    lofreq_params = lofreq_config[1]
    """
    lofreq call-parallel $lofreq_params \
	--pp-threads ${task.cpus} \
	-f $ref \
	-o ${sample_id}_lofreq_${name}.vcf \
	$preprocessed_bam
    """
}

process cliquesnv{
    publishDir "${params.out}/cliquesnv", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam) from cliquesnv_ch

    output:
    file("${sample_id}_cliquesnv.vcf") into cliquesnv_bzip_tabix_vcf_ch
    set val(sample_id),
        file("${sample_id}_cliquesnv.vcf") \
        into cliquesnv_vcf_ch
    file '*' into cliquesnv_out_ch

    when:
    false

    script:
    """
    java  -Xmx58G -jar \$CLIQUESNV_JAR \
	-m snv-illumina-vc \
	-in $preprocessed_bam \
	-outDir vcf_out/ \
	-t 1 \
	-tf 0.005
    mv vcf_out/${preprocessed_bam.baseName}.vcf ${sample_id}_cliquesnv.vcf
    
    # CliqueSNV names contig as 'ref' for some reason, rename it to 'SARS-CoV2'
    # here for downstream analysis (ex: bcftools isec). Todo: Make this param
    sed -i 's/ref/SARS-CoV2/g' ${sample_id}_cliquesnv.vcf
    """
}

process freebayes{
    publishDir "${params.out}/freebayes", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam),
        file(preprocessed_bam_index) from freebayes_ch
    each freebayes_config from params.freebayes_configs

    output:
    file("${sample_id}_freebayes_${name}.vcf") \
	into freebayes_bzip_tabix_vcf_ch, freebayes_reps
    set val(sample_id),
        file("${sample_id}_freebayes_${name}.vcf") \
        into freebayes_vcf_ch, freebayes_bcftools_stats_ch
    file '*' into freebayes_out_ch
    //file("${sample_id}_freebayes_${name}.vcf") into freebayes_reps
    set val("${sample_id}"), 
	val("${sample_id}_freebayes_${name}") into freebayes_rep_ids

    script:
    name = freebayes_config[0]
    fb_params = freebayes_config[1]
    """
    freebayes $fb_params -f $ref $preprocessed_bam > ${sample_id}_freebayes_${name}.vcf
    """
}

process mutect2{
    publishDir "${params.out}/mutect2", mode:'copy'

    input:
    file genome from ref
    file index from index_ch
    set val(sample_id),
        file(preprocessed_bam) from mutect2_ch
    each m2_config from params.m2_configs

    output:
    set val(sample_id), 
	file("${sample_id}_mutect2_${name}_filtered.vcf") \
	into m2_vcf_ch, m2_bcftools_stats_ch
    set val(sample_id),
        file("${sample_id}_mutect2_${name}_unfiltered.vcf") \
        into m2_unfiltered_vcf_ch, m2_unfiltered_bcftools_stats_ch
    file("${sample_id}_mutect2_${name}_unfiltered.vcf") into mutect2_bzip_tabix_vcf_ch
    file '*' into mutect2_out_ch
    file("${sample_id}_mutect2_${name}_filtered.vcf") into m2_reps
    set val("${sample_id}"),
	val("${sample_id}_mutect2_${name}_filtered") into m2_rep_ids
    file("${sample_id}_mutect2_${name}_unfiltered.vcf") into m2_unfiltered_reps
    set val("${sample_id}"),
	val("${sample_id}_mutect2_${name}_unfiltered") into m2_unfiltered_rep_ids

    script:
    name = m2_config[0]
    m2_params = m2_config[1]
    """
    gatk Mutect2 \
	-R $genome \
	$m2_params \
	--max-reads-per-alignment-start 0 \
	-I $preprocessed_bam \
	-O ${sample_id}_mutect2_${name}_unfiltered.vcf
    gatk FilterMutectCalls \
	-R $genome \
	-V ${sample_id}_mutect2_${name}_unfiltered.vcf \
	-O tmp.vcf
    gatk SelectVariants \
	-R $genome \
	-V tmp.vcf \
	--exclude-filtered \
	-O ${sample_id}_mutect2_${name}_filtered.vcf
    """
}

process pilon{
    publishDir "${params.out}/pilon", mode:'copy'

    input:
    set val(sample_id),
	file(preprocessed_bam) from pilon_ch

    output:
    file("${sample_id}_pilon.vcf") into pilon_bzip_tabix_vcf_ch
    file '*' into pilon_out_ch

    when:
    false

    script:
    """
    java -Xmx16G -jar \$PILON_JAR \
	--genome $ref \
	--bam $preprocessed_bam \
	--fix bases \
	--changes \
	--vcf \
	--threads ${task.cpus} \
	--mindepth 10 \
	--output ${sample_id}_pilon_g
    
    gatk SelectVariants \
	-V ${sample_id}_pilon_g.vcf \
	-O ${sample_id}_pilon.vcf \
	--exclude-non-variants \
	--exclude-filtered
    """
}

process bcftools{
    publishDir "${params.out}/bcftools", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam) from bcftools_ch

    output:
    file("${sample_id}_bcftools.vcf") into bcftools_bzip_tabix_vcf_ch

    when:
    false

    script:
    """
    bcftools mpileup \
	--redo-BAQ \
	--adjust-MQ 50 \
	--gap-frac 0.05 \
	--max-depth 10000 \
	--max-idepth 200000 \
	--fasta-ref $ref \
	$preprocessed_bam | bcftools call \
	--ploidy 1 \
	--keep-alts \
	--multiallelic-caller \
	--variants-only \
	--output ${sample_id}_bcftools.vcf
    """
}

process haplotypeCaller {
    publishDir "${params.out}/haplotypecaller", mode:'copy'

    input:
    file genome from ref
    file index from index_ch
    set val(sample_id), 
	file(preprocessed_bam) from sorted_dedup_bam_ch
    each hc_config from params.hc_configs

    output:
    set val(sample_id),
        file("${sample_id}_hc_${name}.vcf") \
        into raw_snps_ch, raw_snps_qc_ch, hc_vcf_ch, \
	hc_bcftools_stats_ch
    file("${sample_id}_hc_${name}.vcf") into hc_reps
    set val("${sample_id}"),
	val("${sample_id}_hc_${name}") into hc_rep_ids

    script:
    name = hc_config[0]
    hc_params = hc_config[1]
    """
    gatk HaplotypeCaller $hc_params \
	--max-reads-per-alignment-start 0 \
	-R $genome \
	-I $preprocessed_bam \
	-O ${sample_id}_raw_variants.vcf
    gatk SelectVariants \
        -R $genome \
        -V ${sample_id}_raw_variants.vcf \
        -select-type SNP \
        -O ${sample_id}_hc_${name}.vcf
    """
}

process filterSnps {
    publishDir "${params.out}/filtered_snps", mode:'copy'
    
    input:
    set val(sample_id), 
	file(raw_snps) from raw_snps_ch

    output:
    set val(sample_id),
        file("${sample_id}_filtered_snps.vcf") \
        into filtered_snps_qc_ch

    // we're not filtering haplotypecaller snps
    when:
    false

    script:
    """
    gatk VariantFiltration \
	-R $ref \
	-V $raw_snps \
	-O ${sample_id}_filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

    # This script generates the _consensus_snps.vcf
    # and _eaf.vcf (empirical AF) files
    filter_variants.py ${sample_id}
    """
}

process snpEff {
    publishDir "${params.out}/snpeff", mode:'copy'

    input:
    //set val(pair_id), \
        //val(round), \
        //file(filtered_snps), \
        //file(filtered_snps_index) \
        //from filtered_snps_for_snpeff

    output:
    //file '*' into snpeff_out
    //file ("${pair_id}_filtered_snps.ann.vcf") into snpeff_bzip_tabix_vcf
    //file ("${pair_id}.csv") into multiqc_snpeff_csv_ch

    script:
    """
    #java -jar SNPEFF_JAR -v \
    #	-dataDir $params.snpeff_data \
    #	$params.snpeff_db \
    #   csvStats {pair_id}.csv \
    #   {filtered_snps} > {pair_id}_filtered_snps.ann.vcf
    """
}

process make_bw{
    publishDir "${params.out}/bigwig", mode:'copy'

    input:
    set val(id), 
	file(bam),
	file(bam_index) \
	from bw_ch

    output:
    file("${id}_coverage.bam.bw") into jbrowse_bw_ch 

    when:
    false

    script:
    """
    bamCoverage \
        -p max  \
        --bam $bam \
	--binSize 1 \
	--ignoreDuplicates \
	--minMappingQuality 20 \
        -o ${id}_coverage.bam.bw
    """
}

process downsample_bam{
    input:
    set val(sample_id), file(bam) from downsample_bam_ch

    output:
    set file("${sample_id}_downsampled.bam"),
        file("${sample_id}_downsampled.bam.bai") into jbrowse_bam_ch

    when:
    false

    script:
    """
    java -jar \$SORTSAMREFNAME_JAR \
        --samoutputformat BAM \
        $bam |\
        java -jar \$BIOSTAR_JAR \
        -n 75 \
        --samoutputformat BAM |\
        samtools sort -o ${sample_id}_downsampled.bam
    samtools index ${sample_id}_downsampled.bam
    """
}

process bzip_tabix_vcf{
    input:
    file(vcf) from pilon_bzip_tabix_vcf_ch
	.mix(mutect2_bzip_tabix_vcf_ch)
	.mix(freebayes_bzip_tabix_vcf_ch)
	.mix(tims_bzip_tabix_vcf_ch)
	.mix(tims_bzip_tabix_vcf_ch_2)
	.mix(varscan_bzip_tabix_vcf_ch)
	.mix(ivar_bzip_tabix_vcf_ch)
        .mix(lofreq_bzip_tabix_vcf_ch)
	.mix(cliquesnv_bzip_tabix_vcf_ch)
	//.mix(downsample_bzip_tabix_vcf_ch)

    output:
    file("*.vcf.gz*") into jbrowse_vcf_ch

    when:
    false

    script:
    """
    bgzip -c ${vcf} > ${vcf}.gz
    tabix -p vcf ${vcf}.gz
    """
}

/* Collect all the pair_ids and send them to pair_id_list_ch for use in jbrowse.
 */
pair_id_ch.collectFile(storeDir: "${params.out}/reports") { item ->
       [ "sample_ids.txt", item + '\n' ]
}.tap{pair_id_list_ch}

process jbrowse{
    cache false

    publishDir "${params.out}/trackList", mode:'copy'

    input:
    file '*' from jbrowse_bw_ch.collect()
    file '*' from jbrowse_bam_ch.collect()
    file '*' from jbrowse_vcf_ch.collect()
    file(pair_id_list) from pair_id_list_ch

    when:
    false

    script:
    """
    cgsb_upload2jbrowse -p MAD -d $params.fcid -s $pair_id_list -f . $ref
    """
}

process compare_replicates{
    cache false

    publishDir "${params.out}/compare_replicates", mode:'copy'

    input:
    file(vcf) from timo_reps
	.mix(timo_reps_2)
	.mix(vs_reps)
	.mix(freebayes_reps)
	.mix(m2_reps)
	.mix(m2_unfiltered_reps)
	.mix(hc_reps)
	.mix(lofreq_reps)
	.mix(ivar_reps)
	.collect()
    set val(id), val(filename) from timo_rep_ids
	.mix(timo_rep_ids_2)
	.mix(vs_rep_ids)
	.mix(freebayes_rep_ids)
        .mix(m2_rep_ids)
        .mix(m2_unfiltered_rep_ids)
        .mix(hc_rep_ids)
        .mix(lofreq_rep_ids)
	.mix(ivar_rep_ids)

    output:
    set val(id), 
	file("${mx}.vcf") \
	into compare_replicates_ch, \
	compare_replicates_bcftools_stats_ch

    when:
    id.contains("_m1_")

    script:
    m2 = filename.replace('_m1_', '_m2_')
    mx = filename.replace('_m1_', '_mx_')
    id = id.replace('_m1_', '_mx_') 
    """
    bgzip -c ${filename}.vcf > ${filename}.vcf.gz
    tabix -p vcf ${filename}.vcf.gz
    bgzip -c ${m2}.vcf > ${m2}.vcf.gz
    tabix -p vcf ${m2}.vcf.gz

    bcftools isec -c none ${filename}.vcf.gz ${m2}.vcf.gz -p ${mx}_isec
    # give it the 0002.vcf (m1-intersection) + m2.vcf
    merge_rep_vcfs.py ${mx}_isec/0002.vcf ${m2}.vcf
    cp out.vcf ${mx}.vcf
    """
}

process compare_afs{
    cache false

    publishDir "${params.out}/compare_afs", mode:'copy'

    input:
    set val(pair_id), 
	file(workflow_vcf) \
	from ivar_vcf_ch
	.mix(freebayes_vcf_ch)
	.mix(m2_vcf_ch)
	.mix(m2_unfiltered_vcf_ch)
	.mix(tims_vcf_ch)
	.mix(tims_vcf_ch_2)
	.mix(hc_vcf_ch)
	.mix(varscan_vcf_ch)
        .mix(lofreq_vcf_ch)
        .mix(cliquesnv_vcf_ch)
	.mix(compare_replicates_ch)
    file(vcf) from golden_vcf_comp_ch.collect()


    output:
    file("${workflow_vcf.baseName}_af_report.csv") into make_af_csv_output
    file("${workflow_vcf.baseName}_af_report_snp_count.txt") into af_report_snp_counts

    script:
    golden_vcf = pair_id.replaceFirst(/_mx_/, '_m1_')
	.replaceFirst(/_frac_\d\.\d+/, '')
	.replaceFirst(/_BWA/, '')
	.replaceFirst(/_STAR/, '') + '_golden.vcf'

    """
    # need to bzip and index the vcf for bcftools
    bgzip -c ${golden_vcf} > ${golden_vcf}.gz
    tabix -p vcf ${golden_vcf}.gz
    bgzip -c ${workflow_vcf} > ${workflow_vcf}.gz
    tabix -p vcf ${workflow_vcf}.gz

    bcftools isec -c none ${golden_vcf}.gz ${workflow_vcf}.gz -p isec

    compare_afs.py --golden_vcf ${golden_vcf}.gz --workflow_vcf ${workflow_vcf}.gz --out ${workflow_vcf.baseName}_af_report
    """
}

process qc {
    input:
    set val(sample_id),
	file("${sample_id}_alignment_metrics.txt"),
	file("${sample_id}_insert_metrics.txt"),
	file("${sample_id}_insert_size_histogram.pdf"),
	file("${sample_id}_depth_out.txt"),
	file("${sample_id}_dedup_metrics.txt"),
	file("${sample_id}_raw_snps.vcf"),
        file("${sample_id}_filtered_snps.vcf") \
	from metrics_output
	.join(dedup_qc_ch)
	.join(raw_snps_qc_ch)
	.join(filtered_snps_qc_ch)

    output:
    file("${sample_id}_report.csv") into parse_metrics_output

    when:
    false

    script:
    """
    parse_metrics.sh ${sample_id} > ${sample_id}_report.csv 
    """
}

/* Process qc above creates a report for each sample.
 * Below we compile these into a single report.
 * Same for af_csv files
 */
//parse_metrics_output.collectFile(name: "${workflow.runName}_report.csv", keepHeader: true, storeDir: "${params.out}/reports")
make_af_csv_output.collectFile(name: "${params.fcid}_${workflow.runName}_af_data.csv", keepHeader: true, storeDir: "${params.out}/reports").tap{af_report_in}
af_report_snp_counts.collectFile(name: "${params.fcid}_${workflow.runName}_af_report_snp_counts.txt", keepHeader: false, storeDir: "${params.out}/reports")
failed_ch.collectFile(name: "failed_varscan.txt", keepHeader: false, storeDir: "${params.out}/reports")
no_mapped_reads_ch.collectFile(name: "no_mapped_reads.txt", keepHeader: false, storeDir: "${params.out}/reports")

process analyze_af_report {
    publishDir "${params.out}/reports", mode:'copy'	

    input:
    file notebook from Channel.fromPath("bin/analyze_af_report.Rmd")
    file(af_report) from af_report_in
    file(golden_vcf) from analyze_af_report_vcf.take(1)

    output:
    file("*.html") into analyze_af_report_out

    script:
    """
    cp ${notebook} ${params.fcid}_${workflow.runName}_MAD_report.Rmd
    R -e "af_report = '${af_report}'; golden_vcf = '${golden_vcf}'; rmarkdown::render('${params.fcid}_${workflow.runName}_MAD_report.Rmd')"
    """
}

process qualimap{
    publishDir "${params.out}/qualimap", mode:'copy'
 
    container 'docker://pegi3s/qualimap'

    input:
    set val(pair_id),
        file(bam),
        file(bam_index) from qualimap_ch

    output:
    file('*') into multiqc_qualimap_ch

    script:
    """
    qualimap BamQC -bam $bam -outdir ${pair_id} -outformat HTML
    """
}

process bcftools_stats{
    input:
    set val(pair_id),
        file(vcf) from golden_vcf_bcftools_stats_ch
	.mix(ivar_bcftools_stats_ch)
        .mix(freebayes_bcftools_stats_ch)
        .mix(m2_bcftools_stats_ch)
        .mix(m2_unfiltered_bcftools_stats_ch)
        .mix(tims_bcftools_stats_ch)
	.mix(tims_bcftools_stats_ch_2)
        .mix(hc_bcftools_stats_ch)
        .mix(varscan_bcftools_stats_ch)
        .mix(lofreq_bcftools_stats_ch)
        .mix(compare_replicates_bcftools_stats_ch)

    output:
    file ("${vcf.baseName}_vcf_stats.vchk") into multiqc_bcftools_stats_ch

    script:
    """
    bcftools stats $vcf > ${vcf.baseName}_vcf_stats.vchk
    """
}

process multiqc{
    container 'docker://ewels/multiqc'

    publishDir "${params.out}/reports", mode:'copy'

    input:
    //file(snpeff_csv) from multiqc_snpeff_csv_ch.collect()
    file(bcftools_stats) from multiqc_bcftools_stats_ch.collect()
    file('*') from multiqc_qualimap_ch.collect()
    file('*') from metrics_multiqc_ch.collect()

    output:
    file("*.html") into multiqc_out

    script:
    """
    for f in \$(ls *metrics.txt);do sed -i 's/_sorted_dedup//g' \$f; done
    for f in \$(ls */genome_results.txt);do sed -i 's/_sorted_dedup//g' \$f; done
    for f in \$(ls *.vchk);do sed -i 's/_filtered_snps_2//g' \$f; done
    #for f in \$(ls *_snpeff_stats.csv);do sed -i 's/_snpeff_stats//g' \$f; done

    # make the config file for multiqc
    echo "remove_sections:" > multiqc_conf.yml
    echo "    - bcftools-stats_indel_plot" >> multiqc_conf.yml
    echo "    - qualimap-insert-size-histogram" >> multiqc_conf.yml
    echo "    - bcftools-stats_indel_plot" >> multiqc_conf.yml

    multiqc -f -c multiqc_conf.yml --ignore *.run .
    sed -i 's/>Bcftools</>VCF Stats</g' multiqc_report.html
    """
}

