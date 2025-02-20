// Adapter trim with Trimmomatic
process trimReads {
    container "quay.io/biocontainers/trimmomatic:0.35--6"

	// Retry on fail at most three times
    errorStrategy 'retry'
    maxRetries 1

    input:
        tuple val(base), file(R1), file(R2)// from input_read_ch
        file ADAPTERS
    output:
        tuple val(base), file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz")// into Trim_out_ch

    publishDir "${params.OUTDIR}trimmed_fastqs", mode: 'copy', pattern: '*'

    script:
    """
    #!/bin/bash

    ls -latr
    trimmomatic PE -threads ${task.cpus} ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}.R2.paired.fastq.gz ${base}.R2.unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 -trimlog ./logFile.txt

    """
}

// Use bbduk to filter reads that match rRNA more stringently
process filterTp {
    container "staphb/bbtools:39.00"

    // Retry on fail at most three times
    errorStrategy 'retry'
    maxRetries 2

    input:
        tuple val(base), file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz")// from Trim_out_ch
        file REF_FASTAS
    output:
        tuple val(base),file("${base}_matched_rRNA_r1.fastq.gz"), file("${base}_matched_rRNA_r2.fastq.gz")// into Trimmed_filtered_reads_ch1
        tuple val(base),file("${base}_unmatched_rRNA_r1.fastq.gz"), file("${base}_unmatched_rRNA_r2.fastq.gz")// into Trimmed_unmatched_reads

    publishDir "${params.OUTDIR}rRNA_filtered_fastqs", mode: 'copy', pattern: '*matched_rRNA*.fastq.gz'

    script:
    """
    #!/bin/bash

    ls -latr
    echo "Filtering trimmed reads of ${base} against TPA rRNA reference..."
    bbduk.sh in1='${base}.R1.paired.fastq.gz' in2='${base}.R2.paired.fastq.gz' out1='${base}_unmatched_rRNA_r1.fastq.gz' out2='${base}_unmatched_rRNA_r2.fastq.gz' outm1='${base}_matched_rRNA_r1.fastq.gz' outm2='${base}_matched_rRNA_r2.fastq.gz' ref=${REF_FASTAS} k=31 hdist=1 mcf=0.98 rieb=f stats='${base}_stats_tp.txt' overwrite=TRUE t=16 -Xmx50g

    """
}

// Take unmatched reads from rRNA filter and map to TPA genome less stringently
process mapUnmatchedReads {
    container "staphb/bbtools:39.00"

    // Retry on fail at most three times
    errorStrategy 'retry'
    maxRetries 2

    input:
        tuple val(base),file("${base}_unmatched_rRNA_r1.fastq.gz"), file("${base}_unmatched_rRNA_r2.fastq.gz")// from Trimmed_unmatched_reads
        file REF_FASTAS_MASKED
    output:
        tuple val(base),file("${base}_matched_tpa_r1.fastq.gz"), file("${base}_matched_tpa_r2.fastq.gz")// into TPA_matched_reads
        tuple val(base),file("${base}_matched_tpa_r1.fastq.gz"), file("${base}_matched_tpa_r2.fastq.gz")// into TPA_matched_reads2

    publishDir "${params.OUTDIR}TPA_filtered_fastqs", mode: 'copy', pattern: '*matched_tpa*.fastq.gz'

    script:
    """
    #!/bin/bash

    ls -latr
    bbduk.sh in1='${base}_unmatched_rRNA_r1.fastq.gz' in2='${base}_unmatched_rRNA_r2.fastq.gz' out1='${base}_unmatched_tpa_r1.fastq.gz' out2='${base}_unmatched2_tpa_r2.fastq.gz' outm1='${base}_matched_tpa_r1.fastq.gz' outm2='${base}_matched_tpa_r2.fastq.gz' ref=${REF_FASTAS_MASKED} k=31 hdist=2 stats='${base}_stats_tp.txt' overwrite=TRUE t=14 -Xmx105g

    """
}

// More filtering of reads for de novo assembly
process moreFiltering {
    container "staphb/bbtools:39.00"

    // Retry on fail at most three times
    errorStrategy 'retry'
    maxRetries 2

    input:
        tuple val(base),file("${base}_matched_tpa_r1.fastq.gz"), file("${base}_matched_tpa_r2.fastq.gz")// from TPA_matched_reads2
        file REPEAT_FILTER_FASTA
    output:
        tuple val(base),file("${base}_no_repeat_genes_r1.fastq.gz"), file("${base}_no_repeat_genes_r2.fastq.gz")// into extra_reads

    publishDir "${params.OUTDIR}extra_filtered_fastqs_for_denovo", mode: 'copy', pattern: '*no_repeat_genes*.fastq.gz'

    script:
    """
    #!/bin/bash

    ls -latr
    bbduk.sh in1='${base}_matched_tpa_r1.fastq.gz' in2='${base}_matched_tpa_r2.fastq.gz' out1='${base}_no_repeat_genes_r1.fastq.gz' out2='${base}_no_repeat_genes_r2.fastq.gz' ref=${REPEAT_FILTER_FASTA} k=45 hdist=2 stats='${base}_stats_tp.txt' overwrite=TRUE t=14 -Xmx105g

    """
}

// Take matched reads and map to NC_021508 reference
process mapReads {
    container "quay.io/biocontainers/bowtie2:2.4.1--py37h8270d21_3"

    publishDir "${params.OUTDIR}mapSams", mode: 'copy', pattern: '*'

    input:
        tuple val(base),file("${base}_matched_tpa_r1.fastq.gz"), file("${base}_matched_tpa_r2.fastq.gz"),file("${base}_matched_rRNA_r1.fastq.gz"), file("${base}_matched_rRNA_r2.fastq.gz")
        
        file(NC_021508)
        file(NC_021508_1)
        file(NC_021508_2)
        file(NC_021508_3)
        file(NC_021508_4)
        file(NC_021508_5)
        file(NC_021508_6)

    output:
        tuple val(base),file("${base}.sam")// into Aligned_sam_ch

    script:
    """
    #!/bin/bash

    ls -latr
    echo "Concatenating TPA and rRNA reads for ${base}..."
    cat ${base}_matched_tpa_r1.fastq.gz ${base}_matched_rRNA_r1.fastq.gz > ${base}_matched_r1.fastq.gz
    cat ${base}_matched_tpa_r2.fastq.gz ${base}_matched_rRNA_r2.fastq.gz > ${base}_matched_r2.fastq.gz
    bowtie2 -x ${params.REFERENCE} -1 '${base}_matched_r1.fastq.gz' -2 '${base}_matched_r2.fastq.gz' -p ${task.cpus} > ${base}.sam

    """
}

// Convert sam to bam
process samToBam {
    container "staphb/samtools:1.19"

    input:
        tuple val(base),file("${base}.sam")// from Aligned_sam_ch
    output:
        tuple val(base),file("${base}_firstmap.sorted.bam")// into Sorted_bam_ch
        file("*.fasta") optional true

    publishDir "${params.OUTDIR}firstmap_sorted_bam", mode: 'copy', pattern: '*_firstmap.sorted.bam'

    script:
    """
    #!/bin/bash

    ls -latr
    /usr/local/bin/samtools view -bS ${base}.sam > ${base}.bam
    /usr/local/bin/samtools sort -o ${base}_firstmap.sorted.bam ${base}.bam

    if [[ ${params.SKIP_DENOVO} == true ]]

    then

        samtools consensus ${base}_firstmap.sorted.bam  > ${base}_from_bam.fasta

    fi

    """

}

// Use Picard to remove duplicates and convert bam back to fastq for downstream remapping
process removeDuplicates{
    container "quay.io/biocontainers/picard:2.23.3--0"

    input:
        tuple val(base),file("${base}_firstmap.sorted.bam")// from Sorted_bam_ch
    output:
        tuple val(base),file("${base}_firstmap_dedup.bam")// into Sorted_dedup_bam_ch
        tuple val(base),file("${base}_firstmap_dedup.bam")// into Sorted_dedup_bam_ch2
        tuple val(base),file("${base}_firstmap_dedup.bam")// into Sorted_dedup_bam_ch3
        tuple val(base),file("${base}_firstmap_dedup.bam")// into Sorted_dedup_bam_ch4

        tuple val(base),file("${base}_deduped_r1.fastq"),file("${base}_deduped_r2.fastq")// into Deduped_reads
        tuple val(base),file("${base}_deduped_r1.fastq"),file("${base}_deduped_r2.fastq")// into Deduped_reads_ch2
        tuple val(base),file("${base}_deduped_r1.fastq"),file("${base}_deduped_r2.fastq")// into Deduped_reads_ch3


    publishDir "${params.OUTDIR}deduped_bams", mode: 'copy', pattern: '*_firstmap_dedup.bam'
    publishDir ("${params.OUTDIR}deduped_fastqs", mode: 'copy', pattern: '*.fastq')

    script:
    """
    #!/bin/bash

    /usr/local/bin/picard MarkDuplicates INPUT=${base}_firstmap.sorted.bam OUTPUT=${base}_firstmap_dedup.bam METRICS_FILE=${base}_metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT
    /usr/local/bin/picard SamToFastq I=${base}_firstmap_dedup.bam FASTQ=${base}_deduped_r1.fastq SECOND_END_FASTQ=${base}_deduped_r2.fastq VALIDATION_STRINGENCY=SILENT

    """
}

// Call variants with Freebayes
process callVariants {
    container "quay.io/biocontainers/freebayes:1.3.2--py37h26878c9_2"

    input:
        tuple val(base),file("${base}_firstmap_dedup.bam")// from Sorted_dedup_bam_ch3
        file(NC_021508)
    output:
        tuple val(base), file("${base}.vcf")// into Freebayes_vcf

    publishDir "${params.OUTDIR}vcfs", mode: 'copy', pattern: '*.vcf'

    script:
    """
    #!/bin/bash

    ls -latr
    /usr/local/bin/freebayes -f ${NC_021508} -p 1 ${base}_firstmap_dedup.bam -v ${base}.vcf
    """
}

process downSample {
    container "staphb/seqtk:1.4"

    input:
        tuple val(base),file("${base}_deduped_r1.fastq"),file("${base}_deduped_r2.fastq")
    output:
        tuple val(base),file("${base}_DownSample_r1.fastq"),file("${base}_DownSample_r2.fastq")

    publishDir "${params.OUTDIR}downsample/${base}/", mode: 'copy', pattern: '*'

    script:
    """
    #!/bin/bash

    seqtk sample -s100 ${base}_deduped_r1.fastq ${params.DOWN_SAMPLE} > ${base}_DownSample_r1.fastq 
    seqtk sample -s100 ${base}_deduped_r2.fastq ${params.DOWN_SAMPLE} > ${base}_DownSample_r2.fastq 

    """
}

//Skip process if Skip Denovo Selected
// De novo assemble matched reads with Unicycler
process deNovoAssembly {
    container "quay.io/biocontainers/unicycler:0.4.4--py37h13b99d1_3"
    
    errorStrategy 'ignore'

    input:
        tuple val(base),file("${base}_DownSample_r1.fastq"),file("${base}_DownSample_r2.fastq")
    output:
        tuple val(base),file("${base}_assembly.gfa"),file("${base}_assembly.fasta")// into Unicycler_ch
        file("*")// into Unicycler_dump_ch

    publishDir "${params.OUTDIR}unicycler_output/${base}/", mode: 'copy', pattern: '*'

    script:
    """
    #!/bin/bash

    unicycler -1 ${base}_DownSample_r1.fastq -2 ${base}_DownSample_r2.fastq -o ./ -t ${task.cpus}

    cp assembly.gfa ${base}_assembly.gfa
    cp assembly.fasta ${base}_assembly.fasta

    """
}

//Skip process if Skip Denovo Selected
// Merges assembly and mapping to make consensus sequence
process mergeAssemblyMapping {
    container "quay.io/michellejlin/tpallidum_wgs"

    errorStrategy 'ignore'

    input:
        tuple val(base),file("${base}_assembly.gfa"),file("${base}_assembly.fasta"),file("${base}_firstmap_dedup.bam")// from Unicycler_ch
        file(NC_021508)
        file(NC_021508_BWA1)
        file(NC_021508_BWA2)
        file(NC_021508_BWA3)
        file(NC_021508_BWA4)
        file(NC_021508_BWA5)
        file(TP_MAKE_SEQ)
    output:
        tuple val(base),file("${base}_consensus.fasta") optional true // into Consensus_ch
        tuple val(base),file("${base}_consensus.fasta") optional true // into Consensus_ch2
        tuple val(base),file("*.sorted.bam") optional true// into Scaffold_bams_ch
        tuple val(base),file("*assembly*") optional true// into Assembly_ch

    publishDir "${params.OUTDIR}merged_assembly_mapping_consensus", mode: 'copy', pattern: '*_consensus.fasta'
    publishDir "${params.OUTDIR}scaffold_bams", mode: 'copy', pattern: '*.sorted.bam'

    script:
    """

    #!/bin/bash
    echo ${base}
    ls -latr

    cp ${base}_assembly.gfa assembly.gfa
    cp ${base}_assembly.fasta assembly.fasta

    Rscript --vanilla ${TP_MAKE_SEQ} \'${base}\' \'${params.REFERENCE}\'


    """
}

process remapReads {
    container "quay.io/michellejlin/tpallidum_wgs"

    input:
        tuple val(base),file("${base}_consensus.fasta"),file("${base}_deduped_r1.fastq"),file("${base}_deduped_r2.fastq")// from Consensus_ch
    output:
        tuple val(base),file("${base}_remapped.sorted.bam")// into Remapped_bam_ch
        tuple val(base),file("${base}_remapped.sorted.bam"),file("${base}_remapped.sorted.bam.bai"),file("${base}_consensus.fasta")// into Pilon_ch
        tuple val(base),file("*")// into Remap_reads_all_ch

    publishDir "${params.OUTDIR}remapped_bams", mode: 'copy', pattern: '*_remapped.sorted.bam'
    publishDir "${params.OUTDIR}remapped_bams/${base}/", mode: 'copy', pattern: '*'


    script:
    """
    ls -latr

    /usr/local/miniconda/bin/bowtie2-build -q ${base}_consensus.fasta ${base}_aligned_scaffolds_${params.REFERENCE}
    /usr/local/miniconda/bin/bowtie2 -x ${base}_aligned_scaffolds_${params.REFERENCE} -1 ${base}_deduped_r1.fastq -2 ${base}_deduped_r2.fastq -p ${task.cpus} | /usr/src/samtools-1.9/samtools view -bS - > ${base}_remapped.bam
    
    /usr/src/samtools-1.9/samtools sort -o ${base}_remapped.sorted.bam ${base}_remapped.bam

    /usr/src/samtools-1.9/samtools index -b ${base}_remapped.sorted.bam ${base}_remapped.sorted.bam.bai
    """
}

process remapReads_2 {
    container "quay.io/michellejlin/tpallidum_wgs"

    input:
        file("${base}_consensus.fasta")
        tuple val(base),file("${base}_deduped_r1.fastq"),file("${base}_deduped_r2.fastq")
    output:
        tuple val(base),file("${base}_remapped.sorted.bam")// into Remapped_bam_ch
        tuple val(base),file("${base}_remapped.sorted.bam"),file("${base}_remapped.sorted.bam.bai"),file("${base}_consensus.fasta")// into Pilon_ch
        tuple val(base),file("*")// into Remap_reads_all_ch

    publishDir "${params.OUTDIR}remapped_bams", mode: 'copy', pattern: '*_remapped.sorted.bam'
    publishDir "${params.OUTDIR}remapped_bams/${base}/", mode: 'copy', pattern: '*'


    script:
    """
    ls -latr

    /usr/local/miniconda/bin/bowtie2-build -q ${base}_consensus.fasta ${base}_aligned_scaffolds_${params.REFERENCE}
    /usr/local/miniconda/bin/bowtie2 -x ${base}_aligned_scaffolds_${params.REFERENCE} -1 ${base}_deduped_r1.fastq -2 ${base}_deduped_r2.fastq -p ${task.cpus} | /usr/src/samtools-1.9/samtools view -bS - > ${base}_remapped.bam
    
    /usr/src/samtools-1.9/samtools sort -o ${base}_remapped.sorted.bam ${base}_remapped.bam

    /usr/src/samtools-1.9/samtools index -b ${base}_remapped.sorted.bam ${base}_remapped.sorted.bam.bai
    """
}

process pilonPolishing {
    container "staphb/pilon:1.23.0"

    input:
        tuple val(base),file("${base}_remapped.sorted.bam"),file("${base}_remapped.sorted.bam.bai"),file("${base}_consensus.fasta")// from Pilon_ch
    output:
        tuple val(base),file("${base}_pilon.fasta")// into Pilon_fasta_ch
        tuple val(base),file("${base}_pilon.fasta")// into Pilon_fasta_ch2

    publishDir "${params.OUTDIR}pilon", mode: 'copy', pattern: '*'

    script:
    """
    ls -latr

    /pilon/pilon --genome ${base}_consensus.fasta --frags ${base}_remapped.sorted.bam --output ${base}_pilon --vcf --changes
    """
}

process remapPilon {
    container "quay.io/michellejlin/tpallidum_wgs"

    input:
        tuple val(base),file("${base}_pilon.fasta"),file("${base}_deduped_r1.fastq"),file("${base}_deduped_r2.fastq")// from Pilon_fasta_ch
    output:
        tuple val(base),file("${base}_pilon_remapped.sorted.bam")// into Pilon_bam_ch
        tuple val(base),file("*")// into Remap_pilon_all_ch

    publishDir "${params.OUTDIR}remapped_pilon_bams", mode: 'copy', pattern: '*pilon_remapped.sorted.bam'
    publishDir "${params.OUTDIR}remapped_pilon_bams/${base}", mode: 'copy', pattern: '*'


    script:
    """
    ls -latr

    /usr/local/miniconda/bin/bowtie2-build -q ${base}_pilon.fasta ${base}_pilon_aligned_scaffolds_${params.REFERENCE}
    /usr/local/miniconda/bin/bowtie2 -x ${base}_pilon_aligned_scaffolds_${params.REFERENCE} -1 ${base}_deduped_r1.fastq -2 ${base}_deduped_r2.fastq -p ${task.cpus} | /usr/src/samtools-1.9/samtools view -bS - > ${base}_pilon_remapped.bam
    /usr/src/samtools-1.9/samtools sort -o ${base}_pilon_remapped.sorted.bam ${base}_pilon_remapped.bam

    /usr/src/samtools-1.9/samtools index -b ${base}_pilon_remapped.sorted.bam ${base}_pilon_remapped.sorted.bam.bai
    """
}

process generatePilonConsensus {
    container "quay.io/michellejlin/tpallidum_wgs"

    input:
        tuple val(base),file("${base}_pilon_remapped.sorted.bam"),file("${base}_pilon.fasta"),file("${base}_firstmap_dedup.bam")// from Pilon_bam_ch
        file(TP_GENERATE_CONSENSUS)
    output:
        tuple val(base),file("${base}_pilon_finalconsensusv2.fasta"),file("${base}_pilon_mappingstats.csv")// into Prokka_pilon_consensus_ch

    publishDir "${params.OUTDIR}finalconsensus_v2", mode: 'copy', pattern: '*_finalconsensusv2.fasta'

    script:
    """
    cp ${base}_firstmap_dedup.bam ${base}_pilon_firstmap_dedup.bam

    ls -latr
    Rscript --vanilla ${TP_GENERATE_CONSENSUS} \'${base}_pilon' \'${params.REFERENCE}\'
    """
}

// Annotate final consensus with prokka
process annotatePilonConsensus {
    container "quay.io/biocontainers/prokka:1.14.6--pl526_0"

    input:
    tuple val(base),file("${base}_pilon_finalconsensusv2.fasta"),file("${base}_pilon_mappingstats.csv")// from Prokka_pilon_consensus_ch
    file(REF_GB)

    output:
        file("*")//// into Annotated_pilon_ch

    publishDir "${params.OUTDIR}finalconsensus_pilon_prokka_annotations", mode: 'copy', pattern: '*'

    script:
    """

    ls -latr

    prokka --proteins ${REF_GB} --outdir ./ --force --prefix ${base} ${base}_pilon_finalconsensusv2.fasta

    """
}


process annotateVCFs {
    container "cave42/snippy_tp_wgs"

    input:

        tuple val(base), file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz")
        tuple val(base),file("${base}_firstmap_dedup.bam")
        tuple val(base),file("${base}_firstmap.sorted.bam")
        file(REF_GB)
        file(NC_021508)
        tuple val(base),file("${base}_deduped_r1.fastq"),file("${base}_deduped_r2.fastq")
        

    output:

        file("*")

    publishDir "${params.OUTDIR}VCF_Annotations", mode: 'copy', pattern: '*'

    script:

    """
    #!/bin/bash

    snippy --outdir ${base} --ref ${REF_GB} --R1 ${base}_deduped_r1.fastq --R2 ${base}_deduped_r2.fastq --cpus 16  --report
    
    """

}

process stats {

    container "staphb/samtools:1.19"

    input:

        tuple val(sample), path(raw_1), path(raw_2), path(trim_1), path(trim_2), path(firstmap_bam), path(dedup_bam), path(pilon_bam), path(final_fasta), path(final_fasta_csv) 

    output:

        file("*.csv")

    script:

    """
    #!/bin/bash

    raw=\$(gzip -dc ${raw_1} | wc -l)
    raw=\$((raw / 2))

    trimmed=\$(gzip -dc ${trim_1} | wc -l)
    trimmed=\$((trimmed / 2))

    firstmap_mapped=\$(samtools view -c -F 4 ${sample}_firstmap.sorted.bam)

    dedup_mapped=\$(samtools view -c -F 4 ${sample}_firstmap_dedup.bam)

    final=\$(samtools view -c -F 4 ${sample}_pilon_remapped.sorted.bam)

    amb=\$(grep -o 'N' ${sample}_pilon_finalconsensusv2.fasta | wc -l)

    echo "${sample},"\$raw","\$trimmed","\$firstmap_mapped","\$on_target_pct","\$dedup_mapped","\$Dedup","\$final","\$amb >> ${sample}_stats.csv

    """

}

process mlst {

    container "staphb/mlst:2.23.0-2024-04-01"

    input:

    tuple val(base),file("${base}_pilon_finalconsensusv2.fasta"),file("${base}_pilon_mappingstats.csv")

    output:

    file '*.csv'

    script:
    """

    #!/bin/bash

    mlst --csv ${base}_pilon_finalconsensusv2.fasta > ${base}_mlst.csv

    """

}

process concatenateCSV {
    
    container "ubuntu:noble-20240407.1"

    input:

    file CSV
    file MLST

    output:

    file 'final_stats.csv'
    file 'mlst_concatenated.csv'

    script:
    """

    #!/bin/bash

    cat *_stats.csv > concatenated.csv

    echo "sample,raw_reads,trimmed_reads,firstmap_bam_reads,on_target_pct,dedup_bam_reads,Dedup Percent,final_mapped_reads,ambiguities" > final_stats.csv
    cat concatenated.csv >> final_stats.csv

    cat *_mlst.csv > mlst_concatenated.csv

    """
}

process finalStats {
    
    container "pandas/pandas:pip-all"

    input:

    file CSV
    file MLST
    file MLST_BD

    output:

    file ('*csv')

    publishDir "${params.OUTDIR}", mode: 'copy', pattern: '*'

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    stats_db = pd.read_csv('final_stats.csv')
    mlst_db = pd.read_csv('mlst_concatenated.csv', header=None)
    TP_mlst_db = pd.read_csv('TP_TypingDB.tsv', sep='\t')

    stats_db['on_target_pct'] = stats_db['firstmap_bam_reads'] / stats_db['trimmed_reads'] * 100

    stats_db['Dedup Percent'] = ( ( stats_db['trimmed_reads'] - stats_db['dedup_bam_reads'] ) * 100 ) / stats_db['trimmed_reads'] 

    mlst_db.iloc[:, 0] = mlst_db.iloc[:, 0].str.replace('_pilon_finalconsensusv2.fasta', '')

    headers = ['sample', 'Organism', 'ST', 'TP0136', 'TP0548', 'TP0705']
    mlst_db.columns = headers 

    mlst_db['ST'] = mlst_db['ST'].astype(str)
    TP_mlst_db['ST'] = TP_mlst_db['ST'].astype(str)
    mlst_db = pd.merge(mlst_db, TP_mlst_db[['ST','clonal_complex']], on='ST', how='left')

    stats_db = pd.merge(stats_db, mlst_db, on='sample', how='left')

    stats_db.to_csv('FINAL_STATS.csv', index=False)

    """
}