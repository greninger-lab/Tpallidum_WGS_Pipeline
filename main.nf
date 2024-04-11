#!/usr/bin/env nextflow

// Using Nextflow DSL-2 to account for logic flow of this workflow
nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    T. pallidum whole genome sequencing
    Usage:
    An example command for running the pipeline is as follows:
    nextflow run michellejlin/Tpallidum_WGS \\
        --INPUT         Input folder where all fastqs are located.
                        ./ can be used for current directory.
                        Fastqs should all be gzipped. This can be done with the command gzip *.fastq. [REQUIRED]
        --OUTDIR        Output directory. [REQUIRED]
        --REFERENCE     Reference used to map samples to, default is SS14 (NC_021508), options are: 
                        Nichols (NC_021490), Endemicum (NZ_CP007548), Pertenue (NC_016842)
        --SKIP_DENOVO   Skips denovo assembly uses assembly mapped to reference, much less memory intensive
        to limit default memory overide, without it is half of memory so left user input, remove all options have user input what they want
    """.stripIndent()
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*          SET UP CONFIGURATION VARIABLES            */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

params.INPUT = false
params.OUTDIR= false
params.SINGLE_END = false
//Choose what reference, default SS14
params.REFERENCE = "NC_021508"
params.SKIP_DENOVO = false

if (params.REFERENCE == "SS14"){
    params.REFERENCE = "NC_021508"
}

if (params.REFERENCE == "Nichols"){
    params.REFERENCE = "NC_021490"
}

if (params.REFERENCE == "Endemicum"){
    params.REFERENCE = "NZ_CP007548"
}

if (params.REFERENCE == "Pertenue"){
    params.REFERENCE = "NC_016842"
}

// if INPUT not set
if (params.INPUT == false) {
    println( "Must provide an input directory with --INPUT")
    exit(1)
}
// Make sure INPUT ends with trailing slash
if (!params.INPUT.endsWith("/")){
    params.INPUT = "${params.INPUT}/"
}
// if OUTDIR not set
if (params.OUTDIR == false) {
    println( "Must provide an output directory with --OUTDIR")
    exit(1)
}
// Make sure OUTDIR ends with trailing slash
if (!params.OUTDIR.endsWith("/")){
    params.OUTDIR = "${params.OUTDIR}/"
}

// Reference files
//s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/
ADAPTERS = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/All_adapters.fa")
//REF_FASTAS = file("${baseDir}/refs/TPA_refseqs.fasta")
REF_FASTAS = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/TPA_rRNA_refs.fasta")
REF_FASTAS_MASKED = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/Tp_refs_rRNA_masked.fasta")
//REF_FASTAS_MASKED = file("${baseDir}/refs/Tp_refs_rRNA_masked.fasta")
REF_FASTAS_TRIM = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/TPA_refseqs_trim.fasta")

NC_021508 = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/${params.REFERENCE}.fasta")
// bowtie2 indexes
NC_021508_1 = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/${params.REFERENCE}.1.bt2")
NC_021508_2 = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/${params.REFERENCE}.2.bt2")
NC_021508_3 = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/${params.REFERENCE}.3.bt2")
NC_021508_4 = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/${params.REFERENCE}.4.bt2")
NC_021508_5 = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/${params.REFERENCE}.rev.1.bt2")
NC_021508_6 = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/${params.REFERENCE}.rev.2.bt2")

// bwa indexes
NC_021508_BWA1 = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/${params.REFERENCE}.fasta.amb")
NC_021508_BWA2 = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/${params.REFERENCE}.fasta.ann")
NC_021508_BWA3 = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/${params.REFERENCE}.fasta.bwt")
NC_021508_BWA4 = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/${params.REFERENCE}.fasta.pac")
NC_021508_BWA5 = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/${params.REFERENCE}.fasta.sa")

REF_GB = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/${params.REFERENCE}.gb")

// Scripts
TP_MAKE_SEQ = file("${baseDir}/tp_make_seq.R")
TP_GENERATE_CONSENSUS = file("${baseDir}/tp_generate_consensus.R")

//REPEAT_FILTER_FASTA = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/repeat_filter_UPDATE.fasta")
REPEAT_FILTER_FASTA = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/repeat_filter_with_TprK.fasta")

CONVERT_TO_ANNOVAR = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/ANNOVAR_Tools/annovar/convert2annovar.pl")
GFF3_TO_GENE_PRED = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/ANNOVAR_Tools/gff3ToGenePred")
RETRIEVE_SEQ_FROM_FASTA = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/ANNOVAR_Tools/annovar/retrieve_seq_from_fasta.pl")
ANNOTATE_VARIATION = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/ANNOVAR_Tools/annovar/annotate_variation.pl")

GFF = file("s3://fh-pi-jerome-k-eco/greninger-lab/clomp-reference-data/tool_specific_data/Tpallidum_WGS/ANNOVAR_Tools/Anotations/${params.REFERENCE}.gff3")

// Read in fastq pairs into input_read_ch
if(params.SINGLE_END == false){
    input_read_ch = Channel
        .fromFilePairs("${params.INPUT}*_R{1,2}*.gz")
        .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT} ending with .gz" }
        .map { it -> [it[0], it[1][0], it[1][1]]}
}

//
// Import processes
//

include {trimReads} from './modules'
include {filterTp} from './modules'
include {mapUnmatchedReads} from './modules'
include {moreFiltering} from './modules'
include {mapReads} from './modules'
include {samToBam} from './modules'
include {removeDuplicates} from './modules'
include {callVariants} from './modules'
include {deNovoAssembly} from './modules'
include {mergeAssemblyMapping} from './modules'
include {remapReads} from './modules'
include {remapReads_2} from './modules'
include {pilonPolishing} from './modules'
include {remapPilon} from './modules'
//include {generateConsensus} from './modules'
include {generatePilonConsensus} from './modules'
//include {annotateConsensus} from './modules'
include {annotatePilonConsensus} from './modules'

include {annotateVCFs} from './modules'

//include {mlst} from './modules'
//include {stats} from './modules'

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    input_read_ch = Channel
        .fromFilePairs("${params.INPUT}*_R{1,2}*.gz")
        .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT} ending with .gz" }
        .map { it -> [it[0], it[1][0], it[1][1]]}

    trimReads (
        input_read_ch,
        ADAPTERS
    )

    filterTp (
        //trimReads.out,
        trimReads.out[0],
        REF_FASTAS
    )

    mapUnmatchedReads (
        filterTp.out[1],
        REF_FASTAS_MASKED
    )

    moreFiltering (
        mapUnmatchedReads.out[1],
        REPEAT_FILTER_FASTA
    )

    mapReads (
        filterTp.out[0].groupTuple()
            .join(
                mapUnmatchedReads.out[0].groupTuple()
        ),

        NC_021508,
        NC_021508_1,
        NC_021508_2,
        NC_021508_3,
        NC_021508_4,
        NC_021508_5,
        NC_021508_6

    )

    samToBam (
        mapReads.out,
    )
    removeDuplicates (
        samToBam.out[0]
    )
    callVariants (
        removeDuplicates.out[2],
        NC_021508
    )

    if(params.SKIP_DENOVO == false){    

    deNovoAssembly (
        moreFiltering.out[0]
    )
    mergeAssemblyMapping (
        //deNovoAssembly.out[0],
        //removeDuplicates.out[1],

        deNovoAssembly.out[0].groupTuple()
            .join(
                removeDuplicates.out[1].groupTuple()
        ),

        file(NC_021508),
        file(NC_021508_BWA1),
        file(NC_021508_BWA2),
        file(NC_021508_BWA3),
        file(NC_021508_BWA4),
        file(NC_021508_BWA5),
        file(TP_MAKE_SEQ)
    )
    remapReads (
        mergeAssemblyMapping.out[0].groupTuple()
            .join(
                removeDuplicates.out[5].groupTuple()
        )
    )
    
    pilonPolishing (
        remapReads.out[1]
    )
    remapPilon (
        pilonPolishing.out[0].groupTuple()
            .join(
                removeDuplicates.out[6].groupTuple()
        )
    )
    generatePilonConsensus (
        remapPilon.out[0].groupTuple(
            ).join(
                pilonPolishing.out[1].groupTuple()
            ).join(
                removeDuplicates.out[3].groupTuple()
        ),
        TP_GENERATE_CONSENSUS
    )

    annotatePilonConsensus (
        generatePilonConsensus.out,
        REF_GB,
    )

    

    annotateVCFs (

        trimReads.out[0], 
        removeDuplicates.out[2],
        samToBam.out[0],
        REF_GB,
        NC_021508,
        removeDuplicates.out[4]
    )

    /*
    annotateVCFs (
        trimReads.out[0], 
        REF_GB
    )
    */

    //mlst (
    //    generatePilonConsensus.out
    //)

    }
    
    if(params.SKIP_DENOVO == true){ 

    remapReads_2 (
        samToBam.out[1],
        removeDuplicates.out[4]
    )
    pilonPolishing (
        remapReads_2.out[1]
    )
    remapPilon (
        pilonPolishing.out[0].groupTuple()
            .join(
                removeDuplicates.out[6].groupTuple()
        )
    )

    generatePilonConsensus (
        remapPilon.out[0].groupTuple(
            ).join(
                pilonPolishing.out[1].groupTuple()
            ).join(
                removeDuplicates.out[3].groupTuple()
        ),
        TP_GENERATE_CONSENSUS
    )

    annotatePilonConsensus (
        generatePilonConsensus.out,
        REF_GB
    )

    annotateVCFs (

        trimReads.out[0], 
        removeDuplicates.out[2],
        samToBam.out[0],
        REF_GB,
        NC_021508,
        removeDuplicates.out[4]
    )

    

    /*

    annotateVCFs (
        file(CONVERT_TO_ANNOVAR),
        file(GFF3_TO_GENE_PRED),
        file(RETRIEVE_SEQ_FROM_FASTA),
        file(ANNOTATE_VARIATION),
        file(GFF),

        callVariants.out[0],

        file(NC_021508)
    )
    */

    //mlst (
    //    generatePilonConsensus.out
    //)

    }
    
    // stats (
    //     input_read_ch,
    //     //annotatePilonConsensus.out[1]
    //
    // )

}