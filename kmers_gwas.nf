#!/usr/bin/env nextflow

/*
 K-MERS GWAS PIPELINE - CHUNK 2 (with resource allocations)
 Handles multiple lanes per sample + k-mer counting
*/

nextflow.enable.dsl=2

// Parameters
params.samples = null          // TSV: Sample, ID, Lane, Read1, Read2
params.outdir = "results"
params.kmer_length = 31
params.kmc_ci = 2              // KMC count threshold for canonical run
params.mac = 2           // Minor allele count - k-mer must appear in at least 5 samples
params.min_strand_pct = 0.2  // 20% - k-mer must appear in each strand form in 20% of samples where it's found
params.maf = 0.05        // Minor allele frequency for kinship calculation
params.help = false
params.phenotype = null              // Phenotype file (required for GWAS)
params.gwas_permutations = 100       // Number of permutations for threshold
params.gwas_maf = 0.05               // MAF for GWAS
params.gwas_mac = 5                  // MAC for GWAS
params.gwas_n_kmers = 10000000       // Number of k-mers to filter in first step
params.gwas_pattern_counter = false  // Count unique k-mer patterns
params.gwas_remove_intermediates = true  // Clean up intermediate files
params.gemma_path = "gemma"          // Path to GEMMA executable
params.run_gwas = false              // Set to true to run GWAS

// if you already ran everything and now just want to run the GWAS
params.precomputed_kmers_table = null   // Path to existing kmers_table (without .table extension)
params.precomputed_kinship = null       // Path to existing kinship matrix
params.skip_kmer_counting = false       // Skip all k-mer counting steps

params.kmersgwas_container = 'docker://bmansfeld/kmersgwas:v0.2-beta'

if (params.help) {
    println """
    Usage: nextflow run kmers_gwas.nf --samples samples.tsv --outdir results
    
    Required:
      --samples       TSV with columns: Sample, ID, Lane, Read1, Read2
    
    Optional:
      --outdir        Output directory (default: results)
      --kmer_length   K-mer length (default: 31, max: 31)
      --kmc_ci        KMC count threshold for canonical (default: 2)
    """.stripIndent()
    exit 0
}

if (!params.samples) error "ERROR: --samples is required"

// Process 1: Merge lanes per sample
process MERGE_LANES {
    tag "${sample}"
    publishDir "${params.outdir}/00_merged/${sample}", mode: 'copy', pattern: "*.txt"
    
    cpus 1
    memory 4.GB
    time 4.h
    
    input:
    tuple val(sample), path(read1_files), path(read2_files)
    
    output:
    tuple val(sample), path("${sample}_R1.merged.fastq.gz"), path("${sample}_R2.merged.fastq.gz"), emit: reads
    path "${sample}_merge_info.txt"
    
    script:
    """
    # Merge R1 files
    cat ${read1_files.join(' ')} > ${sample}_R1.merged.fastq.gz
    
    # Merge R2 files
    cat ${read2_files.join(' ')} > ${sample}_R2.merged.fastq.gz
    
    # Record what was merged
    echo "Sample: ${sample}" > ${sample}_merge_info.txt
    echo "R1 files merged:" >> ${sample}_merge_info.txt
    echo "${read1_files.join('\n')}" >> ${sample}_merge_info.txt
    echo "" >> ${sample}_merge_info.txt
    echo "R2 files merged:" >> ${sample}_merge_info.txt
    echo "${read2_files.join('\n')}" >> ${sample}_merge_info.txt
    """
}

// Process 2: QC with fastp
process FASTP_QC {
    tag "${sample}"
    publishDir "${params.outdir}/01_fastp/${sample}", mode: 'copy'
    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'
    
    cpus 8
    memory 16.GB
    time 8.h
    
    input:
    tuple val(sample), path(read1), path(read2)
    
    output:
    tuple val(sample), path("${sample}_R1.trimmed.fastq.gz"), path("${sample}_R2.trimmed.fastq.gz"), emit: reads
    path "${sample}_fastp.json"
    path "${sample}_fastp.html"
    
    script:
    """
    fastp \
        -i ${read1} \
        -I ${read2} \
        -o ${sample}_R1.trimmed.fastq.gz \
        -O ${sample}_R2.trimmed.fastq.gz \
        --json ${sample}_fastp.json \
        --html ${sample}_fastp.html \
        --thread ${task.cpus} \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50
    """
}

// Process 3: KMC canonical count (with canonization)
process KMC_CANON {
    tag "${sample}"
    publishDir "${params.outdir}/02_kmc/${sample}", mode: 'copy', pattern: "*.{1,2}"
    container 'quay.io/biocontainers/kmc:3.2.1--hf1761c0_2'
      
    cpus 16
    memory 64.GB
    time 24.h
    
    input:
    tuple val(sample), path(read1), path(read2)
    
    output:
    tuple val(sample), path("${sample}_canon.kmc_pre"), path("${sample}_canon.kmc_suf"), emit: kmc_db
    path "kmc_canon.1"
    path "kmc_canon.2"
    
    script:
    """
    # Create input file list
    echo "${read1}" > input_files.txt
    echo "${read2}" >> input_files.txt
    
    # Run KMC with canonization
    kmc -t${task.cpus} \
        -k${params.kmer_length} \
        -ci${params.kmc_ci} \
        @input_files.txt \
        ${sample}_canon \
        . \
        1> kmc_canon.1 2> kmc_canon.2
    """
}

// Process 4: KMC all k-mers count (no canonization)
process KMC_ALL {
    tag "${sample}"
    publishDir "${params.outdir}/02_kmc/${sample}", mode: 'copy', pattern: "*.{1,2}"
    container 'quay.io/biocontainers/kmc:3.2.1--hf1761c0_2'
       
    cpus 16
    memory 64.GB
    time 24.h
    
    input:
    tuple val(sample), path(read1), path(read2)
    
    output:
    tuple val(sample), path("${sample}_all.kmc_pre"), path("${sample}_all.kmc_suf"), emit: kmc_db
    path "kmc_all.1"
    path "kmc_all.2"
    
    script:
    """
    # Create input file list
    echo "${read1}" > input_files.txt
    echo "${read2}" >> input_files.txt
    
    # Run KMC without canonization
    kmc -t${task.cpus} \
        -k${params.kmer_length} \
        -ci0 \
        -b \
        @input_files.txt \
        ${sample}_all \
        . \
        1> kmc_all.1 2> kmc_all.2
    """
}

// Process 5: Combine strand information (UPDATED - now captures log)
process COMBINE_STRAND_INFO {
    tag "${sample}"
    publishDir "${params.outdir}/03_kmers_strand/${sample}", mode: 'copy'
    container params.kmersgwas_container
    
       
    cpus 2
    memory 8.GB
    time 4.h
    
    input:
    tuple val(sample), 
          path(canon_pre), path(canon_suf),
          path(all_pre), path(all_suf)
    
    output:
    tuple val(sample), path("${sample}_kmers_with_strand"), emit: kmers
    path "${sample}_strand_info.log", emit: log
    
    script:
    def canon_prefix = canon_pre.baseName.replaceAll(/\.kmc_pre$/, '')
    def all_prefix = all_pre.baseName.replaceAll(/\.kmc_pre$/, '')
    """
    kmers_add_strand_information \
        -c ${canon_prefix} \
        -n ${all_prefix} \
        -k ${params.kmer_length} \
        -o ${sample}_kmers_with_strand \
        2>&1 | tee ${sample}_strand_info.log
    """
}

// Process 6: Summarize QC stats from ALL samples into ONE file
process SUMMARIZE_QC {
    publishDir "${params.outdir}/04_qc", mode: 'copy'
    
    cpus 1
    memory 4.GB
    time 1.h
    
    input:
    path strand_logs  // This receives ALL logs via .collect()
    
    output:
    path "kmc_qc_summary.tsv"
    
    script:
    """
    #!/bin/bash
    echo -e "Sample\\tCanonized_kmers\\tNonCanon_kmers\\tNonCanon_found\\tFlag_1\\tFlag_2\\tFlag_3" > kmc_qc_summary.tsv
    
    for log in *_strand_info.log; do
        sample=\$(basename \$log _strand_info.log)
        canon=\$(grep "Canonized kmers:" \$log | awk '{print \$NF}')
        noncanon=\$(grep "Non-canon kmers:" \$log | awk '{print \$NF}')
        found=\$(grep "Non-canon kmers found:" \$log | awk '{print \$NF}')
        flag1=\$(grep "flag.*1" \$log | awk '{print \$NF}')
        flag2=\$(grep "flag.*2" \$log | awk '{print \$NF}')
        flag3=\$(grep "flag.*3" \$log | awk '{print \$NF}')
        
        echo -e "\$sample\\t\$canon\\t\$noncanon\\t\$found\\t\$flag1\\t\$flag2\\t\$flag3" >> kmc_qc_summary.tsv
    done
    """
}

// Process 7: List all k-mers found across multiple samples
process LIST_KMERS_MULTI_SAMPLES {
    publishDir "${params.outdir}/05_kmers_list", mode: 'copy'
    container params.kmersgwas_container
    
    cpus 4
    memory 32.GB
    time 12.h
    
    input:
    path kmers_files  // All kmers_with_strand files collected
    
    output:
    path "kmers_to_use", emit: kmers_list
    path "kmers_to_use.shareness"
    path "kmers_to_use.stats.both"
    path "kmers_to_use.stats.only_canonical"
    path "kmers_to_use.stats.only_non_canonical"
    path "kmers_list_paths.txt"
    
    script:
    """
    # Create list of k-mers files with sample names
    # Format: /path/to/file<TAB>sample_name
    for kmers_file in *_kmers_with_strand; do
        sample=\$(basename \$kmers_file _kmers_with_strand)
        echo -e "\${PWD}/\${kmers_file}\\t\${sample}" >> kmers_list_paths.txt
    done
    
    # List k-mers found in multiple samples with filtering
    list_kmers_found_in_multiple_samples \
        -l kmers_list_paths.txt \
        -k ${params.kmer_length} \
        --mac ${params.mac} \
        -p ${params.min_strand_pct} \
        -o kmers_to_use
    """
}


// Process 8: Build k-mers presence/absence table
process BUILD_KMERS_TABLE {
    publishDir "${params.outdir}/06_kmers_table", mode: 'copy'
    container params.kmersgwas_container
    
    cpus 8
    memory 64.GB
    time 24.h
    
    input:
    path kmers_files      // All kmers_with_strand files
    path kmers_to_use     // Filtered k-mer list from LIST_KMERS_MULTI_SAMPLES
    
    output:
    path "kmers_table.table", emit: table
    path "kmers_table.names", emit: names
    path "kmers_list_paths.txt"
    
    script:
    """
    # Create list of k-mers files with sample names
    # Format: /path/to/file<TAB>sample_name
    for kmers_file in *_kmers_with_strand; do
        sample=\$(basename \$kmers_file _kmers_with_strand)
        echo -e "\${PWD}/\${kmers_file}\\t\${sample}" >> kmers_list_paths.txt
    done
    
    # Build the k-mers presence/absence table
    build_kmers_table \
        -l kmers_list_paths.txt \
        -k ${params.kmer_length} \
        -a ${kmers_to_use} \
        -o kmers_table
    """
}

// Process 9: Calculate kinship matrix
process CALCULATE_KINSHIP {
    publishDir "${params.outdir}/07_kinship", mode: 'copy'
    container params.kmersgwas_container
    
    cpus 16
    memory 128.GB
    time '5.d'  // 5 days - this can be slow for large datasets
    
    input:
    path kmers_table
    path kmers_names
    
    output:
    path "kmers_table.kinship", emit: kinship
    
    script:
    """
    # emma_kinship_kmers outputs to stdout, redirect to file
    emma_kinship_kmers \
        -t kmers_table \
        -k ${params.kmer_length} \
        --maf ${params.maf} \
        > kmers_table.kinship
    """
}

// Process 10: Run k-mers GWAS
process RUN_KMERS_GWAS {
    publishDir "${params.outdir}/08_gwas", mode: 'copy'
    container params.kmersgwas_container
    
    cpus 16
    memory 128.GB
    time '7.d'
    
    input:
    path kmers_table
    path kmers_names
    path kinship
    path phenotype
    
    output:
    path "gwas_output/kmers/output/phenotype_value.assoc.txt.gz", emit: results
    path "gwas_output/kmers/pass_threshold_5per", emit: sig_5per, optional: true
    path "gwas_output/kmers/pass_threshold_10per", emit: sig_10per, optional: true
    path "gwas_output/kmers/threshold_5per", emit: thresh_5per
    path "gwas_output/kmers/threshold_10per", emit: thresh_10per
    path "gwas_output/log_file"
    path "gwas_output/kmers/output/*", optional: true
    
    script:
    def pattern_flag = params.gwas_pattern_counter ? "--pattern_counter" : ""
    def remove_flag = params.gwas_remove_intermediates ? "" : "--dont_remove_intermediates"
    """
    # Ensure kinship matrix has the expected name relative to kmers_table
    cp ${kinship} kmers_table.kinship
    
    # Get the kmersGWAS installation from the container
    KMERS_GWAS_BASE=\$(dirname \$(which kmers_add_strand_information))/../
    
    # Run the GWAS pipeline
    python2.7 \${KMERS_GWAS_BASE}/kmers_gwas.py \
        --pheno ${phenotype} \
        --kmers_table kmers_table \
        -l ${params.kmer_length} \
        -p ${task.cpus} \
        --outdir gwas_output \
        --maf ${params.gwas_maf} \
        --mac ${params.gwas_mac} \
        --permutations ${params.gwas_permutations} \
        --kmers_table kmers_table \
        -k ${params.gwas_n_kmers} \
        --gemma_path ${params.gemma_path} \
        ${pattern_flag} \
        ${remove_flag}
    
    # Ensure output files exist (even if empty)
    touch gwas_output/kmers/threshold_5per
    touch gwas_output/kmers/threshold_10per
    """
}

workflow {
    
    // ==============================================================================
    // BRANCH 1: Full pipeline from raw reads
    // ==============================================================================
    if (!params.skip_kmer_counting && params.precomputed_kmers_table == null) {
        
        // Read sample sheet and group by Sample
        Channel
            .fromPath(params.samples)
            .splitCsv(header: true, sep: '\t')
            .map { row -> 
                tuple(row.Sample, file(row.Read1), file(row.Read2))
            }
            .groupTuple(by: 0)
            .map { sample, read1s, read2s ->
                tuple(sample, read1s, read2s)
            }
            .set { samples_grouped }
        
        // Merge lanes per sample
        MERGE_LANES(samples_grouped)
        
        // QC merged reads
        FASTP_QC(MERGE_LANES.out.reads)
        
        // Run both KMC processes in parallel
        KMC_CANON(FASTP_QC.out.reads)
        KMC_ALL(FASTP_QC.out.reads)
        
        // Combine the outputs from both KMC runs
        KMC_CANON.out.kmc_db
            .join(KMC_ALL.out.kmc_db)
            .map { sample, canon_pre, canon_suf, all_pre, all_suf ->
                tuple(sample, canon_pre, canon_suf, all_pre, all_suf)
            }
            .set { combined_kmc }
        
        // Combine strand information
        COMBINE_STRAND_INFO(combined_kmc)
        
        // Collect QC stats
        SUMMARIZE_QC(COMBINE_STRAND_INFO.out.log.collect())
        
        // List all k-mers across all samples
        LIST_KMERS_MULTI_SAMPLES(
            COMBINE_STRAND_INFO.out.kmers.map { sample, kmers -> kmers }.collect()
        )
        
        // Build k-mers table
        BUILD_KMERS_TABLE(
            COMBINE_STRAND_INFO.out.kmers.map { sample, kmers -> kmers }.collect(),
            LIST_KMERS_MULTI_SAMPLES.out.kmers_list
        )
        
        // Set channels for downstream
        kmers_table_ch = BUILD_KMERS_TABLE.out.table
        kmers_names_ch = BUILD_KMERS_TABLE.out.names
        
        // Calculate kinship if not provided
        if (params.precomputed_kinship == null) {
            CALCULATE_KINSHIP(kmers_table_ch, kmers_names_ch)
            kinship_ch = CALCULATE_KINSHIP.out.kinship
        } else {
            kinship_ch = Channel.fromPath(params.precomputed_kinship)
        }
    }
    
    // ==============================================================================
    // BRANCH 2: Use pre-computed k-mers table
    // ==============================================================================
    else if (params.precomputed_kmers_table != null) {
        
        // Load pre-computed files
        kmers_table_ch = Channel.fromPath("${params.precomputed_kmers_table}.table")
        kmers_names_ch = Channel.fromPath("${params.precomputed_kmers_table}.names")
        
        // Use pre-computed kinship or calculate new one
        if (params.precomputed_kinship != null) {
            kinship_ch = Channel.fromPath(params.precomputed_kinship)
        } else {
            CALCULATE_KINSHIP(kmers_table_ch, kmers_names_ch)
            kinship_ch = CALCULATE_KINSHIP.out.kinship
        }
    }
    
    // ==============================================================================
    // GWAS STEP: Run if phenotype provided
    // ==============================================================================
    if (params.run_gwas && params.phenotype != null) {
        RUN_KMERS_GWAS(
            kmers_table_ch,
            kmers_names_ch,
            kinship_ch,
            file(params.phenotype)
        )
    }
}

workflow.onComplete {
    println "==================================="
    println "PIPELINE COMPLETE"
    println "Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    println ""
    println "Outputs:"
    println "  - K-mers table: ${params.outdir}/06_kmers_table/kmers_table.table"
    println "  - Sample names: ${params.outdir}/06_kmers_table/kmers_table.names"
    println "  - Kinship matrix: ${params.outdir}/07_kinship/kmers_table.kinship"
    if (params.run_gwas && params.phenotype != null) {
        println "  - GWAS results: ${params.outdir}/08_gwas/kmers/output/"
        println "  - Significant k-mers (5%): ${params.outdir}/08_gwas/kmers/pass_threshold_5per"
        println "  - Significant k-mers (10%): ${params.outdir}/08_gwas/kmers/pass_threshold_10per"
    }
    println "==================================="
}
