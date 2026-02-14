# K-mers GWAS Pipeline - Lab Guide

## Overview

This pipeline performs k-mer based genome-wide association studies (GWAS) using the kmersGWAS method (Voichek & Weigel, 2020). It takes raw sequencing reads through to significant k-mer associations with your phenotype.

**Pipeline Steps:**
1. Merge sequencing lanes per sample
2. Quality control with fastp
3. K-mer counting with KMC (canonical and non-canonical)
4. Combine strand information
5. List k-mers across all samples with filtering
6. Build k-mers presence/absence table
7. Calculate kinship matrix
8. Run GWAS analysis (optional)

**Expected Runtime (based on actual 161 Malus samples, 10-20X WGS):**
- **Full pipeline through kinship (161 samples): 6h 33m**
  - Total processes: 809 (161 samples × 5 processes + 4 global processes)
  - CPU hours: 243.3
  - Parallelization efficiency: ~37× speedup
- GWAS only (pre-computed tables): 12-48 hours (depending on permutations)
- Test run (3 samples): <30 minutes

**Note:** These times are significantly faster than original manual estimates due to modern hardware, parallelization, and efficient resource allocation.

---

## Setup on Compute2 Cluster

### 1. Load Required Modules

**Every time you run the pipeline, you must load these modules first:**

```bash
ml ris
ml openjdk/17.0.11_9
ml apptainer
```

### 2. Clone the Pipeline Repository

```bash
cd /storage1/fs1/<PI_NAME>/Active/<YOUR_PROJECT>
git clone <repository_url> kmers_gwas_pipeline
cd kmers_gwas_pipeline
```

### 3. Verify Files

Make sure you have:
- `kmers_gwas.nf` - Main pipeline file
- `nextflow.config` - Configuration file

---

## Input Files Required

### Sample Sheet (TSV format)

Create a tab-separated file with sequencing read information:

**Format for Paired-End (PE) reads:** `samples.tsv`
```
Sample	ID	Lane	Read1	Read2
M1	M1_1	1	/path/to/M1_L1_R1.fastq.gz	/path/to/M1_L1_R2.fastq.gz
M1	M1_2	2	/path/to/M1_L2_R1.fastq.gz	/path/to/M1_L2_R2.fastq.gz
M2	M2_1	1	/path/to/M2_L1_R1.fastq.gz	/path/to/M2_L1_R2.fastq.gz
M2	M2_2	2	/path/to/M2_L2_R1.fastq.gz	/path/to/M2_L2_R2.fastq.gz
```

**Format for Single-End (SE) reads:** `samples_se.tsv`
```
Sample	ID	Lane	Read1
M1	M1_1	1	/path/to/M1_L1.fastq.gz
M1	M1_2	2	/path/to/M1_L2.fastq.gz
M2	M2_1	1	/path/to/M2_L1.fastq.gz
M2	M2_2	2	/path/to/M2_L2.fastq.gz
```

**Columns:**
- `Sample` - Sample identifier (biological replicate)
- `ID` - Unique ID for each sequencing run/lane
- `Lane` - Lane number
- `Read1` - Full path to reads (R1 for PE, or single file for SE)
- `Read2` - Full path to reverse reads (PE only, omit this column for SE)

**Note:** 
- The pipeline automatically detects PE vs SE based on whether the `Read2` column exists
- Samples with multiple lanes will be automatically merged before k-mer counting

### Phenotype File (TSV format)

Create a tab-separated file with phenotype values:

**Format:** `phenotype.tsv`
```
accession_id	phenotype_value
M1	12.5
M2	15.3
M3	11.8
```

**Columns:**
- `accession_id` - Must match `Sample` column in sample sheet
- `phenotype_value` - Numeric phenotype value

---

## Running the Pipeline

### First Time: Full Pipeline

Run the complete pipeline from raw reads to k-mer table and kinship matrix:

```bash
# Load modules
ml ris openjdk/17.0.11_9 apptainer

# Run pipeline
nextflow run kmers_gwas.nf \
  --samples samples.tsv \
  --outdir results_full \
  --mac 5 \
  --maf 0.05 \
  --kmer_length 31
```

**This creates:**
- K-mer counts for each sample
- Filtered k-mer list across all samples
- K-mers presence/absence table
- Kinship matrix

**Expected time:** **6h 33m for 161 samples** (10-20X WGS coverage, default settings)

### Running GWAS (First Time)

Add GWAS analysis to the full pipeline:

```bash
nextflow run kmers_gwas.nf \
  --samples samples.tsv \
  --phenotype phenotype.tsv \
  --run_gwas true \
  --outdir results_full \
  --mac 5 \
  --maf 0.05 \
  --gwas_permutations 100
```

### Subsequent Runs: Using Pre-computed Tables

**If you already have k-mer tables and kinship matrix**, skip the expensive k-mer counting:

```bash
nextflow run kmers_gwas.nf \
  --precomputed_kmers_table results_full/06_kmers_table/kmers_table \
  --precomputed_kinship results_full/07_kinship/kmers_table.kinship \
  --phenotype phenotype_v2.tsv \
  --run_gwas true \
  --outdir results_phenotype_v2
```

**Expected time:** 12-48 hours (GWAS only, no k-mer counting)

### Testing Multiple Phenotypes

Once you have computed tables, test different phenotypes quickly:

```bash
# Phenotype 1
nextflow run kmers_gwas.nf \
  --precomputed_kmers_table results_full/06_kmers_table/kmers_table \
  --precomputed_kinship results_full/07_kinship/kmers_table.kinship \
  --phenotype flowering_time.tsv \
  --run_gwas true \
  --outdir results_flowering

# Phenotype 2
nextflow run kmers_gwas.nf \
  --precomputed_kmers_table results_full/06_kmers_table/kmers_table \
  --precomputed_kinship results_full/07_kinship/kmers_table.kinship \
  --phenotype disease_resistance.tsv \
  --run_gwas true \
  --outdir results_disease
```

### Recalculating Kinship with Different Parameters

If you want to test different MAF thresholds for the kinship matrix:

```bash
nextflow run kmers_gwas.nf \
  --precomputed_kmers_table results_full/06_kmers_table/kmers_table \
  --maf 0.01 \
  --phenotype phenotype.tsv \
  --run_gwas true \
  --outdir results_maf01
```

---

## Key Parameters

### K-mer Counting Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--kmer_length` | 31 | Length of k-mers (max 31) |
| `--kmc_ci` | 2 | Minimum k-mer count threshold for canonical run |
| `--mac` | 5 | Minor allele count - k-mer must appear in N samples |
| `--min_strand_pct` | 0.2 | Minimum % of samples where k-mer appears in each strand |

**Recommendations:**
- `--mac`: Use 2-3 for small test datasets (<10 samples), 5+ for real analysis (>100 samples)
- `--kmer_length`: 31 is standard and recommended
- Adjust `--mac` based on your expected minor allele frequency and sample size

### Kinship Matrix Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--maf` | 0.05 | Minor allele frequency for kinship calculation |

**Notes:**
- Lower MAF (e.g., 0.01) = more variants used, slower runtime, potentially better for closely related samples
- Higher MAF (e.g., 0.05) = fewer variants, faster runtime, good for diverse populations

### GWAS Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_gwas` | false | Set to `true` to run GWAS |
| `--phenotype` | null | Path to phenotype file (required for GWAS) |
| `--gwas_permutations` | 100 | Number of permutations for threshold calculation |
| `--gwas_maf` | 0.05 | MAF filter for GWAS |
| `--gwas_mac` | 5 | MAC filter for GWAS |
| `--gwas_n_kmers` | 10000000 | Number of k-mers in first filtering step |
| `--gwas_pattern_counter` | false | Count unique k-mer patterns (slower) |

### Pre-computed Data Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--precomputed_kmers_table` | null | Path to existing kmers_table (without .table extension) |
| `--precomputed_kinship` | null | Path to existing kinship matrix file |

---

## Output Files

### Directory Structure

```
results/
├── 00_merged/                    # Merged FASTQ files per sample
├── 01_fastp/                     # QC reports
├── 02_kmc/                       # KMC k-mer counts (large, auto-cleaned)
├── 03_kmers_strand/              # K-mers with strand info per sample
├── 04_qc/
│   └── kmc_qc_summary.tsv       # QC statistics for all samples
├── 05_kmers_list/
│   ├── kmers_to_use             # Filtered k-mer list (binary)
│   └── kmers_to_use.shareness   # K-mer sharing statistics
├── 06_kmers_table/
│   ├── kmers_table.table        # Presence/absence matrix (binary)
│   └── kmers_table.names        # Sample names in order
├── 07_kinship/
│   └── kmers_table.kinship      # Kinship matrix
└── 08_gwas/
    └── kmers/
        ├── output/
        │   └── phenotype_value.assoc.txt.gz  # Full GWAS results
        ├── pass_threshold_5per     # Significant k-mers (5% FDR)
        ├── pass_threshold_10per    # Significant k-mers (10% FDR)
        ├── threshold_5per          # -log10(p) threshold for 5%
        └── threshold_10per         # -log10(p) threshold for 10%
```

### Important Output Files

**QC Summary:** `04_qc/kmc_qc_summary.tsv`
- Check this to identify problematic samples
- Look for outliers in k-mer counts vs read counts

**K-mers Table:** `06_kmers_table/kmers_table.table`
- Binary presence/absence matrix
- Reuse this for multiple GWAS runs

**Kinship Matrix:** `07_kinship/kmers_table.kinship`
- Tab-separated N×N matrix of genetic relatedness
- Reuse this for multiple phenotypes

**GWAS Results:** `08_gwas/kmers/output/phenotype_value.assoc.txt.gz`
- Full association results for all k-mers
- Columns: chromosome, position, alleles, frequencies, test statistics, p-values

**Significant K-mers:** `08_gwas/kmers/pass_threshold_5per`
- K-mers passing 5% family-wise error rate threshold
- These are your candidate k-mers associated with the phenotype

---

## Monitoring and Troubleshooting

### Monitoring Running Jobs

Check Nextflow progress:
```bash
# View current status
tail -f .nextflow.log

# Check SLURM queue
squeue -u $USER
```

### Resuming Failed Jobs

Nextflow automatically tracks completed tasks. If the pipeline fails, resume from where it stopped:

```bash
nextflow run kmers_gwas.nf \
  --samples samples.tsv \
  --outdir results_full \
  -resume
```

**The `-resume` flag will:**
- Skip completed processes
- Rerun only failed or pending processes
- Save time by not repeating successful steps

### Checking Disk Space

K-mer counting generates large intermediate files. Check usage:

```bash
# Check work directory size
du -sh work/

# Check output directory size
du -sh results_full/
```

**Disk space requirements (approximate):**
- Per sample during KMC: ~28 GB (temporary)
- Final outputs per sample: ~2 GB
- Work directory cleaned automatically with nf-boost

### Common Issues

**Issue:** "No k-mers pass filtering"
- **Cause:** MAC threshold too high for sample size
- **Solution:** Lower `--mac` parameter (e.g., use `--mac 2` for <10 samples)

**Issue:** "Out of memory error" in KMC
- **Cause:** Not enough memory allocated
- **Solution:** KMC uses 64GB by default. For very high coverage, may need manual adjustment

**Issue:** "Kinship calculation taking too long"
- **Cause:** Too many k-mers (low MAF)
- **Solution:** Increase `--maf` to 0.05 or 0.1 for faster computation

**Issue:** "GWAS process fails"
- **Cause:** Phenotype file format issue or missing samples
- **Solution:** Check phenotype file has correct header and sample names match exactly

---

## Example Workflows

### Workflow 1: First Analysis (Complete Pipeline)

```bash
# Load modules
ml ris openjdk/17.0.11_9 apptainer

# Run full pipeline with GWAS
nextflow run kmers_gwas.nf \
  --samples samples.tsv \
  --phenotype flowering_time.tsv \
  --run_gwas true \
  --outdir results_flowering \
  --mac 5 \
  --maf 0.05 \
  --gwas_permutations 100

# Check results
less results_flowering/08_gwas/kmers/pass_threshold_5per
```

### Workflow 2: Testing Different Phenotypes

```bash
# First phenotype (uses pre-computed tables)
nextflow run kmers_gwas.nf \
  --precomputed_kmers_table results_flowering/06_kmers_table/kmers_table \
  --precomputed_kinship results_flowering/07_kinship/kmers_table.kinship \
  --phenotype disease_resistance.tsv \
  --run_gwas true \
  --outdir results_disease

# Second phenotype
nextflow run kmers_gwas.nf \
  --precomputed_kmers_table results_flowering/06_kmers_table/kmers_table \
  --precomputed_kinship results_flowering/07_kinship/kmers_table.kinship \
  --phenotype height.tsv \
  --run_gwas true \
  --outdir results_height
```

### Workflow 3: Small Test Run

Test with a subset of samples first:

```bash
# Create test sample sheet (3 samples)
head -n 4 samples.tsv > samples_test.tsv

# Run test
nextflow run kmers_gwas.nf \
  --samples samples_test.tsv \
  --outdir test_run \
  --mac 2 \
  --maf 0.05

# If successful, run on full dataset
nextflow run kmers_gwas.nf \
  --samples samples.tsv \
  --outdir results_full \
  --mac 5
```

---

## Best Practices

### 1. Start with QC

Always check the QC summary after k-mer counting:

```bash
# View QC summary
column -t results_full/04_qc/kmc_qc_summary.tsv | less -S

# Look for:
# - Correlation between reads and k-mers (should be linear)
# - Outlier samples (very high/low k-mer counts)
# - Flag distribution (flag_0 should be 0)
```

### 2. Use Appropriate MAC

Choose MAC based on sample size and expected allele frequency:

- **10-50 samples:** `--mac 2` or `--mac 3`
- **50-100 samples:** `--mac 3` or `--mac 5`
- **100+ samples:** `--mac 5` or higher

### 3. Preserve Intermediate Results

Keep your k-mer tables and kinship matrices:

```bash
# After first successful run, back up key files
mkdir -p backup
cp results_full/06_kmers_table/kmers_table.* backup/
cp results_full/07_kinship/kmers_table.kinship backup/
```

### 4. Document Your Parameters

Keep a log of parameters used:

```bash
# Save your exact command
echo "nextflow run kmers_gwas.nf ..." > analysis_log.txt
date >> analysis_log.txt
```

### 5. Clean Up Work Directory

After successful completion, clean up work directory to save space:

```bash
# Nextflow automatically cleans with nf-boost
# But you can also manually clean after confirming results
nextflow clean -f
```

---

## Getting Help

### Check Logs

```bash
# Main nextflow log
tail -100 .nextflow.log

# Process-specific logs (in work directory)
cat work/<hash>/<hash>/.command.log
cat work/<hash>/<hash>/.command.err
```

### Useful Commands

```bash
# List all processes that ran
nextflow log

# Get detailed run information
nextflow log <run_name> -f process,status,duration,realtime

# Clean specific run
nextflow clean <run_name> -f
```

---

## Performance Analysis (Real-world Data)

**Test Case:** 161 wild Malus samples, 10-20X WGS coverage, paired-end Illumina

**Results:**
- **Total runtime:** 6h 33m 19s
- **Total processes:** 809 (161 samples × 5 per-sample processes + 4 global processes)
- **CPU hours consumed:** 243.3
- **Parallelization efficiency:** 37× speedup vs sequential execution
- **Peak disk usage:** Managed automatically by nf-boost cleanup

**Breakdown by process type:**
1. **MERGE_LANES** (161 jobs): Fast I/O operation
2. **FASTP_QC** (161 jobs): 8 CPUs each, parallelized well
3. **KMC_CANON** (161 jobs): 16 CPUs each, most compute-intensive per-sample
4. **KMC_ALL** (161 jobs): 16 CPUs each, similar to CANON
5. **COMBINE_STRAND_INFO** (161 jobs): Light processing, quick
6. **SUMMARIZE_QC** (1 job): Aggregates all sample stats
7. **LIST_KMERS_MULTI_SAMPLES** (1 job): Filters k-mers across population
8. **BUILD_KMERS_TABLE** (1 job): Creates presence/absence matrix
9. **CALCULATE_KINSHIP** (1 job): 16 CPUs, completed in <1 hour with MAF 0.05

**Why so fast compared to manual estimates?**
- Modern hardware (faster CPUs, SSDs)
- Efficient parallelization across samples
- Proper resource allocation (16 CPUs for KMC)
- Automatic cleanup preventing I/O bottlenecks
- Good cluster scheduling (queueSize=50)

**Scalability:** The pipeline scales linearly with sample count. Expect ~2.4 minutes per sample for k-mer counting steps.

---

## Citations

If you use this pipeline, please cite:

**K-mers GWAS Method:**
Voichek, Y., & Weigel, D. (2020). Identifying genetic variants underlying phenotypic variation in plants without complete genomes. *Nature Genetics*, 52(5), 534-540.

**Software:**
- KMC: Kokot, M., Długosz, M., & Deorowicz, S. (2017). *Bioinformatics*, 33(17), 2759-2761.
- GEMMA: Zhou, X., & Stephens, M. (2012). *Nature Genetics*, 44(7), 821-824.
- Nextflow: Di Tommaso, P., et al. (2017). *Nature Biotechnology*, 35(4), 316-319.

---

## Pipeline Versions

**Current Version:** 1.0

**Container:** `docker://bmansfeld/kmersgwas:v0.2-beta`

**Last Updated:** February 2026

---

## Contact

For questions or issues:
- Check the troubleshooting section above
- Review pipeline logs in `.nextflow.log`
- Contact your bioinformatics support team
