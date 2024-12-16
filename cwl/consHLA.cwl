cwlVersion: v1.2
class: Workflow
label: consensus_hla
doc: |-
  ## About this Workflow

  Consensus HLA typing from NGS data using HLA-HD (v1.4).

  ### File Inputs
  - Tumour WGS paired-end sequencing reads
  - Germline WGS paired-end sequencing reads
  - Tumour RNAseq paired-end sequencing reads (optional)
  - IPD-IMGT/HLA reference index generated using bowtie2 (*.tar)


  ### Parameters
  - **patient_id**: for naming output files and report generation.
  - **bowtie2_index_prefix**: the tar file name (without the tar suffix).
  - **to_subsample**: Tumour WGS hay have high sequencing depth which will increase runtime. Use this option to subsample to the desired depth and reduce runtime.
  - **number_of_subsample_reads**: Number of reads to be subsampled from tumour paired-end WGS fastq files.


  ### Outputs
  - **sample1 json**: Tumour WGS HLA result
  - **sample2 json**: Germline WGS HLA result
  - **sample3 json**: Tumour RNAseq HLA result
  - **consensus txt and json**: Consensus alleles for all HLA genes
  - **clin_sig txt and json**: Consensus alleles for clinically significant genes
  - **hla report**: Report all consensus alleles in a PDF format


requirements:
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: patient_id
  label: Patient ID
  type: string

- id: tumour_DNA_read2_sequences
  label: Tumour DNA Read 2 Sequences
  type: File

- id: tumour_DNA_read1_sequences
  label: Tumour DNA Read 1 Sequences
  type: File

- id: Normal_DNA_read2_sequences
  label: Normal DNA Read 2 Sequences
  type: File

- id: Normal_DNA_read1_sequences
  label: Normal DNA Read 1 Sequences
  type: File

- id: bowtie2_index
  label: HLA Bowtie2 Index Archive
  doc: Bowtie2 Index Archive for an HLA reference.
  type: File

- id: RNA_read2_sequences
  label: RNA Read 2 Sequences
  doc: Read 2 sequences in FASTA or FASTQ format (may be bgzipped).
  type: File?

- id: RNA_read1_sequences
  label: RNA Read 1 Sequences
  type: File?

- id: alignment_threads
  label: alignment-threads
  type: int?

- id: bowtie2_index_prefix
  label: Bowtie2 Index Prefix
  doc: |-
    The prefix of index files contained in the Bowtie2 index TAR. Note that all Bowtie2 nidex files in the TAR should have this prefix.
  type: string


outputs:
- id: hla_report
  label: HLA Report
  doc: A PDF report containing the HLA consensus results.
  type: File
  outputSource:
  - hla_reports/hla_report

- id: sample3_json
  label: Tumour RNA HLA-HD Results JSON
  type: File?
  outputSource:
  - three_sample_hlatyping_simplified/sample3_json

- id: sample2_json
  label: Normal DNA HLA-HD Results JSON
  type: File
  outputSource:
  - three_sample_hlatyping_simplified/sample2_json

- id: sample1_json
  label: Tumour DNA HLA-HD Results JSON
  type: File
  outputSource:
  - three_sample_hlatyping_simplified/sample1_json

- id: consensus_txt
  label: HLA Consensus Text File
  type: File
  outputSource:
  - three_sample_hlatyping_simplified/consensus_txt

- id: consensus_json
  label: HLA Consensus JSON
  type: File
  outputSource:
  - three_sample_hlatyping_simplified/consensus_json

- id: clin_sig_txt
  label: Clinically Significant HLA Consensus Text File
  type: File
  outputSource:
  - three_sample_hlatyping_simplified/clin_sig_txt

- id: clin_sig_json
  label: Clinically Significant HLA Consensus JSON
  type: File
  outputSource:
  - three_sample_hlatyping_simplified/clin_sig_json

steps:
- id: hla_reports
  label: hla-report
  in:
  - id: clin_sig_hla
    source: three_sample_hlatyping_simplified/clin_sig_json
  - id: patient_id
    source: patient_id
  - id: germline_hla
    source: three_sample_hlatyping_simplified/sample2_json
  - id: tumour_hla
    source: three_sample_hlatyping_simplified/sample1_json
  - id: rna_hla
    source: three_sample_hlatyping_simplified/sample3_json
  run: consHLA.cwl.steps/hla_reports.cwl
  out:
  - id: hla_report

- id: three_sample_hlatyping_simplified
  label: three-sample-hlatyping-simplified
  in:
  - id: RNA_read1_sequences
    source: RNA_read1_sequences
  - id: bowtie2_index
    source: bowtie2_index
  - id: tumour_DNA_read2_sequences
    source: tumour_DNA_read2_sequences
  - id: tumour_DNA_read1_sequences
    source: tumour_DNA_read1_sequences
  - id: Normal_DNA_read2_sequences
    source: Normal_DNA_read2_sequences
  - id: Normal_DNA_read1_sequences
    source: Normal_DNA_read1_sequences
  - id: patient_id
    source: patient_id
  - id: RNA_read2_sequences
    source: RNA_read2_sequences
  - id: alignment_threads
    source: alignment_threads
  - id: bowtie2_index_prefix
    source: bowtie2_index_prefix
  run: consHLA.cwl.steps/three_sample_hlatyping_simplified.cwl
  out:
  - id: clin_sig_json
  - id: consensus_txt
  - id: sample1_json
  - id: consensus_json
  - id: sample2_json
  - id: sample3_json
  - id: clin_sig_txt
