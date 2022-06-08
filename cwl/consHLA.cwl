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
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: tumour_DNA_read2_sequences
  label: Tumour DNA Read 2 Sequences
  type: File
  sbg:fileTypes: FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ
  sbg:x: -740
  sbg:y: -529.8776245117188
- id: tumour_DNA_read1_sequences
  label: Tumour DNA Read 1 Sequences
  type: File
  sbg:fileTypes: FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ
  sbg:x: -738
  sbg:y: -398.8865966796875
- id: Germline_DNA_read2_sequences
  label: Germline DNA Read 2 Sequences
  type: File
  sbg:fileTypes: FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ
  sbg:x: -721.3968505859375
  sbg:y: -55.5
- id: Germline_DNA_read1_sequences
  label: Germline DNA Read 1 Sequences
  type: File
  sbg:fileTypes: FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ
  sbg:x: -722.3968505859375
  sbg:y: 69.5
- id: bowtie2_index
  label: HLA Bowtie2 Index Archive
  doc: Bowtie2 Index Archive for an HLA reference.
  type: File
  sbg:fileTypes: TAR
  sbg:x: -711.3968505859375
  sbg:y: 229.5
- id: RNA_read1_sequences
  label: RNA Read 1 Sequences
  type: File?
  sbg:fileTypes: FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ
  sbg:x: -722.4848022460938
  sbg:y: -163.3471221923828
- id: RNA_read2_sequences
  label: RNA Read 2 Sequences
  doc: Read 2 sequences in FASTA or FASTQ format (may be bgzipped).
  type: File?
  sbg:fileTypes: FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ
  sbg:x: -731.5330810546875
  sbg:y: -274.7503356933594
- id: bowtie2_index_prefix
  label: HLA Bowtie2 Index Prefix
  type: string
  sbg:exposed: true
- id: to_subsample
  label: subsample tumour DNA
  doc: |-
    Tumour DNA fastq may be large due to high sequencing depth. Subsample it to reduce runtime.
  type: boolean
  sbg:exposed: true
- id: number_of_subsample_reads
  doc: The number of reads to subsample for read2
  type: int?
  sbg:exposed: true
- id: patient_id
  label: Patient ID
  type: string
  sbg:x: -543.841796875
  sbg:y: 281.01336669921875

outputs:
- id: hla_report
  label: HLA Report
  doc: A PDF report containing the HLA consensus results.
  type: File
  outputSource:
  - hla_reports/hla_report
  sbg:fileTypes: PDF
  sbg:x: 231.20985412597656
  sbg:y: 274.91827392578125
- id: sample3_json
  label: Tumour RNA HLA-HD Results JSON
  type: File?
  outputSource:
  - three_sample_hlatyping/sample3_json
  sbg:fileTypes: JSON
  sbg:x: -30.096677780151367
  sbg:y: -306.7400207519531
- id: sample2_json
  label: Normal DNA HLA-HD Results JSON
  type: File
  outputSource:
  - three_sample_hlatyping/sample2_json
  sbg:fileTypes: JSON
  sbg:x: -25.071552276611328
  sbg:y: -182.36814880371094
- id: sample1_json
  label: Tumour DNA HLA-HD Results JSON
  type: File
  outputSource:
  - three_sample_hlatyping/sample1_json
  sbg:fileTypes: JSON
  sbg:x: -19.28140640258789
  sbg:y: -60.256282806396484
- id: consensus_txt
  label: HLA Consensus Text File
  type: File
  outputSource:
  - three_sample_hlatyping/consensus_txt
  sbg:fileTypes: TXT
  sbg:x: 302.8392028808594
  sbg:y: 27.28643226623535
- id: clin_sig_txt
  label: Clinically Significant HLA Consensus Text File
  type: File
  outputSource:
  - three_sample_hlatyping/clin_sig_txt
  sbg:fileTypes: TXT
  sbg:x: 305.53265380859375
  sbg:y: 142.72361755371094
- id: consensus_json
  label: HLA Consensus JSON
  type: File
  outputSource:
  - three_sample_hlatyping/consensus_json
  sbg:fileTypes: JSON
  sbg:x: 121.91336822509766
  sbg:y: 81.4509506225586
- id: clin_sig_json
  label: Clinically Significant HLA Consensus JSON
  type: File
  outputSource:
  - three_sample_hlatyping/clin_sig_json
  sbg:fileTypes: JSON
  sbg:x: -18
  sbg:y: 156.28140258789062

steps:
- id: three_sample_hlatyping
  label: three-sample-hlatyping
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
    source: Germline_DNA_read2_sequences
  - id: Normal_DNA_read1_sequences
    source: Germline_DNA_read1_sequences
  - id: patient_id
    source: patient_id
  - id: bowtie2_index_prefix
    source: bowtie2_index_prefix
  - id: RNA_read2_sequences
    source: RNA_read2_sequences
  - id: to_subsample
    source: to_subsample
  - id: number_of_subsample_reads
    source: number_of_subsample_reads
  run: consHLA.steps/three_sample_hlatyping.cwl
  out:
  - id: clin_sig_json
  - id: consensus_txt
  - id: sample1_json
  - id: consensus_json
  - id: sample2_json
  - id: sample3_json
  - id: clin_sig_txt
  sbg:x: -345
  sbg:y: 23
- id: hla_reports
  label: hla-report
  in:
  - id: full_hla
    source: three_sample_hlatyping/consensus_json
  - id: clin_sig_hla
    source: three_sample_hlatyping/clin_sig_json
  - id: patient_id
    source: patient_id
  run: consHLA.steps/hla_reports.cwl
  out:
  - id: hla_report
  sbg:x: -69.48744201660156
  sbg:y: 273.517578125
sbg:appVersion:
- v1.2
- v1.0
sbg:contributors:
- alanwu
sbg:createdBy: alanwu
sbg:sbgMaintained: false
sbg:toolAuthor: Weilin Wu <wwu@ccia.org.au>
sbg:validationErrors: []
sbg:workflowLanguage: CWL
