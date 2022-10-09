cwlVersion: v1.0
class: Workflow
label: three-sample-hlatyping-simplified
doc: |-
  # About this workflow
  This workflow runs HLA-HD (with restricted reads for speed) on three sets of paired-end FASTQs (tumour DNA, tumour RNA, and normal DNA) and generates a consensus type for each HLA allele.  See the `hlahd-consensus-parser` tool documentation for more information on how a consensus HLA type is determined.

  ## Outputs
  There are 4 main outputs from this tool:
  - **HLA Consensus JSON:** a JSON file containing consensus alleles for **all typed HLA genes at the highest possible accuracy level** (two- or three-field).
  - **HLA Consensus Text File:** a text file containing all **unique consensus alleles** from all HLA genes which were successfully typed by HLA-HD for which a consensus could be determined (i.e. 'Not typed' and 'No consensus' alleles are removed). The alleles are truncated to two-field accuracy. This file can be used as input to pVACtools.
  - **Clinically Significant HLA Consensus JSON:** a JSON file containing consensus alleles for **a subset of clinically significant HLA genes (HLA-A, -B, -C, -DRA, -DRB1, -DRB3, -DRB4, -DRB5, -DQA1, -DQB1, -DPA1, -DPB1) at the highest possible accuracy level** (two- or three-field).
  - **Clinically Significant HLA Consensus Text File:** a text file containing all **unique consensus alleles** from the clinically significant subset of HLA genes(defined above) which were successfully typed by HLA-HD for which a consensus could be determined (i.e. 'Not typed' and 'No consensus' alleles are removed). The alleles are truncated to two-field accuracy. This file can be used as input to pVACtools.

  In addition, a JSON file is generated containing the two/three sets of input HLA-HD results in a format consistent with that of the consensus allele JSONs.

  ## Documentation
  - [HLA-HD](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/)

requirements:
- class: LoadListingRequirement
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: RNA_read1_sequences
  label: RNA Read 1 Sequences
  type: File?
- id: bowtie2_index
  label: HLA Bowtie2 Index Archive
  doc: Bowtie2 Index Archive for an HLA reference.
  type: File
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
- id: patient_id
  label: Patient ID
  type: string
- id: bowtie2_index_prefix
  label: HLA Bowtie2 Index Prefix
  type: string
- id: RNA_read2_sequences
  label: RNA Read 2 Sequences
  doc: Read 2 sequences in FASTA or FASTQ format (may be bgzipped).
  type: File?
- id: alignment_threads
  label: alignment-threads
  type: int?

outputs:
- id: clin_sig_json
  label: Clinically Significant HLA Consensus JSON
  type: File
  outputSource:
  - hla_consensus_parser/clin_sig_json
- id: consensus_txt
  label: HLA Consensus Text File
  type: File
  outputSource:
  - hla_consensus_parser/filtered_txt
- id: sample1_json
  label: Tumour DNA HLA-HD Results JSON
  type: File
  outputSource:
  - hla_consensus_parser/tumour_dna_json
- id: consensus_json
  label: HLA Consensus JSON
  type: File
  outputSource:
  - hla_consensus_parser/consensus_json
- id: sample2_json
  label: Normal DNA HLA-HD Results JSON
  type: File
  outputSource:
  - hla_consensus_parser/normal_dna_json
- id: sample3_json
  label: Tumour RNA HLA-HD Results JSON
  type: File?
  outputSource:
  - hla_consensus_parser/rna_json
- id: clin_sig_txt
  label: Clinically Significant HLA Consensus Text File
  type: File
  outputSource:
  - hla_consensus_parser/clin_sig_txt

steps:
- id: hla_consensus_parser
  label: hlahd-consensus-parser
  in:
  - id: tumour_dna
    source: tumour_DNA_hlahd/hlahd_final
  - id: tumour_rna
    source: tumour_RNA_hlahd/hlahd_final
  - id: normal_dna
    source: normal_DNA_hlahd/hlahd_final
  - id: sample_name
    source: patient_id
  run: three_sample_hlatyping.cwl.steps/hla_consensus_parser.cwl
  out:
  - id: clin_sig_json
  - id: filtered_txt
  - id: tumour_dna_json
  - id: consensus_json
  - id: normal_dna_json
  - id: rna_json
  - id: clin_sig_txt
- id: tumour_RNA_hlahd
  label: tumour-RNA-hlahd
  in:
  - id: sample_name
    source: patient_id
  - id: output_prefix
    default: tumour_rna
  - id: bowtie2_index_prefix_1
    source: bowtie2_index_prefix
  - id: read2_sequences
    source: RNA_read2_sequences
  - id: read1_sequences
    source: RNA_read1_sequences
  - id: bowtie2_index
    source: bowtie2_index
  - id: threads
    source: alignment_threads
  run: three_sample_hlatyping.cwl.steps/tumour_RNA_hlahd.cwl
  out:
  - id: hlahd_output
  - id: hlahd_final
- id: tumour_DNA_hlahd
  label: tumour-DNA-hlahd
  in:
  - id: sample_name
    source: patient_id
  - id: output_prefix
    default: tumour_dna
  - id: bowtie2_index_prefix_1
    source: bowtie2_index_prefix
  - id: read2_sequences
    source: tumour_DNA_read2_sequences
  - id: read1_sequences
    source: tumour_DNA_read1_sequences
  - id: bowtie2_index
    source: bowtie2_index
  - id: threads
    source: alignment_threads
  run: three_sample_hlatyping.cwl.steps/tumour_DNA_hlahd.cwl
  out:
  - id: hlahd_output
  - id: hlahd_final
- id: normal_DNA_hlahd
  label: normal-DNA-hlahd
  in:
  - id: sample_name
    source: patient_id
  - id: output_prefix
    default: normal_dna
  - id: bowtie2_index_prefix_1
    source: bowtie2_index_prefix
  - id: read2_sequences
    source: Normal_DNA_read2_sequences
  - id: read1_sequences
    source: Normal_DNA_read1_sequences
  - id: bowtie2_index
    source: bowtie2_index
  - id: threads
    source: alignment_threads
  run: three_sample_hlatyping.cwl.steps/normal_DNA_hlahd.cwl
  out:
  - id: hlahd_output
  - id: hlahd_final
