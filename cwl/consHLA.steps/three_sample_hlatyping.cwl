cwlVersion: v1.2
class: Workflow
label: three-sample-hlatyping
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
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: LoadListingRequirement
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: RNA_read1_sequences
  label: RNA Read 1 Sequences
  type: File?
  sbg:fileTypes: FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ
  sbg:x: 0
  sbg:y: 103.44762420654297
- id: bowtie2_index
  label: HLA Bowtie2 Index Archive
  doc: Bowtie2 Index Archive for an HLA reference.
  type: File
  sbg:fileTypes: TAR
  sbg:x: 0
  sbg:y: 852.4296875
- id: tumour_DNA_read2_sequences
  label: Tumour DNA Read 2 Sequences
  type: File
  sbg:fileTypes: FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ
  sbg:x: 0
  sbg:y: 213.15625
- id: tumour_DNA_read1_sequences
  label: Tumour DNA Read 1 Sequences
  type: File
  sbg:fileTypes: FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ
  sbg:x: 0
  sbg:y: 319.734375
- id: Normal_DNA_read2_sequences
  label: Normal DNA Read 2 Sequences
  type: File
  sbg:fileTypes: FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ
  sbg:x: -3.9518749713897705
  sbg:y: 538.58642578125
- id: Normal_DNA_read1_sequences
  label: Normal DNA Read 1 Sequences
  type: File
  sbg:fileTypes: FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ
  sbg:x: -1.3172919750213623
  sbg:y: 652.1729125976562
- id: patient_id
  label: Patient ID
  type: string
  sbg:x: -1.5587761402130127
  sbg:y: 430.67633056640625
- id: bowtie2_index_prefix
  label: HLA Bowtie2 Index Prefix
  type: string
  sbg:x: 0
  sbg:y: 762.1248168945312
- id: RNA_read2_sequences
  label: RNA Read 2 Sequences
  doc: Read 2 sequences in FASTA or FASTQ format (may be bgzipped).
  type: File?
  sbg:fileTypes: FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ
  sbg:x: -0.7971814274787903
  sbg:y: -9.009395599365234
- id: to_subsample
  label: subsample tumour DNA
  doc: |-
    Tumour DNA fastq may be large due to high sequencing depth. Subsample it to reduce runtime.
  type: boolean
  sbg:x: -226.47021484375
  sbg:y: 343.0579833984375
- id: number_of_subsample_reads
  doc: The number of reads to subsample for read2
  type: int?
  sbg:exposed: true

outputs:
- id: clin_sig_json
  label: Clinically Significant HLA Consensus JSON
  type: File
  outputSource:
  - hlahd_consensus_parser/clin_sig_json
  sbg:fileTypes: JSON
  sbg:x: 1265.2135009765625
  sbg:y: 746.2421875
- id: consensus_txt
  label: HLA Consensus Text File
  type: File
  outputSource:
  - hlahd_consensus_parser/filtered_txt
  sbg:fileTypes: TXT
  sbg:x: 1265.2135009765625
  sbg:y: 426.1171875
- id: sample1_json
  label: Tumour DNA HLA-HD Results JSON
  type: File
  outputSource:
  - hlahd_consensus_parser/sample1_json
  sbg:fileTypes: JSON
  sbg:x: 1265.2135009765625
  sbg:y: 319.4609375
- id: consensus_json
  label: HLA Consensus JSON
  type: File
  outputSource:
  - hlahd_consensus_parser/consensus_json
  sbg:fileTypes: JSON
  sbg:x: 1265.2135009765625
  sbg:y: 532.7734375
- id: sample2_json
  label: Normal DNA HLA-HD Results JSON
  type: File
  outputSource:
  - hlahd_consensus_parser/sample2_json
  sbg:fileTypes: JSON
  sbg:x: 1265.2135009765625
  sbg:y: 212.7265625
- id: sample3_json
  label: Tumour RNA HLA-HD Results JSON
  type: File?
  outputSource:
  - hlahd_consensus_parser/sample3_json
  sbg:fileTypes: JSON
  sbg:x: 1265.2135009765625
  sbg:y: 105.9921875
- id: clin_sig_txt
  label: Clinically Significant HLA Consensus Text File
  type: File
  outputSource:
  - hlahd_consensus_parser/clin_sig_txt
  sbg:fileTypes: TXT
  sbg:x: 1265.2135009765625
  sbg:y: 639.5078125

steps:
- id: tumour_rna_hlahd
  label: RNA-hlahd
  doc: Run restricted reads HLA-HD typing from tumour RNA.
  in:
  - id: sample_name
    source: patient_id
  - id: read2_sequences
    source: RNA_read2_sequences
  - id: read1_sequences
    source: RNA_read1_sequences
  - id: bowtie2_index_prefix
    source: bowtie2_index_prefix
  - id: bowtie2_index
    source: bowtie2_index
  - id: output_prefix
    default: tumour_rna
  run: three_sample_hlatyping.cwl.steps/tumour_rna_hlahd.cwl
  when: |-
    ${
        if (inputs.read1_sequences === null || inputs.read2_sequences === null) {
            return false
        } else {
            return true
        }
    }
  out:
  - id: hla_fastq_reads1
  - id: hla_fastq_reads2
  - id: hlahd_output
  - id: hlahd_final
  sbg:x: 306.1981201171875
  sbg:y: 231.8018798828125
- id: normal_dna_hlahd
  label: normal-DNA-hlahd
  doc: Run restricted reads HLA-HD typing from normal DNA.
  in:
  - id: sample_name
    source: patient_id
  - id: read2_sequences
    source: Normal_DNA_read2_sequences
  - id: read1_sequences
    source: Normal_DNA_read1_sequences
  - id: bowtie2_index_prefix
    source: bowtie2_index_prefix
  - id: bowtie2_index
    source: bowtie2_index
  - id: output_prefix
    default: normal_dna
  run: three_sample_hlatyping.cwl.steps/normal_dna_hlahd.cwl
  out:
  - id: hla_fastq_reads1
  - id: hla_fastq_reads2
  - id: hlahd_output
  - id: hlahd_final
  sbg:x: 303.71875
  sbg:y: 560.734375
- id: hlahd_consensus_parser
  label: hlahd-consensus-parser
  in:
  - id: tumour_dna
    source: tumour_dna_hlahd/hlahd_final
  - id: tumour_rna
    source: tumour_rna_hlahd/hlahd_final
  - id: normal_dna
    source: normal_dna_hlahd/hlahd_final
  - id: sample_name
    source: patient_id
  run: three_sample_hlatyping.cwl.steps/hlahd_consensus_parser.cwl
  out:
  - id: clin_sig_json
  - id: filtered_txt
  - id: sample1_json
  - id: consensus_json
  - id: sample2_json
  - id: sample3_json
  - id: clin_sig_txt
  sbg:x: 740.4396362304688
  sbg:y: 408.8583068847656
- id: tumour_dna_hlahd
  label: tumout-DNA-hlahd
  in:
  - id: sample_name
    source: patient_id
  - id: read1_sequences
    source: tumour_DNA_read1_sequences
  - id: bowtie2_index_prefix
    source: bowtie2_index_prefix
  - id: bowtie2_index
    source: bowtie2_index
  - id: output_prefix
    default: tumour_dna
  - id: read2_sequences
    source: tumour_DNA_read2_sequences
  - id: to_subsample
    default: false
    source: to_subsample
  - id: number_of_subsample_reads
    source: number_of_subsample_reads
  run: three_sample_hlatyping.cwl.steps/tumour_dna_hlahd.cwl
  out:
  - id: hlahd_output
  - id: hlahd_final
  sbg:x: 301.8226623535156
  sbg:y: 399.18768310546875
sbg:appVersion:
- v1.2
- v1.0
sbg:content_hash: a492bbf10043e57be7ae62b411e4e5f0bbc06a84b52273ec4d22c19fa7812f870
sbg:contributors:
- alanwu
sbg:createdBy: alanwu
sbg:createdOn: 1644192825
sbg:sbgMaintained: false
sbg:toolAuthor: Rachel Bowen-James <rbowen-james@ccia.org.au>, Weilin Wu <wwu@ccia.org.au>
sbg:validationErrors: []
sbg:workflowLanguage: CWL
sbg:wrapperAuthor: Rachel Bowen-James <rbowen-james@ccia.org.au>, Weilin Wu <wwu@ccia.org.au>
