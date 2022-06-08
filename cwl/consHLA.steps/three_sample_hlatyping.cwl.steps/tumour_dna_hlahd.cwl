cwlVersion: v1.2
class: Workflow
label: restricted-reads-hlahd
doc: |-
  # About this workflow
  This workflow runs fast HLA typing by restricting reads to those aligning with the HLA region prior to running HLA-HD.  
  Aligning the input FASTQs to a reference containing only the HLA region, then converting the mapped reads back to FASTQs allows us to restrict the reads used by HLA-HD to only those which are critically relevant to HLA typing. While it is optimal to provide all reads to HLA-HD, the runtime is prohibitively long. In our testing, this full read restriction method (including HLA-HD) has been shown to take approximately two thirds of the time taken to run HLA-HD alone on unrestricted reads.

  ## Before running this workflow
  Prior to running this workflow, you must create/download a Bowtie2 Index (generated using `bowtie2 build`) for a reference file **which only includes the HLA region on chromosome 6**. We recommend using the HLA region reference `hla_gen.fasta` provided by IMGT. It can be downloaded from the IMGT GitHub repo [here](https://github.com/ANHIG/IMGTHLA/tree/Latest/fasta) or using this download link [ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta](ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta). The Bowtie2 Index files must then be archived using the `tar` command before being used as the Bowtie2 Index Archive input to this workflow.  
  **A Bowtie2 index archive for `hla_gen.fasta` can be downloaded from the repo to be used as an input for this workflow. The corresponding Bowtie2 Index Prefix parameter should be `hla_gen`.**

  ## Steps
  This workflow follows the steps recommended by the HLA-HD authors [here](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/) under the subheadings Tips > Filtering of reads (March 6, 2019).
  1. Use `bowtie2` to map the input paired-end FASTQs to the HLA reference sequence (provided as the Bowtie2 Index Archive).
  2. Use `samtools view` to extract mapped reads (using option `-F 4`).
  3. Use `samtools fastq` to convert the BAM of aligned HLA reads to paired-end FASTQs.
  4. Run HLA-HD using the new, smaller FASTQs (containing only those reads which aligned to the HLA region) as input.

  ## Documentation
  - [HLA-HD docs](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/)
  - [HLA-HD publication](https://pubmed.ncbi.nlm.nih.gov/28419628/)
  - [Bowtie2](https://github.com/BenLangmead/bowtie2)
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: LoadListingRequirement
- class: MultipleInputFeatureRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: sample_name
  label: Patient ID
  doc: Patient ID to be used for naming the output SAM.
  type: string
  sbg:x: 0
  sbg:y: 0
- id: read1_sequences
  label: Read 1 Sequences
  doc: Read 1 sequences in FASTA or FASTQ format (may be bgzipped).
  type: File
  sbg:fileTypes: FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ
  sbg:x: -275.3482971191406
  sbg:y: 339.40435791015625
- id: bowtie2_index_prefix
  label: Bowtie2 Index Prefix
  doc: |-
    The prefix of index files contained in the Bowtie2 index TAR. Note that all Bowtie2 nidex files in the TAR should have this prefix.
  type: string
  sbg:x: 0
  sbg:y: 426.203125
- id: bowtie2_index
  label: Bowtie2 Index Archive
  doc: |-
    A TAR archive containing Bowtie2 index files. For the purposes of speeding up HLA-HD, this should be an archive of Bowtie2 index files for an HLA region reference, such as `hla_gen` provided by the IMGT.
  type: File
  sbg:fileTypes: TAR
  sbg:x: 0
  sbg:y: 532.5625
- id: output_prefix
  label: HLA-HD Output Prefix
  doc: Optional prefix for HLA-HD output files and directory.
  type: string?
  sbg:x: 0
  sbg:y: 319.734375
- id: read2_sequences
  label: Read 2 Sequences
  doc: Read 2 sequences in FASTA or FASTQ format (may be bgzipped).
  type: File
  sbg:fileTypes: FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ
  sbg:x: -256.0113830566406
  sbg:y: 27.584135055541992
- id: to_subsample
  label: to subsample?
  type: boolean
  sbg:x: -267
  sbg:y: 207
- id: number_of_subsample_reads
  doc: The number of reads to subsample for read2
  type: int?
  sbg:x: -394.8772888183594
  sbg:y: 130.20858764648438

outputs:
- id: hlahd_output
  label: HLA-HD Output
  doc: Directory containing all HLA-HD output files.
  type: Directory
  outputSource:
  - hla_hd_1/hlahd_results
  sbg:x: 1407.3223876953125
  sbg:y: 213.046875
- id: hlahd_final
  label: HLA-HD Final Results File
  doc: The final results text file produced by HLA-HD.
  type: File
  outputSource:
  - hla_hd_1/hlahd_final_results
  sbg:fileTypes: TXT
  sbg:x: 1407.3223876953125
  sbg:y: 319.515625

steps:
- id: bowtie2
  label: bowtie2
  doc: Run Bowtie2 alignment of input FASTQs to an HLA region reference.
  in:
  - id: bowtie2_index
    source: bowtie2_index
  - id: read1_sequences
    source:
    - read1_seqtk/subsampled_fastq
    - read1_sequences
    pickValue: first_non_null
  - id: read2_sequences
    source:
    - read2_seqtk/subsampled_fastq
    - read2_sequences
    pickValue: first_non_null
  - id: no_unaligned
    default: true
  - id: sample_name
    source: sample_name
  - id: bowtie2_index_prefix
    source: bowtie2_index_prefix
  run: tumour_dna_hlahd.cwl.steps/bowtie2.cwl
  out:
  - id: aligned_sam
  sbg:x: 216.125
  sbg:y: 238.3359375
- id: samtools_view
  label: samtools-view
  in:
  - id: output_format
    default: BAM
  - id: input_alignment
    source: bowtie2/aligned_sam
  - id: fast_bam_compression
    default: true
  - id: exclude_reads_any
    default: '4'
  run: tumour_dna_hlahd.cwl.steps/samtools_view.cwl
  out:
  - id: output_alignment
  sbg:x: 533.77783203125
  sbg:y: 266.3359375
- id: samtools_fastq_1
  label: samtools-fastq
  in:
  - id: input_alignment
    source: samtools_view/output_alignment
  run: tumour_dna_hlahd.cwl.steps/samtools_fastq_1.cwl
  out:
  - id: output_fastq_1
  - id: output_fastq_2
  sbg:x: 799.46533203125
  sbg:y: 259.2265625
- id: hla_hd_1
  label: hla-hd
  in:
  - id: threads
    default: 2
  - id: minimum_read_length
    default: 0
  - id: fastq_reads1
    source: samtools_fastq_1/output_fastq_1
  - id: fastq_reads2
    source: samtools_fastq_1/output_fastq_2
  - id: sample_id
    source: sample_name
  - id: output_prefix
    source: output_prefix
  run: tumour_dna_hlahd.cwl.steps/hla_hd_1.cwl
  out:
  - id: hlahd_results
  - id: hlahd_final_results
  sbg:x: 1095.0543212890625
  sbg:y: 138.9765625
- id: read2_seqtk
  label: read2_seqtk
  in:
  - id: seed
    default: 123
  - id: input_fastq
    source: read2_sequences
  - id: num_reads
    source: number_of_subsample_reads
  - id: custom_input
    source: to_subsample
  run: tumour_dna_hlahd.cwl.steps/read2_seqtk.cwl
  when: $(inputs.custom_input)
  out:
  - id: subsampled_fastq
  sbg:x: -105
  sbg:y: 117
- id: read1_seqtk
  label: read1_seqtk
  in:
  - id: seed
    default: 123
  - id: input_fastq
    source: read1_sequences
  - id: num_reads
    source: number_of_subsample_reads
  - id: custom_input
    source: to_subsample
  run: tumour_dna_hlahd.cwl.steps/read1_seqtk.cwl
  when: $(inputs.custom_input)
  out:
  - id: subsampled_fastq
  sbg:x: -102.25846099853516
  sbg:y: 249.2247314453125
sbg:appVersion:
- v1.2
- v1.0
sbg:contributors:
- alanwu
- mwonge
sbg:createdBy: mwonge
sbg:createdOn: 1642571972
sbg:sbgMaintained: false
sbg:toolAuthor: Rachel Bowen-James <rbowen-james@ccia.org.au>, Weilin Wu <wwu@ccia.org.au>
sbg:validationErrors: []
sbg:workflowLanguage: CWL
sbg:wrapperAuthor: Rachel Bowen-James <rbowen-james@ccia.org.au>, Weilin Wu <wwu@ccia.org.au>
