cwlVersion: v1.0
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

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: sample_name
  label: Patient ID
  doc: Patient ID to be used for naming the output SAM.
  type: string
- id: read2_sequences
  label: Read 2 Sequences
  doc: Read 2 sequences in FASTA or FASTQ format (may be bgzipped).
  type: File
- id: read1_sequences
  label: Read 1 Sequences
  doc: Read 1 sequences in FASTA or FASTQ format (may be bgzipped).
  type: File
- id: bowtie2_index_prefix
  label: Bowtie2 Index Prefix
  doc: |-
    The prefix of index files contained in the Bowtie2 index TAR. Note that all Bowtie2 nidex files in the TAR should have this prefix.
  type: string
- id: bowtie2_index
  label: Bowtie2 Index Archive
  doc: |-
    A TAR archive containing Bowtie2 index files. For the purposes of speeding up HLA-HD, this should be an archive of Bowtie2 index files for an HLA region reference, such as `hla_gen` provided by the IMGT.
  type: File
- id: output_prefix
  label: HLA-HD Output Prefix
  doc: Optional prefix for HLA-HD output files and directory.
  type: string?

outputs:
- id: hla_fastq_reads1
  label: HLA Reads 1
  doc: A FASTQ file containing end 1 reads which were found to align to the HLA region.
  type: File
  outputSource:
  - samtools_fastq_1/output_fastq_2
- id: hla_fastq_reads2
  label: HLA Reads 2
  doc: A FASTQ file containing end 2 reads which were found to align to the HLA region.
  type: File
  outputSource:
  - samtools_fastq_1/output_fastq_1
- id: hlahd_output
  label: HLA-HD Output
  doc: Directory containing all HLA-HD output files.
  type: Directory
  outputSource:
  - hla_hd_1/hlahd_results
- id: hlahd_final
  label: HLA-HD Final Results File
  doc: The final results text file produced by HLA-HD.
  type: File
  outputSource:
  - hla_hd_1/hlahd_final_results

steps:
- id: bowtie2
  label: bowtie2
  doc: Run Bowtie2 alignment of input FASTQs to an HLA region reference.
  in:
  - id: bowtie2_index
    source: bowtie2_index
  - id: read1_sequences
    source: read1_sequences
  - id: read2_sequences
    source: read2_sequences
  - id: no_unaligned
    default: true
  - id: sample_name
    source: sample_name
  - id: bowtie2_index_prefix
    source: bowtie2_index_prefix
  run: tumour_rna_hlahd.cwl.steps/bowtie2.cwl
  out:
  - id: aligned_sam
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
  run: tumour_rna_hlahd.cwl.steps/samtools_view.cwl
  out:
  - id: output_alignment
- id: samtools_fastq_1
  label: samtools-fastq
  in:
  - id: input_alignment
    source: samtools_view/output_alignment
  run: tumour_rna_hlahd.cwl.steps/samtools_fastq_1.cwl
  out:
  - id: output_fastq_1
  - id: output_fastq_2
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
  run: tumour_rna_hlahd.cwl.steps/hla_hd_1.cwl
  out:
  - id: hlahd_results
  - id: hlahd_final_results

