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
  run:
    cwlVersion: v1.2
    class: CommandLineTool
    label: hla-report
    doc: |-
      # About this tool

      This tool generates an HLA report from the output of disTIL HLA consensus HLA typing.
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: alanwucci/conshla_report:2.0.0
    - class: InlineJavascriptRequirement

    inputs:
    - id: clin_sig_hla
      label: Clinically Significant HLA Consensus JSON
      doc: |-
        A JSON file containing consensus alleles for the clinically significant classical HLA genes. Only classical Class-I and Class-II genes are included and all alleles are truncated to two-field accuracy.
      type: File
      
    - id: patient_id
      label: Patient ID
      doc: ID for the patient.
      type: string
    - id: germline_hla
      label: germline_hla
      type: File
      
    - id: tumour_hla
      label: tumour_hla
      type: File
      
    - id: rna_hla
      label: rna_hla
      type: File

    outputs:
    - id: hla_report
      label: HLA Report
      doc: A PDF report containing the HLA consensus results.
      type: File
      outputBinding:
        glob: "${\n    return inputs.patient_id + '_hlaReport.pdf'\n}"
      

    baseCommand: []
    arguments:
    - prefix: ''
      position: 2
      valueFrom: |-
        ${
            var cmd = 'Rscript -e "rmarkdown::render(\'hla_report_generator.Rmd\',params=list(germline_hla_json=\'' + inputs.germline_hla.path + '\', tumour_hla_json=\'' + inputs.tumour_hla.path + '\', rna_hla_json=\'' + inputs.rna_hla.path + '\', clin_hla_json=\'' + inputs.clin_sig_hla.path + '\', pid=\'' + inputs.patient_id + '\'), output_file=paste(\'' + inputs.patient_id + '\', \'_hlaReport.pdf\', sep=\'\'))\"'
            return cmd
        }
      shellQuote: false
    - prefix: ''
      position: 10
      valueFrom: 1>&2
      shellQuote: false
    - prefix: ''
      position: 0
      valueFrom: ' cp /hla_report_generator.Rmd . && ls 1>&2 &&'
      shellQuote: false
    - prefix: ''
      position: 1
      valueFrom: echo "CLIN" 1>&2 && head $(inputs.clin_sig_hla.path) 1>&2 &&
      shellQuote: false
    id: hla-reports/5
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
  run:
    cwlVersion: v1.2
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
      
      
    - id: RNA_read2_sequences
      label: RNA Read 2 Sequences
      doc: Read 2 sequences in FASTA or FASTQ format (may be bgzipped).
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
      run:
        cwlVersion: v1.0
        class: CommandLineTool
        label: hlahd-consensus-parser
        doc: |-
          # About this tool
          This tool produces a consensus HLA type based on two or three sets of candidate alleles produced by running HLA-HD on any combination of tumour/normal WGS/WES/RNA-seq (e.g. tumour RNA-seq, normal WGS, tumour WGS) samples for a single patient.

          ## Before running this tool
          Before running this tool, you need to have obtained the HLA-HD results for two/three samples. You should use the 'final results' text files output by HLA-HD as inputs to this tool.

          ## Logic
          This tool compares the HLA types determined from three samples, and defines a 'consensus' type for each allele according to the following rules:
          - An allele must have at least 2-field (4-digit) accuracy in order to be considered, otherwise it is considered 'Not typed'.
          - If a 3-field allele (e.g. HLA-A*01:01:01) is an exact match across at least two samples, it is taken as the consensus type.
          - If there is no 3-field match across at least two samples, then the allele is truncated to 2-fields (e.g. HLA-A*01:01) and again compared between samples. If this 2-field allele is an exact match across at least two samples, then it is taken as the consensus type.
          - If there is inadequate support for an allele, it is set as 'No consensus'.

          Using the above logic, this tool only provides a consensus type for an allele if there is adequate support for it (i.e. exact match across at least two samples). These consensus types may be alleles with 2- or 3-field accuracy.  
          Note that if an allele is set as 'Not typed', this means HLA-HD was unsuccessful in typing the allele. If an allele is set as 'No consensus', it was typed by HLA-HD but there was insufficient evidence to call a consensus allele.

          ## Outputs
          There are 4 main outputs from this tool:
          - **HLA Consensus JSON:** a JSON file containing consensus alleles for **all typed HLA genes at the highest possible accuracy level** (two- or three-field).
          - **HLA Consensus Text File:** a text file containing all **unique consensus alleles** from all HLA genes which were successfully typed by HLA-HD for which a consensus could be determined (i.e. 'Not typed' and 'No consensus' alleles are removed). The alleles are truncated to two-field accuracy. This file can be used as input to pVACtools.
          - **Clinically Significant HLA Consensus JSON:** a JSON file containing consensus alleles for **a subset of clinically significant HLA genes (HLA-A, -B, -C, -DRA, -DRB1, -DRB3, -DRB4, -DRB5, -DQA1, -DQB1, -DPA1, -DPB1) at the highest possible accuracy level** (two- or three-field).
          - **Clinically Significant HLA Consensus Text File:** a text file containing all **unique consensus alleles** from the clinically significant subset of HLA genes(defined above) which were successfully typed by HLA-HD for which a consensus could be determined (i.e. 'Not typed' and 'No consensus' alleles are removed). The alleles are truncated to two-field accuracy. This file can be used as input to pVACtools.

          In addition, a JSON file is generated containing the two/three sets of input HLA-HD results in a format consistent with that of the consensus allele JSONs.

          ## Docker
          This tool uses the Docker image: `alanwucci/consensus_parser:1.0.0`.

          ## Documentation
          - [HLA-HD](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/)
        $namespaces:
          sbg: https://sevenbridges.com

        requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: alanwucci/consensus_parser:1.0.0

        inputs:
        - id: tumour_dna
          label: Sample 1 HLA-HD Results
          doc: The final results text file produced by running HLA-HD on the tumour
            DNA sample.
          type: File
          inputBinding:
            prefix: -tumour
            position: 1
            shellQuote: false
          
        - id: tumour_rna
          label: Sample 3 HLA-HD Results
          doc: The final results text file produced by running HLA-HD on the tumour
            RNA sample.
          type: File?
          inputBinding:
            prefix: -rnaseq
            position: 3
            shellQuote: false
          
        - id: normal_dna
          label: Sample 2 HLA-HD Results
          doc: The final results text file produced by running HLA-HD on the normal
            DNA sample.
          type: File
          inputBinding:
            prefix: -germline
            position: 2
            shellQuote: false
          
        - id: sample_name
          label: Sample Name
          doc: The sample name used to name the output files.
          type: string
          inputBinding:
            position: 1
            shellQuote: false

        outputs:
        - id: clin_sig_json
          label: Clinically Significant HLA Genes Consensus JSON
          doc: |-
            A JSON file containing consensus alleles for the clinically significant classical HLA genes: A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1.
          type: File
          outputBinding:
            glob: '*Sample_hla.consensus.clinSig.json'
          
        - id: filtered_txt
          label: HLA Consensus Text File
          doc: |-
            A text file containing a comma-separated list of consensus alleles for all typed HLA genes. Only includes alleles for which a consensus could be determined. All alleles are truncated to two-field accuracy. This file should be used as input for pVACseq.
          type: File
          outputBinding:
            glob: '*Sample_hla.consensus.trunc.txt'
          
        - id: tumour_dna_json
          label: Sample 1 HLA-HD Results JSON
          doc: HLA-HD final results obtained from sample 1 and represented in JSON
            format.
          type: File
          outputBinding:
            glob: '*tumour_hla.json'
          
        - id: consensus_json
          label: HLA Consensus JSON
          doc: A JSON file containing cosnensus alleles for all typed HLA genes.
          type: File
          outputBinding:
            glob: '*Sample_hla.consensus.json'
          
        - id: normal_dna_json
          label: Sample 2 HLA-HD Results JSON
          doc: HLA-HD final results obtained from sample 2 and represented in JSON
            format.
          type: File
          outputBinding:
            glob: '*germline_hla.json'
          
        - id: rna_json
          label: Sample 3 HLA-HD Results JSON
          doc: HLA-HD final results obtained from sample 3 and represented in JSON
            format.
          type: File?
          outputBinding:
            glob: '*rnaseq_hla.json'
          
        - id: clin_sig_txt
          label: Clinically Significant HLA Genes Consensus Text File
          doc: |-
            A text file containing a comma-separated list of consensus alleles for the clinically significant classical HLA genes: A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1.. Only includes alleles for which a consensus could be determined. All alleles are truncated to two-field accuracy. This file should be used as input for pVACseq.
          type: File
          outputBinding:
            glob: '*Sample_hla.consensus.clinSig.trunc.txt'
          

        baseCommand:
        - python3
        - /app/hlahd_consensus_parser_v2.py
        arguments:
        - prefix: ''
          position: 4
          valueFrom: '> hla'
          separate: false
          shellQuote: false
        id: hla-consensus-parser/7
        
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
      - id: read2_sequences
        source: RNA_read2_sequences
      - id: read1_sequences
        source: RNA_read1_sequences
      - id: bowtie2_index
        source: bowtie2_index
      - id: threads
        source: alignment_threads
      - id: bowtie2_index_prefix
        source: bowtie2_index_prefix
      run:
        cwlVersion: v1.2
        class: Workflow
        label: restricted-reads-hlahd-simplified
        doc: |-
          # About this workflow
          This workflow runs fast HLA typing by restricting reads to those aligning with the HLA region prior to running HLA-HD.  
          Aligning the input FASTQs to a reference containing only the HLA region, then converting the mapped reads back to FASTQs allows us to restrict the reads used by HLA-HD to only those which are critically relevant to HLA typing. While it is optimal to provide all reads to HLA-HD, the runtime is prohibitively long. In our testing, this full read restriction method (including HLA-HD) has been shown to take approximately two thirds of the time taken to run HLA-HD alone on unrestricted reads.

          ## Before running this workflow
          Prior to running this workflow, you must create/download a Bowtie2 Index (generated using `bowtie2 build`) for a reference file **which only includes the HLA region on chromosome 6**. We recommend using the HLA region reference `hla_gen.fasta` provided by IMGT. It can be downloaded from the IMGT GitHub repo [here](https://github.com/ANHIG/IMGTHLA/tree/Latest/fasta) or using this download link [ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta](ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta). The Bowtie2 Index files must then be archived using the `tar` command before being used as the Bowtie2 Index Archive input to this workflow.  
          **A Bowtie2 index archive for `hla_gen.fasta` can be downloaded from the disTIL repo to be used as an input for this workflow. The corresponding Bowtie2 Index Prefix parameter should be `hla_gen`.**

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
        - class: InlineJavascriptRequirement
        - class: StepInputExpressionRequirement

        inputs:
        - id: sample_name
          label: Patient ID
          doc: Patient ID to be used for naming the output SAM.
          type: string
          
          
        - id: output_prefix
          label: HLA-HD Output Prefix
          doc: Optional prefix for HLA-HD output files and directory.
          type: string?
          
          
        - id: read2_sequences
          label: Read 2 Sequences
          doc: Read 2 sequences in FASTA or FASTQ format (may be bgzipped).
          type: File
          
          
          
        - id: read1_sequences
          label: Read 1 Sequences
          doc: Read 1 sequences in FASTA or FASTQ format (may be bgzipped).
          type: File
          
          
          
        - id: bowtie2_index
          label: Bowtie2 Index Archive
          doc: A TAR archive containing Bowtie2 index files.
          type: File
          
          
          
        - id: threads
          type: int?
          
          
        - id: bowtie2_index_prefix
          label: Bowtie2 Index Prefix
          doc: |-
            The prefix of index files contained in the Bowtie2 index TAR. Note that all Bowtie2 nidex files in the TAR should have this prefix.
          type: string
          

        outputs:
        - id: hlahd_final
          label: HLA-HD Final Results File
          doc: The final results text file produced by HLA-HD.
          type: File
          outputSource:
          - hla_hd/hlahd_final_results
          
          
          

        steps:
        - id: samtools_view
          label: samtools-view
          in:
          - id: output_format
            default: BAM
          - id: input_alignment
            source: bowtie3/aligned_sam
          - id: fast_bam_compression
            default: true
          - id: exclude_reads_any
            default: '4'
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            label: samtools-view
            doc: |-
              # About this tool
              This app runs samtools view on an input alignment.

              ## Docker
              This CWL tool uses the Docker image `rachelbj/samtools:1.10.0` which contains:
              - htslib v1.10.2
              - bcftools v1.10.2
              - samtools v1.10

              ## Documentation
              - [samtools](http://www.htslib.org/doc/samtools.html)
            $namespaces:
              sbg: https://www.sevenbridges.com/

            requirements:
            - class: ShellCommandRequirement
            - class: DockerRequirement
              dockerPull: rachelbj/samtools:1.10.0
            - class: InlineJavascriptRequirement
              expressionLib:
              - |2-

                var setMetadata = function(file, metadata) {
                    if (!('metadata' in file))
                        file['metadata'] = metadata;
                    else {
                        for (var key in metadata) {
                            file['metadata'][key] = metadata[key];
                        }
                    }
                    return file
                };

                var inheritMetadata = function(o1, o2) {
                    var commonMetadata = {};
                    if (!Array.isArray(o2)) {
                        o2 = [o2]
                    }
                    for (var i = 0; i < o2.length; i++) {
                        var example = o2[i]['metadata'];
                        for (var key in example) {
                            if (i == 0)
                                commonMetadata[key] = example[key];
                            else {
                                if (!(commonMetadata[key] == example[key])) {
                                    delete commonMetadata[key]
                                }
                            }
                        }
                    }
                    if (!Array.isArray(o1)) {
                        o1 = setMetadata(o1, commonMetadata)
                    } else {
                        for (var i = 0; i < o1.length; i++) {
                            o1[i] = setMetadata(o1[i], commonMetadata)
                        }
                    }
                    return o1;
                };

            inputs:
            - id: output_format
              label: Output Format
              doc: Specifies BAM, SAM or CRAM output format.
              type:
              - 'null'
              - name: output_format
                type: enum
                symbols:
                - BAM
                - SAM
                - CRAM
              inputBinding:
                prefix: --output-fmt
                position: 0
                shellQuote: false
            - id: input_alignment
              label: Input Alignment
              doc: An input BAM, SAM or CRAM file.
              type: File
              inputBinding:
                position: 10
                shellQuote: false
              
            - id: fast_bam_compression
              label: Fast BAM Compression
              doc: Whether the output file (which must be BAM format) should be compressed.
              type: boolean?
              inputBinding:
                prefix: '-1'
                position: 0
                shellQuote: false
            - id: include_header
              label: Include Header
              doc: Whether the input alignment header should be included in the output
                file.
              type: boolean?
              inputBinding:
                prefix: -h
                position: 0
                shellQuote: false
            - id: header_only
              label: Output Header Only
              doc: |-
                When this option is selected, the output file will only contain the header from the input alignment.
              type: boolean?
              inputBinding:
                prefix: -H
                position: 0
                shellQuote: false
            - id: include_reads
              label: Include reads with all of these flags
              doc: |-
                Only output alignments with all bits set in this integer present in the FLAG field.
              type: string?
              inputBinding:
                prefix: -f
                position: 0
                shellQuote: false
            - id: exclude_reads_any
              label: Exclude reads with any of these flags
              doc: |-
                Do not output alignments with any bits set in this integer present in the FLAG field.
              type: string?
              inputBinding:
                prefix: -F
                position: 0
                shellQuote: false
            - id: exclude_reads_all
              label: Exclude reads with all of these flags
              doc: |-
                Only exclude reads with all of the bits set in this integer present in the FLAG field.
              type: string?
              inputBinding:
                prefix: -G
                position: 0
                shellQuote: false
            - id: bed_file
              label: BED File
              doc: Only output alignments overlapping the regions specified in this
                BED file.
              type: File?
              inputBinding:
                prefix: -L
                position: 0
                shellQuote: false
              

            outputs:
            - id: output_alignment
              label: Output Alignment
              doc: Output from samtools view.
              type: File
              outputBinding:
                glob: |-
                  ${
                      //Find the input basename (without the file extension)
                      var input_split = inputs.input_alignment.path.split('/')
                      var input_base = input_split[input_split.length - 1].split('.').slice(0,-1). join('.')
                      
                      var ext = ""
                      //Determine the output file extension
                      if (inputs.fast_bam_compression || inputs.output_format == "BAM") {
                          ext = ".bam"
                      } else if (inputs.output_format == "SAM") {
                          ext = ".sam"
                      } else if (inputs.output_format == "CRAM") {
                          ext = ".cram"
                      } else {
                          ext = "." + input_split[input_split.length - 1].split('.').slice(-1)
                      }
                      
                      //If only output header then add '.header'
                      if (inputs.header_only) {
                          return input_base + '.header' + ext
                      }
                      
                      //If filtered on flags/bed file then add '.filtered'
                      if (inputs.include_reads || inputs.exclude_reads_any || inputs.exclude_reads_all || inputs.bed_file) {
                          return input_base + '.filtered' + ext
                      }
                      
                      return input_base + ext
                  }
                outputEval: $(inheritMetadata(self, inputs.input_alignment))
              

            baseCommand:
            - samtools
            - view
            arguments:
            - prefix: -o
              position: 0
              valueFrom: |-
                ${
                    //Find the input basename (without the file extension)
                    var input_split = inputs.input_alignment.path.split('/')
                    var input_base = input_split[input_split.length - 1].split('.').slice(0,-1). join('.')
                    
                    var ext = ""
                    //Determine the output file extension
                    if (inputs.fast_bam_compression || inputs.output_format == "BAM") {
                        ext = ".bam"
                    } else if (inputs.output_format == "SAM") {
                        ext = ".sam"
                    } else if (inputs.output_format == "CRAM") {
                        ext = ".cram"
                    } else {
                        ext = "." + input_split[input_split.length - 1].split('.').slice(-1)
                    }
                    
                    //If only output header then add '.header'
                    if (inputs.header_only) {
                        return input_base + '.header' + ext
                    }
                    
                    //If filtered on flags/bed file then add '.filtered'
                    if (inputs.include_reads || inputs.exclude_reads_any || inputs.exclude_reads_all || inputs.bed_file) {
                        return input_base + '.filtered' + ext
                    }
                    
                    return input_base + ext
                }
              shellQuote: false
            id: samtools_view
            
          out:
          - id: output_alignment
          
          
        - id: samtools_fastq_1
          label: samtools-fastq
          in:
          - id: input_alignment
            source: samtools_view/output_alignment
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            label: samtools-fastq
            doc: |-
              # About this tool
              This app runs samtools fastq on an input alignment.  
              Note that the output paired-end FASTQs have the 'paired_end' metadata field populated with '1' or '2'.

              ## Docker
              This CWL tool uses the Docker image `rachelbj/samtools:1.10.0` which contains:
              - htslib v1.10.2
              - bcftools v1.10.2
              - samtools v1.10

              ## Documentation
              - [samtools](http://www.htslib.org/doc/samtools.html)
            $namespaces:
              sbg: https://www.sevenbridges.com/

            requirements:
            - class: ShellCommandRequirement
            - class: DockerRequirement
              dockerPull: rachelbj/samtools:1.10.0
            - class: InlineJavascriptRequirement

            inputs:
            - id: input_alignment
              label: Input Alignment
              doc: An input BAM, SAM or CRAM file.
              type: File
              inputBinding:
                position: 10
                shellQuote: false
              

            outputs:
            - id: output_fastq_1
              label: Paired-End FASTQ 1
              doc: Paired-end FASTQ 1 output by samtools fastq.
              type: File
              outputBinding:
                glob: |-
                  ${
                      var input_split = inputs.input_alignment.path.split('/')
                      var input_base = input_split[input_split.length - 1]

                      return input_base + '.pe_1.fastq'
                  }
                outputEval: |-
                  ${
                    var out = self[0];
                    out.metadata = {'paired_end' : ''}
                    out.metadata['paired_end'] = 1;
                    
                    return out
                  }
              
            - id: output_fastq_2
              label: Paired-End FASTQ 2
              doc: Paired-end FASTQ 2 output by samtools fastq.
              type: File
              outputBinding:
                glob: |-
                  ${
                      var input_split = inputs.input_alignment.path.split('/')
                      var input_base = input_split[input_split.length - 1]

                      return input_base + '.pe_2.fastq'
                  }
                outputEval: |-
                  ${
                    var out = self[0];
                    out.metadata = {'paired_end' : ''}
                    out.metadata['paired_end'] = 1;
                    
                    return out
                  }
              

            baseCommand:
            - samtools
            - fastq
            arguments:
            - prefix: '-1'
              position: 0
              valueFrom: |-
                ${
                    var input_split = inputs.input_alignment.path.split('/')
                    var input_base = input_split[input_split.length - 1]

                    return input_base + '.pe_1.fastq'
                }
              shellQuote: false
            - prefix: '-2'
              position: 0
              valueFrom: |-
                ${
                    var input_split = inputs.input_alignment.path.split('/')
                    var input_base = input_split[input_split.length - 1]

                    return input_base + '.pe_2.fastq'
                }
              shellQuote: false
            id: samtools_fastq
            
          out:
          - id: output_fastq_1
          - id: output_fastq_2
          
          
        - id: bowtie3
          label: bowtie2
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
          - id: threads
            source: threads
          run:
            cwlVersion: v1.2
            class: CommandLineTool
            label: bowtie2
            doc: |-
              # About this tool
              This tool runs Bowtie2 (**v2.4.1**) alignment of reads to a reference sequence.

              ## Inputs and parameters
              - **Bowtie2 Index Archive:** an archive (TAR) of index files generated from a reference genome (using Bowtie2 build command).
              - **Bowtie2 Index Prefix:** the basename of the index files in the Bowtie2 Index Archive. The basename should be shared by all index files in the TAR, and is the name of the index files up to but not including the final `.1.bt2` etc.
              - **Read 1 Sequences:** file containing mate 1 reads (for paired-end sequencing).
              - **Read 2 Sequences:** file containing mate 2 reads (for paired-end sequencing).
              - **Sample Name:** the name of the sample being analysed (used as the name for the output SAM file).

              ## Docker
              This tool uses the Docker image: `biocontainers/bowtie2:v2.4.1_cv1`.

              ## Documentation
              - [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
            $namespaces:
              sbg: https://sevenbridges.com

            requirements:
            - class: ShellCommandRequirement
            - class: ResourceRequirement
              coresMin: $(inputs.threads)
            - class: DockerRequirement
              dockerPull: biocontainers/bowtie2:v2.4.1_cv1
            - class: InlineJavascriptRequirement
              expressionLib:
              - |2-

                var setMetadata = function(file, metadata) {
                    if (!('metadata' in file))
                        file['metadata'] = metadata;
                    else {
                        for (var key in metadata) {
                            file['metadata'][key] = metadata[key];
                        }
                    }
                    return file
                };

                var inheritMetadata = function(o1, o2) {
                    var commonMetadata = {};
                    if (!Array.isArray(o2)) {
                        o2 = [o2]
                    }
                    for (var i = 0; i < o2.length; i++) {
                        var example = o2[i]['metadata'];
                        for (var key in example) {
                            if (i == 0)
                                commonMetadata[key] = example[key];
                            else {
                                if (!(commonMetadata[key] == example[key])) {
                                    delete commonMetadata[key]
                                }
                            }
                        }
                    }
                    if (!Array.isArray(o1)) {
                        o1 = setMetadata(o1, commonMetadata)
                    } else {
                        for (var i = 0; i < o1.length; i++) {
                            o1[i] = setMetadata(o1[i], commonMetadata)
                        }
                    }
                    return o1;
                };

            inputs:
            - id: bowtie2_index
              label: Bowtie2 Index Archive
              doc: A TAR archive containing Bowtie2 index files.
              type: File
              
            - id: read1_sequences
              label: Read 1 Sequences
              doc: Read 1 sequences in FASTA or FASTQ format (may be bgzipped).
              type: File
              inputBinding:
                prefix: '-1'
                position: 2
                shellQuote: false
              
            - id: read2_sequences
              label: Read 2 Sequences
              doc: Read 2 sequences in FASTA or FASTQ format (may be bgzipped).
              type: File
              inputBinding:
                prefix: '-2'
                position: 3
                shellQuote: false
              
            - id: no_unaligned
              label: Suppress Unaligned Reads
              doc: Suppres SAM records for unaligned reads.
              type: boolean?
              inputBinding:
                position: 1
                valueFrom: |-
                  ${
                      if (self == true) {
                          return "--no-unal"
                      } else {
                          return ""
                      }
                  }
                shellQuote: false
            - id: sample_name
              label: Sample Name
              doc: Sample name used to name the output SAM file.
              type: string
            - id: bowtie2_index_prefix
              label: Bowtie2 Index Prefix
              doc: |-
                The prefix of index files contained in the Bowtie2 index TAR. Note that all Bowtie2 nidex files in the TAR should have this prefix.
              type: string
              inputBinding:
                prefix: -x
                position: 1
                shellQuote: false
            - id: threads
              type: int?
              default: 8
              inputBinding:
                prefix: -p
                position: 0
                shellQuote: false

            outputs:
            - id: aligned_sam
              type: File
              outputBinding:
                glob: '*.sam'
                outputEval: $(inheritMetadata(self, inputs.read1_sequences))

            baseCommand: []
            arguments:
            - prefix: -S
              position: 4
              valueFrom: $(inputs.sample_name).sam
              shellQuote: false
            - prefix: ''
              position: 0
              valueFrom: |-
                ${
                    var ref = inputs.bowtie2_index.path
                    return "tar -xvf " + ref + " &&"
                }
              shellQuote: false
            - prefix: ''
              position: 0
              valueFrom: bowtie2
              shellQuote: false
            id: bowtie2/6
            
          out:
          - id: aligned_sam
          
          
        - id: hla_hd
          label: hla-hd
          in:
          - id: threads
            source: threads
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
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            label: hla-hd
            doc: |-
              ## About HLA-HD
              HLA-HD documentation and release notes can be found [here](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/).
              HLA-HD (HLA typing from High-quality Dictionary) can accurately determine HLA alleles with 6-digit precision from NGS data (FASTQ format). RNA-Seq data can also be applied.

              ## About this CWL tool
              This CWL tool runs HLA-HD version 1.4.0 on paired-end FASTQ files to determine HLA type.

              ### Inputs and parameters
              - The input paired-end read files can be from **WGS/WES or RNA-seq**.
              - The input paired-end read files must be in FASTQ format (**not zipped**).
              - The default minimum read length is 100, however this is often too strict. Choose a lower threshold to include more reads.

              ### Output
              - HLA-HD results are output to a directory named using the input sample id.
              - The final summary of HLA typing results can be found at the following path: `<output_dir_name>/result/<sample_id>_final.result.txt`.

              ### Other notes
              - This tool uses the HLA dictionary created from release 3.15.0 of the [IPD-IMGT/HLA](https://github.com/ANHIG/IMGTHLA) database.
              - This tool by default uses HLA allele frequency data included with the HLA-HD release 1.4.0.
            $namespaces:
              sbg: https://sevenbridges.com

            requirements:
            - class: ShellCommandRequirement
            - class: LoadListingRequirement
            - class: DockerRequirement
              dockerPull: pgc-images.sbgenomics.com/rbowen_james/hla-hd
            - class: InlineJavascriptRequirement

            inputs:
            - id: threads
              label: Threads
              doc: Number of cores used to execute the program.
              type: int?
              inputBinding:
                prefix: -t
                position: 1
                shellQuote: false
            - id: minimum_read_length
              label: Minimum read length
              doc: A read whose length is shorter than this parameter is ignored.
              type: int?
              inputBinding:
                prefix: -m
                position: 1
                shellQuote: false
              
            - id: fastq_reads1
              label: FASTQ Reads 1
              doc: Paired-end reads 1 in FASTQ format.
              type: File
              inputBinding:
                position: 2
                shellQuote: false
              
            - id: fastq_reads2
              label: FASTQ Reads 2
              doc: Paired-end reads 2 in FASTQ format.
              type: File
              inputBinding:
                position: 2
                shellQuote: false
              
            - id: sample_id
              label: Sample ID
              doc: |-
                Sample ID for the input FASTQs. This will be used as the name of the output directory.
              type: string
            - id: output_prefix
              label: Output Prefix
              doc: Optional prefix for output directory and files.
              type: string?

            outputs:
            - id: hlahd_final_results
              type: File
              outputBinding:
                glob: |-
                  ${
                      if (!inputs.output_prefix) {
                          return inputs.sample_id + "_hlahd/" + inputs.sample_id + "/result/" + inputs.sample_id + "_final.result.txt"
                      } else {
                          return inputs.sample_id + "_" + inputs.output_prefix + "_hlahd" + "/" + inputs.sample_id + "_" + inputs.output_prefix + "/result/" + inputs.sample_id + "_" + inputs.output_prefix + "_final.result.txt"
                      }
                  }

            baseCommand: []
            arguments:
            - prefix: -f
              position: 2
              valueFrom: /app/hlahd.1.4.0/freq_data
              shellQuote: false
            - prefix: ''
              position: 3
              valueFrom: /app/hlahd.1.4.0/HLA_gene.split.txt
              separate: false
              shellQuote: false
            - prefix: ''
              position: 3
              valueFrom: /app/hlahd.1.4.0/dictionary
              separate: false
              shellQuote: false
            - prefix: ''
              position: 1
              valueFrom: hlahd.sh
              separate: false
              shellQuote: false
            - prefix: ''
              position: 5
              valueFrom: ./
              shellQuote: false
            - prefix: ''
              position: 6
              valueFrom: |-
                ${
                    if (!inputs.output_prefix) {
                        return "&& mv ./" + inputs.sample_id + " ./" + inputs.sample_id + "_hlahd"
                    } else {
                        return "&& mv ./" + inputs.sample_id + "_" + inputs.output_prefix + " ./" + inputs.sample_id + "_" + inputs.output_prefix + "_hlahd"
                    }
                }
              separate: false
              shellQuote: false
            - prefix: ''
              position: 0
              valueFrom: |-
                ${
                    if (!inputs.output_prefix) {
                        return "mkdir ./" + inputs.sample_id + "_hlahd &&"
                    } else {
                        return "mkdir ./" + inputs.sample_id + "_" + inputs.output_prefix + "_hlahd" + " &&"
                    }
                }
              shellQuote: false
            - prefix: ''
              position: 4
              valueFrom: |-
                ${
                    if (!inputs.output_prefix) {
                        return inputs.sample_id
                    } else {
                        return inputs.sample_id + "_" + inputs.output_prefix
                    }
                }
              shellQuote: false
            id: hla-hd/3
            
          out:
          - id: hlahd_final_results
          
          
        id: restricted-reads-hlahd-simplified/6
        
      when: |-
        ${
            if (inputs.read2_sequences == null || inputs.read1_sequences == null) {
                return false;
            } else {
                return true;
            }
        }
      out:
      - id: hlahd_final
      
      
    - id: tumour_DNA_hlahd
      label: tumour-DNA-hlahd
      in:
      - id: sample_name
        source: patient_id
      - id: output_prefix
        default: tumour_dna
      - id: read2_sequences
        source: tumour_DNA_read2_sequences
      - id: read1_sequences
        source: tumour_DNA_read1_sequences
      - id: bowtie2_index
        source: bowtie2_index
      - id: threads
        source: alignment_threads
      - id: bowtie2_index_prefix
        source: bowtie2_index_prefix
      run:
        cwlVersion: v1.2
        class: Workflow
        label: restricted-reads-hlahd-simplified
        doc: |-
          # About this workflow
          This workflow runs fast HLA typing by restricting reads to those aligning with the HLA region prior to running HLA-HD.  
          Aligning the input FASTQs to a reference containing only the HLA region, then converting the mapped reads back to FASTQs allows us to restrict the reads used by HLA-HD to only those which are critically relevant to HLA typing. While it is optimal to provide all reads to HLA-HD, the runtime is prohibitively long. In our testing, this full read restriction method (including HLA-HD) has been shown to take approximately two thirds of the time taken to run HLA-HD alone on unrestricted reads.

          ## Before running this workflow
          Prior to running this workflow, you must create/download a Bowtie2 Index (generated using `bowtie2 build`) for a reference file **which only includes the HLA region on chromosome 6**. We recommend using the HLA region reference `hla_gen.fasta` provided by IMGT. It can be downloaded from the IMGT GitHub repo [here](https://github.com/ANHIG/IMGTHLA/tree/Latest/fasta) or using this download link [ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta](ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta). The Bowtie2 Index files must then be archived using the `tar` command before being used as the Bowtie2 Index Archive input to this workflow.  
          **A Bowtie2 index archive for `hla_gen.fasta` can be downloaded from the disTIL repo to be used as an input for this workflow. The corresponding Bowtie2 Index Prefix parameter should be `hla_gen`.**

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
        - class: InlineJavascriptRequirement
        - class: StepInputExpressionRequirement

        inputs:
        - id: sample_name
          label: Patient ID
          doc: Patient ID to be used for naming the output SAM.
          type: string
          
          
        - id: output_prefix
          label: HLA-HD Output Prefix
          doc: Optional prefix for HLA-HD output files and directory.
          type: string?
          
          
        - id: read2_sequences
          label: Read 2 Sequences
          doc: Read 2 sequences in FASTA or FASTQ format (may be bgzipped).
          type: File
          
          
          
        - id: read1_sequences
          label: Read 1 Sequences
          doc: Read 1 sequences in FASTA or FASTQ format (may be bgzipped).
          type: File
          
          
          
        - id: bowtie2_index
          label: Bowtie2 Index Archive
          doc: A TAR archive containing Bowtie2 index files.
          type: File
          
          
          
        - id: threads
          type: int?
          
          
        - id: bowtie2_index_prefix
          label: Bowtie2 Index Prefix
          doc: |-
            The prefix of index files contained in the Bowtie2 index TAR. Note that all Bowtie2 nidex files in the TAR should have this prefix.
          type: string
          

        outputs:
        - id: hlahd_final
          label: HLA-HD Final Results File
          doc: The final results text file produced by HLA-HD.
          type: File
          outputSource:
          - hla_hd/hlahd_final_results
          
          
          

        steps:
        - id: samtools_view
          label: samtools-view
          in:
          - id: output_format
            default: BAM
          - id: input_alignment
            source: bowtie3/aligned_sam
          - id: fast_bam_compression
            default: true
          - id: exclude_reads_any
            default: '4'
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            label: samtools-view
            doc: |-
              # About this tool
              This app runs samtools view on an input alignment.

              ## Docker
              This CWL tool uses the Docker image `rachelbj/samtools:1.10.0` which contains:
              - htslib v1.10.2
              - bcftools v1.10.2
              - samtools v1.10

              ## Documentation
              - [samtools](http://www.htslib.org/doc/samtools.html)
            $namespaces:
              sbg: https://www.sevenbridges.com/

            requirements:
            - class: ShellCommandRequirement
            - class: DockerRequirement
              dockerPull: rachelbj/samtools:1.10.0
            - class: InlineJavascriptRequirement
              expressionLib:
              - |2-

                var setMetadata = function(file, metadata) {
                    if (!('metadata' in file))
                        file['metadata'] = metadata;
                    else {
                        for (var key in metadata) {
                            file['metadata'][key] = metadata[key];
                        }
                    }
                    return file
                };

                var inheritMetadata = function(o1, o2) {
                    var commonMetadata = {};
                    if (!Array.isArray(o2)) {
                        o2 = [o2]
                    }
                    for (var i = 0; i < o2.length; i++) {
                        var example = o2[i]['metadata'];
                        for (var key in example) {
                            if (i == 0)
                                commonMetadata[key] = example[key];
                            else {
                                if (!(commonMetadata[key] == example[key])) {
                                    delete commonMetadata[key]
                                }
                            }
                        }
                    }
                    if (!Array.isArray(o1)) {
                        o1 = setMetadata(o1, commonMetadata)
                    } else {
                        for (var i = 0; i < o1.length; i++) {
                            o1[i] = setMetadata(o1[i], commonMetadata)
                        }
                    }
                    return o1;
                };

            inputs:
            - id: output_format
              label: Output Format
              doc: Specifies BAM, SAM or CRAM output format.
              type:
              - 'null'
              - name: output_format
                type: enum
                symbols:
                - BAM
                - SAM
                - CRAM
              inputBinding:
                prefix: --output-fmt
                position: 0
                shellQuote: false
            - id: input_alignment
              label: Input Alignment
              doc: An input BAM, SAM or CRAM file.
              type: File
              inputBinding:
                position: 10
                shellQuote: false
              
            - id: fast_bam_compression
              label: Fast BAM Compression
              doc: Whether the output file (which must be BAM format) should be compressed.
              type: boolean?
              inputBinding:
                prefix: '-1'
                position: 0
                shellQuote: false
            - id: include_header
              label: Include Header
              doc: Whether the input alignment header should be included in the output
                file.
              type: boolean?
              inputBinding:
                prefix: -h
                position: 0
                shellQuote: false
            - id: header_only
              label: Output Header Only
              doc: |-
                When this option is selected, the output file will only contain the header from the input alignment.
              type: boolean?
              inputBinding:
                prefix: -H
                position: 0
                shellQuote: false
            - id: include_reads
              label: Include reads with all of these flags
              doc: |-
                Only output alignments with all bits set in this integer present in the FLAG field.
              type: string?
              inputBinding:
                prefix: -f
                position: 0
                shellQuote: false
            - id: exclude_reads_any
              label: Exclude reads with any of these flags
              doc: |-
                Do not output alignments with any bits set in this integer present in the FLAG field.
              type: string?
              inputBinding:
                prefix: -F
                position: 0
                shellQuote: false
            - id: exclude_reads_all
              label: Exclude reads with all of these flags
              doc: |-
                Only exclude reads with all of the bits set in this integer present in the FLAG field.
              type: string?
              inputBinding:
                prefix: -G
                position: 0
                shellQuote: false
            - id: bed_file
              label: BED File
              doc: Only output alignments overlapping the regions specified in this
                BED file.
              type: File?
              inputBinding:
                prefix: -L
                position: 0
                shellQuote: false
              

            outputs:
            - id: output_alignment
              label: Output Alignment
              doc: Output from samtools view.
              type: File
              outputBinding:
                glob: |-
                  ${
                      //Find the input basename (without the file extension)
                      var input_split = inputs.input_alignment.path.split('/')
                      var input_base = input_split[input_split.length - 1].split('.').slice(0,-1). join('.')
                      
                      var ext = ""
                      //Determine the output file extension
                      if (inputs.fast_bam_compression || inputs.output_format == "BAM") {
                          ext = ".bam"
                      } else if (inputs.output_format == "SAM") {
                          ext = ".sam"
                      } else if (inputs.output_format == "CRAM") {
                          ext = ".cram"
                      } else {
                          ext = "." + input_split[input_split.length - 1].split('.').slice(-1)
                      }
                      
                      //If only output header then add '.header'
                      if (inputs.header_only) {
                          return input_base + '.header' + ext
                      }
                      
                      //If filtered on flags/bed file then add '.filtered'
                      if (inputs.include_reads || inputs.exclude_reads_any || inputs.exclude_reads_all || inputs.bed_file) {
                          return input_base + '.filtered' + ext
                      }
                      
                      return input_base + ext
                  }
                outputEval: $(inheritMetadata(self, inputs.input_alignment))
              

            baseCommand:
            - samtools
            - view
            arguments:
            - prefix: -o
              position: 0
              valueFrom: |-
                ${
                    //Find the input basename (without the file extension)
                    var input_split = inputs.input_alignment.path.split('/')
                    var input_base = input_split[input_split.length - 1].split('.').slice(0,-1). join('.')
                    
                    var ext = ""
                    //Determine the output file extension
                    if (inputs.fast_bam_compression || inputs.output_format == "BAM") {
                        ext = ".bam"
                    } else if (inputs.output_format == "SAM") {
                        ext = ".sam"
                    } else if (inputs.output_format == "CRAM") {
                        ext = ".cram"
                    } else {
                        ext = "." + input_split[input_split.length - 1].split('.').slice(-1)
                    }
                    
                    //If only output header then add '.header'
                    if (inputs.header_only) {
                        return input_base + '.header' + ext
                    }
                    
                    //If filtered on flags/bed file then add '.filtered'
                    if (inputs.include_reads || inputs.exclude_reads_any || inputs.exclude_reads_all || inputs.bed_file) {
                        return input_base + '.filtered' + ext
                    }
                    
                    return input_base + ext
                }
              shellQuote: false
            id: samtools_view
            
          out:
          - id: output_alignment
          
          
        - id: samtools_fastq_1
          label: samtools-fastq
          in:
          - id: input_alignment
            source: samtools_view/output_alignment
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            label: samtools-fastq
            doc: |-
              # About this tool
              This app runs samtools fastq on an input alignment.  
              Note that the output paired-end FASTQs have the 'paired_end' metadata field populated with '1' or '2'.

              ## Docker
              This CWL tool uses the Docker image `rachelbj/samtools:1.10.0` which contains:
              - htslib v1.10.2
              - bcftools v1.10.2
              - samtools v1.10

              ## Documentation
              - [samtools](http://www.htslib.org/doc/samtools.html)
            $namespaces:
              sbg: https://www.sevenbridges.com/

            requirements:
            - class: ShellCommandRequirement
            - class: DockerRequirement
              dockerPull: rachelbj/samtools:1.10.0
            - class: InlineJavascriptRequirement

            inputs:
            - id: input_alignment
              label: Input Alignment
              doc: An input BAM, SAM or CRAM file.
              type: File
              inputBinding:
                position: 10
                shellQuote: false
              

            outputs:
            - id: output_fastq_1
              label: Paired-End FASTQ 1
              doc: Paired-end FASTQ 1 output by samtools fastq.
              type: File
              outputBinding:
                glob: |-
                  ${
                      var input_split = inputs.input_alignment.path.split('/')
                      var input_base = input_split[input_split.length - 1]

                      return input_base + '.pe_1.fastq'
                  }
                outputEval: |-
                  ${
                    var out = self[0];
                    out.metadata = {'paired_end' : ''}
                    out.metadata['paired_end'] = 1;
                    
                    return out
                  }
              
            - id: output_fastq_2
              label: Paired-End FASTQ 2
              doc: Paired-end FASTQ 2 output by samtools fastq.
              type: File
              outputBinding:
                glob: |-
                  ${
                      var input_split = inputs.input_alignment.path.split('/')
                      var input_base = input_split[input_split.length - 1]

                      return input_base + '.pe_2.fastq'
                  }
                outputEval: |-
                  ${
                    var out = self[0];
                    out.metadata = {'paired_end' : ''}
                    out.metadata['paired_end'] = 1;
                    
                    return out
                  }
              

            baseCommand:
            - samtools
            - fastq
            arguments:
            - prefix: '-1'
              position: 0
              valueFrom: |-
                ${
                    var input_split = inputs.input_alignment.path.split('/')
                    var input_base = input_split[input_split.length - 1]

                    return input_base + '.pe_1.fastq'
                }
              shellQuote: false
            - prefix: '-2'
              position: 0
              valueFrom: |-
                ${
                    var input_split = inputs.input_alignment.path.split('/')
                    var input_base = input_split[input_split.length - 1]

                    return input_base + '.pe_2.fastq'
                }
              shellQuote: false
            id: samtools_fastq
            
          out:
          - id: output_fastq_1
          - id: output_fastq_2
          
          
        - id: bowtie3
          label: bowtie2
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
          - id: threads
            source: threads
          run:
            cwlVersion: v1.2
            class: CommandLineTool
            label: bowtie2
            doc: |-
              # About this tool
              This tool runs Bowtie2 (**v2.4.1**) alignment of reads to a reference sequence.

              ## Inputs and parameters
              - **Bowtie2 Index Archive:** an archive (TAR) of index files generated from a reference genome (using Bowtie2 build command).
              - **Bowtie2 Index Prefix:** the basename of the index files in the Bowtie2 Index Archive. The basename should be shared by all index files in the TAR, and is the name of the index files up to but not including the final `.1.bt2` etc.
              - **Read 1 Sequences:** file containing mate 1 reads (for paired-end sequencing).
              - **Read 2 Sequences:** file containing mate 2 reads (for paired-end sequencing).
              - **Sample Name:** the name of the sample being analysed (used as the name for the output SAM file).

              ## Docker
              This tool uses the Docker image: `biocontainers/bowtie2:v2.4.1_cv1`.

              ## Documentation
              - [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
            $namespaces:
              sbg: https://sevenbridges.com

            requirements:
            - class: ShellCommandRequirement
            - class: ResourceRequirement
              coresMin: $(inputs.threads)
            - class: DockerRequirement
              dockerPull: biocontainers/bowtie2:v2.4.1_cv1
            - class: InlineJavascriptRequirement
              expressionLib:
              - |2-

                var setMetadata = function(file, metadata) {
                    if (!('metadata' in file))
                        file['metadata'] = metadata;
                    else {
                        for (var key in metadata) {
                            file['metadata'][key] = metadata[key];
                        }
                    }
                    return file
                };

                var inheritMetadata = function(o1, o2) {
                    var commonMetadata = {};
                    if (!Array.isArray(o2)) {
                        o2 = [o2]
                    }
                    for (var i = 0; i < o2.length; i++) {
                        var example = o2[i]['metadata'];
                        for (var key in example) {
                            if (i == 0)
                                commonMetadata[key] = example[key];
                            else {
                                if (!(commonMetadata[key] == example[key])) {
                                    delete commonMetadata[key]
                                }
                            }
                        }
                    }
                    if (!Array.isArray(o1)) {
                        o1 = setMetadata(o1, commonMetadata)
                    } else {
                        for (var i = 0; i < o1.length; i++) {
                            o1[i] = setMetadata(o1[i], commonMetadata)
                        }
                    }
                    return o1;
                };

            inputs:
            - id: bowtie2_index
              label: Bowtie2 Index Archive
              doc: A TAR archive containing Bowtie2 index files.
              type: File
              
            - id: read1_sequences
              label: Read 1 Sequences
              doc: Read 1 sequences in FASTA or FASTQ format (may be bgzipped).
              type: File
              inputBinding:
                prefix: '-1'
                position: 2
                shellQuote: false
              
            - id: read2_sequences
              label: Read 2 Sequences
              doc: Read 2 sequences in FASTA or FASTQ format (may be bgzipped).
              type: File
              inputBinding:
                prefix: '-2'
                position: 3
                shellQuote: false
              
            - id: no_unaligned
              label: Suppress Unaligned Reads
              doc: Suppres SAM records for unaligned reads.
              type: boolean?
              inputBinding:
                position: 1
                valueFrom: |-
                  ${
                      if (self == true) {
                          return "--no-unal"
                      } else {
                          return ""
                      }
                  }
                shellQuote: false
            - id: sample_name
              label: Sample Name
              doc: Sample name used to name the output SAM file.
              type: string
            - id: bowtie2_index_prefix
              label: Bowtie2 Index Prefix
              doc: |-
                The prefix of index files contained in the Bowtie2 index TAR. Note that all Bowtie2 nidex files in the TAR should have this prefix.
              type: string
              inputBinding:
                prefix: -x
                position: 1
                shellQuote: false
            - id: threads
              type: int?
              default: 8
              inputBinding:
                prefix: -p
                position: 0
                shellQuote: false

            outputs:
            - id: aligned_sam
              type: File
              outputBinding:
                glob: '*.sam'
                outputEval: $(inheritMetadata(self, inputs.read1_sequences))

            baseCommand: []
            arguments:
            - prefix: -S
              position: 4
              valueFrom: $(inputs.sample_name).sam
              shellQuote: false
            - prefix: ''
              position: 0
              valueFrom: |-
                ${
                    var ref = inputs.bowtie2_index.path
                    return "tar -xvf " + ref + " &&"
                }
              shellQuote: false
            - prefix: ''
              position: 0
              valueFrom: bowtie2
              shellQuote: false
            id: bowtie2/6
            
          out:
          - id: aligned_sam
          
          
        - id: hla_hd
          label: hla-hd
          in:
          - id: threads
            source: threads
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
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            label: hla-hd
            doc: |-
              ## About HLA-HD
              HLA-HD documentation and release notes can be found [here](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/).
              HLA-HD (HLA typing from High-quality Dictionary) can accurately determine HLA alleles with 6-digit precision from NGS data (FASTQ format). RNA-Seq data can also be applied.

              ## About this CWL tool
              This CWL tool runs HLA-HD version 1.4.0 on paired-end FASTQ files to determine HLA type.

              ### Inputs and parameters
              - The input paired-end read files can be from **WGS/WES or RNA-seq**.
              - The input paired-end read files must be in FASTQ format (**not zipped**).
              - The default minimum read length is 100, however this is often too strict. Choose a lower threshold to include more reads.

              ### Output
              - HLA-HD results are output to a directory named using the input sample id.
              - The final summary of HLA typing results can be found at the following path: `<output_dir_name>/result/<sample_id>_final.result.txt`.

              ### Other notes
              - This tool uses the HLA dictionary created from release 3.15.0 of the [IPD-IMGT/HLA](https://github.com/ANHIG/IMGTHLA) database.
              - This tool by default uses HLA allele frequency data included with the HLA-HD release 1.4.0.
            $namespaces:
              sbg: https://sevenbridges.com

            requirements:
            - class: ShellCommandRequirement
            - class: LoadListingRequirement
            - class: DockerRequirement
              dockerPull: pgc-images.sbgenomics.com/rbowen_james/hla-hd
            - class: InlineJavascriptRequirement

            inputs:
            - id: threads
              label: Threads
              doc: Number of cores used to execute the program.
              type: int?
              inputBinding:
                prefix: -t
                position: 1
                shellQuote: false
            - id: minimum_read_length
              label: Minimum read length
              doc: A read whose length is shorter than this parameter is ignored.
              type: int?
              inputBinding:
                prefix: -m
                position: 1
                shellQuote: false
            - id: fastq_reads1
              label: FASTQ Reads 1
              doc: Paired-end reads 1 in FASTQ format.
              type: File
              inputBinding:
                position: 2
                shellQuote: false
              
            - id: fastq_reads2
              label: FASTQ Reads 2
              doc: Paired-end reads 2 in FASTQ format.
              type: File
              inputBinding:
                position: 2
                shellQuote: false
              
            - id: sample_id
              label: Sample ID
              doc: |-
                Sample ID for the input FASTQs. This will be used as the name of the output directory.
              type: string
            - id: output_prefix
              label: Output Prefix
              doc: Optional prefix for output directory and files.
              type: string?

            outputs:
            - id: hlahd_final_results
              type: File
              outputBinding:
                glob: |-
                  ${
                      if (!inputs.output_prefix) {
                          return inputs.sample_id + "_hlahd/" + inputs.sample_id + "/result/" + inputs.sample_id + "_final.result.txt"
                      } else {
                          return inputs.sample_id + "_" + inputs.output_prefix + "_hlahd" + "/" + inputs.sample_id + "_" + inputs.output_prefix + "/result/" + inputs.sample_id + "_" + inputs.output_prefix + "_final.result.txt"
                      }
                  }

            baseCommand: []
            arguments:
            - prefix: -f
              position: 2
              valueFrom: /app/hlahd.1.4.0/freq_data
              shellQuote: false
            - prefix: ''
              position: 3
              valueFrom: /app/hlahd.1.4.0/HLA_gene.split.txt
              separate: false
              shellQuote: false
            - prefix: ''
              position: 3
              valueFrom: /app/hlahd.1.4.0/dictionary
              separate: false
              shellQuote: false
            - prefix: ''
              position: 1
              valueFrom: hlahd.sh
              separate: false
              shellQuote: false
            - prefix: ''
              position: 5
              valueFrom: ./
              shellQuote: false
            - prefix: ''
              position: 6
              valueFrom: |-
                ${
                    if (!inputs.output_prefix) {
                        return "&& mv ./" + inputs.sample_id + " ./" + inputs.sample_id + "_hlahd"
                    } else {
                        return "&& mv ./" + inputs.sample_id + "_" + inputs.output_prefix + " ./" + inputs.sample_id + "_" + inputs.output_prefix + "_hlahd"
                    }
                }
              separate: false
              shellQuote: false
            - prefix: ''
              position: 0
              valueFrom: |-
                ${
                    if (!inputs.output_prefix) {
                        return "mkdir ./" + inputs.sample_id + "_hlahd &&"
                    } else {
                        return "mkdir ./" + inputs.sample_id + "_" + inputs.output_prefix + "_hlahd" + " &&"
                    }
                }
              shellQuote: false
            - prefix: ''
              position: 4
              valueFrom: |-
                ${
                    if (!inputs.output_prefix) {
                        return inputs.sample_id
                    } else {
                        return inputs.sample_id + "_" + inputs.output_prefix
                    }
                }
              shellQuote: false
            id: hla-hd/3
            
        id: restricted-reads-hlahd-simplified/6
        
      out:
      - id: hlahd_final
      
      
    - id: normal_DNA_hlahd
      label: normal-DNA-hlahd
      in:
      - id: sample_name
        source: patient_id
      - id: output_prefix
        default: normal_dna
      - id: read2_sequences
        source: Normal_DNA_read2_sequences
      - id: read1_sequences
        source: Normal_DNA_read1_sequences
      - id: bowtie2_index
        source: bowtie2_index
      - id: threads
        source: alignment_threads
      - id: bowtie2_index_prefix
        source: bowtie2_index_prefix
      run:
        cwlVersion: v1.2
        class: Workflow
        label: restricted-reads-hlahd-simplified
        doc: |-
          # About this workflow
          This workflow runs fast HLA typing by restricting reads to those aligning with the HLA region prior to running HLA-HD.  
          Aligning the input FASTQs to a reference containing only the HLA region, then converting the mapped reads back to FASTQs allows us to restrict the reads used by HLA-HD to only those which are critically relevant to HLA typing. While it is optimal to provide all reads to HLA-HD, the runtime is prohibitively long. In our testing, this full read restriction method (including HLA-HD) has been shown to take approximately two thirds of the time taken to run HLA-HD alone on unrestricted reads.

          ## Before running this workflow
          Prior to running this workflow, you must create/download a Bowtie2 Index (generated using `bowtie2 build`) for a reference file **which only includes the HLA region on chromosome 6**. We recommend using the HLA region reference `hla_gen.fasta` provided by IMGT. It can be downloaded from the IMGT GitHub repo [here](https://github.com/ANHIG/IMGTHLA/tree/Latest/fasta) or using this download link [ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta](ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta). The Bowtie2 Index files must then be archived using the `tar` command before being used as the Bowtie2 Index Archive input to this workflow.  
          **A Bowtie2 index archive for `hla_gen.fasta` can be downloaded from the disTIL repo to be used as an input for this workflow. The corresponding Bowtie2 Index Prefix parameter should be `hla_gen`.**

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
        - class: InlineJavascriptRequirement
        - class: StepInputExpressionRequirement

        inputs:
        - id: sample_name
          label: Patient ID
          doc: Patient ID to be used for naming the output SAM.
          type: string
          
          
        - id: output_prefix
          label: HLA-HD Output Prefix
          doc: Optional prefix for HLA-HD output files and directory.
          type: string?
          
          
        - id: read2_sequences
          label: Read 2 Sequences
          doc: Read 2 sequences in FASTA or FASTQ format (may be bgzipped).
          type: File
          
          
          
        - id: read1_sequences
          label: Read 1 Sequences
          doc: Read 1 sequences in FASTA or FASTQ format (may be bgzipped).
          type: File
          
          
          
        - id: bowtie2_index
          label: Bowtie2 Index Archive
          doc: A TAR archive containing Bowtie2 index files.
          type: File
          
          
          
        - id: threads
          type: int?
          
          
        - id: bowtie2_index_prefix
          label: Bowtie2 Index Prefix
          doc: |-
            The prefix of index files contained in the Bowtie2 index TAR. Note that all Bowtie2 nidex files in the TAR should have this prefix.
          type: string
          

        outputs:
        - id: hlahd_final
          label: HLA-HD Final Results File
          doc: The final results text file produced by HLA-HD.
          type: File
          outputSource:
          - hla_hd/hlahd_final_results
          
          
          

        steps:
        - id: samtools_view
          label: samtools-view
          in:
          - id: output_format
            default: BAM
          - id: input_alignment
            source: bowtie3/aligned_sam
          - id: fast_bam_compression
            default: true
          - id: exclude_reads_any
            default: '4'
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            label: samtools-view
            doc: |-
              # About this tool
              This app runs samtools view on an input alignment.

              ## Docker
              This CWL tool uses the Docker image `rachelbj/samtools:1.10.0` which contains:
              - htslib v1.10.2
              - bcftools v1.10.2
              - samtools v1.10

              ## Documentation
              - [samtools](http://www.htslib.org/doc/samtools.html)
            $namespaces:
              sbg: https://www.sevenbridges.com/

            requirements:
            - class: ShellCommandRequirement
            - class: DockerRequirement
              dockerPull: rachelbj/samtools:1.10.0
            - class: InlineJavascriptRequirement
              expressionLib:
              - |2-

                var setMetadata = function(file, metadata) {
                    if (!('metadata' in file))
                        file['metadata'] = metadata;
                    else {
                        for (var key in metadata) {
                            file['metadata'][key] = metadata[key];
                        }
                    }
                    return file
                };

                var inheritMetadata = function(o1, o2) {
                    var commonMetadata = {};
                    if (!Array.isArray(o2)) {
                        o2 = [o2]
                    }
                    for (var i = 0; i < o2.length; i++) {
                        var example = o2[i]['metadata'];
                        for (var key in example) {
                            if (i == 0)
                                commonMetadata[key] = example[key];
                            else {
                                if (!(commonMetadata[key] == example[key])) {
                                    delete commonMetadata[key]
                                }
                            }
                        }
                    }
                    if (!Array.isArray(o1)) {
                        o1 = setMetadata(o1, commonMetadata)
                    } else {
                        for (var i = 0; i < o1.length; i++) {
                            o1[i] = setMetadata(o1[i], commonMetadata)
                        }
                    }
                    return o1;
                };

            inputs:
            - id: output_format
              label: Output Format
              doc: Specifies BAM, SAM or CRAM output format.
              type:
              - 'null'
              - name: output_format
                type: enum
                symbols:
                - BAM
                - SAM
                - CRAM
              inputBinding:
                prefix: --output-fmt
                position: 0
                shellQuote: false
            - id: input_alignment
              label: Input Alignment
              doc: An input BAM, SAM or CRAM file.
              type: File
              inputBinding:
                position: 10
                shellQuote: false
              
            - id: fast_bam_compression
              label: Fast BAM Compression
              doc: Whether the output file (which must be BAM format) should be compressed.
              type: boolean?
              inputBinding:
                prefix: '-1'
                position: 0
                shellQuote: false
            - id: include_header
              label: Include Header
              doc: Whether the input alignment header should be included in the output
                file.
              type: boolean?
              inputBinding:
                prefix: -h
                position: 0
                shellQuote: false
            - id: header_only
              label: Output Header Only
              doc: |-
                When this option is selected, the output file will only contain the header from the input alignment.
              type: boolean?
              inputBinding:
                prefix: -H
                position: 0
                shellQuote: false
            - id: include_reads
              label: Include reads with all of these flags
              doc: |-
                Only output alignments with all bits set in this integer present in the FLAG field.
              type: string?
              inputBinding:
                prefix: -f
                position: 0
                shellQuote: false
            - id: exclude_reads_any
              label: Exclude reads with any of these flags
              doc: |-
                Do not output alignments with any bits set in this integer present in the FLAG field.
              type: string?
              inputBinding:
                prefix: -F
                position: 0
                shellQuote: false
            - id: exclude_reads_all
              label: Exclude reads with all of these flags
              doc: |-
                Only exclude reads with all of the bits set in this integer present in the FLAG field.
              type: string?
              inputBinding:
                prefix: -G
                position: 0
                shellQuote: false
            - id: bed_file
              label: BED File
              doc: Only output alignments overlapping the regions specified in this
                BED file.
              type: File?
              inputBinding:
                prefix: -L
                position: 0
                shellQuote: false
              

            outputs:
            - id: output_alignment
              label: Output Alignment
              doc: Output from samtools view.
              type: File
              outputBinding:
                glob: |-
                  ${
                      //Find the input basename (without the file extension)
                      var input_split = inputs.input_alignment.path.split('/')
                      var input_base = input_split[input_split.length - 1].split('.').slice(0,-1). join('.')
                      
                      var ext = ""
                      //Determine the output file extension
                      if (inputs.fast_bam_compression || inputs.output_format == "BAM") {
                          ext = ".bam"
                      } else if (inputs.output_format == "SAM") {
                          ext = ".sam"
                      } else if (inputs.output_format == "CRAM") {
                          ext = ".cram"
                      } else {
                          ext = "." + input_split[input_split.length - 1].split('.').slice(-1)
                      }
                      
                      //If only output header then add '.header'
                      if (inputs.header_only) {
                          return input_base + '.header' + ext
                      }
                      
                      //If filtered on flags/bed file then add '.filtered'
                      if (inputs.include_reads || inputs.exclude_reads_any || inputs.exclude_reads_all || inputs.bed_file) {
                          return input_base + '.filtered' + ext
                      }
                      
                      return input_base + ext
                  }
                outputEval: $(inheritMetadata(self, inputs.input_alignment))
              

            baseCommand:
            - samtools
            - view
            arguments:
            - prefix: -o
              position: 0
              valueFrom: |-
                ${
                    //Find the input basename (without the file extension)
                    var input_split = inputs.input_alignment.path.split('/')
                    var input_base = input_split[input_split.length - 1].split('.').slice(0,-1). join('.')
                    
                    var ext = ""
                    //Determine the output file extension
                    if (inputs.fast_bam_compression || inputs.output_format == "BAM") {
                        ext = ".bam"
                    } else if (inputs.output_format == "SAM") {
                        ext = ".sam"
                    } else if (inputs.output_format == "CRAM") {
                        ext = ".cram"
                    } else {
                        ext = "." + input_split[input_split.length - 1].split('.').slice(-1)
                    }
                    
                    //If only output header then add '.header'
                    if (inputs.header_only) {
                        return input_base + '.header' + ext
                    }
                    
                    //If filtered on flags/bed file then add '.filtered'
                    if (inputs.include_reads || inputs.exclude_reads_any || inputs.exclude_reads_all || inputs.bed_file) {
                        return input_base + '.filtered' + ext
                    }
                    
                    return input_base + ext
                }
              shellQuote: false
            id: samtools_view
            
          out:
          - id: output_alignment
          
          
        - id: samtools_fastq_1
          label: samtools-fastq
          in:
          - id: input_alignment
            source: samtools_view/output_alignment
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            label: samtools-fastq
            doc: |-
              # About this tool
              This app runs samtools fastq on an input alignment.  
              Note that the output paired-end FASTQs have the 'paired_end' metadata field populated with '1' or '2'.

              ## Docker
              This CWL tool uses the Docker image `rachelbj/samtools:1.10.0` which contains:
              - htslib v1.10.2
              - bcftools v1.10.2
              - samtools v1.10

              ## Documentation
              - [samtools](http://www.htslib.org/doc/samtools.html)
            $namespaces:
              sbg: https://www.sevenbridges.com/

            requirements:
            - class: ShellCommandRequirement
            - class: DockerRequirement
              dockerPull: rachelbj/samtools:1.10.0
            - class: InlineJavascriptRequirement

            inputs:
            - id: input_alignment
              label: Input Alignment
              doc: An input BAM, SAM or CRAM file.
              type: File
              inputBinding:
                position: 10
                shellQuote: false
              

            outputs:
            - id: output_fastq_1
              label: Paired-End FASTQ 1
              doc: Paired-end FASTQ 1 output by samtools fastq.
              type: File
              outputBinding:
                glob: |-
                  ${
                      var input_split = inputs.input_alignment.path.split('/')
                      var input_base = input_split[input_split.length - 1]

                      return input_base + '.pe_1.fastq'
                  }
                outputEval: |-
                  ${
                    var out = self[0];
                    out.metadata = {'paired_end' : ''}
                    out.metadata['paired_end'] = 1;
                    
                    return out
                  }
              
            - id: output_fastq_2
              label: Paired-End FASTQ 2
              doc: Paired-end FASTQ 2 output by samtools fastq.
              type: File
              outputBinding:
                glob: |-
                  ${
                      var input_split = inputs.input_alignment.path.split('/')
                      var input_base = input_split[input_split.length - 1]

                      return input_base + '.pe_2.fastq'
                  }
                outputEval: |-
                  ${
                    var out = self[0];
                    out.metadata = {'paired_end' : ''}
                    out.metadata['paired_end'] = 1;
                    
                    return out
                  }
              

            baseCommand:
            - samtools
            - fastq
            arguments:
            - prefix: '-1'
              position: 0
              valueFrom: |-
                ${
                    var input_split = inputs.input_alignment.path.split('/')
                    var input_base = input_split[input_split.length - 1]

                    return input_base + '.pe_1.fastq'
                }
              shellQuote: false
            - prefix: '-2'
              position: 0
              valueFrom: |-
                ${
                    var input_split = inputs.input_alignment.path.split('/')
                    var input_base = input_split[input_split.length - 1]

                    return input_base + '.pe_2.fastq'
                }
              shellQuote: false
            id: samtools_fastq
            
          out:
          - id: output_fastq_1
          - id: output_fastq_2
          
          
        - id: bowtie3
          label: bowtie2
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
          - id: threads
            source: threads
          run:
            cwlVersion: v1.2
            class: CommandLineTool
            label: bowtie2
            doc: |-
              # About this tool
              This tool runs Bowtie2 (**v2.4.1**) alignment of reads to a reference sequence.

              ## Inputs and parameters
              - **Bowtie2 Index Archive:** an archive (TAR) of index files generated from a reference genome (using Bowtie2 build command).
              - **Bowtie2 Index Prefix:** the basename of the index files in the Bowtie2 Index Archive. The basename should be shared by all index files in the TAR, and is the name of the index files up to but not including the final `.1.bt2` etc.
              - **Read 1 Sequences:** file containing mate 1 reads (for paired-end sequencing).
              - **Read 2 Sequences:** file containing mate 2 reads (for paired-end sequencing).
              - **Sample Name:** the name of the sample being analysed (used as the name for the output SAM file).

              ## Docker
              This tool uses the Docker image: `biocontainers/bowtie2:v2.4.1_cv1`.

              ## Documentation
              - [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
            $namespaces:
              sbg: https://sevenbridges.com

            requirements:
            - class: ShellCommandRequirement
            - class: ResourceRequirement
              coresMin: $(inputs.threads)
            - class: DockerRequirement
              dockerPull: biocontainers/bowtie2:v2.4.1_cv1
            - class: InlineJavascriptRequirement
              expressionLib:
              - |2-

                var setMetadata = function(file, metadata) {
                    if (!('metadata' in file))
                        file['metadata'] = metadata;
                    else {
                        for (var key in metadata) {
                            file['metadata'][key] = metadata[key];
                        }
                    }
                    return file
                };

                var inheritMetadata = function(o1, o2) {
                    var commonMetadata = {};
                    if (!Array.isArray(o2)) {
                        o2 = [o2]
                    }
                    for (var i = 0; i < o2.length; i++) {
                        var example = o2[i]['metadata'];
                        for (var key in example) {
                            if (i == 0)
                                commonMetadata[key] = example[key];
                            else {
                                if (!(commonMetadata[key] == example[key])) {
                                    delete commonMetadata[key]
                                }
                            }
                        }
                    }
                    if (!Array.isArray(o1)) {
                        o1 = setMetadata(o1, commonMetadata)
                    } else {
                        for (var i = 0; i < o1.length; i++) {
                            o1[i] = setMetadata(o1[i], commonMetadata)
                        }
                    }
                    return o1;
                };

            inputs:
            - id: bowtie2_index
              label: Bowtie2 Index Archive
              doc: A TAR archive containing Bowtie2 index files.
              type: File
              
            - id: read1_sequences
              label: Read 1 Sequences
              doc: Read 1 sequences in FASTA or FASTQ format (may be bgzipped).
              type: File
              inputBinding:
                prefix: '-1'
                position: 2
                shellQuote: false
              
            - id: read2_sequences
              label: Read 2 Sequences
              doc: Read 2 sequences in FASTA or FASTQ format (may be bgzipped).
              type: File
              inputBinding:
                prefix: '-2'
                position: 3
                shellQuote: false
              
            - id: no_unaligned
              label: Suppress Unaligned Reads
              doc: Suppres SAM records for unaligned reads.
              type: boolean?
              inputBinding:
                position: 1
                valueFrom: |-
                  ${
                      if (self == true) {
                          return "--no-unal"
                      } else {
                          return ""
                      }
                  }
                shellQuote: false
            - id: sample_name
              label: Sample Name
              doc: Sample name used to name the output SAM file.
              type: string
            - id: bowtie2_index_prefix
              label: Bowtie2 Index Prefix
              doc: |-
                The prefix of index files contained in the Bowtie2 index TAR. Note that all Bowtie2 nidex files in the TAR should have this prefix.
              type: string
              inputBinding:
                prefix: -x
                position: 1
                shellQuote: false
            - id: threads
              type: int?
              default: 8
              inputBinding:
                prefix: -p
                position: 0
                shellQuote: false

            outputs:
            - id: aligned_sam
              type: File
              outputBinding:
                glob: '*.sam'
                outputEval: $(inheritMetadata(self, inputs.read1_sequences))

            baseCommand: []
            arguments:
            - prefix: -S
              position: 4
              valueFrom: $(inputs.sample_name).sam
              shellQuote: false
            - prefix: ''
              position: 0
              valueFrom: |-
                ${
                    var ref = inputs.bowtie2_index.path
                    return "tar -xvf " + ref + " &&"
                }
              shellQuote: false
            - prefix: ''
              position: 0
              valueFrom: bowtie2
              shellQuote: false
            id: bowtie2/6
            
          out:
          - id: aligned_sam
          
          
        - id: hla_hd
          label: hla-hd
          in:
          - id: threads
            source: threads
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
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            label: hla-hd
            doc: |-
              ## About HLA-HD
              HLA-HD documentation and release notes can be found [here](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/).
              HLA-HD (HLA typing from High-quality Dictionary) can accurately determine HLA alleles with 6-digit precision from NGS data (FASTQ format). RNA-Seq data can also be applied.

              ## About this CWL tool
              This CWL tool runs HLA-HD version 1.4.0 on paired-end FASTQ files to determine HLA type.

              ### Inputs and parameters
              - The input paired-end read files can be from **WGS/WES or RNA-seq**.
              - The input paired-end read files must be in FASTQ format (**not zipped**).
              - The default minimum read length is 100, however this is often too strict. Choose a lower threshold to include more reads.

              ### Output
              - HLA-HD results are output to a directory named using the input sample id.
              - The final summary of HLA typing results can be found at the following path: `<output_dir_name>/result/<sample_id>_final.result.txt`.

              ### Other notes
              - This tool uses the HLA dictionary created from release 3.15.0 of the [IPD-IMGT/HLA](https://github.com/ANHIG/IMGTHLA) database.
              - This tool by default uses HLA allele frequency data included with the HLA-HD release 1.4.0.
            $namespaces:
              sbg: https://sevenbridges.com

            requirements:
            - class: ShellCommandRequirement
            - class: LoadListingRequirement
            - class: DockerRequirement
              dockerPull: pgc-images.sbgenomics.com/rbowen_james/hla-hd
            - class: InlineJavascriptRequirement

            inputs:
            - id: threads
              label: Threads
              doc: Number of cores used to execute the program.
              type: int?
              inputBinding:
                prefix: -t
                position: 1
                shellQuote: false
            - id: minimum_read_length
              label: Minimum read length
              doc: A read whose length is shorter than this parameter is ignored.
              type: int?
              inputBinding:
                prefix: -m
                position: 1
                shellQuote: false
            - id: fastq_reads1
              label: FASTQ Reads 1
              doc: Paired-end reads 1 in FASTQ format.
              type: File
              inputBinding:
                position: 2
                shellQuote: false
              
            - id: fastq_reads2
              label: FASTQ Reads 2
              doc: Paired-end reads 2 in FASTQ format.
              type: File
              inputBinding:
                position: 2
                shellQuote: false
              
            - id: sample_id
              label: Sample ID
              doc: |-
                Sample ID for the input FASTQs. This will be used as the name of the output directory.
              type: string
            - id: output_prefix
              label: Output Prefix
              doc: Optional prefix for output directory and files.
              type: string?

            outputs:
            - id: hlahd_final_results
              type: File
              outputBinding:
                glob: |-
                  ${
                      if (!inputs.output_prefix) {
                          return inputs.sample_id + "_hlahd/" + inputs.sample_id + "/result/" + inputs.sample_id + "_final.result.txt"
                      } else {
                          return inputs.sample_id + "_" + inputs.output_prefix + "_hlahd" + "/" + inputs.sample_id + "_" + inputs.output_prefix + "/result/" + inputs.sample_id + "_" + inputs.output_prefix + "_final.result.txt"
                      }
                  }

            baseCommand: []
            arguments:
            - prefix: -f
              position: 2
              valueFrom: /app/hlahd.1.4.0/freq_data
              shellQuote: false
            - prefix: ''
              position: 3
              valueFrom: /app/hlahd.1.4.0/HLA_gene.split.txt
              separate: false
              shellQuote: false
            - prefix: ''
              position: 3
              valueFrom: /app/hlahd.1.4.0/dictionary
              separate: false
              shellQuote: false
            - prefix: ''
              position: 1
              valueFrom: hlahd.sh
              separate: false
              shellQuote: false
            - prefix: ''
              position: 5
              valueFrom: ./
              shellQuote: false
            - prefix: ''
              position: 6
              valueFrom: |-
                ${
                    if (!inputs.output_prefix) {
                        return "&& mv ./" + inputs.sample_id + " ./" + inputs.sample_id + "_hlahd"
                    } else {
                        return "&& mv ./" + inputs.sample_id + "_" + inputs.output_prefix + " ./" + inputs.sample_id + "_" + inputs.output_prefix + "_hlahd"
                    }
                }
              separate: false
              shellQuote: false
            - prefix: ''
              position: 0
              valueFrom: |-
                ${
                    if (!inputs.output_prefix) {
                        return "mkdir ./" + inputs.sample_id + "_hlahd &&"
                    } else {
                        return "mkdir ./" + inputs.sample_id + "_" + inputs.output_prefix + "_hlahd" + " &&"
                    }
                }
              shellQuote: false
            - prefix: ''
              position: 4
              valueFrom: |-
                ${
                    if (!inputs.output_prefix) {
                        return inputs.sample_id
                    } else {
                        return inputs.sample_id + "_" + inputs.output_prefix
                    }
                }
              shellQuote: false
            id: hla-hd/3
            
      out:
      - id: hlahd_final
      
      
    id: three-sample-hlatyping-simplified/9
    
  out:
  - id: clin_sig_json
  - id: consensus_txt
  - id: sample1_json
  - id: consensus_json
  - id: sample2_json
  - id: sample3_json
  - id: clin_sig_txt
  
  
id: |-
  consensus-hla/15
