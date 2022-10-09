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

requirements:
- class: ShellCommandRequirement
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
- id: hlahd_results
  label: Output directory
  doc: Directory containing results of the HLA-HD run.
  type: Directory
  outputBinding:
    glob: |-
      ${
          if (!inputs.output_prefix) {
              return inputs.sample_id + "_hlahd"
          } else {
              return inputs.sample_id + "_" + inputs.output_prefix + "_hlahd"
          }
      }
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
  shellQuote: true
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
id: hla_hd
