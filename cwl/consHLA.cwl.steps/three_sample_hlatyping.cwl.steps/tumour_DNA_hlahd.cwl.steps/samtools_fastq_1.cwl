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
