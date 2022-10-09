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
  doc: Whether the input alignment header should be included in the output file.
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
  doc: Only output alignments overlapping the regions specified in this BED file.
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
