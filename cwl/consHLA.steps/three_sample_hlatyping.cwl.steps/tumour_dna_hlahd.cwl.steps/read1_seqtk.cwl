cwlVersion: v1.2
class: CommandLineTool
label: seqtk
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: nanozoo/seqtk:latest
- class: InlineJavascriptRequirement
  expressionLib:
  - |2-

    var setMetadata = function(file, metadata) {
        if (!('metadata' in file)) {
            file['metadata'] = {}
        }
        for (var key in metadata) {
            file['metadata'][key] = metadata[key];
        }
        return file
    };
    var inheritMetadata = function(o1, o2) {
        var commonMetadata = {};
        if (!o2) {
            return o1;
        };
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
            for (var key in commonMetadata) {
                if (!(key in example)) {
                    delete commonMetadata[key]
                }
            }
        }
        if (!Array.isArray(o1)) {
            o1 = setMetadata(o1, commonMetadata)
            if (o1.secondaryFiles) {
                o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)
            }
        } else {
            for (var i = 0; i < o1.length; i++) {
                o1[i] = setMetadata(o1[i], commonMetadata)
                if (o1[i].secondaryFiles) {
                    o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)
                }
            }
        }
        return o1;
    };

inputs:
- id: seed
  label: seed
  doc: Seed for value for random subsampling
  type: int
  default: 1
  inputBinding:
    prefix: -s
    position: 0
    shellQuote: false
- id: input_fastq
  type: File
  inputBinding:
    position: 1
    shellQuote: false
  sbg:fileTypes: fastq.gz, FASTQ.GZ
- id: num_reads
  type: int
  inputBinding:
    position: 2
    shellQuote: false

outputs:
- id: subsampled_fastq
  label: subsampled compressed fastq
  type: File
  outputBinding:
    glob: '*subsample.fastq.gz'
    outputEval: $(inheritMetadata(self, inputs.input_fastq))
  sbg:fileTypes: fastq.gz, FASTQ.GZ
stdout: |-
  ${
      var elements = inputs.input_fastq.basename.split('.');
      var filename = elements[0];
      
      return filename + '_subsample.fastq.gz'
  }

baseCommand:
- seqtk
- sample
arguments:
- prefix: ''
  position: 3
  valueFrom: '| gzip -c'
  shellQuote: false
- prefix: ''
  position: 0
  valueFrom: '-2'
  shellQuote: false
sbg:wrapperAuthor: Rachel Bowen-James <rbowen-james@ccia.org.au>, Weilin Wu <wwu@ccia.org.au>