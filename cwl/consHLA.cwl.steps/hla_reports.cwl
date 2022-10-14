cwlVersion: v1.0
class: CommandLineTool
label: hla-report
doc: |-
  # About this tool

  This tool generates an HLA report from the output of disTIL HLA consensus HLA typing.

requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: alanwucci/conshla_report:1.0.0
- class: InlineJavascriptRequirement

inputs:
- id: full_hla
  label: HLA Consensus JSON
  doc: |-
    A JSON file containing all HLA consensus alleles. All HLA genes typed by HLA-HD are included at two- or three-field accuracies.
  type: File
- id: clin_sig_hla
  label: Clinically Significant HLA Consensus JSON
  doc: |-
    A JSON file containing consensus alleles for the clinically significant classical HLA genes. Only classical Class-I and Class-II genes are included and all alleles are truncated to two-field accuracy.
  type: File
- id: patient_id
  label: Patient ID
  doc: ID for the patient.
  type: string

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
        var cmd = 'Rscript -e "rmarkdown::render(\'hla_report_generator.Rmd\',params=list(full_hla_json=\'' + inputs.full_hla.path + '\', clin_hla_json=\'' + inputs.clin_sig_hla.path + '\', pid=\'' + inputs.patient_id + '\'), output_file=paste(\'' + inputs.patient_id + '\', \'_hlaReport.pdf\', sep=\'\'))\"'
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
  valueFrom: |-
    echo "RES" 1>&2 && head $(inputs.full_hla.path) 1>&2 && echo "CLIN" 1>&2 && head $(inputs.clin_sig_hla.path) 1>&2 &&
  shellQuote: false

