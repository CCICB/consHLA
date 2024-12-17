cwlVersion: v1.0
class: CommandLineTool
label: hla-report
doc: |-
  # About this tool

  This tool generates an HLA report from the output of disTIL HLA consensus HLA typing.


requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: alanwucci/conshla_report:2.1.0
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
       if (inputs.rna_hla) {
          var cmd = 'Rscript -e "rmarkdown::render(\'hla_report_generator.Rmd\',params=list(germline_hla_json=\'' + inputs.germline_hla.path + '\', tumour_hla_json=\'' + inputs.tumour_hla.path + '\', rna_hla_json=\'' + inputs.rna_hla.path + '\', clin_hla_json=\'' + inputs.clin_sig_hla.path + '\', pid=\'' + inputs.patient_id + '\'), output_file=paste(\'' + inputs.patient_id + '\', \'_hlaReport.pdf\', sep=\'\'))\"'
          return cmd
      } else {
          var cmd = 'Rscript -e "rmarkdown::render(\'hla_report_generator.Rmd\',params=list(germline_hla_json=\'' + inputs.germline_hla.path + '\', tumour_hla_json=\'' + inputs.tumour_hla.path + '\', clin_hla_json=\'' + inputs.clin_sig_hla.path + '\', pid=\'' + inputs.patient_id + '\'), output_file=paste(\'' + inputs.patient_id + '\', \'_hlaReport.pdf\', sep=\'\'))\"'
          return cmd
      }
    }
  shellQuote: false
- prefix: ''
  position: 10
  valueFrom: 1>&2
  shellQuote: false
- prefix: ''
  position: 0
  valueFrom: ' cp /hla_report_generator.Rmd . &&'
  shellQuote: false