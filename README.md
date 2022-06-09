# consHLA
intro here
<br>

## Running the workflow
### Requirements
inputs, RAM, storage, etc
### Local 
Running a `.cwl` workflow requires specific software. Here we pick `cwltool`. Install it following these [instructions](https://github.com/common-workflow-language/cwltool). `cwltool` usage is shown below where `[tool-or-workflow-description]` is the `.cwl` file and `[input-job-settings]` is a `.json` or `.yml` file specifying the input parameters. <br>
`cwltool [tool-or-workflow-description] [input-job-settings]`

### Cloud platform
how to run it on the cloud i.e. CAVATICA
<br>
## Test samples
Publicly available NGS data for two cell lines COLO829 and HCC1954 were used to demonstrate consHLA functionality. Download the files to validate consHLA installation. The expected output is provided in `./sample_output` 
- COLO829 tumour WGS [link](https://trace.ncbi.nlm.nih.gov/Traces/sra?run=DRR260182)
- COLO829 germline WGS [link](https://trace.ncbi.nlm.nih.gov/Traces/sra?run=DRR260183) 
- COLO829 tumour RNAseq [link](https://www.ncbi.nlm.nih.gov/sra/SRX5414783)
- HCC1954 tumour WGS [link](https://trace.ncbi.nlm.nih.gov/Traces/sra?run=DRR260184)
- HCC1954 germline WGS [link](https://trace.ncbi.nlm.nih.gov/Traces/sra?run=DRR260185)

## Funding

## Development

## LICENSE

consHLA is a wrapper on HLA-HD and is protected by MIT open source software license. For commercial use of consHLA, please contact the author of [HLA-HD](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/) to obtain a commercial license.  
