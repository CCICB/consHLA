# consHLA
![overall workflow](/assests/figures/consHLA_overall.png)
<br>

## Running the workflow
### RAM Requirements (indicative only)
RAM depend on input file size and whether subsampling is turned on. <br>
For WGS results with 30x coverage: min RAM = 2Gb <br>
For WGS results with 100x coverage and subsample to 30x: min RAM = 30Gb <br>

### Cloud platforms (recommended)
You can run the workflow on any cloud platform supporting CWL execution (i.e. [Cavatica](https://cavatica.sbgenomics.com/))
<br>

### Local 
Running a `.cwl` workflow requires specific software. Here we pick `cwltool`. Install it following these [instructions](https://github.com/common-workflow-language/cwltool). `cwltool` usage is shown below <br>
```
cwltool cwl/consHLA.cwl sample_input.yml
```
You can run the whole or part of the consHLA workflow by specifing the `.cwl` file and supplying the correct `input.yml`


## Output files 
`*_sample1_hla.json`: HLA alleles typed from tumour WGS <br>
`*_sample2_hla.json`: HLA alleles typed from germline WGS <br>
`*_sample3_hla.json`: HLA alleles typed from tumour RNAseq <br>
`*_threeSample_hla.consensus.clinSig.[json|txt]`: Consensus HLA alleles for clinically significant genes <br>
`*_threeSample_hla.consensus.[json|txt]`: Consensus HLA alleles for all genes <br>


## Test samples
Publicly available NGS data for two cell lines COLO829 and HCC1954 were used to demonstrate consHLA functionality. Download the files to validate consHLA installation. The expected output is provided in `./sample_output` 
- COLO829 tumour WGS [link](https://trace.ncbi.nlm.nih.gov/Traces/sra?run=DRR260182)
- COLO829 germline WGS [link](https://trace.ncbi.nlm.nih.gov/Traces/sra?run=DRR260183) 
- COLO829 tumour RNAseq [link](https://www.ncbi.nlm.nih.gov/sra/SRX5414783)
- HCC1954 tumour WGS [link](https://trace.ncbi.nlm.nih.gov/Traces/sra?run=DRR260184)
- HCC1954 germline WGS [link](https://trace.ncbi.nlm.nih.gov/Traces/sra?run=DRR260185)

## Runtime
Runtime tested with 30x WGS and RNAseq with 180M reads on amazon cloud computing EC2 instance model c5.4xlarge with 16 CPUs, 32Gb of RAM, and 1024Gb of attached storage
![Runtime analysis](/assests/figures/runtime_with_total.png)

## Funding

## Development

## LICENSE

consHLA is a wrapper on HLA-HD and is protected by MIT open source software license. For commercial use of consHLA, please contact the author of [HLA-HD](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/) to obtain a commercial license.  
