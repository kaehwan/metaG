# Nextflow pipeline for metagenomics-based analysis of Group B *Streptococcus*

## System Requirements
- **Disk space**: Approximately 100 GB or above.
- **CPU**: 16 cores or above. 
- **Memory**: 256 GB or above.
- **Dependencies**: Nextflow version 22.04.5 and Docker version 20.10.12.
> [!NOTE]
> Users are encouraged to modify the parameters of this pipeline in the nextflow.config file based on the specs of their computing resources.

## Usage
To start the pipeline, use the following command:
```
nextflow run main.nf -profile docker
```


The essential parameters to be set include:

**--input "/path/to/sampleList.csv"**
- Specify the directory of fastq.gz file(s) in CSV file format.
- CSV headers for Illumina reads: "sample_id,short_read_1,short_read_2",
- CSV headers for Nanopore reads: "sample_id,long_read".

**--outdir "/path/to/output/folder"**
- Specify the output directory. [Default: "./results/"]

**--mode "mapping"**
- Specify the analysis mode: mapping, assembly, nanopore. [Default: "mapping"]

**--metaphlan4_db "/path/to/metaphlan_db"**
- Specify the location of unarchived MetaPhlAn4 database.

**--kraken2_db "/path/to/kraken2_db/"**
- Specify the location of unarchived Kraken2 database.

## Publication
For additional implementation details and guidance on using this pipeline, users can also refer to:
- K. H. Sim, *et. al.*. **A metagenomics-based workflow for detection and genomic characterization of GBS in raw freshwater fish**.