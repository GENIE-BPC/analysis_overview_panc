# analysis_overview_panc

Analysis code for:

> Genomic Actionability and Analysis in Pancreatic Cancer: Insights from the GENIE BPC Multi-Institutional Cohort. *Clinical Cancer Research*, 2026. (under review)

This repository contains the R code used to generate the figures and tables in the manuscript, based on the AACR Project GENIE Biopharma Collaborative (BPC) pancreatic cancer cohort.


## Data Availability

**Data Cohort** 

N = 1,040 patients from the GENIE BPC pancreas cohort (v1.0-public) with eligible tumor histology types were included. Clinical characteristics, survival outcomes, and genomic data from each patient's first next-generation sequencing (NGS) report were analyzed.

The GENIE BPC data are not included in this repository. To obtain the data:

1. Register for a free [Synapse](https://www.synapse.org/) account and accept the AACR Project GENIE terms of access.
2. Download the PANC v1.0-public release programmatically with the [*{genieBPC}*](https://genie-bpc.github.io/genieBPC/) R package:

   ```r
   library(genieBPC)
   set_synapse_credentials()   # uses your Synapse personal access token
   panc_data <- pull_data_synapse(cohort = "PANC", version = "v1.0-public")
   ```

3. Place the downloaded files in the `data/` folder.

**Germline Data**

Germline mutations in the homologous recombination deficiency (HRD) genes *BRCA1*, *BRCA2*, and *PALB2* were identified in two steps: first through OncoKB germline pathogenicity annotation, and then through manual review of the calculated variant allele frequency (VAF) for each pathogenic variant. Genetic mutation data are publicly available through [cBioPortal for GENIE](https://genie.cbioportal.org/) and can be accessed programmatically using the [*{cbioportalR}*](https://www.karissawhiting.com/cbioportalR/) R package.

## Software Requirements

Analyses were run in **R 4.4.1**.

## Contact

Joshua S. Jolissaint, MD, MSc, Division of Surgical Oncology, University of Maryland School of Medicine. [jjolissaint@som.umaryland.edu]

## License

This code is released under the [MIT License](LICENSE). 