# *CSMeD:* Citation Screening Meta-Dataset for systematic review automation evaluation

This package serves as basis for the paper: _"CSMeD: Bridging the Dataset Gap in Automated Citation Screening for Systematic Literature Reviews"_ by Wojciech Kusa, Oscar E. Mendoza, Matthias Samwald, Petr Knoth, Allan Hanbury

[![https://openreview.net/forum?id=ZbmS3MU25p](http://img.shields.io/badge/NeurIPS_2023-CSMeD-7514ae.svg)](https://openreview.net/forum?id=ZbmS3MU25p) 

Table of Contents

- [1. CSMeD: Citation screening datasets for title and abstract screening](#csmed)
- [2. CSMeD-FT: Full-text screening dataset](#csmed-ft)
- [3. Installation](#installation)
- [4. Visualisations](#visualisations)
- [5. Experiments](#experiments)

____

/ add link from table of contents to the sections

## <a name="csmed" /> 1. CSMeD: Citation screening datasets for title and abstract screening


|        | Introduced in                                                                             | # reviews | Domain           | Avg. size | Avg. ratio of included (TA) | Avg. ratio of included (FT) | Additional data | Data URL                                                                                                     | Cochrane | Publicly available | Included in **CSMeD** |
|-------:|:------------------------------------------------------------------------------------------|----------:|:-----------------|----------:|----------------------------:|----------------------------:|-----------------|--------------------------------------------------------------------------------------------------------------|----------|--------------------|-----------------------|
|      1 | [Cohen et al. (2006)](https://doi.org/10.1197/jamia.M1929)                                |        15 | Drug             |     1,249 |                        7.7% |                           — | —               | [Web](https://dmice.ohsu.edu/cohenaa/systematic-drug-class-review-data.html)                                 | —        | ✓                  | ✓                     |
|      2 | [Wallace et al. (2010)](https://doi.org/10.1145/1835804.1835829)                          |         3 | Clinical         |     3,456 |                        7.9% |                           — | —               | [GiitHub](https://github.com/bwallace/citation-screening)                                                    | —        | ✓                  | ✓                     |
|      3 | [Howard et al. (2015)](https://doi.org/10.1186/s13643-016-0263-z)                         |         5 | Mixed            |    19,271 |                        4.6% |                           — | —               | [Supplementary](https://systematicreviewsjournal.biomedcentral.com/articles/10.1186/s13643-016-0263-z#Sec30) | —        | ✓                  | ✓                     |
|      4 | [Miwa et al. (2015)](https://doi.org/10.1016/j.jbi.2014.06.005)                           |         4 | Social science   |     8,933 |                        6.4% |                           — | —               | —                                                                                                            | —        | —                  | —                     |
|      5 | [Scells et al. (2017)](https://dl.acm.org/doi/10.1145/3077136.3080707)                    |        93 | Clinical         |     1,159 |                        1.2% |                           — | Search queries  | [GitHub](https://github.com/ielab/SIGIR2017-SysRev-Collection)                                               | ✓        | ✓                  | ✓                     |
|      6 | [CLEF TAR 2017](https://ceur-ws.org/Vol-1866/invited_paper_12.pdf)                        |        50 | DTA              |     5,339 |                        4.4% |                           — | Review protocol | [GitHub](https://github.com/CLEF-TAR/tar/tree/master/2017-TAR)                                               | ✓        | ✓                  | ✓                     |
|      7 | [CLEF TAR 2018](https://ceur-ws.org/Vol-2125/invited_paper_6.pdf)                         |        30 | DTA              |     7,283 |                        4.7% |                           — | Review protocol | [GitHub](https://github.com/CLEF-TAR/tar/tree/master/2018-TAR)                                               | ✓        | ✓                  | ✓                     |
|      8 | [CLEF TAR 2019](https://ceur-ws.org/Vol-2380/paper_250.pdf)                               |        49 | Mixed**          |     2,659 |                        8.9% |                           — | Review protocol | [GitHub](https://github.com/CLEF-TAR/tar/tree/master/2019-TAR)                                               | ✓        | ✓                  | ✓                     |
|      9 | [Alharbi et al. (2019)](https://dl.acm.org/doi/10.1145/3331184.3331358)                   |        25 | Clinical         |     4,402 |                        0.4% |                           — | Review updates  | [GitHub](https://github.com/Amal-Alharbi/Systematic_Reviews_Update)                                          | ✓        | ✓                  | ✓                     |
|     10 | [Parmar (2021)](https://keep.lib.asu.edu/_flysystem/fedora/c7/Parmar_asu_0010N_21179.pdf) |         6 | Biomedical       |     3,019 |                       21.6% |                        7.3% | —               | —                                                                                                            | —        | —                  | —                     |
|     11 | [Hannousse et al. (2022)](https://doi.org/10.1007/978-3-031-04112-9_15)                   |         7 | Computer Science |       340 |                       11.7% |                           — | Review protocol | [GitHub](https://github.com/hannousse/Semantic-Scholar-Evaluation)                                           | —        | ✓                  | ✓                     |
|     12 | [Wang et al. (2022)](https://doi.org/10.1145/3477495.3531748)                             |        40 | Clinical         |     1,326 |                           — |                           — | Review protocol | [GitHub](https://github.com/ielab/sysrev-seed-collection)                                                    | —        | ✓                  | —                     |

** CLEF TAR 2019 contains 38 reviews of interventions, 8 DTA, 1 Prognosis and 2 Qualitative systematic reviews.

TA stands for Title + Abstract screening phase, FT for Full-text screening phase.
`Avg. size` describes the size of a review in terms of the number records retrieved from the search
query. `Avg. ratio of included (TA)` describes the average ratio of included records in the TA
phase. `Avg. ratio of included (FT)` describes the average ratio of included records in the FT phase.

##  <a name="csmed-ft" /> 2. CSMeD-FT: Full-text screening dataset

| Dataset name        | #reviews | #docs. | #included | %included | Avg. #words in document | Avg. #words in review |
|---------------------|----------|--------|-----------|-----------|-------------------------|-----------------------|
| CSMeD-FT-train      | 148      | 2,053  | 904       | 44.0%     | 4,535                   | 1,493                 |
| CSMeD-FT-dev        | 36       | 644    | 202       | 31.4%     | 4,419                   | 1,402                 |
| CSMeD-FT-test       | 29       | 636    | 278       | 43.7%     | 4,957                   | 2,318                 |
| CSMeD-FT-test-small | 16       | 50     | 22        | 44.0%     | 5,042                   | 2,354                 |

*Column '#docs' refers to the total number of documents included in the dataset and '#included' mentions number of
included documents on the full-text step. CSMeD-test-small is a subset of
CSMeD-test.*

## <a name="installation" /> 3. Installation

### Requirements

Assuming you have `conda` installed, to create environment for loading CSMeD run:

```zsh
$ conda create -n csmed python=3.10
$ conda activate csmed
(csmed)$ pip install -r requirements.txt
```


### Data acquisition prerequisites

To obtain the datasets, you need to configure the following:

- Get a cookie from https://www.cochranelibrary.com/

Furthermore, to obtain full-text PDFs, you need to configure the following:

1. SemanticScholar API key: https://www.semanticscholar.org/product/api
2. CORE API key: https://core.ac.uk/services/api
3. GROBID: https://grobid.readthedocs.io/en/latest/Install-Grobid/

If you have all the prerequisites, run:

```zsh
(csmed)$ python confgure.py
```

And follow the prompts providing API keys, cookies, email address to use PubMed Entrez and paths to GROBID server.
You don't need to provide all the information, the bare minimum to construct the datasets is the cookie from Cochrane
Library and the email address for PubMed Entrez.

### Downloading datasets

First install additional requirements:

```zsh
(csmed)$ pip install -r dev-requirements.txt
```

To download the datasets, run:

```zsh
(csmed)$ python scripts/prepare_prospective_dataset.py
```

## <a name="visualisations" /> 4. Visualisations

To run visualisations first you need to install additional requirements:

```zsh
(csmed)$ pip install -r vis-requirements.txt
```

Then you can run the visualisations using streamlit:

```zsh
(csmed)$ streamlit run visualisation/_🏠_Home.py.py
```


## <a name="experiments" /> 5. Experiments

Baseline experiments from the paper are described in the at: [WojciechKusa/CSMeD-baselines](https://github.com/WojciechKusa/CSMeD-baselines) repository.