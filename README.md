# *CSMeD:* Citation Screening Meta-Dataset for systematic review automation evaluation

This package serves as basis for the paper: _"CSMeD: Bridging the Dataset Gap in Automated Citation Screening for
Systematic Literature Reviews"_ by Wojciech Kusa, Oscar E. Mendoza, Matthias Samwald, Petr Knoth, Allan Hanbury (2023)

[![https://proceedings.neurips.cc/paper_files/paper/2023/hash/4962a23916103301b27bde29a27642e8-Abstract-Datasets_and_Benchmarks.html](http://img.shields.io/badge/NeurIPS_2023-CSMeD-7514ae.svg)](https://proceedings.neurips.cc/paper_files/paper/2023/hash/4962a23916103301b27bde29a27642e8-Abstract-Datasets_and_Benchmarks.html)

---

Table of Contents

1. [CSMeD: Title and abstract screening datasets](#csmed)
2. [CSMeD-FT: Full-text screening dataset](#csmed-ft)
3. [Installation](#installation)
4. [Examples](#examples)
5. [Visualisations](#visualisations)
6. [Experiments](#experiments)

____

## <a name="csmed" /> 1. CSMeD: Citation screening datasets for title and abstract screening

Original datasets used to create CSMeD are described in the table below:

|    | Introduced in                                                           | # reviews | Domain           | Avg. size | Avg. ratio of included (TA) | Avg. ratio of included (FT) | Additional data | Data URL                                                                                                     | Cochrane | Publicly available | Included in **CSMeD**                              |
|---:|:------------------------------------------------------------------------|----------:|:-----------------|----------:|----------------------------:|----------------------------:|-----------------|--------------------------------------------------------------------------------------------------------------|----------|--------------------|----------------------------------------------------|
|  1 | [Cohen et al. (2006)](https://doi.org/10.1197/jamia.M1929)              |        15 | Drug             |     1,249 |                        7.7% |                           — | —               | [Web](https://dmice.ohsu.edu/cohenaa/systematic-drug-class-review-data.html)                                 | —        | ✓                  | [✓](csmed%2Fdatasets%2Fcohen%2Fcohen.py)           |
|  2 | [Wallace et al. (2010)](https://doi.org/10.1145/1835804.1835829)        |         3 | Clinical         |     3,456 |                        7.9% |                           — | —               | [GiitHub](https://github.com/bwallace/citation-screening)                                                    | —        | ✓                  | [✓](csmed%2Fdatasets%2Fsas2010%2Fsas2010.py)       |
|  3 | [Howard et al. (2015)](https://doi.org/10.1186/s13643-016-0263-z)       |         5 | Mixed            |    19,271 |                        4.6% |                           — | —               | [Supplementary](https://systematicreviewsjournal.biomedcentral.com/articles/10.1186/s13643-016-0263-z#Sec30) | —        | ✓                  | [✓](csmed%2Fdatasets%2Fswift%2Fswift.py)           |
|  4 | [Miwa et al. (2015)](https://doi.org/10.1016/j.jbi.2014.06.005)         |         4 | Social science   |     8,933 |                        6.4% |                           — | —               | —                                                                                                            | —        | —                  | —                                                  |
|  5 | [Scells et al. (2017)](https://dl.acm.org/doi/10.1145/3077136.3080707)  |        93 | Clinical         |     1,159 |                        1.2% |                           — | Search queries  | [GitHub](https://github.com/ielab/SIGIR2017-SysRev-Collection)                                               | ✓        | ✓                  | [✓](csmed%2Fdatasets%2Fsigir2017%2Fsigir2017.py)   |
|  6 | [CLEF TAR 2017](https://ceur-ws.org/Vol-1866/invited_paper_12.pdf)      |        50 | DTA              |     5,339 |                        4.4% |                           — | Review protocol | [GitHub](https://github.com/CLEF-TAR/tar/tree/master/2017-TAR)                                               | ✓        | ✓                  | [✓](csmed%2Fdatasets%2Ftar2017%2Ftar2017.py)       |
|  7 | [CLEF TAR 2018](https://ceur-ws.org/Vol-2125/invited_paper_6.pdf)       |        30 | DTA              |     7,283 |                        4.7% |                           — | Review protocol | [GitHub](https://github.com/CLEF-TAR/tar/tree/master/2018-TAR)                                               | ✓        | ✓                  | [✓](csmed%2Fdatasets%2Ftar2018%2Ftar2018.py)       |
|  8 | [CLEF TAR 2019](https://ceur-ws.org/Vol-2380/paper_250.pdf)             |        49 | Mixed**          |     2,659 |                        8.9% |                           — | Review protocol | [GitHub](https://github.com/CLEF-TAR/tar/tree/master/2019-TAR)                                               | ✓        | ✓                  | [✓](csmed%2Fdatasets%2Ftar2019%2Ftar2019.py)       |
|  9 | [Alharbi et al. (2019)](https://dl.acm.org/doi/10.1145/3331184.3331358) |        25 | Clinical         |     4,402 |                        0.4% |                           — | Review updates  | [GitHub](https://github.com/Amal-Alharbi/Systematic_Reviews_Update)                                          | ✓        | ✓                  | [✓](csmed%2Fdatasets%2Fsr_updates%2Fsr_updates.py) |
| 10 | [Hannousse et al. (2022)](https://doi.org/10.1007/978-3-031-04112-9_15) |         7 | Computer Science |       340 |                       11.7% |                           — | Review protocol | [GitHub](https://github.com/hannousse/Semantic-Scholar-Evaluation)                                           | —        | ✓                  | [✓](csmed%2Fdatasets%2Fcs_reviews%2Fcs_reviews.py) |

TA stands for Title + Abstract screening phase, FT for Full-text screening phase.
`Avg. size` describes the size of a review in terms of the number records retrieved from the search
query. `Avg. ratio of included (TA)` describes the average ratio of included records in the TA
phase. `Avg. ratio of included (FT)` describes the average ratio of included records in the FT phase.

### CSMeD datasets

CSMeD beyond offering unified access to the original datasets, provides a unified meta-dataset containing all the
original datasets. Statistics of the CSMeD datasets are presented in the table below.

| Dataset name                                                            | #reviews | #docs   | #included | Avg. #docs | Avg. %included | Avg. #words in document | 
|-------------------------------------------------------------------------|----------|---------|-----------|------------|----------------|-------------------------| 
| [CSMeD-basic](csmed%2Fdatasets%2Fcsmed_basic%2Fcsmed_basic.py)          |          |         |           |            |                |                         |
| CSMeD-basic-train                                                       | 30       | 128,438 | 7,958     | 4,281      | 9.6%           | 229                     | 
|                                                                         |          |         |           |            |                |                         |
| [CSMeD-cochrane](csmed%2Fdatasets%2Fcsmed_cochrane%2Fcsmed_cochrane.py) |          |         |           |            |                |                         |
| CSMeD-cochrane-train                                                    | 195      | 372,422 | 7,589     | 1,910      | 21.9%          | 180                     | 
| CSMeD-cochrane-dev                                                      | 100      | 229,376 | 4,365     | 2,294      | 20.8%          | 201                     |
|                                                                         |          |         |           |            |                |                         |
| CSMeD-all                                                               | 325      | 730,236 | 19,912    | 2,247      | 20.5%          | 195                     |

## <a name="csmed-ft" /> 2. CSMeD-FT: Full-text screening dataset

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

To obtain the metadata for CSMeD-Cochrane datasets, you need to configure the cookie for the [Cochrane Library](https://www.cochranelibrary.com/) website.

Furthermore, to obtain full-text PDFs for CSMeD-FT, you need to configure the following:

1. SemanticScholar API key: https://www.semanticscholar.org/product/api
2. CORE API key: https://core.ac.uk/services/api
3. GROBID: https://grobid.readthedocs.io/en/latest/Install-Grobid/

If you have all the prerequisites, run:

```zsh
(csmed)$ python confgure.py
```

And follow the prompts providing API keys, cookies, email address to use PubMed Entrez APIs and paths to GROBID server.
You don't need to provide all the information, the bare minimum to construct the datasets is the cookie from Cochrane
Library and the email address for PubMed Entrez.

### Downloading raw full-text datasets

First install additional requirements:

```zsh
(csmed)$ pip install -r dev-requirements.txt
```

To download the datasets, run:

```zsh
(csmed)$ python scripts/prepare_full_texts.py
```

## <a name="examples" /> 4. Examples

Examples presenting how to use the datasets are available in the [notebooks/](notebooks) directory.

## <a name="visualisations" /> 5. Visualisations

To run visualisations first you need to install additional requirements:

```zsh
(csmed)$ pip install -r vis-requirements.txt
```

Then you can run the visualisations using streamlit:

```zsh
(csmed)$ streamlit run visualisation/_🏠_Home.py.py
```

## <a name="experiments" /> 6. Experiments

Baseline experiments from the paper are described in the
at: [WojciechKusa/CSMeD-baselines](https://github.com/WojciechKusa/CSMeD-baselines) repository.