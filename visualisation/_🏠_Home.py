import streamlit as st

st.title(
    "This is the visualisation for the CSMeD meta-dataset. Use the sidebar for navigation."
)

st.write(
    "Some of the data is not available in the repository, but can be downloaded from the original source."
    "These pages will not work on the deployed version of the app, but can be run locally."
)

st.write(
    "Pages available in the deployed version of the app are marked with a ✔️ emoji."
)

datasets_table = """
|        | Introduced in                                                                             | # reviews | Domain           | Avg. size | Avg. ratio of included (TA) | Avg. ratio of included (FT) | Additional data | Data URL                                                                                                     | Cochrane | Publicly available | Included in **CSMeD** |
|-------:|:------------------------------------------------------------------------------------------|----------:|:-----------------|----------:|----------------------------:|----------------------------:|-----------------|--------------------------------------------------------------------------------------------------------------|----------|--------------------|-----------------------|
|      1 | [Cohen et al. (2006)](https://doi.org/10.1197/jamia.M1929)                                |        15 | Drug             |     1,249 |                        7.7% |                           — | —               | [Web](https://dmice.ohsu.edu/cohenaa/systematic-drug-class-review-data.html)                                 | —        | ✓                  | ✓                     |
|      2 | [Wallace et al. (2010)](https://doi.org/10.1145/1835804.1835829)                          |         3 | Clinical         |     3,456 |                        7.9% |                           — | —               | [GiitHub](https://github.com/bwallace/citation-screening)                                                    | —        | ✓                  | ✓                     |
|      3 | [Howard et al. (2015)](https://doi.org/10.1186/s13643-016-0263-z)                         |         5 | Mixed            |    19,271 |                        4.6% |                           — | —               | [Supplementary](https://systematicreviewsjournal.biomedcentral.com/articles/10.1186/s13643-016-0263-z#Sec30) | —        | ✓                  | ✓                     |
|      4 | [Miwa et al. (2015)](https://doi.org/10.1016/j.jbi.2014.06.005)                           |         4 | Social science   |     8,933 |                        6.4% |                           — | —               | —                                                                                                            | —        | —                  | —                     |
|      5 | [Scells et al. (2017)](https://dl.acm.org/doi/10.1145/3077136.3080707)                    |        93 | Clinical         |     1,159 |                        1.2% |                           — | Search queries  | [GitHub](https://github.com/ielab/SIGIR2017-SysRev-Collection)                                               | ✓        | ✓                  | ✓                     |
|      6 | [CLEF TAR 2017](https://ceur-ws.org/Vol-1866/invited_paper_12.pdf)                        |        50 | DTA              |     5,339 |                        4.4% |                           — | Review protocol | [GitHub](https://github.com/CLEF-TAR/tar/tree/master/2017-TAR)                                               | ✓        | ✓                  | ✓                     |
|      7 | [CLEF TAR 2018](https://ceur-ws.org/Vol-2125/invited_paper_6.pdf)                         |        30 | DTA              |     7,283 |                        4.7% |                           — | Review protocol | [GitHub](https://github.com/CLEF-TAR/tar/tree/master/2018-TAR)                                               | ✓        | ✓                  | ✓                     |
|      8 | [CLEF TAR 2019](https://ceur-ws.org/Vol-2380/paper_250.pdf)                               |        49 | Mixed          |     2,659 |                        8.9% |                           — | Review protocol | [GitHub](https://github.com/CLEF-TAR/tar/tree/master/2019-TAR)                                               | ✓        | ✓                  | ✓                     |
|      9 | [Alharbi et al. (2019)](https://dl.acm.org/doi/10.1145/3331184.3331358)                   |        25 | Clinical         |     4,402 |                        0.4% |                           — | Review updates  | [GitHub](https://github.com/Amal-Alharbi/Systematic_Reviews_Update)                                          | ✓        | ✓                  | ✓                     |
|     10 | [Parmar (2021)](https://keep.lib.asu.edu/_flysystem/fedora/c7/Parmar_asu_0010N_21179.pdf) |         6 | Biomedical       |     3,019 |                       21.6% |                        7.3% | —               | —                                                                                                            | —        | —                  | —                     |
|     11 | [Hannousse et al. (2022)](https://doi.org/10.1007/978-3-031-04112-9_15)                   |         7 | Computer Science |       340 |                       11.7% |                           — | Review protocol | [GitHub](https://github.com/hannousse/Semantic-Scholar-Evaluation)                                           | —        | ✓                  | ✓                     |
|     12 | [Wang et al. (2022)](https://doi.org/10.1145/3477495.3531748)                             |        40 | Clinical         |     1,326 |                           — |                           — | Review protocol | [GitHub](https://github.com/ielab/sysrev-seed-collection)                                                    | —        | ✓                  | —                     |"""

st.markdown(datasets_table, unsafe_allow_html=True)
