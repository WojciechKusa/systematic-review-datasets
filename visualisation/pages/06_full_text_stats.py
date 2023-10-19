import pandas as pd
import os
import streamlit as st
from xml.etree import ElementTree as ET
from bs4 import BeautifulSoup
import json
import plotly.express as px

data_path = "data/CSMeD/"


def elem_to_text(elem, default=""):
    if elem:
        return elem.getText(separator=" ", strip=True)
    else:
        return default


dataset = st.sidebar.selectbox(
    "Select dataset",
    [x for x in os.listdir(data_path) if x.endswith(".csv")],
    index=0,
)

df = pd.read_csv(f"{data_path}/{dataset}")

st.dataframe(df)

st.write("Number of documents: ", len(df))
st.write("Number of reviews: ", len(df["review_id"].unique()))

st.write(df["decision"].value_counts())

fig = px.histogram(df, x="review_id", color="decision", barmode="group")
st.plotly_chart(fig)

papers = {}
#todo
if dataset == "CSMeD-FT-train.csv":
    ft_path = "data/external/sigir2017/"
elif dataset == "CSMeD-FT-dev.csv":
    ft_path = "data/external/ec2s/"
else:
    ft_path = "data/external/pcs/"



reviews_metadata = {}

for review in df["review_id"].unique():
    review_references_df = pd.read_csv(
        f"{ft_path}/{review}/references.csv"
    )
    with open(f"{ft_path}/{review}/review_details.json") as f:
        review_details = json.load(f)
    review_details['criteria_text'] = " ".join(review_details['criteria'].values())
    reviews_metadata[review] = review_details

    for paper in df[df["review_id"] == review]["document_id"].unique():

        paper_file = f"{ft_path}/{review}/pdfs/{paper}.tei.xml"
        if os.path.exists(paper_file):
            with open(paper_file, "r") as tei:
                soup = BeautifulSoup(tei, "lxml")

            text_tag = soup.find("text", {"xml:lang": "en"})
            main_text = ""
            if text_tag is not None:
                for element in text_tag.find_all(recursive=False):
                    if element.name == "back":
                        break
                    main_text += elem_to_text(element)

            pubmed_id = review_references_df[
                    review_references_df["reference_id"] == paper
                ]["PUBMED_ID"].values[0]
            if pubmed_id is None or pd.isna(pubmed_id):
                pubmed_id = ""
            else:
                pubmed_id = str(int(pubmed_id))

            papers[paper] = {
                "review_id": review,
                "document_id": paper,
                "decision": df[
                    (df["review_id"] == review) & (df["document_id"] == paper)
                ]["decision"].values[0],
                "reason_for_exclusion": df[
                    (df["review_id"] == review) & (df["document_id"] == paper)
                ]["reason_for_exclusion"].values[0],
                "publication_date": elem_to_text(soup.date),
                "doi": review_references_df[
                    review_references_df["reference_id"] == paper
                ]["doi"].values[0],
                'journal': review_references_df[
                    review_references_df["reference_id"] == paper
                ]["journal"].values[0],
                "year": review_references_df[
                    review_references_df["reference_id"] == paper
                ]["year"].values[0],
                "PubMed ID": pubmed_id,
                # pdf links are stored in CORE_PDF_LINK and S2_PDF_LINK columns
                "PDF links": [x for x in [
                    review_references_df[
                        review_references_df["reference_id"] == paper
                    ]["CORE_PDF_LINK"].values[0],
                    review_references_df[
                        review_references_df["reference_id"] == paper
                    ]["S2_PDF_LINK"].values[0],
                ] if x is not None],
                "title": elem_to_text(soup.title),
                "authors": elem_to_text(soup.author),
                "abstract": elem_to_text(soup.abstract),
                "main_text": main_text,
                "citation": review_references_df[
                    review_references_df["reference_id"] == paper
                    ]["citation"].values[0],
            }


papers_df = pd.DataFrame.from_dict(papers, orient="index")

st.write("Number of documents: ", len(papers_df))
st.dataframe(papers_df)

papers_df["main_text_word_count"] = papers_df["main_text"].apply(lambda x: len(x.split()))
papers_df["abstract_word_count"] = papers_df["abstract"].apply(lambda x: len(x.split()))
papers_df["title_word_count"] = papers_df["title"].apply(lambda x: len(x.split()))

# all three in one plot
fig = px.histogram(
    papers_df,
    x=["main_text_word_count", "abstract_word_count", "title_word_count"],
    barmode="overlay",
    color_discrete_sequence=["#636EFA", "#EF553B", "#00CC96"],
    opacity=0.75,
    marginal="box",
)
st.plotly_chart(fig)

if not os.path.exists(f"{data_path}/CSMeD-FT"):
    os.makedirs(f"{data_path}/CSMeD-FT")

st.write("Mean word count in main text: ", papers_df["main_text_word_count"].mean())
st.write("Mean word in paper title+abstract+main text: ", (papers_df["title_word_count"] + papers_df["abstract_word_count"] + papers_df["main_text_word_count"]).mean())
papers_df.to_csv(f"{data_path}/CSMeD-FT/{dataset}", index=False)


reviews_metadata_df = pd.DataFrame.from_dict(reviews_metadata, orient="index")
with open(f"{data_path}/CSMeD-FT/{dataset[:-4]}_reviews_metadata.json", "w") as f:
    json.dump(reviews_metadata, f, indent=4)

# count length of "abstract", "title" and "criteria" columns
reviews_metadata_df["abstract_word_count"] = reviews_metadata_df["abstract"].apply(lambda x: len(x.split()))
reviews_metadata_df["title_word_count"] = reviews_metadata_df["title"].apply(lambda x: len(x.split()))
reviews_metadata_df["criteria_word_count"] = reviews_metadata_df["criteria_text"].apply(lambda x: len(x.split()))
st.dataframe(reviews_metadata_df)

st.write("Mean word count in systematic reviews combined: ", (reviews_metadata_df["abstract_word_count"] + reviews_metadata_df["title_word_count"] + reviews_metadata_df["criteria_word_count"]).mean())

# all three in one plot
fig = px.histogram(
    reviews_metadata_df,
    x=["abstract_word_count", "title_word_count", "criteria_word_count"],
    barmode="overlay",
    color_discrete_sequence=["#636EFA", "#EF553B", "#00CC96"],
    opacity=0.75,
    marginal="box",
)
st.plotly_chart(fig)


combined_df = pd.merge(
    papers_df,
    reviews_metadata_df,
    left_on="review_id",
    right_index=True,
    how="left",
    suffixes=("", "_review"),
)

# measure word count of "main_text", "abstract" and "title", "criteria_text_review", "title_review", "abstract_review" all together
combined_df['word_count'] = combined_df['main_text_word_count'] + combined_df['abstract_word_count'] + combined_df['title_word_count'] + combined_df['criteria_word_count'] + combined_df['title_word_count_review'] + combined_df['abstract_word_count_review']
st.dataframe(combined_df)

# all three in one plot
fig = px.histogram(
    combined_df,
    x=["word_count"],
    barmode="overlay",
    color_discrete_sequence=["#636EFA", "#EF553B", "#00CC96"],
    opacity=0.75,
    marginal="box",
)
st.plotly_chart(fig)

st.write("Mean word count in main text: ", combined_df["word_count"].mean())

