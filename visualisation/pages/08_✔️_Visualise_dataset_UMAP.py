import pandas as pd
import streamlit as st
import umap.plot
from datasets import load_dataset
from sklearn.feature_extraction.text import TfidfVectorizer

from visualisation.create_data_card import load_helper


def parse_metrics(metadata, st=None):
    for k, m in metadata.items():
        mattrs = m.__dict__
        for m, attr in mattrs.items():
            if type(attr) == int and attr > 0:
                st.metric(label=f"{k}-{m}", value=attr)


if __name__ == "__main__":

    conhelps = load_helper()
    configs_set = set()

    for conhelper in conhelps:
        configs_set.add(conhelper.dataset_name)

    st.set_page_config(layout="wide")
    s = st.session_state

    if not s:
        s.pressed_first_button = False
    data_name = st.sidebar.selectbox("dataset", configs_set)
    st.sidebar.write("you selected:", data_name)

    # setup data configs
    data_helpers = conhelps.for_dataset(data_name)
    data_configs = [d.config for d in data_helpers]
    data_config_names = [d.config.name for d in data_helpers]
    data_config_name = st.sidebar.selectbox("config", set(data_config_names))

    if st.sidebar.button("fetch") or s.pressed_first_button:
        s.pressed_first_button = True
        helper = conhelps.for_config_name(data_config_name)
        metadata_helper = helper.get_metadata()
        st.header(f"Dataset stats for {data_config_name}")

        parse_metrics(metadata_helper, st.sidebar)

        # load HF dataset
        data_idx = data_config_names.index(data_config_name)
        data_config = data_configs[data_idx]
        dataset = load_dataset(
            f"csmed/datasets/{data_name}/{data_name}.py",
            name=data_config_name,
        )

        category_labels = []
        for split in dataset:
            for item in dataset[split]:
                category_labels.append(f"{split}_{item['labels'][0]}")

        hover_df = pd.DataFrame(category_labels, columns=["category"])
        text = []
        for split in dataset:
            for item in dataset[split]:
                text.append(item["text"])

        hover_df["text"] = text
        st.write(f"Converting {len(text)} examples using TF-IDF vectoriser.")

        tfidf_vectorizer = TfidfVectorizer(
            min_df=0.01, stop_words="english", ngram_range=(1, 1), max_features=10000
        )
        tfidf_word_doc_matrix = tfidf_vectorizer.fit_transform(text)
        st.write(
            f"Created TF-IDF matrix shape: {tfidf_word_doc_matrix.shape}. Now reducing to 2 dimensions."
        )

        embedding = umap.UMAP(n_components=2, metric="hellinger").fit(
            tfidf_word_doc_matrix
        )

        fig = umap.plot.interactive(
            embedding,
            labels=hover_df["category"],
            hover_data=hover_df,
            point_size=5,
            interactive_text_search=True,
            cmap="BrBG",
        )
        st.bokeh_chart(fig)
