import pandas as pd
import os
import streamlit as st

import plotly.express as px
data_path = "data/LivSB/"


dataset = st.sidebar.selectbox(
    "Select dataset",
    [x for x in os.listdir(data_path) if x.endswith(".csv")],
    index=0,
)

df = pd.read_csv(f"{data_path}/{dataset}")

st.dataframe(df)

st.write(df["decision"].value_counts())

# add plot of decision by review_id
# not stacked
# fig = px.histogram(df, x="review_id", color="decision")
# st.plotly_chart(fig)
fig = px.histogram(df, x="review_id", color="decision", barmode="group")
st.plotly_chart(fig)


