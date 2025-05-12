import streamlit as st
import pandas as pd
import requests
import networkx as nx
from pyvis.network import Network
import streamlit.components.v1 as components
import py3Dmol
import matplotlib.pyplot as plt
import numpy as np
import tempfile
import plotly.express as px
import seaborn as sns

# --- PAGE CONFIG ---
st.set_page_config(page_title="NEUROGEN PPI", layout="wide")

# --- Inject MADEVoyager Font ---
st.markdown("""
    <style>
    @font-face {
        font-family: 'MADEVoyager';
        src: url('https://raw.githubusercontent.com/MeghanaVaddella/Neurodegenerative_database_codes/main/MADEVoyagerPERSONAL_USE-Bold.otf') format('opentype');
    }
    </style>
""", unsafe_allow_html=True)

# --- Theme Colors (Dark Mode Palette) ---
body_bg = "#C4D8E2"
header_bg = "#3B5875"
header_text_color = "#C4AEAD"
general_text_color = "#001C3D"
button_bg = "#36454F"
button_text_color = "#DBE9F4"
table_bg = "#36454F"
table_border_color = "#5D8AA8"
search_bg = "#A7C7E7"
remaining_bg = "#8BA8B7"

# --- Inject Custom Unified Styling ---
st.markdown(f"""
    <style>
    body, .stApp {{
        background-color: {body_bg};
        color: {general_text_color};
    }}
    .block-container {{
        background-color: {body_bg} !important;
    }}
    .header-text {{
        background-color: {header_bg};
        padding: 1.2rem;
        font-family: 'MADEVoyager', sans-serif;
        font-size: 58px;
        text-align: center;
        margin-top: 0.5em;
        margin-bottom: 0.3em;
        color: {header_text_color};
        letter-spacing: 2px;
        border-radius: 12px;
    }}
    div[data-baseweb="tab-list"] {{
        background-color: {body_bg} !important;
        border-bottom: none !important;
        padding-left: 1rem;
    }}
    button[data-baseweb="tab"] {{
        background-color: {body_bg} !important;
        color: {general_text_color} !important;
        border: none !important;
        font-weight: bold;
        font-size: 18px;
        margin-right: 1.5rem;
    }}
    button[data-baseweb="tab"]:hover {{
        color: #002B5B !important;
        background-color: {body_bg} !important;
        border: none !important;
    }}
    button[data-baseweb="tab"][aria-selected="true"] {{
        border-bottom: 2px solid {general_text_color} !important;
        font-weight: bold;
    }}
    .stButton button, button {{
        background-color: {button_bg} !important;
        color: {button_text_color} !important;
        font-weight: bold;
        border: none;
        border-radius: 8px;
        padding: 0.5rem 1rem;
    }}
    .stTextInput > div > input,
    .stSelectbox > div,
    .stMultiSelect > div,
    .stSlider > div,
    .stNumberInput > div,
    .stTextArea > div > textarea {{
        background-color: {search_bg} !important;
        color: {general_text_color} !important;
        border: 1px solid {general_text_color} !important;
        border-radius: 6px;
    }}
    .stDataFrame, .data-box {{
        background-color: {table_bg} !important;
        color: {general_text_color} !important;
        border-radius: 10px;
        padding: 1rem;
        border: 1px solid {table_border_color};
    }}
    table {{
        color: {general_text_color} !important;
        background-color: {table_bg} !important;
        border: 1px solid {table_border_color} !important;
    }}
    h1, h2, h3, h4, h5, h6 {{
        color: {general_text_color} !important;
        font-family: 'MADEVoyager', sans-serif;
    }}
    </style>
""", unsafe_allow_html=True)

# --- Header ---
st.markdown("<div class='header-text'>NEUROGEN PPI</div>", unsafe_allow_html=True)

# ---- LOAD DATA FUNCTIONS ----
@st.cache_data(show_spinner=False)
def load_ppi_data():
    url = "https://raw.githubusercontent.com/MeghanaVaddella/my-cv-dataset/refs/heads/main/my-cv-data.csv"
    try:
        return pd.read_csv(url)
    except Exception as e:
        st.error(f"Error loading PPI data: {e}")
        return pd.DataFrame()

@st.cache_data(show_spinner=False)
def load_3d_data():
    urls = [
        "https://raw.githubusercontent.com/MeghanaVaddella/Neurodegenerative_Database/refs/heads/main/3D%20Structure-1.csv",
        "https://raw.githubusercontent.com/MeghanaVaddella/Neurodegenerative_Database/refs/heads/main/3D%20Structure-2.csv"
    ]
    try:
        dfs = [pd.read_csv(url) for url in urls]
        return pd.concat(dfs, ignore_index=True)
    except Exception as e:
        st.error(f"Error loading 3D structure data: {e}")
        return pd.DataFrame()

@st.cache_data(show_spinner=False)
def load_no_3d_data():
    url = "https://raw.githubusercontent.com/MeghanaVaddella/Neurodegenerative_Database/refs/heads/main/No%203D%20Structure.csv"
    try:
        return pd.read_csv(url)
    except Exception as e:
        st.error(f"Error loading No 3D structure data: {e}")
        return pd.DataFrame()

@st.cache_data(show_spinner=False)
def load_disease_data():
    url = "https://raw.githubusercontent.com/MeghanaVaddella/Neurodegenerative_database_codes/refs/heads/main/disease%20data.txt"
    try:
        response = requests.get(url)
        lines = response.text.splitlines()
        return [line.strip() for line in lines if line.strip()][1:]
    except Exception as e:
        st.error(f"Error loading disease text: {e}")
        return []

# ---- LOAD ALL ----
ppi_df = load_ppi_data()
df_3d = load_3d_data()
no_structure_df = load_no_3d_data()
disease_text = load_disease_data()

# ---- TABS ----
tabs = st.tabs([
    "Home", 
    "Data", 
    "3D Structure Data", 
    "3D Visualizer", 
    "Data Visualizer",
    "GitHub Edit"
])

# ---- HOME TAB ----
with tabs[0]:
    keywords = [
        "Alzheimer's Disease", 
        "Parkinson's Disease", 
        "Amyotrophic Lateral Sclerosis (ALS)", 
        "Multiple Sclerosis (MS)", 
        "Friedreich’s Ataxia (FA)"
    ]

    for paragraph in disease_text:
        for keyword in keywords:
            if keyword in paragraph:
                paragraph = paragraph.replace(
                    keyword, 
                    f"<span style='color:#d62728; font-weight:bold; font-size:25px;'>{keyword}</span>"
                )

        replacements = {
            "PROTEIN-PROTEIN INTERACTIONS (PPI)": "<span style='font-weight:bold; color:#800020; font-size:30px;'>PROTEIN-PROTEIN INTERACTIONS (PPI)</span>",
            "What are Protein-Protein Interactions (PPIs) ?": "<span style='font-weight:bold; color:#8B0000; font-size:30px;'>What are Protein-Protein Interactions (PPIs) ?</span>",
            "Applications of the Protein-Protein Interactions": "<span style='font-weight:bold; color:#800020; font-size:30px;'>Applications of the Protein-Protein Interactions</span>",
            "1.Drug Discovery:": "<span style='font-weight:bold; color:#8B4513; font-size:25px;'>1.Drug Discovery:</span>",
            "2.Diagnostics:": "<span style='font-weight:bold; color:#8B4513; font-size:25px;'>2.Diagnostics:</span>",
            "3.Synthetic Biology:": "<span style='font-weight:bold; color:#8B4513; font-size:25px;'>3.Synthetic Biology:</span>",
            "4.Functional Genomics:": "<span style='font-weight:bold; color:#8B4513; font-size:25px;'>4.Functional Genomics:</span>",
            "5.Structural Biology:": "<span style='font-weight:bold; color:#8B4513; font-size:25px;'>5.Structural Biology:</span>",
            "6.Systems Biology and Target Validation:": "<span style='font-weight:bold; color:#8B4513; font-size:25px;'>6.Systems Biology and Target Validation:</span>",
            "What does the database Have?": "<span style='font-weight:bold; color:#800020; font-size:30px;'>What does the database Have?</span>",
        }

        for old, new in replacements.items():
            paragraph = paragraph.replace(old, new)

        paragraph = f"<span style='font-size:20px;'>{paragraph}</span>"
        st.markdown(paragraph, unsafe_allow_html=True)

# ---- DATA TAB ----
with tabs[1]:
    st.header("Protein-Protein Interaction Data")
    st.dataframe(ppi_df, use_container_width=True, hide_index=True)
    st.download_button("Download PPI CSV", ppi_df.to_csv(index=False), "PPI_data.csv", "text/csv")

    st.subheader("Visualize Interactions")
    selected_protein = st.selectbox("Choose Protein", pd.unique(ppi_df[['Protein A', 'Protein B']].values.ravel('K')))

    def build_ppi_graph(protein, df):
        G = nx.Graph()
        edges = df[(df['Protein A'] == protein) | (df['Protein B'] == protein)]
        for _, row in edges.iterrows():
            G.add_edge(row['Protein A'], row['Protein B'])
        return G

    if selected_protein:
        G = build_ppi_graph(selected_protein, ppi_df)
        net = Network(height="600px", width="100%", bgcolor="#FFFFFF", font_color="#000000")
        net.from_nx(G)
        path = tempfile.NamedTemporaryFile(delete=False, suffix=".html").name
        net.show(path)
        components.html(open(path, 'r').read(), height=600, scrolling=True)
# ---- 3D STRUCTURE DATA TAB ----
with tabs[2]:
    st.header("3D Structure Available Proteins")
    st.dataframe(df_3d, use_container_width=True, hide_index=True)
    st.download_button("Download 3D Structure Data", df_3d.to_csv(index=False), "3D_structure_data.csv", "text/csv")

    st.header("Proteins Without 3D Structure")
    st.dataframe(no_structure_df, use_container_width=True, hide_index=True)
    st.download_button("Download No 3D Structure Data", no_structure_df.to_csv(index=False), "No_3D_structure_data.csv", "text/csv")

# ---- 3D VISUALIZER TAB ----
with tabs[3]:
    st.header("Protein 3D Structure Viewer")

    protein_ids = df_3d['UniProt_ID'].dropna().unique()
    selected_protein_id = st.selectbox("Select Protein (UniProt ID)", protein_ids)

    if selected_protein_id:
        pdb_id = df_3d[df_3d['UniProt_ID'] == selected_protein_id]['PDB_ID'].values[0]
        st.subheader(f"3D Structure for {selected_protein_id} (PDB: {pdb_id})")

        viewer = py3Dmol.view(query='pdb:' + pdb_id)
        viewer.setStyle({'cartoon': {'color': 'spectrum'}})
        viewer.setBackgroundColor('white')
        viewer.zoomTo()
        viewer_html = viewer.render().replace("\n", "")
        components.html(viewer_html, height=600)

# ---- DATA VISUALIZER TAB ----
with tabs[4]:
    st.header("Data Visualizer")

    st.subheader("Protein Interaction Count")
    all_proteins = pd.concat([ppi_df['Protein A'], ppi_df['Protein B']])
    protein_counts = all_proteins.value_counts().reset_index()
    protein_counts.columns = ['Protein', 'Interaction Count']

    fig_bar = px.bar(protein_counts.head(20), x='Protein', y='Interaction Count',
                     title="Top 20 Proteins by Interaction Count", color='Interaction Count')
    st.plotly_chart(fig_bar, use_container_width=True)

    st.subheader("Interaction Matrix Heatmap")
    pivot = pd.crosstab(ppi_df['Protein A'], ppi_df['Protein B'])
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(pivot, cmap="coolwarm", linewidths=0.5, ax=ax)
    st.pyplot(fig)

# ---- GITHUB EDIT TAB ----
with tabs[5]:
    st.header("🛠️ GitHub Edit Zone")

    st.markdown("""
      USE THIS SECTION TO ACCESS AND EDIT THE DATASETS DIRECTLY FROM THE GITHUB
    """) 

    github_links = {
        "PPI Data (CSV)": "https://github.com/MeghanaVaddella/my-cv-dataset/blob/main/my-cv-data.csv",
        "3D Structure Data 1": "https://github.com/MeghanaVaddella/Neurodegenerative_Database/blob/main/3D%20Structure-1.csv",
        "3D Structure Data 2": "https://github.com/MeghanaVaddella/Neurodegenerative_Database/blob/main/3D%20Structure-2.csv",
        "No 3D Structure (CSV)": "https://github.com/MeghanaVaddella/Neurodegenerative_Database/blob/main/No%203D%20Structure.csv"
    }

    for label, url in github_links.items():
        st.markdown(f"- 🔗 **[{label}]({url})**")

    st.markdown("""
    📢 **CHANGES IN THE GITHUB WILL BE REFLECTED IN THE APP WHEN THE PAGE IS RELOADED!!**
    """)   
