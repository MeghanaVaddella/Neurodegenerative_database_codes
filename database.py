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
from stmol import showmol

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
        net = Network(height="600px", width="100%", directed=False)
        net.from_nx(G)
        net.save_graph("ppi_graph.html")
        return "ppi_graph.html"

    if st.button("Show PPI Network"):
        if not ppi_df.empty:
            file_path = build_ppi_graph(selected_protein, ppi_df)
            components.html(open(file_path, 'r').read(), height=600)
            with open(file_path, "rb") as f:
                st.download_button("Download Network HTML", f, "ppi_network.html", "text/html")
        else:
            st.warning("PPI data is empty.")

# ---- 3D STRUCTURE DATA TAB ----
with tabs[2]:
    st.header("3D Structure Available Proteins")
    st.dataframe(df_3d, use_container_width=True, hide_index=True)
    st.download_button("Download 3D Structure Data", df_3d.to_csv(index=False), "3D_structure_data.csv", "text/csv")

    st.header("Proteins Without 3D Structure")
    st.dataframe(no_structure_df, use_container_width=True, hide_index=True)
    st.download_button("Download No 3D Structure Data", no_structure_df.to_csv(index=False), "No_3D_structure_data.csv", "text/csv")

# ---- 3D VISUALIZER TAB ----
with tabs[3]:  # 3D Visualizer tab
    st.write("### 3D Protein Structure Visualizer")

    # MolStar Viewer using PDB IDs from dataset
    if not df_3d.empty:
        col1, col2 = st.columns(2)

        with col1:
            protein_a_options = df_3d['Protein A'].dropna().unique().tolist()
            selected_protein_a = st.selectbox("🔍 Select Protein A", options=[""] + protein_a_options, key="select_protein_a")

        with col2:
            protein_b_options = df_3d['Protein B'].dropna().unique().tolist()
            selected_protein_b = st.selectbox("🔍 Select Protein B", options=[""] + protein_b_options, key="select_protein_b")

        result_col1, result_col2 = st.columns(2)
        pdb_ids = []

        # Visualize Protein A
        if selected_protein_a:
            protein_a_data = df_3d[df_3d['Protein A'] == selected_protein_a]
            if not protein_a_data.empty:
                row = protein_a_data.iloc[0]
                with result_col1:
                    st.write(f"**🧬 Protein A:** {row['Protein A']}")
                    st.write(f"**UniProt ID A:** {row['UniProtID A']}")
                    pdb_ids_a = row['PDB ID A'].split(", ")
                    if pdb_ids_a[0] != "NA":
                        pdb_links_a = " | ".join([f"[{pdb}](https://www.rcsb.org/structure/{pdb})" for pdb in pdb_ids_a])
                        st.markdown(f"🔗 **PDB IDs A:** {pdb_links_a}", unsafe_allow_html=True)
                        pdb_ids.extend(pdb_ids_a)
            else:
                with result_col1:
                    st.warning("No matching Protein A found.")

        # Visualize Protein B
        if selected_protein_b:
            protein_b_data = df_3d[df_3d['Protein B'] == selected_protein_b]
            if not protein_b_data.empty:
                row = protein_b_data.iloc[0]
                with result_col2:
                    st.write(f"**🧬 Protein B:** {row['Protein B']}")
                    st.write(f"**UniProt ID B:** {row['UniProtID B']}")
                    pdb_ids_b = row['PDB ID B'].split(", ")
                    if pdb_ids_b[0] != "NA":
                        pdb_links_b = " | ".join([f"[{pdb}](https://www.rcsb.org/structure/{pdb})" for pdb in pdb_ids_b])
                        st.markdown(f"🔗 **PDB IDs B:** {pdb_links_b}", unsafe_allow_html=True)
                        pdb_ids.extend(pdb_ids_b)
            else:
                with result_col2:
                    st.warning("No matching Protein B found.")

        # Mol* Viewer
        st.write("### 🧬 Mol* (MolStar) Viewer")
        pdb_ids = list(filter(lambda x: x != "NA", pdb_ids))

        if pdb_ids:
            molstar_url = "https://molstar.org/viewer/?url=" + ",".join([f"https://files.rcsb.org/download/{pdb}.pdb" for pdb in pdb_ids])
            st.components.v1.iframe(molstar_url, width=1000, height=600)
        else:
            st.warning("No valid PDB IDs found for visualization.")

    st.markdown("---")

    # ---- AlphaFold 3D Viewer ----
    st.write("### 💾 Download AlphaFold Predicted Structures")

# Function to fetch AlphaFold PDB
def fetch_alphafold_pdb(uniprot_id):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    response = requests.get(url)
    return response.text if response.status_code == 200 else None

# Protein selection
col3, col4 = st.columns(2)

with col3:
    uniprot_a_options = df_3d['UniProtID A'].dropna().unique().tolist()
    selected_uniprot_a = st.selectbox("🔍 Select UniProt ID A", options=[""] + uniprot_a_options, key="select_uniprot_a")

with col4:
    uniprot_b_options = df_3d['UniProtID B'].dropna().unique().tolist()
    selected_uniprot_b = st.selectbox("🔍 Select UniProt ID B", options=[""] + uniprot_b_options, key="select_uniprot_b")

# Download combined PDB file
if selected_uniprot_a and selected_uniprot_b:
    pdb_a = fetch_alphafold_pdb(selected_uniprot_a)
    pdb_b = fetch_alphafold_pdb(selected_uniprot_b)

    if pdb_a and pdb_b:
        combined_pdb = (
            f"REMARK   Protein A: {selected_uniprot_a}\n{pdb_a}\n"
            f"REMARK   Protein B: {selected_uniprot_b}\n{pdb_b}"
        )

        st.subheader("⬇️ Download Combined AlphaFold Structure")
        st.download_button(
            label="Download PDB File",
            data=combined_pdb,
            file_name=f"{selected_uniprot_a}_{selected_uniprot_b}_combined.pdb",
            mime="chemical/x-pdb"
        )
    else:
        st.error("❌ Failed to fetch one or both AlphaFold PDB files.")

st.markdown("---")

    # ---- AlphaFold-Multimer FASTA Generator ----
st.write("### 🧬 Predict Interactions using AlphaFold-Multimer")

def fetch_sequence(uniprot_id):
    """Fetch protein sequence from UniProt"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    return response.text if response.ok else None

# Select UniProt IDs from dropdowns
fasta_col1, fasta_col2 = st.columns(2)

with fasta_col1:
    fasta_uid1 = st.selectbox(
        "🔍 Select UniProt ID for Protein A (FASTA)",
        options=[""] + uniprot_a_options,
        key="select_fasta_a"
    )

with fasta_col2:
    fasta_uid2 = st.selectbox(
        "🔍 Select UniProt ID for Protein B (FASTA)",
        options=[""] + uniprot_b_options,
        key="select_fasta_b"
    )

# Generate FASTA on button click
if st.button("Generate AlphaFold-Multimer Input (FASTA)"):
    if fasta_uid1 and fasta_uid2:
        seq1 = fetch_sequence(fasta_uid1)
        seq2 = fetch_sequence(fasta_uid2)

        if seq1 and seq2:
            combined_fasta = f"{seq1.strip()}\n{seq2.strip()}"
            st.success("✅ FASTA file generated successfully.")
            st.download_button(
                "⬇️ Download FASTA",
                data=combined_fasta,
                file_name="multimer_input.fasta",
                mime="text/plain"
            )
            st.code(combined_fasta)

            # ColabFold link
            colab_link = "https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb"
            st.markdown(
                f"🔗 **[Open AlphaFold-Multimer in Google Colab]({colab_link})**",
                unsafe_allow_html=True
            )
        else:
            st.error("❌ Error fetching sequences. Check UniProt IDs.")
    else:
        st.warning("⚠️ Please select both UniProt IDs for FASTA generation.")

st.markdown("---")


     # ---- Upload PDB file ----
     st.subheader("📦 Upload Predicted PDB File from AlphaFold")
     pdb_file = st.file_uploader("Upload PDB file", type=["pdb"], key="upload_pdb")

    if pdb_file:
        pdb_str = pdb_file.read().decode("utf-8", errors="replace")
        st.success("✅ PDB uploaded successfully!")

        col1, col2 = st.columns([2, 1])
        with col1:
            st.markdown("### 📄 PDB File Preview")
            st.text_area("PDB File Content", pdb_str, height=500)
            st.download_button(
                label="📥 Download PDB",
                data=pdb_str,
                file_name="uploaded_structure.pdb",
                mime="chemical/x-pdb"
            )

        with col2:
            viewer = py3Dmol.view(width=400, height=300)
            viewer.addModel(pdb_str, "pdb")
            viewer.setStyle({'cartoon': {'color': 'spectrum'}})
            viewer.setBackgroundColor("white")
            viewer.zoomTo()
            components.html(viewer._make_html(), height=350)

# ---- DATA VISUALIZER TAB ----
with tabs[4]:
    st.header("Data Visualizer")

    # Load Data
    github_url = "https://raw.githubusercontent.com/MeghanaVaddella/Neurodegenerative_Database/refs/heads/main/my-cv-data.csv"
    df = pd.read_csv(github_url)

    # Filter for high-confidence interactions
    high_conf = df[df['Combined Score'] > 0.85].copy()
    high_conf = high_conf[['Protein A', 'Protein B', 'Combined Score']].dropna()

    # Unique proteins list
    all_proteins = pd.unique(pd.concat([high_conf['Protein A'], high_conf['Protein B']]))

    # Select protein for heatmap
    selected_protein = st.selectbox("Select a Protein A to visualize its interactions", sorted(all_proteins))

    # Function to create horizontal heatmap
    def create_horizontal_heatmap(protein):
        interactions = high_conf[
            (high_conf['Protein A'] == protein) | 
            (high_conf['Protein B'] == protein)
        ]
        data = []
        for _, row in interactions.iterrows():
            partner = row['Protein B'] if row['Protein A'] == protein else row['Protein A']
            data.append((protein, partner, row['Combined Score']))
        interaction_df = pd.DataFrame(data, columns=['Protein A', 'Protein B', 'Score'])
        pivot_df = interaction_df.pivot_table(
            index='Protein A',
            columns='Protein B',
            values='Score',
            aggfunc='max'
        )
        fig, ax = plt.subplots(figsize=(max(8, len(pivot_df.columns) * 0.6), 3))
        sns.heatmap(
            pivot_df,
            annot=True,
            cmap='rocket_r',
            vmin=0.85,
            vmax=1,
            linewidths=0.5,
            cbar_kws={'label': 'Interaction Confidence'},
            annot_kws={'size': 9},
            ax=ax
        )
        ax.set_title(f"High-Confidence Interactions for {protein}", pad=15)
        ax.set_xlabel("Interaction Partners", labelpad=12)
        ax.set_ylabel("Selected Protein", labelpad=12)
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        st.pyplot(fig)

    # Show heatmap
    create_horizontal_heatmap(selected_protein)

    # ---- Top Interacting Proteins ----
    st.subheader("Top Proteins by Interaction Count & Combined Score")
    protein_stats = (
        pd.concat([df['Protein A'], df['Protein B']])
        .value_counts()
        .reset_index(name='Interaction Count')
        .merge(
            df.groupby('Protein A')['Combined Score'].mean().reset_index(),
            left_on='index', right_on='Protein A', how='left'
        )
        .rename(columns={'Combined Score': 'Avg Combined Score'})
        .head(20)
    )
    fig_bar = px.bar(
        protein_stats,
        x='Interaction Count',
        y='index',
        orientation='h',
        color='Avg Combined Score',
        color_continuous_scale='Viridis',
        title="Top 20 Proteins: Interaction Count vs. Combined Score"
    )
    st.plotly_chart(fig_bar, use_container_width=True)

    # ---- Disease-Specific Analysis ----
    if 'Disease Associated' in df.columns:
        st.subheader("Disease-Specific Analysis")
        disease_counts = df['Disease Associated'].value_counts().reset_index()
        disease_counts.columns = ['Disease Associated', 'count']
        fig_pie = px.pie(
            disease_counts,
            names='Disease Associated',
            values='count',
            title="Distribution by Disease"
        )
        st.plotly_chart(fig_pie, use_container_width=True)

        selected_disease = st.selectbox(
            "Select Disease for Detailed Analysis",
            df['Disease Associated'].dropna().unique()
        )
        filtered_data = df[df['Disease Associated'] == selected_disease]
        st.write(f"Interactions associated with {selected_disease}:")
        st.dataframe(filtered_data[['Protein A', 'Protein B', 'Experimental System', 'Pubmed ID']])

    # ---- Experimental System Distribution ----
    st.subheader("Experimental System Distribution")
    exp_systems = df['Experimental System'].value_counts().reset_index()
    exp_systems.columns = ['Experimental System', 'count']
    fig_exp = px.bar(
        exp_systems,
        x='count',
        y='Experimental System',
        orientation='h',
        title="Types of Experimental Evidence"
    )
    st.plotly_chart(fig_exp, use_container_width=True)

    # ---- Combined Score Analysis ----
    st.subheader("Combined Score Analysis")
    fig_score = px.histogram(
        df,
        x='Combined Score',
        nbins=50,
        title="Distribution of Combined Confidence Scores"
    )
    st.plotly_chart(fig_score, use_container_width=True)

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
