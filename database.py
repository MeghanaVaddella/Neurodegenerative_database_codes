import streamlit as st 
import pandas as pd
import requests
import networkx as nx
from pyvis.network import Network
import streamlit.components.v1 as components
import py3Dmol
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import tempfile

# --- PAGE CONFIG ---
st.set_page_config(
    page_title="NEUROGEN PPI",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# --- THEME CONSTANTS ---
BODY_BG = "#C4D8E2"
HEADER_BG = "#3B5875"
HEADER_TEXT = "#C4AEAD"
TEXT_COLOR = "#001C3D"

# --- CUSTOM CSS (updated for responsive tab font, 2-space gap, and responsive components) ---
st.markdown(f"""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&family=Georgia:wght@400;700&display=swap');

    /* Global Styles */
    .stApp {{
        background-color: {BODY_BG};
        color: {TEXT_COLOR};
        font-family: 'Inter', sans-serif;
    }}

    /* Header Styling */
    .main-header {{
        background-color: {HEADER_BG};
        padding: 3rem;
        text-align: center;
        color: #e2e8f0;
        font-family: 'Georgia', serif;
        font-size: 4rem;
        font-weight: bold;
        letter-spacing: 0.1em;
        text-transform: uppercase;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
        margin-bottom: 2rem;
        border-radius: 12px;
        text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
    }}

    /* Tab Styling - larger & responsive; gap set to two spaces using '2ch' */
    .stTabs [data-baseweb="tab-list"] {{
        gap: 2ch; /* two character-spaces between tabs */
        background-color: {HEADER_BG};
        padding: 12px 16px;
        border-radius: 12px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        align-items: center;
    }}
    .stTabs [data-baseweb="tab"] {{
        height: auto;
        white-space: nowrap;
        background-color: transparent;
        border-radius: 8px;
        color: #cbd5e1; 
        font-weight: 700;
        font-size: 1.6rem; /* base size for laptops/desktops */
        padding: 8px 18px;
        margin-right: 0.4rem;
        line-height: 1.1;
    }}
    /* Selected tab appearance */
    .stTabs [aria-selected="true"] {{
        background-color: {BODY_BG};
        color: {TEXT_COLOR};
        border-bottom: 4px solid {TEXT_COLOR};
        transform: translateY(-2px);
    }}

    /* Responsive font scaling: mobile, tablet, large screens */
    @media (max-width: 600px) {{
        .stTabs [data-baseweb="tab"] {{
            font-size: 1.05rem; /* mobile */
            padding: 6px 10px;
        }}
        .main-header {{
            font-size: 2.2rem;
            padding: 1.5rem;
        }}
        .bubble-card p {{ font-size: 16px !important; line-height:1.5 !important; }}
    }}
    @media (min-width: 601px) and (max-width: 1199px) {{
        .stTabs [data-baseweb="tab"] {{
            font-size: 1.35rem; /* tablets / small laptops */
        }}
        .main-header {{ font-size: 3rem; padding:2rem; }}
    }}
    @media (min-width: 1200px) {{
        .stTabs [data-baseweb="tab"] {{
            font-size: 2rem; /* large desktop */
        }}
    }}

    /* Container/Card Styling */
    .white-card {{
        background-color: white;
        padding: 2rem;
        border-radius: 16px;
        box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1);
        border: 1px solid #e2e8f0;
        margin-bottom: 2rem;
    }}

    /* Home Bubble Styling */
    .bubble-card {{
        background-color: white;
        padding: 2.5rem;
        border-radius: 24px;
        box-shadow: 0 20px 25px -5px rgba(0, 0, 0, 0.1);
        border: 1px solid #f1f5f9;
        margin-bottom: 2rem;
        transition: transform 0.2s;
    }}
    .bubble-card:hover {{
        transform: translateY(-4px);
    }}

    /* Headings */
    h1, h2, h3, h4 {{
        font-family: 'Georgia', serif;
        color: #1e3a8a !important;
    }}

    /* Button Styling */
    .stButton > button {{
        background-color: #3B5875;
        color: white;
        border-radius: 8px;
        font-weight: bold;
        border: none;
        padding: 0.6rem 1.2rem;
        transition: all 0.2s;
        font-size: 1rem;
    }}
    .stButton > button:hover {{
        background-color: #2c435a;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }}

    /* Helper for Colab Button */
    .colab-btn {{
        background-color: white; 
        color: #3B5875; 
        border: 2px solid #3B5875; 
        padding: 10px 20px; 
        border-radius: 8px; 
        font-weight: bold; 
        text-decoration: none; 
        display: inline-block;
        margin-top: 10px;
    }}
    .colab-btn:hover {{
        background-color: #f0f9ff;
    }}

    /* Responsive container for embedded viewers (MolStar / py3Dmol) */
    .viewer-embed {{
        width: 100% !important;
        max-width: 1400px;
        margin: 0 auto;
    }}
    /* In-case the generated py3Dmol HTML has fixed px values, we will also ensure iframe/div shrink */
    .viewer-embed iframe, .viewer-embed div {{
        width: 100% !important;
        height: 60vh !important;
        min-height: 350px;
    }}
    @media (max-width: 600px) {{
        .viewer-embed iframe, .viewer-embed div {{
            height: 48vh !important; /* mobile-friendly height */
        }}
    }}
    </style>
""", unsafe_allow_html=True)

# --- HEADER ---
st.markdown('<div class="main-header">NEUROGEN PPI</div>', unsafe_allow_html=True)

# --- RAW TEXT DATA ---
HOME_TEXT_RAW = """
... (your existing HOME_TEXT_RAW content unchanged) ...
"""

# --- DATA LOADING ---
@st.cache_data(show_spinner=False)
def load_data():
    try:
        ppi_df = pd.read_csv("https://raw.githubusercontent.com/MeghanaVaddella/Neurodegenerative_Database/refs/heads/main/my-cv-data.csv")
        df_3d_1 = pd.read_csv("https://raw.githubusercontent.com/MeghanaVaddella/Neurodegenerative_Database/refs/heads/main/3D%20Structure-1.csv")
        df_3d_2 = pd.read_csv("https://raw.githubusercontent.com/MeghanaVaddella/Neurodegenerative_Database/refs/heads/main/3D%20Structure-2.csv")
        df_3d = pd.concat([df_3d_1, df_3d_2], ignore_index=True)
        no_struct_df = pd.read_csv("https://raw.githubusercontent.com/MeghanaVaddella/Neurodegenerative_Database/refs/heads/main/No%203D%20Structure.csv")
        return ppi_df, df_3d, no_struct_df
    except Exception as e:
        st.error(f"Error loading data: {e}")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

ppi_df, df_3d, no_structure_df = load_data()

# --- TABS ---
tabs = st.tabs([
    "🏠 Home", 
    "📄 Data", 
    "🧬 3D Structure Data", 
    "🔬 3D Visualizer", 
    "📊 Data Visualizer",
    "🛠 GitHub Edit"
])

# ================= HOME TAB =================
with tabs[0]:
    # (Keep your existing HOME rendering code — omitted here for brevity in this snippet)
    st.markdown("<div class='bubble-card'><h2 style='color:#881337; font-size:36px;'>Welcome to NEUROGEN PPI</h2><p style='font-size:18px;'>(Home content renders here — unchanged)</p></div>", unsafe_allow_html=True)

# ================= DATA TAB =================
with tabs[1]:
    with st.container():
        st.header("Protein-Protein Interaction Data")
        search_term = st.text_input("Search Data", placeholder="Type to search...", key="ppi_search")
        display_df = ppi_df
        if search_term:
            display_df = ppi_df[ppi_df.astype(str).apply(lambda x: x.str.contains(search_term, case=False)).any(axis=1)]
        st.dataframe(display_df, use_container_width=True, hide_index=True)
        st.download_button("Download PPI CSV", ppi_df.to_csv(index=False), "PPI_data.csv", "text/csv")
        st.markdown("---")
        st.subheader("Visualize Interactions")
        all_proteins = pd.concat([ppi_df['Protein A'], ppi_df['Protein B']]).unique()
        selected_protein = st.selectbox("Choose Protein", sorted(all_proteins), key="ppi_viz_select")
        if st.button("Generate Network Graph", key="gen_graph_btn"):
            subset = ppi_df[(ppi_df['Protein A'] == selected_protein) | (ppi_df['Protein B'] == selected_protein)]
            if not subset.empty:
                net = Network(height="600px", width="100%", bgcolor="#ffffff", font_color="black")
                net.add_node(selected_protein, color="#dc2626", size=30, title=f"CENTER: {selected_protein}")
                for _, row in subset.iterrows():
                    partner = row['Protein B'] if row['Protein A'] == selected_protein else row['Protein A']
                    score = row['Combined Score']
                    net.add_node(partner, color="#3b82f6", size=20, title=f"{partner}\nScore: {score}")
                    net.add_edge(selected_protein, partner, value=float(score) if isinstance(score, (int, float)) else 1)
                with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp:
                    net.save_graph(tmp.name)
                    with open(tmp.name, 'r', encoding='utf-8') as f:
                        html_bytes = f.read()
                # embed with responsive wrapper class
                components.html(f"<div class='viewer-embed'>{html_bytes}</div>", height=620)
                st.download_button("Download Network HTML", html_bytes, f"ppi_network_{selected_protein}.html", "text/html")
            else:
                st.warning("No interactions found for this protein.")

# ================= 3D STRUCTURE DATA TAB =================
with tabs[2]:
    with st.container():
        st.header("3D Structure Data")
        search_3d = st.text_input("Search 3D Data", placeholder="Type to search...", key="3d_search")
        display_3d = df_3d
        if search_3d:
            display_3d = df_3d[df_3d.astype(str).apply(lambda x: x.str.contains(search_3d, case=False)).any(axis=1)]
        st.dataframe(display_3d, use_container_width=True, hide_index=True)
        st.download_button("Download 3D Structure Data", df_3d.to_csv(index=False), "3D_structure_data.csv", "text/csv")
        st.markdown("---")
        st.subheader("Proteins Without 3D Structure")
        st.dataframe(no_structure_df, use_container_width=True, hide_index=True)
        st.download_button("Download No 3D Structure Data", no_structure_df.to_csv(index=False), "No_3D_structure_data.csv", "text/csv")

# ================= 3D VISUALIZER TAB =================
with tabs[3]:
    st.header("3D Protein Structure Visualizer")
    top_col1, top_col2 = st.columns([1, 1], gap="large")

    with top_col1:
        st.markdown("#### ⚙ Controls")
        unique_prot_a = sorted(df_3d['Protein A'].dropna().unique()) if not df_3d.empty else []
        sel_prot_a = st.selectbox("Select Protein A", [""] + unique_prot_a, key="viz_a")
        avail_prot_b = []
        if sel_prot_a:
            avail_prot_b = sorted(df_3d[df_3d['Protein A'] == sel_prot_a]['Protein B'].unique())
        sel_prot_b = st.selectbox("Select Protein B", [""] + avail_prot_b, key="viz_b", disabled=not sel_prot_a)

    with top_col2:
        st.markdown("#### ℹ Interaction Info")
        if sel_prot_a and not df_3d.empty:
            row_a = df_3d[df_3d['Protein A'] == sel_prot_a].iloc[0]
            st.info(f"*Protein A: {row_a['Protein A']} | **UniProt: {row_a['UniProtID A']} | **PDBs*: {row_a['PDB ID A']}")
        if sel_prot_b and not df_3d.empty:
            row_b = df_3d[(df_3d['Protein A'] == sel_prot_a) & (df_3d['Protein B'] == sel_prot_b)].iloc[0]
            st.success(f"*Protein B: {row_b['Protein B']} | **UniProt: {row_b['UniProtID B']} | **PDBs*: {row_b['PDB ID B']}")

    st.markdown("---")

    # --- Mol* Viewer (Iframe) ---
    st.markdown("### 🧬 Mol* (MolStar) Viewer")
    pdb_ids = []
    if sel_prot_a and not df_3d.empty:
        row_a = df_3d[df_3d['Protein A'] == sel_prot_a].iloc[0]
        ids_a = str(row_a['PDB ID A']).split(',')
        pdb_ids.extend([x.strip() for x in ids_a if x.strip() != 'NA'])
    if sel_prot_b and not df_3d.empty:
        row_b = df_3d[(df_3d['Protein A'] == sel_prot_a) & (df_3d['Protein B'] == sel_prot_b)].iloc[0]
        ids_b = str(row_b['PDB ID B']).split(',')
        pdb_ids.extend([x.strip() for x in ids_b if x.strip() != 'NA'])
    pdb_ids = list(set(pdb_ids))
    if pdb_ids:
        molstar_url = "https://molstar.org/viewer/?url=" + ",".join([f"https://files.rcsb.org/download/{pdb}.pdb" for pdb in pdb_ids])
        # embed inside responsive wrapper
        st.components.v1.iframe(molstar_url, width="100%", height=600)
    else:
        st.warning("No valid PDB IDs found for visualization or no proteins selected.")

    st.markdown("---")

    # --- AlphaFold Download (unchanged) ---
    st.subheader("💾 Download AlphaFold Predicted Structures")
    af_col1, af_col2 = st.columns(2)
    with af_col1:
        uniprot_opts_a = sorted(df_3d['UniProtID A'].dropna().unique()) if not df_3d.empty else []
        af_uni_a = st.selectbox("Select UniProt ID A (AF)", [""] + uniprot_opts_a, key="af_a")
    with af_col2:
        af_uni_b_opts = []
        if af_uni_a and not df_3d.empty:
            af_uni_b_opts = sorted(df_3d[df_3d['UniProtID A'] == af_uni_a]['UniProtID B'].unique())
        af_uni_b = st.selectbox("Select UniProt ID B (AF)", [""] + af_uni_b_opts, key="af_b", disabled=not af_uni_a)

    if st.button("Fetch & Download AlphaFold PDBs"):
        if af_uni_a and af_uni_b:
            def get_af_pdb(uid):
                for v in [4, 3, 2, 1]:
                    url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v{v}.pdb"
                    try:
                        r = requests.get(url, timeout=10)
                        if r.status_code == 200: return r.text
                    except: pass
                return None
            with st.spinner("Fetching AlphaFold structures..."):
                pdb_a_txt = get_af_pdb(af_uni_a)
                pdb_b_txt = get_af_pdb(af_uni_b)
                if pdb_a_txt and pdb_b_txt:
                    combined = f"REMARK Protein A: {af_uni_a}\n{pdb_a_txt}\nTER\nREMARK Protein B: {af_uni_b}\n{pdb_b_txt}\nEND"
                    st.download_button(label="Download Combined PDB", data=combined, file_name=f"AF_{af_uni_a}_{af_uni_b}.pdb", mime="chemical/x-pdb")
                    st.success("Ready for download!")
                elif pdb_a_txt:
                    st.download_button(label=f"Download {af_uni_a}", data=pdb_a_txt, file_name=f"AF_{af_uni_a}.pdb", mime="chemical/x-pdb")
                elif pdb_b_txt:
                    st.download_button(label=f"Download {af_uni_b}", data=pdb_b_txt, file_name=f"AF_{af_uni_b}.pdb", mime="chemical/x-pdb")
                else:
                    pass
        else:
            st.warning("Select both UniProt IDs first.")

    st.markdown("---")

    # --- FASTA Generator (unchanged) ---
    st.subheader("🧬 Predict Interactions using AlphaFold-Multimer")
    fa_col1, fa_col2 = st.columns(2)
    with fa_col1:
        fa_uni_a = st.selectbox("UniProt ID A (FASTA)", [""] + uniprot_opts_a, key="fa_a")
    with fa_col2:
        fa_uni_b_opts = []
        if fa_uni_a and not df_3d.empty:
            fa_uni_b_opts = sorted(df_3d[df_3d['UniProtID A'] == fa_uni_a]['UniProtID B'].unique())
        fa_uni_b = st.selectbox("UniProt ID B (FASTA)", [""] + fa_uni_b_opts, key="fa_b", disabled=not fa_uni_a)

    if st.button("Generate FASTA"):
        if fa_uni_a and fa_uni_b:
            def get_fasta(uid):
                try:
                    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}.fasta", timeout=10)
                    return r.text if r.ok else ""
                except: return ""
            f1 = get_fasta(fa_uni_a)
            f2 = get_fasta(fa_uni_b)
            if f1 and f2:
                combined_fasta = f"{f1}\n{f2}"
                st.text_area("Generated FASTA", value=combined_fasta, height=200)
                st.download_button("Download FASTA", combined_fasta, "multimer.fasta", "text/plain")
                st.markdown("""
                <a href="https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb" target="_blank" class="colab-btn">
                    Open ColabFold 🔗
                </a>
                """, unsafe_allow_html=True)
            else:
                st.error("Failed to fetch sequences.")

    st.markdown("---")

    # --- PDB Upload (Local Viewer using py3Dmol) - UPDATED for responsive display ---
    st.subheader("📦 Upload Predicted PDB File")
    uploaded_file = st.file_uploader("Upload PDB", type=['pdb'])
    if uploaded_file:
        try:
            pdb_content = uploaded_file.read().decode("utf-8")
        except:
            pdb_content = uploaded_file.read().decode("latin-1")

        st.text_area("File Content Preview", pdb_content, height=150)
        st.download_button("Download Uploaded Structure", pdb_content, "uploaded_structure.pdb", "chemical/x-pdb")

        # Create py3Dmol view (we'll extract the HTML and make it responsive)
        view_up = py3Dmol.view(width=800, height=500)
        view_up.addModel(pdb_content, "pdb")
        view_up.setStyle({'cartoon': {'color': 'spectrum'}})
        view_up.zoomTo()
        view_up.setBackgroundColor('white')

        # Generate HTML and then patch widths/heights to use responsive units.
        raw_html = view_up._make_html()

        # Replace common fixed px occurrences (this is a pragmatic approach - py3Dmol HTML varies)
        # We'll wrap in a responsive container and override inline styles where possible.
        responsive_html = f"""
        <div class='viewer-embed'>
            <style>
                .viewer-embed > div, .viewer-embed iframe {{
                    width: 100% !important;
                    height: 60vh !important;
                    min-height: 350px;
                }}
                @media (max-width:600px) {{
                    .viewer-embed > div, .viewer-embed iframe {{
                        height: 48vh !important;
                    }}
                }}
            </style>
            {raw_html}
        </div>
        """

        # Some py3Dmol outputs use hardcoded pixel sizes inside style attributes. Attempt to replace them.
        responsive_html = responsive_html.replace("width: 800px", "width: 100%").replace("height: 500px", "height: 60vh")
        # Embed the final HTML. The 'height' argument is a fallback container height; the internal CSS uses vh for responsiveness.
        components.html(responsive_html, height=700, scrolling=True)

# ================= DATA VISUALIZER TAB =================
with tabs[4]:
    
    # Top Section: Heatmaps (Left) | Network (Right)
    with st.container():
        top_left, top_right = st.columns([1, 2], gap="large")
        
        # --- LEFT: Heatmaps ---
        with top_left:
            # REMOVED white-card
            st.markdown("### 🧬 Interactive Matrix")
            
            all_prots = sorted(pd.concat([ppi_df['Protein A'], ppi_df['Protein B']]).unique())
            default_prots = all_prots[:5] if len(all_prots) > 5 else all_prots
            selected_heatmap_prots = st.multiselect("Add Proteins", all_prots, default=default_prots, key="hm_multi")
            
            if selected_heatmap_prots:
                # Adjacency Matrix
                hm_data = pd.DataFrame(index=selected_heatmap_prots, columns=selected_heatmap_prots).fillna(0)
                for p1 in selected_heatmap_prots:
                    for p2 in selected_heatmap_prots:
                        if p1 == p2:
                            hm_data.loc[p1, p2] = 1000
                        else:
                            match = ppi_df[((ppi_df['Protein A'] == p1) & (ppi_df['Protein B'] == p2)) | 
                                           ((ppi_df['Protein B'] == p1) & (ppi_df['Protein A'] == p2))]
                            if not match.empty:
                                val = match.iloc[0]['Combined Score']
                                hm_data.loc[p1, p2] = float(val) * 1000 if float(val) <= 1 else float(val)

                fig_hm, ax_hm = plt.subplots(figsize=(5, 4))
                sns.heatmap(hm_data.astype(float), cmap="rocket_r", ax=ax_hm, cbar=False, annot=False)
                st.pyplot(fig_hm)
            
            st.markdown("---")
            st.markdown("### 🔥 High-Conf. (>0.85)")
            
            # Helper for High Conf
            def normalize_score(x):
                try:
                    v = float(x)
                    return v * 1000 if v <= 1 else v
                except: return 0

            ppi_df['norm_score'] = ppi_df['Combined Score'].apply(normalize_score)
            high_conf_df = ppi_df[ppi_df['norm_score'] > 850]
            
            hc_prots = sorted(pd.concat([high_conf_df['Protein A'], high_conf_df['Protein B']]).unique())
            
            if hc_prots:
                sel_hc_prot = st.selectbox("Select Protein", hc_prots, key="hc_sel")
                
                subset_hc = high_conf_df[(high_conf_df['Protein A'] == sel_hc_prot) | (high_conf_df['Protein B'] == sel_hc_prot)]
                
                if not subset_hc.empty:
                    partners = []
                    scores = []
                    for _, r in subset_hc.iterrows():
                        p = r['Protein B'] if r['Protein A'] == sel_hc_prot else r['Protein A']
                        partners.append(p)
                        scores.append(r['norm_score'])
                    
                    hc_viz_df = pd.DataFrame({'Score': scores}, index=partners).sort_values('Score', ascending=False)
                    
                    fig_hc, ax_hc = plt.subplots(figsize=(5, len(partners)*0.4 + 1))
                    sns.heatmap(hc_viz_df, cmap="rocket_r", annot=True, cbar=False, ax=ax_hc)
                    st.pyplot(fig_hc)
                else:
                    st.write("No interactions > 0.85")
            else:
                st.write("No high confidence data available.")
                

        # --- RIGHT: Disease Network ---
        with top_right:
            # REMOVED white-card
            st.markdown("### 🕸 Disease Interaction Network")
            
            diseases = [
                "Alzheimer's Disease", "Parkinson's Disease", 
                "Amyotrophic Lateral Sclerosis (ALS)", "Multiple Sclerosis (MS)", 
                "Friedreich’s Ataxia (FA)"
            ]
            
            sel_disease = st.selectbox("Select Disease", diseases, key="net_dis_sel")
            
            if st.button("Simulate Network", key="sim_net"):
                dis_df = ppi_df[ppi_df['Disease Associated'].astype(str).str.contains(sel_disease, regex=False, na=False)]
                
                if not dis_df.empty:
                    # Filter top hubs to prevent overcrowding
                    G = nx.from_pandas_edgelist(dis_df, 'Protein A', 'Protein B', ['Combined Score'])
                    degrees = dict(G.degree)
                    top_nodes = sorted(degrees, key=degrees.get, reverse=True)[:30]
                    H = G.subgraph(top_nodes)
                    
                    net_viz = Network(height="600px", width="100%", bgcolor="#ffffff", font_color="black")
                    net_viz.from_nx(H)
                    
                    # Physics settings for force layout
                    net_viz.set_options("""
                    var options = {
                      "nodes": {
                        "color": {
                          "background": "#3b82f6",
                          "border": "white"
                        },
                        "font": {
                          "size": 16,
                          "face": "tahoma"
                        }
                      },
                      "physics": {
                        "forceAtlas2Based": {
                          "gravitationalConstant": -50,
                          "centralGravity": 0.01,
                          "springLength": 100,
                          "springConstant": 0.08
                        },
                        "minVelocity": 0.75,
                        "solver": "forceAtlas2Based"
                      }
                    }
                    """)
                    
                    with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_net:
                        net_viz.save_graph(tmp_net.name)
                        with open(tmp_net.name, 'r', encoding='utf-8') as f:
                            net_html = f.read()
                    
                    components.html(net_html, height=620)
                else:
                    st.warning("No interactions found for this disease.")
            else:
                st.info("Select a disease and click 'Simulate Network' to generate the graph.")
            

    # Bottom Charts
    with st.container():
        chart_c1, chart_c2 = st.columns(2)
        
        with chart_c1:
            # REMOVED white-card
            st.markdown("#### Distribution by Disease")
            if 'Disease Associated' in ppi_df.columns:
                dis_counts = ppi_df['Disease Associated'].value_counts().reset_index()
                dis_counts.columns = ['Disease', 'Count']
                fig_pie = px.pie(dis_counts, names='Disease', values='Count', color_discrete_sequence=px.colors.qualitative.Pastel)
                st.plotly_chart(fig_pie, use_container_width=True)
            
            # REMOVED white-card
            st.markdown("#### Top 15 Interacting Proteins")
            protein_counts = pd.concat([ppi_df['Protein A'], ppi_df['Protein B']]).value_counts().head(15).reset_index()
            protein_counts.columns = ['Protein', 'Count']
            fig_top = px.bar(protein_counts, x='Protein', y='Count', color='Count')
            st.plotly_chart(fig_top, use_container_width=True)
            
        with chart_c2:
            # REMOVED white-card
            st.markdown("#### Experimental Systems")
            if 'Experimental System' in ppi_df.columns:
                exp_counts = ppi_df['Experimental System'].value_counts().reset_index()
                exp_counts.columns = ['System', 'Count']
                fig_bar = px.bar(exp_counts, x='Count', y='System', orientation='h', color='Count', color_continuous_scale='Blues')
                st.plotly_chart(fig_bar, use_container_width=True)

            # REMOVED white-card
            st.markdown("#### Confidence Score Distribution")
            fig_hist = px.histogram(ppi_df, x='norm_score', nbins=20, title="Score Distribution", color_discrete_sequence=['#8884d8'])
            st.plotly_chart(fig_hist, use_container_width=True)


# ================= GITHUB EDIT TAB =================
with tabs[5]:
    # REMOVED white-card
    st.header("🛠 GitHub Repository")
    st.write("Access and edit the datasets directly on GitHub. Changes made to the CSV files in the repository will be automatically reflected here upon reloading the application.")
    
    links = {
        "📄 PPI Data (CSV)": "https://github.com/MeghanaVaddella/my-cv-dataset/blob/main/my-cv-data.csv",
        "🧬 3D Structure Part 1": "https://github.com/MeghanaVaddella/Neurodegenerative_Database/blob/main/3D%20Structure-1.csv",
        "🧬 3D Structure Part 2": "https://github.com/MeghanaVaddella/Neurodegenerative_Database/blob/main/3D%20Structure-2.csv",
        "🚫 No 3D Structure Data": "https://github.com/MeghanaVaddella/Neurodegenerative_Database/blob/main/No%203D%20Structure.csv"
    }
    
    for label, url in links.items():
        st.markdown(f"### [{label}]({url})")
