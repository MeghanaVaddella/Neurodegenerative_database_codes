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

# --- CUSTOM CSS ---
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

    /* Tab Styling */
    .stTabs [data-baseweb="tab-list"] {{
        gap: 30px; /* Increased space between tabs */
        background-color: {HEADER_BG};
        padding: 15px 20px;
        border-radius: 12px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }}
    .stTabs [data-baseweb="tab"] {{
        height: 80px;
        white-space: pre-wrap;
        background-color: transparent;
        border-radius: 8px;
        color: #cbd5e1; 
        font-weight: 700;
        font-size: 4.0rem;
        padding: 0 25px;
        margin-right: 15px;
    }}
    .stTabs [aria-selected="true"] {{
        background-color: {BODY_BG};
        color: {TEXT_COLOR};
        border-bottom: 4px solid {TEXT_COLOR};
        transform: translateY(-2px);
    }}

    /* Container/Card Styling - Kept for Home only if needed */
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
        padding: 4rem;
        border-radius: 24px;
        box-shadow: 0 20px 25px -5px rgba(0, 0, 0, 0.1);
        border: 1px solid #f1f5f9;
        margin-bottom: 3rem;
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
    </style>
""", unsafe_allow_html=True)

# --- HEADER ---
st.markdown('<div class="main-header">NEUROGEN PPI</div>', unsafe_allow_html=True)

# --- RAW TEXT DATA ---
HOME_TEXT_RAW = """
PROTEIN-PROTEIN INTERACTIONS (PPI)
What are Protein-Protein Interactions (PPIs) ?
Protein-Protein Interactions (PPIs) refer to the physical contacts between two or more protein molecules. These interactions can be temporary (transient) or permanent (stable). They are fundamental for almost every biological process like forming protein complexes, signaling inside cells, immune responses, metabolism, and more. Proteins interact using specific binding sites through hydrogen bonds, ionic bonds, hydrophobic interactions, and van der Waals forces. PPIs are essential because they govern many biological processes including DNA replication, gene expression, metabolic pathways, immune responses, and programmed cell death (apoptosis). Disruption in these interactions can lead to diseases such as cancer, autoimmune disorders, and neurodegenerative conditions. Thus, understanding PPIs not only helps in understanding how life functions at the molecular level but also in identifying the causes of various diseases.

Applications of the Protein-Protein Interactions
1.Drug Discovery:
Targeting PPIs enables the development of novel therapies, especially for diseases with previously undruggable targets such as cancer and viral infections. 
2.Diagnostics:
Abnormal PPI patterns serve as early biomarkers for detecting diseases before symptoms arise, improving diagnosis and treatment planning.
3.Synthetic Biology:
Engineered PPI networks are used to design programmable cells and synthetic pathways for use in medicine, agriculture, and bioengineering.
4.Functional Genomics:
PPI data helps predict the function of unknown proteins by analyzing their interaction partners within biological networks.
5.Structural Biology:
Studying PPIs reveals the architecture and dynamics of large protein complexes, offering insights into their mechanisms at the atomic level.
6.Systems Biology and Target Validation:
Mapping PPI networks allows researchers to model complex cellular systems and understand how proteins coordinate cellular processes and PPI analysis helps validate protein targets in drug development by confirming their role in disease-related pathways.
    
The database is about Alzheimer's Disease, Parkinson's Disease, Amyotrophic Lateral Sclerosis (ALS), Multiple Sclerosis (MS), Friedreich’s Ataxia (FA). It consists of about 3034 Protein-Protein Interactions involving the listed diseases. 
  
Here is some Information about the Diseases for which the PPIs are present in the database

Alzheimer's Disease  
Alzheimer's is a neurodegenerative disorder characterized by progressive memory loss, cognitive decline, and behavioral changes. It is caused by genetic mutations, aging, beta-amyloid plaques, and Tau protein tangles. Symptoms include confusion, difficulty in problem-solving, and mood swings. The disease primarily affects the hippocampus and spreads to other brain regions, involving neurotransmitters like acetylcholine and proteins such as beta-amyloid and Tau. Diagnosis involves cognitive assessments, MRI, PET scans, and cerebrospinal fluid tests. Treatment includes cholinesterase inhibitors, memantine, cognitive therapy, and lifestyle adjustments. Preventive strategies involve physical exercise, a brain-healthy diet, social engagement, and cognitive training, but the disease ultimately leads to progressive cognitive decline, loss of independence, and is eventually fatal.  

Parkinson's Disease  
Parkinson's primarily affects movement due to the degeneration of dopamine-producing neurons. It is associated with genetic predisposition, environmental toxins, and aging. Symptoms include tremors, rigidity, bradykinesia, and postural instability, as the disease targets the basal ganglia and substantia nigra. Dopamine and alpha-synuclein are key molecules involved in its pathology. Diagnosis is based on clinical symptoms and dopamine transporter scans. Treatment includes medications like levodopa, dopamine agonists, and deep brain stimulation, while lifestyle modifications such as regular exercise and a healthy diet help in management. The disease progressively impairs motor function, leading to difficulties in swallowing and, in advanced stages, dementia.  

Amyotrophic Lateral Sclerosis (ALS)  
ALS, also known as Lou Gehrig’s Disease, is a neurodegenerative condition affecting motor neurons, leading to muscle weakness and eventual paralysis. It is caused by genetic mutations, such as those in SOD1 and C9orf72, and environmental factors. Symptoms include muscle weakness, difficulty speaking and swallowing, and respiratory failure, as motor neurons in the brain and spinal cord deteriorate. Glutamate toxicity and SOD1 protein abnormalities play a role in the disease. Diagnosis is confirmed through electromyography (EMG) and genetic testing. Treatments like riluzole and edaravone, along with respiratory support, help manage symptoms, while physical therapy and assistive devices aid in mobility. ALS leads to progressive paralysis and significantly reduces life expectancy.  

Multiple Sclerosis (MS)  
MS is an autoimmune neurodegenerative disorder where the immune system attacks the myelin sheath of neurons, leading to symptoms such as fatigue, vision problems, muscle weakness, and coordination issues. It is associated with genetic factors, viral infections, and immune dysfunction, primarily affecting the brain and spinal cord. Key molecules involved include myelin proteins and inflammatory cytokines. Diagnosis involves MRI, lumbar puncture, and evoked potentials. Treatment consists of immunomodulatory drugs, corticosteroids, and physiotherapy. Lifestyle modifications, including a healthy diet and regular exercise, can help manage the condition. The disease progression varies, with periods of relapse and remission, and it can significantly impact daily life.  

Friedreich’s Ataxia (FA)  
Friedreich’s Ataxia is a rare inherited neurodegenerative disorder caused by mutations in the FXN gene, leading to reduced frataxin protein levels. It results in progressive damage to the nervous system, affecting movement and coordination. Symptoms include muscle weakness, ataxia, vision and hearing impairment, and scoliosis, as the disease primarily targets the spinal cord, peripheral nerves, and cerebellum. Diagnosis involves genetic testing, nerve conduction studies, and MRI. There is no cure, but treatments such as physical and speech therapy, along with supportive care, can help manage symptoms, while physical therapy and assistive devices aid in mobility. While there is no known prevention, maintaining mobility and heart health may slow progression. FA leads to progressive disability, an increased risk of heart disease, and a shortened lifespan.

What does the database Have?

1. The database presents the Protein-Protein Interactions Co-Expression, Experimentally Determined Interactions, Automated Textmining, Combined Score, Diseases Associated, BioGRID Interaction ID, Enterz Gene Interactor, BioGRID Interactor ID for Protein A and Protein B, Experimental System, Pubmed ID and the Author collected from STRING DATABASE, BIOGRID and IntACT.

2. The 3D-Visualization of the Protein-Protein Interactions is done using the Uniprot ID's, PDB ID's of both Protein A and Protein B and are viewed in MolStar Viewer. The structure of the Protein-Protein Interactions are visualized by using the Uniprot ID's in Pymol-3D Viewer.

3. The Protein Structure can be predicted using the AlphaFold-Multimer by generating the FASTA Sequences which help in generating the Protein Foldings of the Interaction using Google Colab having the AlphaFold2 in which the templates are generated using MMseq2. The Structures of the Protein Interactions can be be viewed in Chimera by downloading the PDB file from the 3D Visualizer: AlphaFold-based 3D Viewer (py3Dmol).

Additionally, the data present can be downloaded and new Data can be added using the GitHub links present in the Github Edit Tab.
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
# Tabs with emojis and names beside them
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
    lines = HOME_TEXT_RAW.strip().split('\n')
    
    headers = [
        "PROTEIN-PROTEIN INTERACTIONS (PPI)",
        "Applications of the Protein-Protein Interactions",
        "What does the database Have?",
        "Here is some Information about the Diseases"
    ]
    
    sub_headers = [
        "What are Protein-Protein Interactions (PPIs) ?",
        "Alzheimer's Disease", "Parkinson's Disease", 
        "Amyotrophic Lateral Sclerosis (ALS)", "Multiple Sclerosis (MS)", 
        "Friedreich’s Ataxia (FA)"
    ]
    
    def render_bubble(content_lines, is_last_card=False):
        if not content_lines: return
        html_content = ""
        
        for line in content_lines:
            line = line.strip()
            if not line: continue
            
            # Check for Main Headers
            is_header = any(line.startswith(h) for h in headers)
            if is_header:
                html_content += f"<h2 style='color:#881337; font-size:46px; border-bottom:3px solid #88133720; padding-bottom:15px; margin-bottom:25px; font-family:Georgia, serif;'>{line}</h2>"
                continue
                
            # Check for Sub Headers (Diseases, PPI Question)
            is_sub = any(line.startswith(s) for s in sub_headers)
            if is_sub:
                html_content += f"<div style='color:#991b1b; font-size:36px; font-weight:bold; margin-top:30px; margin-bottom:15px; font-family:Georgia, serif;'>{line}</div>"
                continue
            
            # Check for Lists
            if line[0].isdigit() and ("." in line[:3]):
                 html_content += f"<div style='color:#92400e; font-size:30px; font-weight:bold; margin-top:20px; margin-bottom:10px; font-family:Georgia, serif;'>{line}</div>"
                 continue
                 
            # Regular Text
            p_style = "color:#000000; font-size:26px; line-height:1.7; margin-bottom:15px;"
            if is_last_card and "Additionally" in line:
                # UPDATED: Black, non-bold
                p_style = "color:#000000; font-size:28px; line-height:1.7; margin-bottom:15px; font-weight:normal;"
            
            html_content += f"<p style='{p_style}'>{line}</p>"

        st.markdown(f"<div class='bubble-card'>{html_content}</div>", unsafe_allow_html=True)

    block_buffer = []
    for line in lines:
        if any(line.strip().startswith(h) for h in headers):
            if block_buffer:
                is_last = "Additionally" in "".join(block_buffer)
                render_bubble(block_buffer, is_last)
                block_buffer = []
        block_buffer.append(line)
    
    if block_buffer:
        is_last = "Additionally" in "".join(block_buffer)
        render_bubble(block_buffer, is_last)


# ================= DATA TAB =================
with tabs[1]:
    with st.container():
        # REMOVED white-card
        st.header("Protein-Protein Interaction Data")
        
        search_term = st.text_input("Search Data", placeholder="Type to search...", key="ppi_search")
        
        display_df = ppi_df
        if search_term:
            display_df = ppi_df[ppi_df.astype(str).apply(lambda x: x.str.contains(search_term, case=False)).any(axis=1)]
            
        st.dataframe(display_df, use_container_width=True, hide_index=True)
        st.download_button("Download PPI CSV", ppi_df.to_csv(index=False), "PPI_data.csv", "text/csv")
        
        st.markdown("---")
        st.subheader("Visualize Interactions")
        
        # Interaction Graph Visualizer (Star Network)
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
                        
                components.html(html_bytes, height=620)
                st.download_button("Download Network HTML", html_bytes, f"ppi_network_{selected_protein}.html", "text/html")
            else:
                st.warning("No interactions found for this protein.")


# ================= 3D STRUCTURE DATA TAB =================
with tabs[2]:
    with st.container():
        # REMOVED white-card
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
    # REMOVED white-card
    st.header("3D Protein Structure Visualizer")
    
    # --- Layout: Controls (Left) | Info (Right) ---
    top_col1, top_col2 = st.columns([1, 1], gap="large")
    
    with top_col1:
        st.markdown("#### ⚙ Controls")
        
        unique_prot_a = sorted(df_3d['Protein A'].dropna().unique())
        sel_prot_a = st.selectbox("Select Protein A", [""] + unique_prot_a, key="viz_a")
        
        # Linked Search Logic: Filter B based on A
        avail_prot_b = []
        if sel_prot_a:
            avail_prot_b = sorted(df_3d[df_3d['Protein A'] == sel_prot_a]['Protein B'].unique())
        
        sel_prot_b = st.selectbox(
            "Select Protein B", 
            [""] + avail_prot_b, 
            key="viz_b", 
            disabled=not sel_prot_a
        )
    
    with top_col2:
        st.markdown("#### ℹ Interaction Info")
        if sel_prot_a:
            row_a = df_3d[df_3d['Protein A'] == sel_prot_a].iloc[0]
            st.info(f"*Protein A: {row_a['Protein A']} | **UniProt: {row_a['UniProtID A']} | **PDBs*: {row_a['PDB ID A']}")
            
        if sel_prot_b:
            row_b = df_3d[(df_3d['Protein A'] == sel_prot_a) & (df_3d['Protein B'] == sel_prot_b)].iloc[0]
            st.success(f"*Protein B: {row_b['Protein B']} | **UniProt: {row_b['UniProtID B']} | **PDBs*: {row_b['PDB ID B']}")

    st.markdown("---")
    
    # --- Mol* Viewer (Iframe) ---
    st.markdown("### 🧬 Mol* (MolStar) Viewer")
    pdb_ids = []
    
    if sel_prot_a:
        row_a = df_3d[df_3d['Protein A'] == sel_prot_a].iloc[0]
        ids_a = str(row_a['PDB ID A']).split(',')
        pdb_ids.extend([x.strip() for x in ids_a if x.strip() != 'NA'])
        
    if sel_prot_b:
        row_b = df_3d[(df_3d['Protein A'] == sel_prot_a) & (df_3d['Protein B'] == sel_prot_b)].iloc[0]
        ids_b = str(row_b['PDB ID B']).split(',')
        pdb_ids.extend([x.strip() for x in ids_b if x.strip() != 'NA'])
        
    pdb_ids = list(set(pdb_ids)) # Unique IDs

    if pdb_ids:
        # Construct dynamic Mol* URL
        molstar_url = "https://molstar.org/viewer/?url=" + ",".join([f"https://files.rcsb.org/download/{pdb}.pdb" for pdb in pdb_ids])
        st.components.v1.iframe(molstar_url, width=1000, height=600)
    else:
        st.warning("No valid PDB IDs found for visualization or no proteins selected.")

    st.markdown("---")
    
    # --- AlphaFold Download ---
    st.subheader("💾 Download AlphaFold Predicted Structures")
    
    af_col1, af_col2 = st.columns(2)
    with af_col1:
        uniprot_opts_a = sorted(df_3d['UniProtID A'].dropna().unique())
        af_uni_a = st.selectbox("Select UniProt ID A (AF)", [""] + uniprot_opts_a, key="af_a")
    
    with af_col2:
        af_uni_b_opts = []
        if af_uni_a:
            af_uni_b_opts = sorted(df_3d[df_3d['UniProtID A'] == af_uni_a]['UniProtID B'].unique())
        af_uni_b = st.selectbox("Select UniProt ID B (AF)", [""] + af_uni_b_opts, key="af_b", disabled=not af_uni_a)
        
    if st.button("Fetch & Download AlphaFold PDBs"):
        if af_uni_a and af_uni_b:
            def get_af_pdb(uid):
                # Try v4, v3, v2
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
                
                # UPDATED Logic: Download whatever is found, suppress "Could not find" error
                if pdb_a_txt and pdb_b_txt:
                    combined = f"REMARK Protein A: {af_uni_a}\n{pdb_a_txt}\nTER\nREMARK Protein B: {af_uni_b}\n{pdb_b_txt}\nEND"
                    st.download_button(
                        label="Download Combined PDB",
                        data=combined,
                        file_name=f"AF_{af_uni_a}_{af_uni_b}.pdb",
                        mime="chemical/x-pdb"
                    )
                    st.success("Ready for download!")
                elif pdb_a_txt:
                    st.download_button(label=f"Download {af_uni_a}", data=pdb_a_txt, file_name=f"AF_{af_uni_a}.pdb", mime="chemical/x-pdb")
                elif pdb_b_txt:
                    st.download_button(label=f"Download {af_uni_b}", data=pdb_b_txt, file_name=f"AF_{af_uni_b}.pdb", mime="chemical/x-pdb")
                else:
                    # Do nothing / suppress error as requested
                    pass
        else:
            st.warning("Select both UniProt IDs first.")
            
    st.markdown("---")

    # --- FASTA Generator ---
    st.subheader("🧬 Predict Interactions using AlphaFold-Multimer")
    
    fa_col1, fa_col2 = st.columns(2)
    with fa_col1:
        fa_uni_a = st.selectbox("UniProt ID A (FASTA)", [""] + uniprot_opts_a, key="fa_a")
    with fa_col2:
        fa_uni_b_opts = []
        if fa_uni_a:
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
                
                # Styled Link Button
                st.markdown("""
                <a href="https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb" target="_blank" class="colab-btn">
                    Open ColabFold 🔗
                </a>
                """, unsafe_allow_html=True)
            else:
                st.error("Failed to fetch sequences.")

    st.markdown("---")

    # --- PDB Upload (Local Viewer using py3Dmol) ---
    st.subheader("📦 Upload Predicted PDB File")
    uploaded_file = st.file_uploader("Upload PDB", type=['pdb'])
    if uploaded_file:
        pdb_content = uploaded_file.read().decode("utf-8")
        st.text_area("File Content Preview", pdb_content, height=150)
        
        # Add Download Option for the uploaded file (or visualize content)
        st.download_button("Download Uploaded Structure", pdb_content, "uploaded_structure.pdb", "chemical/x-pdb")
        
        view_up = py3Dmol.view(width=800, height=500)
        view_up.addModel(pdb_content, "pdb")
        # Ensure rainbow (spectrum)
        view_up.setStyle({'cartoon': {'color': 'spectrum'}})
        view_up.zoomTo()
        view_up.setBackgroundColor('white')
        components.html(view_up._make_html(), height=500)
        


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

