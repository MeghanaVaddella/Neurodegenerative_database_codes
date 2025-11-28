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
import plotly.graph_objects as go
import tempfile
import time

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
    /* Global Styles */
    .stApp {{
        background-color: {BODY_BG};
        color: {TEXT_COLOR};
        font-family: 'Inter', sans-serif;
    }}
    
    /* Header Styling */
    .main-header {{
        background-color: {HEADER_BG};
        padding: 2rem;
        text-align: center;
        color: #e2e8f0;
        font-family: 'Georgia', serif;
        font-size: 3rem;
        font-weight: bold;
        letter-spacing: 0.1em;
        text-transform: uppercase;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
        margin-bottom: 2rem;
        border-radius: 8px;
    }}

    /* Tab Styling */
    .stTabs [data-baseweb="tab-list"] {{
        gap: 8px;
        background-color: {HEADER_BG};
        padding: 10px;
        border-radius: 8px;
    }}
    .stTabs [data-baseweb="tab"] {{
        height: 50px;
        white-space: pre-wrap;
        background-color: transparent;
        border-radius: 4px;
        color: #e2e8f0;
        font-weight: 600;
    }}
    .stTabs [aria-selected="true"] {{
        background-color: {BODY_BG};
        color: {TEXT_COLOR};
    }}

    /* Container/Card Styling */
    .css-1r6slb0, .css-12oz5g7 {{
        background-color: white;
        padding: 2rem;
        border-radius: 12px;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
        border: 1px solid #e2e8f0;
    }}
    
    /* Home Bubble Styling */
    .bubble-card {{
        background-color: white;
        padding: 2.5rem;
        border-radius: 12px;
        box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1);
        border: 1px solid #f1f5f9;
        margin-bottom: 2rem;
        transition: transform 0.2s;
    }}
    .bubble-card:hover {{
        transform: translateY(-2px);
    }}
    
    /* Headings */
    h1, h2, h3 {{
        color: #1e3a8a !important;
        font-family: 'Georgia', serif;
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
Alzheimer's Disease is a neurodegenerative disorder characterized by progressive memory loss, cognitive decline, and behavioral changes. It is caused by genetic mutations, aging, beta-amyloid plaques, and Tau protein tangles. Symptoms include confusion, difficulty in problem-solving, and mood swings. The disease primarily affects the hippocampus and spreads to other brain regions, involving neurotransmitters like acetylcholine and proteins such as beta-amyloid and Tau. Diagnosis involves cognitive assessments, MRI, PET scans, and cerebrospinal fluid tests. Treatment includes cholinesterase inhibitors, memantine, cognitive therapy, and lifestyle adjustments. Preventive strategies involve physical exercise, a brain-healthy diet, social engagement, and cognitive training, but the disease ultimately leads to progressive cognitive decline, loss of independence, and is eventually fatal.  

Parkinson's Disease  
Parkinson's Disease primarily affects movement due to the degeneration of dopamine-producing neurons. It is associated with genetic predisposition, environmental toxins, and aging. Symptoms include tremors, rigidity, bradykinesia, and postural instability, as the disease targets the basal ganglia and substantia nigra. Dopamine and alpha-synuclein are key molecules involved in its pathology. Diagnosis is based on clinical symptoms and dopamine transporter scans. Treatment includes medications like levodopa, dopamine agonists, and deep brain stimulation, while lifestyle modifications such as regular exercise and a healthy diet help in management. The disease progressively impairs motor function, leading to difficulties in swallowing and, in advanced stages, dementia.  

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
tabs = st.tabs([
    "Home", 
    "Data", 
    "3D Structure Data", 
    "3D Visualizer", 
    "Data Visualizer",
    "GitHub Edit"
])

# ================= HOME TAB =================
with tabs[0]:
    # Custom text parser to create bubbles
    lines = HOME_TEXT_RAW.strip().split('\n')
    
    current_bubble_content = []
    
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
    
    def render_bubble(content_lines):
        if not content_lines: return
        html_content = ""
        for line in content_lines:
            line = line.strip()
            if not line: continue
            
            # Check for Main Headers
            is_header = any(line.startswith(h) for h in headers)
            if is_header:
                html_content += f"<h2 style='color:#881337; font-size:32px; border-bottom:2px solid #88133720; padding-bottom:10px; margin-bottom:20px;'>{line}</h2>"
                continue
                
            # Check for Sub Headers
            is_sub = any(line.startswith(s) for s in sub_headers)
            if is_sub:
                html_content += f"<div style='color:#991b1b; font-size:24px; font-weight:bold; margin-top:20px; margin-bottom:10px; font-family:Georgia, serif;'>{line}</div>"
                continue
            
            # Check for Lists
            if line[0].isdigit() and ("." in line[:3]):
                 html_content += f"<div style='color:#92400e; font-size:20px; font-weight:bold; margin-top:15px; margin-bottom:5px; font-family:serif;'>{line}</div>"
                 continue
                 
            # Regular Text
            html_content += f"<p style='color:#1f2937; font-size:18px; line-height:1.6; margin-bottom:10px;'>{line}</p>"

        st.markdown(f"<div class='bubble-card'>{html_content}</div>", unsafe_allow_html=True)

    # Logic to split text into blocks based on Main Headers
    block_buffer = []
    for line in lines:
        if any(line.strip().startswith(h) for h in headers):
            if block_buffer:
                render_bubble(block_buffer)
                block_buffer = []
        block_buffer.append(line)
    
    if block_buffer:
        render_bubble(block_buffer)


# ================= DATA TAB =================
with tabs[1]:
    with st.container():
        st.markdown("<div style='background-color:white; padding:20px; border-radius:10px;'>", unsafe_allow_html=True)
        st.header("Protein-Protein Interaction Data")
        
        # Search functionality
        search_term = st.text_input("Search Data", placeholder="Type to search...", key="ppi_search")
        
        display_df = ppi_df
        if search_term:
            display_df = ppi_df[ppi_df.astype(str).apply(lambda x: x.str.contains(search_term, case=False)).any(axis=1)]
            
        st.dataframe(display_df, use_container_width=True, hide_index=True)
        st.download_button("Download PPI CSV", ppi_df.to_csv(index=False), "PPI_data.csv", "text/csv")
        
        st.markdown("---")
        st.subheader("Visualize Interactions")
        
        # Interaction Graph Visualizer
        all_proteins = pd.concat([ppi_df['Protein A'], ppi_df['Protein B']]).unique()
        selected_protein = st.selectbox("Choose Protein", sorted(all_proteins), key="ppi_viz_select")
        
        if st.button("Generate Network Graph", key="gen_graph_btn"):
            # Filter data for selected protein (Star Network)
            subset = ppi_df[(ppi_df['Protein A'] == selected_protein) | (ppi_df['Protein B'] == selected_protein)]
            
            if not subset.empty:
                net = Network(height="600px", width="100%", bgcolor="#ffffff", font_color="black")
                
                # Add central node
                net.add_node(selected_protein, color="#dc2626", size=25, title=selected_protein)
                
                # Add neighbors
                for _, row in subset.iterrows():
                    partner = row['Protein B'] if row['Protein A'] == selected_protein else row['Protein A']
                    score = row['Combined Score']
                    net.add_node(partner, color="#3b82f6", size=15, title=f"{partner}\nScore: {score}")
                    net.add_edge(selected_protein, partner, value=float(score) if isinstance(score, (int, float)) else 1)
                
                # Save and read
                with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp:
                    net.save_graph(tmp.name)
                    with open(tmp.name, 'r', encoding='utf-8') as f:
                        html_bytes = f.read()
                        
                components.html(html_bytes, height=620)
                st.download_button("Download Network HTML", html_bytes, f"ppi_network_{selected_protein}.html", "text/html")
            else:
                st.warning("No interactions found for this protein.")
        
        st.markdown("</div>", unsafe_allow_html=True)


# ================= 3D STRUCTURE DATA TAB =================
with tabs[2]:
    with st.container():
        st.markdown("<div style='background-color:white; padding:20px; border-radius:10px;'>", unsafe_allow_html=True)
        
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
        
        st.markdown("</div>", unsafe_allow_html=True)


# ================= 3D VISUALIZER TAB =================
with tabs[3]:
    st.markdown("<div style='background-color:white; padding:20px; border-radius:10px;'>", unsafe_allow_html=True)
    st.header("3D Protein Structure Visualizer")
    
    # --- Layout: Controls Left, Viewer Right ---
    viz_col1, viz_col2 = st.columns([1, 2], gap="large")
    
    with viz_col1:
        st.markdown("#### ⚙️ Controls")
        
        # Protein A
        unique_prot_a = sorted(df_3d['Protein A'].dropna().unique())
        sel_prot_a = st.selectbox("Select Protein A", [""] + unique_prot_a, key="viz_a")
        
        # Filter Protein B based on A
        avail_prot_b = []
        if sel_prot_a:
            avail_prot_b = sorted(df_3d[df_3d['Protein A'] == sel_prot_a]['Protein B'].unique())
        
        sel_prot_b = st.selectbox(
            "Select Protein B", 
            [""] + avail_prot_b, 
            key="viz_b", 
            disabled=not sel_prot_a
        )
        
        # Metadata Display
        if sel_prot_a:
            row_a = df_3d[df_3d['Protein A'] == sel_prot_a].iloc[0]
            st.info(f"**Protein A**: {row_a['Protein A']}\n\n**UniProt**: {row_a['UniProtID A']}\n\n**PDBs**: {row_a['PDB ID A']}")
            
        if sel_prot_b:
            row_b = df_3d[(df_3d['Protein A'] == sel_prot_a) & (df_3d['Protein B'] == sel_prot_b)].iloc[0]
            st.warning(f"**Protein B**: {row_b['Protein B']}\n\n**UniProt**: {row_b['UniProtID B']}\n\n**PDBs**: {row_b['PDB ID B']}")

    with viz_col2:
        # 3D Viewer Logic
        viewer_height = 500
        view = py3Dmol.view(height=viewer_height, width="100%")
        view.setBackgroundColor('white')
        
        has_model = False
        
        if sel_prot_a:
            row_a = df_3d[df_3d['Protein A'] == sel_prot_a].iloc[0]
            pdbs_a = [p.strip() for p in str(row_a['PDB ID A']).split(',') if p.strip() != 'NA']
            if pdbs_a:
                view.addModel(f"pdb:{pdbs_a[0]}", 'pdb') # Load from RCSB
                view.setStyle({'model': -1}, {'cartoon': {'color': 'blue'}})
                has_model = True
        
        if sel_prot_b:
            row_b = df_3d[(df_3d['Protein A'] == sel_prot_a) & (df_3d['Protein B'] == sel_prot_b)].iloc[0]
            pdbs_b = [p.strip() for p in str(row_b['PDB ID B']).split(',') if p.strip() != 'NA']
            if pdbs_b:
                view.addModel(f"pdb:{pdbs_b[0]}", 'pdb') # Load from RCSB
                view.setStyle({'model': -1}, {'cartoon': {'color': 'orange'}})
                has_model = True

        if has_model:
            view.zoomTo()
        else:
            # Placeholder text via style since py3Dmol is canvas
            pass 
            
        html_view = view._make_html()
        components.html(html_view, height=viewer_height + 20)

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
                # Try v4, then v3
                for v in [4, 3, 2, 1]:
                    url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v{v}.pdb"
                    r = requests.get(url)
                    if r.status_code == 200: return r.text
                return None
            
            with st.spinner("Fetching AlphaFold structures..."):
                pdb_a_txt = get_af_pdb(af_uni_a)
                pdb_b_txt = get_af_pdb(af_uni_b)
                
                if pdb_a_txt and pdb_b_txt:
                    combined = f"REMARK Protein A: {af_uni_a}\n{pdb_a_txt}\nTER\nREMARK Protein B: {af_uni_b}\n{pdb_b_txt}\nEND"
                    st.download_button(
                        label="Download Combined PDB",
                        data=combined,
                        file_name=f"AF_{af_uni_a}_{af_uni_b}.pdb",
                        mime="chemical/x-pdb"
                    )
                    st.success("Ready for download!")
                else:
                    st.error("Could not find AlphaFold entries for one or both proteins.")
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
                r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}.fasta")
                return r.text if r.ok else ""
            
            f1 = get_fasta(fa_uni_a)
            f2 = get_fasta(fa_uni_b)
            
            if f1 and f2:
                combined_fasta = f"{f1}\n{f2}"
                st.download_button("Download FASTA", combined_fasta, "multimer.fasta", "text/plain")
                st.markdown("[Open ColabFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)")
            else:
                st.error("Failed to fetch sequences.")

    st.markdown("---")

    # --- PDB Upload ---
    st.subheader("📦 Upload Predicted PDB File")
    uploaded_file = st.file_uploader("Upload PDB", type=['pdb'])
    if uploaded_file:
        pdb_content = uploaded_file.read().decode("utf-8")
        view_up = py3Dmol.view(width=800, height=500)
        view_up.addModel(pdb_content, "pdb")
        view_up.setStyle({'cartoon': {'color': 'spectrum'}})
        view_up.zoomTo()
        view_up.setBackgroundColor('white')
        components.html(view_up._make_html(), height=500)
        
    st.markdown("</div>", unsafe_allow_html=True)


# ================= DATA VISUALIZER TAB =================
with tabs[4]:
    # Layout: Top Row (Heatmaps Left, Network Right)
    
    # 1. Top Section Container
    with st.container():
        # Using columns to create the layout: 1/3 width for Heatmaps, 2/3 width for Network
        top_left, top_right = st.columns([1, 2], gap="medium")
        
        # --- LEFT COLUMN: HEATMAPS ---
        with top_left:
            st.markdown("<div style='background-color:white; padding:15px; border-radius:10px; margin-bottom:20px; border:1px solid #e5e7eb;'>", unsafe_allow_html=True)
            st.markdown("### 🧬 Interactive Matrix")
            
            # Interactive Matrix Logic (Simulated with multiselect)
            all_prots = sorted(pd.concat([ppi_df['Protein A'], ppi_df['Protein B']]).unique())
            default_prots = all_prots[:5]
            selected_heatmap_prots = st.multiselect("Add Proteins", all_prots, default=default_prots, key="hm_multi")
            
            if selected_heatmap_prots:
                # Build adjacency matrix
                hm_data = pd.DataFrame(index=selected_heatmap_prots, columns=selected_heatmap_prots).fillna(0)
                for p1 in selected_heatmap_prots:
                    for p2 in selected_heatmap_prots:
                        if p1 == p2:
                            hm_data.loc[p1, p2] = 1000 # Max score for self
                        else:
                            match = ppi_df[((ppi_df['Protein A'] == p1) & (ppi_df['Protein B'] == p2)) | 
                                           ((ppi_df['Protein B'] == p1) & (ppi_df['Protein A'] == p2))]
                            if not match.empty:
                                val = match.iloc[0]['Combined Score']
                                hm_data.loc[p1, p2] = float(val) * 1000 if float(val) <= 1 else float(val)

                fig_hm, ax_hm = plt.subplots(figsize=(5, 4))
                sns.heatmap(hm_data.astype(float), cmap="rocket_r", ax=ax_hm, cbar=False, annot=False)
                st.pyplot(fig_hm)
            st.markdown("</div>", unsafe_allow_html=True)

            st.markdown("<div style='background-color:white; padding:15px; border-radius:10px; border:1px solid #e5e7eb;'>", unsafe_allow_html=True)
            st.markdown("### 🔥 High-Conf. (>0.85)")
            
            # Filter High Confidence > 0.85 (assuming 0-1 scale or 0-1000 scale)
            # Normalize score column first
            temp_df = ppi_df.copy()
            # Handle mixed types or scaling
            def clean_score(x):
                try:
                    v = float(x)
                    return v if v > 1 else v * 1000
                except:
                    return 0
            
            temp_df['score_norm'] = temp_df['Combined Score'].apply(clean_score)
            high_conf_df = temp_df[temp_df['score_norm'] > 850]
            
            hc_prots = sorted(pd.concat([high_conf_df['Protein A'], high_conf_df['Protein B']]).unique())
            
            if hc_prots:
                sel_hc_prot = st.selectbox("Select Protein", hc_prots, key="hc_sel")
                
                # Show neighbors heatmap
                subset_hc = high_conf_df[(high_conf_df['Protein A'] == sel_hc_prot) | (high_conf_df['Protein B'] == sel_hc_prot)]
                
                if not subset_hc.empty:
                    partners = []
                    scores = []
                    for _, r in subset_hc.iterrows():
                        p = r['Protein B'] if r['Protein A'] == sel_hc_prot else r['Protein A']
                        partners.append(p)
                        scores.append(r['score_norm'])
                    
                    hc_viz_df = pd.DataFrame({'Partner': partners, 'Score': scores}).set_index('Partner').sort_values('Score', ascending=False)
                    
                    fig_hc, ax_hc = plt.subplots(figsize=(5, len(partners)*0.3 + 1))
                    sns.heatmap(hc_viz_df, cmap="rocket_r", annot=True, cbar=False, ax=ax_hc)
                    st.pyplot(fig_hc)
                else:
                    st.write("No interactions > 0.85")
            else:
                st.write("No high confidence data available.")
            st.markdown("</div>", unsafe_allow_html=True)

        # --- RIGHT COLUMN: DISEASE NETWORK ---
        with top_right:
            st.markdown("<div style='background-color:white; padding:20px; border-radius:10px; height:100%; border:1px solid #e5e7eb;'>", unsafe_allow_html=True)
            st.markdown("### 🕸️ Disease Interaction Network")
            
            diseases = [
                "Alzheimer's Disease", "Parkinson's Disease", 
                "Amyotrophic Lateral Sclerosis (ALS)", "Multiple Sclerosis (MS)", 
                "Friedreich’s Ataxia (FA)"
            ]
            
            sel_disease = st.selectbox("Select Disease", diseases, key="net_dis_sel")
            
            if st.button("Simulate Network", key="sim_net"):
                # Filter by disease
                dis_df = ppi_df[ppi_df['Disease Associated'].astype(str).str.contains(sel_disease, regex=False, na=False)]
                
                if not dis_df.empty:
                    # Filter for hubs (Top 30 by degree)
                    G = nx.from_pandas_edgelist(dis_df, 'Protein A', 'Protein B', ['Combined Score'])
                    degrees = dict(G.degree)
                    top_nodes = sorted(degrees, key=degrees.get, reverse=True)[:30]
                    H = G.subgraph(top_nodes)
                    
                    net_viz = Network(height="500px", width="100%", bgcolor="#ffffff", font_color="black")
                    net_viz.from_nx(H)
                    
                    # Custom physics/options
                    net_viz.set_options("""
                    var options = {
                      "nodes": {
                        "color": {
                          "background": "#3b82f6",
                          "border": "white"
                        },
                        "font": {
                          "size": 14,
                          "face": "tahoma"
                        }
                      },
                      "edges": {
                        "color": {
                          "inherit": true
                        },
                        "smooth": false
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
                    
                    components.html(net_html, height=520)
                else:
                    st.warning("No interactions found for this disease.")
            else:
                st.info("Click 'Simulate Network' to generate the graph.")
            
            st.markdown("</div>", unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)

    # 2. Bottom Section: Charts
    with st.container():
        chart_c1, chart_c2 = st.columns(2)
        
        with chart_c1:
            st.markdown("<div style='background-color:white; padding:15px; border-radius:10px;'>", unsafe_allow_html=True)
            st.markdown("#### Distribution by Disease")
            if 'Disease Associated' in ppi_df.columns:
                dis_counts = ppi_df['Disease Associated'].value_counts().reset_index()
                dis_counts.columns = ['Disease', 'Count']
                fig_pie = px.pie(dis_counts, names='Disease', values='Count', color_discrete_sequence=px.colors.qualitative.Pastel)
                st.plotly_chart(fig_pie, use_container_width=True)
            st.markdown("</div>", unsafe_allow_html=True)
            
        with chart_c2:
            st.markdown("<div style='background-color:white; padding:15px; border-radius:10px;'>", unsafe_allow_html=True)
            st.markdown("#### Experimental Systems")
            if 'Experimental System' in ppi_df.columns:
                exp_counts = ppi_df['Experimental System'].value_counts().reset_index()
                exp_counts.columns = ['System', 'Count']
                fig_bar = px.bar(exp_counts, x='Count', y='System', orientation='h', color='Count', color_continuous_scale='Blues')
                st.plotly_chart(fig_bar, use_container_width=True)
            st.markdown("</div>", unsafe_allow_html=True)


# ================= GITHUB EDIT TAB =================
with tabs[5]:
    st.markdown("<div style='background-color:white; padding:40px; border-radius:10px; text-align:center;'>", unsafe_allow_html=True)
    st.header("GitHub Repository")
    st.write("Access and edit the datasets directly on GitHub. Changes made to the CSV files in the repository will be automatically reflected here upon reloading the application.")
    
    links = {
        "PPI Data (CSV)": "https://github.com/MeghanaVaddella/my-cv-dataset/blob/main/my-cv-data.csv",
        "3D Structure Part 1": "https://github.com/MeghanaVaddella/Neurodegenerative_Database/blob/main/3D%20Structure-1.csv",
        "3D Structure Part 2": "https://github.com/MeghanaVaddella/Neurodegenerative_Database/blob/main/3D%20Structure-2.csv",
        "No 3D Structure Data": "https://github.com/MeghanaVaddella/Neurodegenerative_Database/blob/main/No%203D%20Structure.csv"
    }
    
    for label, url in links.items():
        st.markdown(f"### [{label}]({url})")
        
    st.markdown("</div>", unsafe_allow_html=True)
