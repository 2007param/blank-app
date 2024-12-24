import streamlit as st
import pandas as pd
from collections import Counter
from PIL import Image

def footer():
st.markdown(
    """
    <footer style='text-align: center; padding: 10px; font-size: 14px; color: #777;'>
        <p>© 2024 Bioinformatics Tool. All rights reserved.</p>
    </footer>
    """, unsafe_allow_html=True
)

def is_valid_sequence(sequence):
return all(base in 'ACTG' for base in sequence)

def apply_background_color(color):
st.markdown(
    f"""
    <style>
    .stApp {{
        background-color: {color};
    }}
    </style>
    """,
    unsafe_allow_html=True
)

def get_nucleotide_count(sequence):
if is_valid_sequence(sequence):
    apply_background_color("#FFC0CB")  
    return {'A': sequence.count('A'), 'T': sequence.count('T'),
            'C': sequence.count('C'), 'G': sequence.count('G')}
else:
    return "Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G."

# Function for k-mer analysis
def kmer_analysis(sequence, k):
if is_valid_sequence(sequence):
    kmers = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        kmers[kmer] = kmers.get(kmer, 0) + 1
    apply_background_color("#008000")  
    return kmers
else:
    return "Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G."

# Function for Hamming distance
def hamming_distance(seq1, seq2):
if len(seq1) != len(seq2):
    return "Error: Sequences must have the same length"
if is_valid_sequence(seq1) and is_valid_sequence(seq2):
    apply_background_color("#D3D3D3")  
    return sum(el1 != el2 for el1, el2 in zip(seq1, seq2))
else:
    return "Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G."

# Function for Gene finding
def find_genes(sequence):
if is_valid_sequence(sequence):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    genes = []
    for i in range(len(sequence)):
        if sequence[i:i + 3] == start_codon:
            for j in range(i + 3, len(sequence), 3):
                if sequence[j:j + 3] in stop_codons:
                    genes.append(sequence[i:j + 3])
                    break
    apply_background_color("#00FF00")  
    return genes
else:
    return "Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G."

# Function for Reverse complement
def reverse_complement(sequence):
if is_valid_sequence(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    apply_background_color("#FFD700")  
    return ''.join(complement[base] for base in reversed(sequence))
else:
    return "Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G."

# Function for GC Content
def gc_content(sequence):
if is_valid_sequence(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    apply_background_color("#A52A2A")  
    return (gc_count / len(sequence)) * 100
else:
    return "Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G."

# Function for Transcription
def transcription(sequence):
if is_valid_sequence(sequence):
    apply_background_color("#FFA500")  
    return sequence.replace('T', 'U')
else:
    return "Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G."

# Function for Translation (simple case, codon-to-amino acid)
def translation(sequence):
if is_valid_sequence(sequence):
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }
    protein = ""
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        protein += codon_table.get(codon, '-')
    apply_background_color("#FAD7A0") 
    return protein
else:
    return "Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G."

# Function for Sequence alignment (simple version)
def sequence_alignment(seq1, seq2):
if is_valid_sequence(seq1) and is_valid_sequence(seq2):
    from difflib import SequenceMatcher
    matcher = SequenceMatcher(None, seq1, seq2)
    match_ratio = matcher.ratio()
    apply_background_color("#808000")  
    return f"Sequence similarity ratio: {match_ratio*100:.2f}%"
else:
    return "Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G."
# Function to display the title page
def title_page():    
#display_image()
apply_background_color("#D6EAF8")  
st.markdown("<h1 style='text-align: center;'>🧬 Bioinformatics Tool</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center;'>Welcome to the Bioinformatics Tool. This app provides several bioinformatics analysis tools such as Nucleotide count, K-mer analysis, Gene Finding, and more. Select an option from the menu to get started!</p>", unsafe_allow_html=True)
footer()
def about_us_page():
#display_image1()
apply_background_color("#F4ECF7")  
st.markdown("<h1 style='text-align: center;'>About Us</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center;'>We are a team of passionate bioinformaticians, biologists, and developers working to provide tools for DNA sequence analysis. Our goal is to make bioinformatics more accessible and provide a range of useful functionalities for researchers and students alike.</p>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center;'>If you have any questions or feedback, please don't hesitate to contact us at:</p>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center;'>Email: ai-tutorial@example.com</p>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center;'>Phone: (123)-456-7890</p>", unsafe_allow_html=True)
footer()

# Function for Nucleotide count page
def nucleotide_count_page():
#display_image2()
apply_background_color("#FFC0CB")  
st.markdown("<h2 style='text-align: center;'>Nucleotide Count</h2>", unsafe_allow_html=True)
sequence = st.text_area("Enter a DNA sequence:").upper()
if sequence:
    if is_valid_sequence(sequence):
        counts = get_nucleotide_count(sequence)
        st.markdown(f"<p style='text-align: center;'>A: {counts['A']}, T: {counts['T']}, C: {counts['C']}, G: {counts['G']}</p>", unsafe_allow_html=True)
    else:
        st.error("Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G.")
footer()
# Function for K-mer analysis page
def kmer_analysis_page():
#display_image3()
apply_background_color("#008000")  
st.markdown("<h2 style='text-align: center;'>K-mer Analysis</h2>", unsafe_allow_html=True)
sequence = st.text_area("Enter a DNA sequence:")
k = st.number_input("Enter the value of k:", min_value=1, step=1)
if sequence and k:
    if is_valid_sequence(sequence):
        kmers = kmer_analysis(sequence, int(k))
        st.markdown("<p style='text-align: center;'>K-mer analysis results:</p>", unsafe_allow_html=True)
        st.write(kmers)
    else:
        st.error("Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G.")
footer()
# Function for Gene Finding page
def gene_finding_page():
#display_image4()
apply_background_color("#00FF00")  
st.markdown("<h2 style='text-align: center;'>Gene Finding</h2>", unsafe_allow_html=True)
sequence = st.text_area("Enter a DNA sequence:")
if sequence:
    if is_valid_sequence(sequence):
        genes = find_genes(sequence.upper())
        st.markdown(f"<p style='text-align: center;'>Found genes: {genes}</p>", unsafe_allow_html=True)
    else:
        st.error("Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G.")
footer()
# Function for Hamming Distance page
def hamming_distance_page():
#display_image5()
apply_background_color("#D3D3D3")  
st.markdown("<h2 style='text-align: center;'>Hamming Distance</h2>", unsafe_allow_html=True)
seq1 = st.text_area("Enter the first DNA sequence:")
seq2 = st.text_area("Enter the second DNA sequence:")
if seq1 and seq2:
    if is_valid_sequence(sequence):
        result = hamming_distance(seq1.upper(), seq2.upper())
        st.markdown(f"<p style='text-align: center;'>Hamming Distance: {result}</p>", unsafe_allow_html=True)
    else:
        st.error("Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G.")
footer()
# Function for Reverse Complement page
def reverse_complement_page():
#display_image6()
apply_background_color("#FFD700")  
st.markdown("<h2 style='text-align: center;'>Reverse Complement</h2>", unsafe_allow_html=True)
sequence = st.text_area("Enter a DNA sequence:")
if sequence:
    if is_valid_sequence(sequence):
        complement = reverse_complement(sequence.upper())
        st.markdown(f"<p style='text-align: center;'>Reverse complement: {complement}</p>", unsafe_allow_html=True)
    else:
        st.error("Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G.")
footer()
# Function for GC Content page
def gc_content_page():
#display_image7()
apply_background_color("#A52A2A")  
st.markdown("<h2 style='text-align: center;'>GC Content</h2>", unsafe_allow_html=True)
sequence = st.text_area("Enter a DNA sequence:")
if sequence:
    if is_valid_sequence(sequence):
        gc_percent = gc_content(sequence.upper())
        st.markdown(f"<p style='text-align: center;'>GC Content: {gc_percent:.2f}%</p>", unsafe_allow_html=True)
    else:
        st.error("Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G.")
footer()
# Function for Transcription page
def transcription_page():
#display_image8()
apply_background_color("#FFA500")  
st.markdown("<h2 style='text-align: center;'>Transcription</h2>", unsafe_allow_html=True)
sequence = st.text_area("Enter a DNA sequence:")
if sequence:
    if is_valid_sequence(sequence):
        transcribed = transcription(sequence.upper())
        st.markdown(f"<p style='text-align: center;'>Transcribed RNA sequence: {transcribed}</p>", unsafe_allow_html=True)
    else:
        st.error("Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G.")
footer()
# Function for Translation page
def translation_page():
#display_image9()
apply_background_color("#FAD7A0")  
st.markdown("<h2 style='text-align: center;'>Translation</h2>", unsafe_allow_html=True)
sequence = st.text_area("Enter a DNA sequence:")
if sequence:
    if is_valid_sequence(sequence):
        protein = translation(sequence.upper())
        st.markdown(f"<p style='text-align: center;'>Translated protein sequence: {protein}</p>", unsafe_allow_html=True)
    else:
        st.error("Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G.")
footer()
# Function for Sequence Alignment page
def sequence_alignment_page():
#display_image10()
apply_background_color("#808000")  
st.markdown("<h2 style='text-align: center;'>Sequence Alignment</h2>", unsafe_allow_html=True)
seq1 = st.text_area("Enter the first DNA sequence:")
seq2 = st.text_area("Enter the second DNA sequence:")
if seq1 and seq2:
    if is_valid_sequence(sequence):
        result = sequence_alignment(seq1.upper(), seq2.upper())
        st.markdown(f"<p style='text-align: center;'>Sequence alignment result: {result}</p>", unsafe_allow_html=True)
    else:
        st.error("Invalid sequence! Please enter a DNA sequence containing only A, C, T, or G.")
footer()
# Function to control page navigation
def main():
# Sidebar navigation to select pages
page = st.sidebar.radio(
    "Select a page:",
    ("Title Page", "Nucleotide Count", "K-mer Analysis", "Gene Finding", "Hamming Distance", 
     "Reverse Complement", "GC Content", "Transcription", "Translation", 
     "Sequence Alignment", "About Us")
)
footer()

# Display selected page
if page == "Title Page":
    title_page()
elif page == "Nucleotide Count":
    nucleotide_count_page()
elif page == "K-mer Analysis":
    kmer_analysis_page()
elif page == "Gene Finding":
    gene_finding_page()
elif page == "Hamming Distance":
    hamming_distance_page()
elif page == "Reverse Complement":
    reverse_complement_page()
elif page == "GC Content":
    gc_content_page()
elif page == "Transcription":
    transcription_page()
elif page == "Translation":
    translation_page()
elif page == "Sequence Alignment":
    sequence_alignment_page()
elif page == "About Us":
    about_us_page()

if __name__ == "__main__":
    main()
