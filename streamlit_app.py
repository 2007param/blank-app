import streamlit as st

# Function for Nucleotide count
def get_nucleotide_count(sequence):
    return {'A': sequence.count('A'), 'T': sequence.count('T'),
            'C': sequence.count('C'), 'G': sequence.count('G')}

# Function for k-mer analysis
def kmer_analysis(sequence, k):
    kmers = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        kmers[kmer] = kmers.get(kmer, 0) + 1
    return kmers

# Function for Hamming distance
def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        return "Error: Sequences must have the same length"
    return sum(el1 != el2 for el1, el2 in zip(seq1, seq2))

# Function for Gene finding
def find_genes(sequence):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    genes = []
    for i in range(len(sequence)):
        if sequence[i:i + 3] == start_codon:
            for j in range(i + 3, len(sequence), 3):
                if sequence[j:j + 3] in stop_codons:
                    genes.append(sequence[i:j + 3])
                    break
    return genes

# Function for Reverse complement
def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))

# Function for GC Content
def gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

# Function for Transcription
def transcription(sequence):
    return sequence.replace('T', 'U')

# Function for Translation (simple case, codon-to-amino acid)
def translation(sequence):
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
    return protein

# Function for Sequence alignment (simple version)
def sequence_alignment(seq1, seq2):
    from difflib import SequenceMatcher
    matcher = SequenceMatcher(None, seq1, seq2)
    match_ratio = matcher.ratio()
    return f"Sequence similarity ratio: {match_ratio*100:.2f}%"

# Function to display the title page
def title_page():
    st.markdown("<h1 style='text-align: center;'>🎈 Bioinformatics Tool</h1>", unsafe_allow_html=True)
    st.markdown("<p style='text-align: center;'>Welcome to the Bioinformatics Tool. This app provides several bioinformatics analysis tools such as Nucleotide count, K-mer analysis, Gene Finding, and more. Select an option from the menu to get started!</p>", unsafe_allow_html=True)

# Function for About Us page
def about_us_page():
    st.markdown("<h2 style='text-align: center;'>About Us</h2>", unsafe_allow_html=True)
    st.markdown("<p style='text-align: center;'>We are a team of passionate bioinformaticians, biologists, and developers working to provide tools for DNA sequence analysis. Our goal is to make bioinformatics more accessible and provide a range of useful functionalities for researchers and students alike.</p>", unsafe_allow_html=True)
    st.markdown("<p style='text-align: center;'>This tool includes features such as nucleotide counting, k-mer analysis, gene finding, and more. We hope these tools will be helpful in your research and studies.</p>", unsafe_allow_html=True)

# Function for Nucleotide count page
def nucleotide_count_page():
    st.markdown("<h2 style='text-align: center;'>Nucleotide Count</h2>", unsafe_allow_html=True)
    sequence = st.text_area("Enter a DNA sequence:")
    if sequence:
        counts = get_nucleotide_count(sequence.upper())
        st.markdown(f"<p style='text-align: center;'>A: {counts['A']}, T: {counts['T']}, C: {counts['C']}, G: {counts['G']}</p>", unsafe_allow_html=True)

# Function for K-mer analysis page
def kmer_analysis_page():
    st.markdown("<h2 style='text-align: center;'>K-mer Analysis</h2>", unsafe_allow_html=True)
    sequence = st.text_area("Enter a DNA sequence:")
    k = st.number_input("Enter the value of k:", min_value=1, step=1)
    if sequence and k:
        kmers = kmer_analysis(sequence.upper(), k)
        st.markdown("<p style='text-align: center;'>K-mer analysis results:</p>", unsafe_allow_html=True)
        st.write(kmers)

# Function for Gene Finding page
def gene_finding_page():
    st.markdown("<h2 style='text-align: center;'>Gene Finding</h2>", unsafe_allow_html=True)
    sequence = st.text_area("Enter a DNA sequence:")
    if sequence:
        genes = find_genes(sequence.upper())
        st.markdown(f"<p style='text-align: center;'>Found genes: {genes}</p>", unsafe_allow_html=True)

# Function for Hamming Distance page
def hamming_distance_page():
    st.markdown("<h2 style='text-align: center;'>Hamming Distance</h2>", unsafe_allow_html=True)
    seq1 = st.text_area("Enter the first DNA sequence:")
    seq2 = st.text_area("Enter the second DNA sequence:")
    if seq1 and seq2:
        result = hamming_distance(seq1.upper(), seq2.upper())
        st.markdown(f"<p style='text-align: center;'>Hamming Distance: {result}</p>", unsafe_allow_html=True)

# Function for Reverse Complement page
def reverse_complement_page():
    st.markdown("<h2 style='text-align: center;'>Reverse Complement</h2>", unsafe_allow_html=True)
    sequence = st.text_area("Enter a DNA sequence:")
    if sequence:
        complement = reverse_complement(sequence.upper())
        st.markdown(f"<p style='text-align: center;'>Reverse complement: {complement}</p>", unsafe_allow_html=True)

# Function for GC Content page
def gc_content_page():
    st.markdown("<h2 style='text-align: center;'>GC Content</h2>", unsafe_allow_html=True)
    sequence = st.text_area("Enter a DNA sequence:")
    if sequence:
        gc_percent = gc_content(sequence.upper())
        st.markdown(f"<p style='text-align: center;'>GC Content: {gc_percent:.2f}%</p>", unsafe_allow_html=True)

# Function for Transcription page
def transcription_page():
    st.markdown("<h2 style='text-align: center;'>Transcription</h2>", unsafe_allow_html=True)
    sequence = st.text_area("Enter a DNA sequence:")
    if sequence:
        transcribed = transcription(sequence.upper())
        st.markdown(f"<p style='text-align: center;'>Transcribed RNA sequence: {transcribed}</p>", unsafe_allow_html=True)

# Function for Translation page
def translation_page():
    st.markdown("<h2 style='text-align: center;'>Translation</h2>", unsafe_allow_html=True)
    sequence = st.text_area("Enter a DNA sequence:")
    if sequence:
        protein = translation(sequence.upper())
        st.markdown(f"<p style='text-align: center;'>Translated protein sequence: {protein}</p>", unsafe_allow_html=True)

# Function for Sequence Alignment page
def sequence_alignment_page():
    st.markdown("<h2 style='text-align: center;'>Sequence Alignment</h2>", unsafe_allow_html=True)
    seq1 = st.text_area("Enter the first DNA sequence:")
    seq2 = st.text_area("Enter the second DNA sequence:")
    if seq1 and seq2:
        result = sequence_alignment(seq1.upper(), seq2.upper())
        st.markdown(f"<p style='text-align: center;'>Sequence alignment result: {result}</p>", unsafe_allow_html=True)

# Function to control page navigation
def main():
    # Sidebar navigation to select pages
    page = st.sidebar.radio(
        "Select a page:",
        ("Title Page", "About Us", "Nucleotide Count", "K-mer Analysis", "Gene Finding", "Hamming Distance", 
         "Reverse Complement", "GC Content", "Transcription", "Translation", 
         "Sequence Alignment")
    )

    # Display selected page
    if page == "Title Page":
        title_page()
    elif page == "About Us":
        about_us_page()
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

if __name__ == "__main__":
    main()
