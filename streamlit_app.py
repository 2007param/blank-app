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
