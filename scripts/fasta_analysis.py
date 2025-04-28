from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.Align import PairwiseAligner

def read_fasta(file_path):
    return {record.id: record.seq for record in SeqIO.parse(file_path, "fasta")}

def basic_stats(sequence):
    return {
        "length": len(sequence),
        "gc_content": gc_fraction(sequence),
        "molecular_weight": molecular_weight(sequence, "DNA"),
        "nucleotide_counts": {nt: sequence.count(nt) for nt in "ATCG"}
    }

def translate_sequence(dna_sequence, to_stop=True):
    return dna_sequence.translate(to_stop=to_stop, table=1)

def find_motifs(sequence, motif_list):
    return {
        motif: [i for i in range(len(sequence)) if sequence.startswith(Seq(motif), i)]
        for motif in motif_list
    }

def align_sequences(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    return aligner.align(seq1, seq2)

if __name__ == "__main__":
    sequences = read_fasta("data/ls_orchid.fasta")
    
    # Ejemplo de uso
    seq1, seq2 = list(sequences.values())[:2]
    print(basic_stats(seq1))
    
    # Alineamiento
    alignments = align_sequences(seq1, seq2)
    print("\nMejor alineamiento:")
    print(alignments[0])