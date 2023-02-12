# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")  

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    # Initialize NW
    NW = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat",
                         gap_open = -10,
                         gap_extend = -1)
    # Calculate alignment scores
    gg_score, gg_hs, gg = NW.align(hs_seq, gg_seq)
    mm_score, mm_hs, ms = NW.align(hs_seq, mm_seq)
    br_score, br_hs, bs = NW.align(hs_seq, br_seq)
    tt_score, tt_hs, tt = NW.align(hs_seq, tt_seq)

    speices_scores = {"GG" : gg_score,
                      "MM" : mm_score,
                      "BR" : br_score,
                      "TT" : tt_score}  
    # Sort the species scores
    species_scores_sorted = sorted(speices_scores.items(),
                                   key=lambda x:x[1])
    
    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print("Alignment scores of species to human BRD2:")
    for species, scores in species_scores_sorted:
        print(species + ": " + str(scores))
    
if __name__ == "__main__":
    main()
