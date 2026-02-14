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
    
    alignpls = NeedlemanWunsch("./data/BLOSUM62.txt", gap_open=-10, gap_extend=-1)
    
    # animal vs human
    gg_score, gg_align, hs_align_gg = alignpls.align(hs_seq, gg_seq)
    mm_score, mm_align, hs_align_mm = alignpls.align(hs_seq, mm_seq)
    br_score, br_align, hs_align_br = alignpls.align(hs_seq, br_seq)
    tt_score, tt_align, hs_align_tt = alignpls.align(hs_seq, tt_seq)
    
    # tuple of name, score
    species_scores = [
        ("Gallus gallus", gg_score),
        ("Mus musculus", mm_score),
        ("Balaeniceps rex", br_score),
        ("Tursiops truncatus", tt_score)
    ]
    
    # descending since highest = most similar 
    species_scores.sort(key=lambda x: x[1], reverse=True)

    print("Species by similarity to human in descending order:")
    for i in range(len(species_scores)):  # ← FIXED: changed all_species to species_scores
        species_name = species_scores[i][0]  # ← FIXED
        score = species_scores[i][1]  # ← FIXED
        print(str(i+1) + ". " + species_name + ": " + str(score))

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    print("Alignment scores between each species BRD2 and human BRD2:")
    print("Gallus gallus vs Human: " + str(gg_score))
    print("Mus musculus vs Human: " + str(mm_score))
    print("Balaeniceps rex vs Human: " + str(br_score))
    print("Tursiops truncatus vs Human: " + str(tt_score))

    

if __name__ == "__main__":
    main()