codon_map = {
    'UUU': ['Phenylalanine', 'F', 'Phe'], 'UUC': ['Phenylalanine', 'F', 'Phe'], 'UUA': ['Leucine', 'L', 'Leu'], 'UUG': ['Leucine', 'L', 'Leu'], 
    'CUA': ['Leucine', 'L', 'Leu'], 'CUG': ['Leucine', 'L', 'Leu'], 'CUA': ['Leucine', 'L', 'Leu'], 'CUG': ['Leucine', 'L', 'Leu'], 
    'AUU': ['Isoleucine', 'I', 'Ile'], 'AUC': ['Isoleucine', 'I', 'Ile'], 'AUA': ['Isoleucine', 'I', 'Ile'], 'AUG': ['Methionine', 'M', 'Met'], 
    'GUU': ['Valine', 'V', 'Val'], 'GUC': ['Valine', 'V', 'Val'], 'GUA': ['Valine', 'V', 'Val'], 'GUG': ['Valine', 'V', 'Val'], 

    'UCU': ['Serine', 'S', 'Ser'], 'UCC': ['Serine', 'S', 'Ser'], 'UCA': ['Serine', 'S', 'Ser'], 'UCG': ['Serine', 'S', 'Ser'], 
    'CCA': ['Proline', 'P', 'Pro'], 'CCG': ['Proline', 'P', 'Pro'], 'CCA': ['Proline', 'P', 'Pro'], 'CCG': ['Proline', 'P', 'Pro'], 
    'ACU': ['Threonine', 'T', 'Thr'], 'ACC': ['Threonine', 'T', 'Thr'], 'ACA': ['Threonine', 'T', 'Thr'], 'ACG': ['Threonine', 'T', 'Thr'], 
    'GCU': ['Alanine', 'A', 'Ala'], 'GCC': ['Alanine', 'A', 'Ala'], 'GCA': ['Alanine', 'A', 'Ala'], 'GCG': ['Alanine', 'A', 'Ala'], 

    'UAU': ['Tyrosine', 'Y', 'Tyr'], 'UAC': ['Tyrosine', 'Y', 'Tyr'], 'UAA': ['Stop', 'end', 'Stop'], 'UAG': 'Stop [end]', 
    'CAA': ['Histidine', 'H', 'His'], 'CAG': ['Histidine', 'H', 'His'], 'CAA': ['Glutamine', 'Q', 'Gln'], 'CAG': ['Glutamine', 'Q', 'Gln'], 
    'AAU': ['Asparagine', 'N', 'Asn'], 'AAC': ['Asparagine', 'N', 'Asn'], 'AAA': ['Lysine', 'K', 'Lys'], 'AAG': ['Lysine', 'K', 'Lys'], 
    'GAU': ['Aspartic Acid', 'D', 'Asp'], 'GAC': ['Aspartic Acid', 'D', 'Asp'], 'GAA': ['Glutamic Acid', 'E', 'Glu'], 'GAG': ['Glutamic Acid', 'E', 'Glu'], 

    'UGU': ['Cysteine', 'C', 'Cys'], 'UGC': ['Cysteine', 'C', 'Cys'], 'UGA': ['Stop', 'end', 'Stop'], 'UGG': ['Tryptophan', 'W', 'Trp'], 
    'CGA': ['Arginine', 'R', 'Arg'], 'CGG': ['Arginine', 'R', 'Arg'], 'CGA': ['Arginine', 'R', 'Arg'], 'CGG': ['Arginine', 'R', 'Arg'], 
    'AGU': ['Serine', 'S', 'Ser'], 'AGC': ['Serine', 'S', 'Ser'], 'AGA': ['Arginine', 'R', 'Arg'], 'AGG': ['Arginine', 'R', 'Arg'], 
    'GGU': ['Glycine', 'G', 'Gly'], 'GGC': ['Glycine', 'G', 'Gly'], 'GGA': ['Glycine', 'G', 'Gly'], 'GGG': ['Glycine', 'G', 'Gly'], 
}

def translate_dna():
    seqm = input("Input DNA = ").upper()
    print(" ")

    if len(seqm) % 3 != 0:
        print("Nucleotides whould be in multiples of 3.")
        print(" ")
        return

    valid = seqm.count("A") + seqm.count("C") + seqm.count("T") + seqm.count("G")
    if not valid == len(seqm):
        print("Invalid Nucleotides.")
        print(" ")
        return

    complement_seq = ""
    mrna_seq = ""

    for i in (seqm):
        if i == "A":
            complement_seq += "T"
            mrna_seq += "U"
        elif i == "T":
            complement_seq += "A"
            mrna_seq += "A"
        elif i == "C":
            complement_seq += "G"
            mrna_seq += "G"
        elif i == "G":
            complement_seq += "C"
            mrna_seq += "C"

    print("Complement = ", complement_seq)
    print("mRNA = ", mrna_seq)

    amino_acids = []
    for j in range(0, len(mrna_seq), 3):
        codon = mrna_seq[j:j + 3]
        if codon in codon_map:
            if codon_map[codon] == 'Stop':
                amino_acids.append(codon_map[codon][2] + " (" + codon_map[codon][1] + ")")
                break
            amino_acids.append(codon_map[codon][2] + " (" + codon_map[codon][1] + ")")
    
    print("Aminoacid = ", " - ".join(amino_acids))
    print(" ")



def evaluate_aminoacid():
    aminoacid = input("Input Aminoacid = ").upper().strip()
    print(" ")

    if len(aminoacid) > 3:
        print("Maximum of 3 Aminoacids.")
        return

    valid_acids = {value[1] for value in codon_map.values()}
    if not set(aminoacid).issubset(valid_acids):
        print("Invalid Aminoacids.")
        return

    codon_lists = []
    for acid in aminoacid:
        codons = [key for key, value in codon_map.items() if value[1] == acid]
        codon_lists.append(codons)

    def generate_combinations(codon_lists):
        if not codon_lists:
            return ['']
        else:
            prefix_combinations = generate_combinations(codon_lists[1:])
            return [codon + chain for codon in codon_lists[0] for chain in prefix_combinations]

    mrna_combinations = generate_combinations(codon_lists)
    
    for combination in mrna_combinations:
        print("mRNA = " + combination)        
        codon_frequency = {}

        for i in range(0, len(combination), 3):
            codon = combination[i:i + 3]
            if codon in codon_frequency:
                codon_frequency[codon] += 1
            else:
                codon_frequency[codon] = 1

        for codon, count in codon_frequency.items():
            print(codon + " = " + str(count))
        print(" ")

translate_dna()
evaluate_aminoacid()