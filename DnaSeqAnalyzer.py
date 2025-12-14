# GC content calculator
sequence = input("Enter your Genomic sequence: ")
sequence = "".join(sequence.split()) # remove any spaces if available
sequence = sequence.upper() # handle the input if it were in lower case

# define the variables to be caluclated in the loop
count_G = 0
count_T = 0
count_A = 0
count_C = 0
total_length = 0

for i in sequence:
    if (i == "U"):
        print("Enter The DNA sequence NOT RNA")
        exit()
    elif (i!= "A" and i!= "T" and i!= "G" and i!= "C"):
        print("ENTER A VALID FUCKING SEQUENCE CONTAINING ONLY ATGC")
        exit()
    #elif (i == " "):
    #    print("there is a fucking empty space i told you not to include space")
    #    exit()
    elif (i == "A"):
        count_A +=1
    elif (i == "T"):
        count_T += 1
    elif (i == "G"):
        count_G += 1
    elif (i == "C"):
        count_C += 1
total_length = int(count_G + count_T + count_A + count_C)
print("\n{}, {}, {}, {} of A, T, G, C bases ,respectively,. Total of {} bases.".format(count_A, count_T, count_G, count_C, total_length))
percentage_of_GC = ((count_G + count_C)/(total_length)) * 100
percentage_of_GC = round(percentage_of_GC, 2)
print("\nThe GC content is: {}%.".format(percentage_of_GC))

# AT/GC skew
pG = (count_G/total_length)*100
pC = (count_C/total_length)*100
pA = (count_A/total_length)*100
pT = (count_T/total_length)*100

GC_skew = round((pG-pC)/(pG+pC),2)
AT_skew = round((pA-pT)/(pA+pT),2)

print("GC skew is: {}% and AT skew is: {}%".format(GC_skew,AT_skew))

# create the complementray strand of the sequence
replacements = str.maketrans({"A": "T", "T": "A", "G": "C", "C": "G"})
res = sequence.translate(replacements)
print("\nThe complementary strand of your sequence is: {}.".format(res))


# Represent the original sequence with its complementary sequence with a | in between 
original_strand = " ".join(sequence) # add a spcace for better visualization
linkage = " ".join("|" for _ in sequence)
complementary_strand = " ".join(res)

print("\nO: " + original_strand) # \n adds a space before the line of the top strand
print("   " + linkage) # 3 spaces because of "O: " and "C: " since they also have 3 spcaces
print("C: " + complementary_strand)

# Central Dogma

# TRANSCRIPTION
DNA_to_mRNA_replacement = str.maketrans({"T": "U"})
mRNA = sequence.translate(DNA_to_mRNA_replacement)
print("\nYour mRNA sequence is: {}".format(mRNA))


# TRANSLATION
codon_table = {
    # U
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y", "UAA": "_", "UAG": "_",
    "UGU": "C", "UGC": "C", "UGA": "_", "UGG": "W",

    # C
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    # A
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",

    # G
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
} # Codon table for Amino Acid mapping

codons = [mRNA[i:i+3] for i in range(0, len(mRNA), 3)] # extract the codons in triplets
if (len(codons[-1]) < 3):
    codons.pop() # remove the remaining bases if not a codon
    
print("\nYour Codons are: {}".format(codons))
print("\nYour Amino Acids Chain: ", end ="")

Amino_Acids_list = [] # initiate an empty Amino Acid list

# the mapping from RNA to Amino Acids
for i in codons:
    i = str(i)
    Amino_Acid = codon_table[i]
    if (Amino_Acid == "_"): # this is the stop codon
        print("\nTranslation terminated at a Stop codon ({}) found on the {}'th codon".format(i, codons.index(i) + 1),end ="")
        break
    else:
        print(Amino_Acid, end="")
        Amino_Acids_list.append(Amino_Acid)
      
Amino_Acid_weights ={
    "A" :89.1, "R" :174.1, "N" :132.1, "D" :133.1, "C" :121.2,
    "E" :147.1, "Q" :146.2, "G" :75.1, "H" :155.2, "I" :131.2,
    "L" :131.2, "K" :146.2, "M" :149.2, "F" :165.2, "P" :115.1,
    "S" :105.1, "T" :119.1, "W" :204.2, "Y" :181.2, "V" :117.1} #dictionary for amino acids weight

weights = 0 # initiate the weight list
print("\nThe weight of your Amino Acids Chain is: ", end="")
for i in Amino_Acids_list:
    i = str(i)
    weights += Amino_Acid_weights[i]
print(str(round(weights,2)) + " g/mol", end="")


# Since water is removed due to the link happenening between each amino acid,
# the weight of water should be considered and removed respectively
corrected_weight = round(weights - ((len(Amino_Acids_list)-1)*18.015), 2)
print("\nThe weight after removing water molecules due to Hydrolysis is: {} g/mol".format(corrected_weight))

Hydrophobic_list = ["I","V","L","F","C","M","A","W"]
Neutral_list = ["G","T","S","Y","P","H"]
Hydrophilic_list = ["N","D","Q","E","K","R"]


# Hydrophobicity of the produced Amino Acid chain 
hydrophobic = 0
hydrophilic = 0
neutral = 0
for i in Amino_Acids_list:
    if i in Hydrophobic_list:
        hydrophobic += 1
    elif i in Hydrophilic_list:
        hydrophilic += 1
    elif i in Neutral_list:
        neutral += 1
if hydrophobic > hydrophilic and hydrophobic > neutral:
    print("Your chain has more hydrophobic Amino Acids")
elif hydrophilic > hydrophobic and hydrophilic > neutral:
    print("Your chain is has more hydrophilic Amino Acids")
elif neutral > hydrophilic and neutral > hydrophobic:
    print("Your chain has more neutral Amino Acids in terms of hydrophobicity ") 

print ("\nYou have {} hydrphobic, {} hydrophilic and {} neutral Amino Acids in the chain".format(hydrophobic,hydrophilic,neutral))

# Kyte-Doolittle score for hydrophibicity
hydrophibicity = {
    'A': 1.8,
    'R': -4.5,
    'N': -3.5,
    'D': -3.5,
    'C': 2.5,
    'Q': -3.5,
    'E': -3.5,
    'G': -0.4,
    'H': -3.2,
    'I': 4.5,
    'L': 3.8,
    'K': -3.9,
    'M': 1.9,
    'F': 2.8,
    'P': -1.6,
    'S': -0.8,
    'T': -0.7,
    'W': -0.9,
    'Y': -1.3,
    'V': 4.2
}

hydrophibicity_score = 0
for i in Amino_Acids_list:
    hydrophibicity_score += hydrophibicity[i]
hydrophibicity_score = round(hydrophibicity_score, 2)
if hydrophibicity_score < -1:
    print("The chain is Hydrophilic with a K&D score of: {}".format(hydrophibicity_score))
elif hydrophibicity_score > 1:
    print("The chain is Hydrophobic with a K&D score of: {}".format(hydrophibicity_score))
else:
    print("The chain is Neutral with a K&D score of: {}".format(hydrophibicity_score))
