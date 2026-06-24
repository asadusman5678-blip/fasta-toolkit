import sys

if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = "dna.example.fasta"

try:
    f = open(filename)
except FileNotFoundError:
    print(f"File {filename} does not exist!!")

seqs = {}
for line in f:
    line = line.rstrip()
    if line[0] == '>':
        words = line.split()
        name = words[0][1:]
        seqs[name] = ''
    else:
        seqs[name] = seqs[name] + line

f.close()
print("-- FASTA file read successfully --")
print("-- list of records in the file --")
print("----------------------------------------------------------------------------------------")

for name, seq in seqs.items():
    print(name, seq)
print("----------------------------------------------------------------------------------------")
"""
answering part 1 How many records are in the file? 
"""
print(f"Number of records in the file: {len(seqs)}")
print("----------------------------------------------------------------------------------------")
"""answering part 2 What is the length of each sequence?
    What is the longest sequence and what is the shortest sequence?
    Is there more than one longest or shortest sequence?
    What are their identifiers?
    """
seqs_lengths ={}
for name, seq in seqs.items():
    seqs_lengths[name] = len(seq)
    print(f"Length of sequence {name}: {len(seq)}")
print("----------------------------------------------------------------------------------------")
longest_length = max(len(seq) for seq in seqs.values())
shortest_length = min(len(seq) for seq in seqs.values())



longest_sequences = [name for name, length in seqs_lengths.items() if length == longest_length]
shortest_sequences = [name for name, length in seqs_lengths.items() if length == shortest_length]

print("are there more than one longest or shortest sequence? What are their identifiers?")
if len(longest_sequences) > 1:
    print(f"There are {len(longest_sequences)} longest sequences")
else:
    print("There is only one longest sequence")
if len(shortest_sequences) > 1:
    print(f"There are {len(shortest_sequences)} shortest sequences")
else:
    print("There is only one shortest sequence")

print(f"Longest sequence(s): {', '.join(longest_sequences)} (Length: {longest_length})")
print(f"Shortest sequence(s): {', '.join(shortest_sequences)} (Length: {shortest_length})")
print("----------------------------------------------------------------------------------------")
"""
answering  part 3
In molecular biology, a reading frame is a way of dividing the DNA sequence of nucleotides into a set of consecutive, non-overlapping triplets (or codons). Depending on where we start, there are six possible reading frames: three in the forward (5' to 3') direction and three in the reverse (3' to 5'). For instance, the three possible forward reading frames for the sequence AGGTGACACCGCAAGCCTTATATTAGC are: 

AGG TGA CAC CGC AAG CCT TAT ATT AGC

A GGT GAC ACC GCA AGC CTT ATA TTA GC

AG GTG ACA CCG CAA GCC TTA TAT TAG C 

These are called reading frames 1, 2, and 3 respectively. An open reading frame (ORF) is the part of a reading frame that has the potential to encode a protein. It starts with a start codon (ATG), and ends with a stop codon (TAA, TAG or TGA). For instance, ATGAAATAG is an ORF of length 9.

Given an input reading frame on the forward strand (1, 2, or 3) your program should be able to identify all ORFs present in each sequence of the FASTA file, and answer the following questions: what is the length of the longest ORF in the file? What is the identifier of the sequence containing the longest ORF? For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier? What is the starting position of the longest ORF in the sequence that contains it? The position should indicate the character number in the sequence. For instance, the following ORF in reading frame 1:

>sequence1

ATGCCCTAG

starts at position 1.

Note that because the following sequence:

>sequence2

ATGAAAAAA

does not have any stop codon in reading frame 1, we do not consider it to be an ORF in reading frame 1
    """
# --- 1. Define the ORF Finder Function ---
def get_orfs(sequence, frame):
    """
    Scans a sequence in a specific forward frame (1, 2, or 3).
    Returns a list of tuples: (orf_sequence, start_position_1_indexed, length)
    """
    orfs = []
    start_index = frame - 1  # Convert frame 1,2,3 to Python index 0,1,2
    
    # Extract codons and track their 1-based character positions
    codons = []
    positions = []
    for i in range(start_index, len(sequence) - 2, 3):
        codons.append(sequence[i:i+3])
        positions.append(i + 1)
        
    in_orf = False
    orf_start = None
    current_orf = []
    
    for codon, pos in zip(codons, positions):
        if not in_orf:
            if codon == "ATG":
                in_orf = True
                orf_start = pos
                current_orf.append(codon)
        else:
            current_orf.append(codon)
            if codon in {"TAA", "TAG", "TGA"}:
                orf_seq = "".join(current_orf)
                orfs.append((orf_seq, orf_start, len(orf_seq)))
                in_orf = False
                current_orf = []
                orf_start = None
                
    return orfs


# ==========================================
# 2. Get and Validate User Inputs
# ==========================================

# --- Input 1: Validate Reading Frame ---
while True:
    frame_input = input("Enter a forward reading frame (1, 2, or 3): ").strip()
    if frame_input in {"1", "2", "3"}:
        chosen_frame = int(frame_input)
        break
    else:
        print("Invalid input! Please enter exactly 1, 2, or 3.\n")

# --- Input 2: Validate Sequence Identifier ---
while True:
    target_id = input("Enter a specific sequence identifier to look up: ").strip()
    if target_id in seqs:
        break
    else:
        print(f"Error: '{target_id}' was not found in your FASTA file.")
        print(f"Available identifiers are: {', '.join(list(seqs.keys())[:5])}... (showing first 5)\n")


# --- 3. Process Your 'seqs' Dictionary ---
all_found_orfs = {}

for name, seq in seqs.items():
    found = get_orfs(seq, chosen_frame)
    if found:
        all_found_orfs[name] = found


# --- 4. Extract Answers to Your Questions ---
if not all_found_orfs:
    print(f"\nNo valid ORFs found anywhere in the file for reading frame {chosen_frame}.")
else:
    max_length = 0
    max_id = None
    max_start = None
    max_seq = ""

    for name, orf_list in all_found_orfs.items():
        for orf_seq, start, length in orf_list:
            if length > max_length:
                max_length = length
                max_id = name
                max_start = start
                max_seq = orf_seq

    # Display Answers
    print(f"\n=== ANALYSIS FOR READING FRAME {chosen_frame} ===")
    print(f"Q1: What is the length of the longest ORF in the file?\n--> {max_length} base pairs")
    print(f"Q2: What is the identifier of the sequence containing the longest ORF?\n--> {max_id}")
    print(f"Q4: What is the starting position of the longest ORF in that sequence?\n--> Position {max_start}")
    
    print("\n--------------------------------------------------")
    print(f"Q3: Specific lookup for identifier '{target_id}':")
    if target_id in all_found_orfs:
        longest_local = max(all_found_orfs[target_id], key=lambda x: x[2])
        print(f"--> Longest ORF sequence: {longest_local[0]}")
        print(f"--> Length: {longest_local[2]}")
        print(f"--> Starts at character position: {longest_local[1]}")
    else:
        print(f"--> No valid ORFs found inside sequence '{target_id}' for frame {chosen_frame}.")

print("----------------------------------------------------------------------------------------")
# --- Helper Function to Get the Reverse Complement ---
def get_reverse_complement(seq):
    """
    Translates a DNA sequence to its complement and reverses it.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    # Replace each character with its complement, then reverse the string
    complement_seq = [complement.get(base, base) for base in seq]
    return "".join(complement_seq)[::-1]


# ==========================================
# 5. Get and Validate Repeat Length (n)
# ==========================================
while True:
    n_input = input("Enter the length of the repeats to look for (n): ").strip()
    if n_input.isdigit() and int(n_input) > 0:
        n = int(n_input)
        break
    else:
        print("Invalid input! Please enter a positive integer greater than 0.\n")


# ==========================================
# 6. Process and Count Repeats on Both Strands
# ==========================================
# Dictionaries to track occurrences separately and combined
fwd_counts = {}
rev_counts = {}
combined_counts = {}

for name, fwd_seq in seqs.items():
    rev_seq = get_reverse_complement(fwd_seq)
    
    # 1. Scan Forward Strand
    for i in range(len(fwd_seq) - n + 1):
        fwd_sub = fwd_seq[i:i+n]
        fwd_counts[fwd_sub] = fwd_counts.get(fwd_sub, 0) + 1
        combined_counts[fwd_sub] = combined_counts.get(fwd_sub, 0) + 1
        
    # 2. Scan Reverse Strand
    for i in range(len(rev_seq) - n + 1):
        rev_sub = rev_seq[i:i+n]
        rev_counts[rev_sub] = rev_counts.get(rev_sub, 0) + 1
        combined_counts[rev_sub] = combined_counts.get(rev_sub, 0) + 1


# ==========================================
# 7. Analyze and Print Results
# ==========================================
def analyze_strand_results(counts_dict, title):
    """
    Helper function to filter repeats (>1 occurrence) and find the maximums.
    """
    print(f"\n=== {title} ===")
    
    # Repeats must occur more than once
    repeats = {rep: count for rep, count in counts_dict.items() if count > 1}
    
    if not repeats:
        print(f"No repeats of length {n} found.")
        return
        
    print("Detected repeats and frequencies:")
    for rep_seq, count in repeats.items():
        print(f"--> {rep_seq}: {count} times")
        
    # Find most frequent
    max_freq = max(repeats.values())
    most_freq = [rep for rep, count in repeats.items() if count == max_freq]
    
    print(f"Most frequent pattern(s): {most_freq} (occurring {max_freq} times)")


# Print the 3 separate reports requested
analyze_strand_results(fwd_counts, f"FORWARD STRAND ANALYSIS (n={n})")
analyze_strand_results(rev_counts, f"REVERSE STRAND ANALYSIS (n={n})")
analyze_strand_results(combined_counts, f"COMBINED STRANDS ANALYSIS (n={n})")
