

# PDB format has each ATOM line corresponding to its associated amino acid as a 3-letter code.
# dictionary is used to convert this 3-letter code into easier 1-letter code for future use.
dic = {
    'CYS': 'C',
    'MET': 'M',
    'ILE': 'I',
    'VAL': 'V',
    'GLN': 'Q',
    'LYS': 'K',
    'PRO': 'P',
    'THR': 'T',
    'PHE': 'F',
    'ALA': 'A',
    'HIS': 'H',
    'GLY': 'G',
    'ASP': 'D',
    'GLU': 'E',
    'LEU': 'L',
    'ARG': 'R',
    'TRP': 'W',
    'ASN': 'N',
    'TYR': 'Y',
    'SER': 'S',
}


#####
# Takes in a '*.pdb' file and returns information to reconstruct the different sequences in a protein and the different connection information 
#####
def read_pdb(file):
    # Opens the file
    with open(file, 'r') as f:
        # Sentinel value until it comes to the first line with 'ATOM'
        start = False
        # Array to hold the different sequences
        sequences = []
        # Array to hold the sequence of amino acids
        sequence = []
        # Array to hold the information regarding which atom connects with which other atom
        conect = []
        # For each line in the pdb file...
        for line in f:
            # Strip away the new line token
            line = line.strip('\n')
            # Post-modify the line into an array where the space is the separator
            post_line = line.split(' ')
            # Checks to see if we have hit the first line with 'ATOM' beginning the transcription process of the protein
            if post_line[0] == 'ATOM' and start == False:
                # Set sentinel to True for future use
                start = True
                continue
            # Now that we are in the important information of the protein...
            if post_line[0] == 'ATOM' and start == True:
                # Take the essential information
                # TUPLE FORMAT: (unique ATOM ID, Amino Acid Residue, residue ID for this particular sequence)
                line_tuple = (int(line[6:11].strip(' ')),dic[line[17:20]], int(line[22:26]))
                # Add tuple to the sequence array
                sequence.append(line_tuple)
                continue
            # When you hit the TERMINAL state of the sequence
            if post_line[0] == 'TER':
                # Sequence is finished, then append sequence to the list of sequences
                sequences.append(sequence)
                # Reset sequence array to empty array for next sequence in the protein
                sequence = []
                continue
            # When you hit the CONECT information of the pdb file
            if post_line[0] == 'CONECT':
                # Get the information from the line that is not an empty string and not the token 'CONECT'
                # CONECT FORMAT: [ATOM1, ATOM2]
                c = [int(x) for x in post_line if x != '' and x != 'CONECT']
                # Append that information to the conect array
                conect.append(c)
    return sequences, conect

# Takes in an array of sequence information in the form of (unique ATOM ID, Amino Acid Residue, residue ID)
# Returns a string respresentation (1-letter code) of each sequence
def reconstruct_sequences(sequences):
    # Array to hold the reconstructed sequences
    reconstructed_sequences = []
    # For each sequence..
    for sequence in sequences:
        # Take the last tuple from the sequence
        last = sequence[-1]
        # Get the residue ID of that last tuple to get the length of the sequence
        _,_,number = last
        # Have a pointer in order to search for the next residue ID
        pointer = 0
        # Empty string to reconstruct the sequence
        sequence_string = ''
        # NOTE: There is a unique ATOM ID, but each amino acid will have a different amount of atoms. Thus each residue ID will have multiple ATOMs
        # As such, in order to differentiate a new residue between the next, we need to traverse to the next unique residue ID
        # We traverse through each tuple that has a unique ATOM ID, and see if a new unique residue ID has been seen => if sequence[pointer][2] == index
        # If that is the case, we take the tuple's second element (residue 1-letter code) and append it to the sequnece string => sequence_string += sequence[pointer][1]
        for index in xrange(1, number + 1):
            while True:
                if sequence[pointer][2] == index:
                    sequence_string += sequence[pointer][1]
                    break
                pointer = pointer + 1
        reconstructed_sequences.append(sequence_string)
    return reconstructed_sequences