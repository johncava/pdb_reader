from pdb_reader import *

sequences, conect = read_pdb('3HFM.pdb')
r_sequences = reconstruct_sequences(sequences)
for seq in r_sequences:
    print seq
    print "------------------------"