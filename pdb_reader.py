
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

def read_pdb(file):
    with open(file, 'r') as f:
        start = False
        sequences = []
        sequence = []
        conect = []
        for line in f:
            line = line.strip('\n')
            post_line = line.split(' ')
            if post_line[0] == 'ATOM' and start == False:
                start = True
                continue
            if post_line[0] == 'ATOM' and start == True:
                line_tuple = (int(line[6:11].strip(' ')),dic[line[17:20]], int(line[22:26]))
                #print line_tuple
                sequence.append(line_tuple)
                continue
            if post_line[0] == 'TER':
                sequences.append(sequence)
                sequence = []
                continue
            if post_line[0] == 'CONECT':
                c = [int(x) for x in post_line if x != '' and x != 'CONECT']
                conect.append(c)
            #print post_line
    #print sequences[0]
    #print sequences[1]
    #print sequences[2]
    for sequence in sequences:
        last = sequence[-1]
        _,_,number = last
        pointer = 0
        sequence_string = ''
        for index in xrange(1, number + 1):
            while True:
                if sequence[pointer][2] == index:
                    sequence_string += sequence[pointer][1]
                    break
                pointer = pointer + 1
        print sequence_string
        print '-------------------------------' 
read_pdb('3HFM.pdb')