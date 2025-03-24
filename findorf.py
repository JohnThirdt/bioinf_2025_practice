START = 'ATG'
STOP  = ['TAA', 'TAG', 'TGA']

# ----------------------------------------------------translation-------------------------------------------------------------------
table = [
    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG',
    'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG',
    'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
    ]

tr_={'T':0, 'C': 1,  'A': 2,  'G': 3}
def translate(triple):

    number = 16*tr_.get(triple[0]) + 4*tr_.get(triple[1]) + tr_.get(triple[2])
    return table[0][number]
# -------------------------------------------------------------------------------------------------------------------------------------
# find start codon index
def find_start(seq, j):
    i=j
    while seq[i:i+3]!=START:
        i+=1
        if i>=len(seq)-3:
            return "Hidden message"
    return i

# get sequence for ORF + translated one
def get_orf(seq, i): # i -  START index
    j = i 
    orf=''
    orf_tr=''
    while not(seq[j:j+3] in STOP):
        if j>=len(seq)-3:
            return [None , None]
        triple =seq[j:j+3]
        orf+=triple
        orf_tr+=translate(triple)
        j+=3
    # triple =seq[j:j+3]
    # orf+=triple
    # orf_tr+=translate(triple)
    return [orf, orf_tr, j] # ORF + Translated ORF + Index of first base in STOP


path = 'ENSG00000107018.txt'
seqs = []
i=-1
F = False
B = False
with open(path) as f:
    while True:
        line=f.readline().rstrip()
        if line=='':  # Somwtimes void lines act strange
            if B:
                break
            B = True
            continue
        if line[0]=='>':
            if not('pep' in line):
                F = True
                seqs = seqs + ['']
                i+=1
            else:
                F =False
        else:
            if F:
                seqs[i]= seqs[i] + line
            

for i in range(len(seqs)):
    seqs[i] = ''.join(seqs[i])
    seqs[i].replace('\n', '')
f = open("OUTPUT.txt", "w")
count_seq = 0
for seq_fasta in seqs:
    count_seq+=1
    i=0
    counter = 1
    while i < len(seq_fasta)-3:
        r = find_start(seq_fasta, i)
        if r == 'Hidden message':break
        o = get_orf(seq_fasta, r)
        if o[1]!= None: # and len(o[0])>=30:
            result = (f"In sequence {count_seq} found res_num {counter}.ORF found at base {r+1} to base {o[2]}: {o[0]} \n \t Translated gene on ptoteome: {o[1]} \n")
            print(result)
            f.write(result)
            counter+=1

        i=r+1 # Next for  ATG...... >>> TG......
        # i=o[2] #It can solve problem with ORFFinder (to remove nested sequences) START1.....STOP1... >>> ... START2....STOP2>>>
