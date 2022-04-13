import sys

def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'N': 'N'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

if __name__ == "__main__":
    db_file = sys.argv[1]
    print("id,sequence,pam,chromosome,position,sense")
    with open(db_file, 'r') as f:
        for l in f:
            if l.startswith('@'):
                continue
    
            chrm = l.split('\t')[2].split('_')[0]
            pos = l.split('\t')[3]
            sense = '-' if l.split('\t')[1] == '16' else '+'
            #chrm, pos, sense = identifier.split(':')
            kmer = l.split('\t')[9]

            if sense == '-':
                kmer = revcom(kmer)
            
            pam = kmer[-3:]
            kmer = kmer[:-3]
            identifier = f'{chrm}:{pos}:{sense}'
            print(identifier + "," + kmer + "," + pam + "," + chrm + "," + pos + "," + sense)

