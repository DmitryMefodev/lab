with open('rosalind_rna.txt') as file:
    s = file.readline()
    print(s.replace('T', 'U'))