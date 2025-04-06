with open('rosalind_dna (1).txt', 'r') as file:
    s = file.readline()
    print(s.count('A'), s.count('C'), s.count('G'), s.count('T'))
