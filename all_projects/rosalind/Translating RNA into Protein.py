with open('RNA codon table.txt', 'r') as file:
    read_list = file.readlines()
    codon_table = {}
    for s in read_list:
        string = s.split()
        for i in range(0, 7, 2):
            codon_table[string[i]] = string[i + 1]

with open('rosalind_prot.txt', 'r') as file:
    s = file.readline()
    stack = ''
    ans = ''
    for i in s:
        stack += i
        if stack in codon_table.keys():
            if codon_table[stack] == 'Stop':
                break
            ans += codon_table[stack]
            stack = ''

print(ans)
