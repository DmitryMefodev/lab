with open('rosalind_splc.txt', 'r') as file:
    read_list = file.readlines()
    key = ''
    string_dict = {}
    for string in read_list:
        if '>' in string:
            key = string.strip()
            string_dict[key] = ''
            continue
        string_dict[key] += string.strip()

string, introns = list(string_dict.keys())[0], list(string_dict.keys())[1::]
string_DNA = string_dict[string]
for i in introns:
    while string_DNA.find(string_dict[i]) != -1:
        ind = string_DNA.find(string_dict[i])
        string_DNA = string_DNA[:ind] + string_DNA[ind + len(string_dict[i]):]

string_RNA = string_DNA.replace('T', 'U')

with open('RNA codon table.txt', 'r') as file:
    read_list = file.readlines()
    codon_table = {}
    for s in read_list:
        string = s.split()
        for i in range(0, 7, 2):
            codon_table[string[i]] = string[i + 1]

stack = ''
ans = ''
for i in string_RNA:
    stack += i
    if stack in codon_table.keys():
        if codon_table[stack] == 'Stop':
            break
        ans += codon_table[stack]
        stack = ''

print(ans)