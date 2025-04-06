with open('rosalind_gc (1).txt', 'r') as file:
    read_list = file.readlines()
    key = ''
    string_dict = {}
    for string in read_list:
        if '>' in string:
            key = string
            string_dict[key] = ''
            continue
        string_dict[key] += string.strip()

ans_key = 'a',
ans_perc = 0
for key, string in string_dict.items():
    if (string.count('C') + string.count('G'))/len(string) > ans_perc:
        ans_perc = (string.count('C') + string.count('G'))/len(string)
        ans_key = key
print(ans_key)
print('%.6f'%round(ans_perc * 100, 6))