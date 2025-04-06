with open('rosalind_subs (4).txt', 'r') as file:
    s = file.readlines()
    s1 = s[0].strip()
    s2 = s[1].strip()
    ans = []
    for i in range(len(s1)):
        if s1[i:len(s2) + i] == s2:
            ans.append(i + 1)

    print(' '.join(map(str, ans)))
