months = 80
dead = 20

pairs = {dead - i: 0 for i in range(dead)}
pairs[dead] = 1
for i in range(months):
    last_pairs = pairs[1]
    for j in range(1, dead + 1):
        if j <= dead - 1:
            pairs[j] = pairs[j + 1]
        else:
            new_pairs = sum([pairs[k] for k in range(1, dead - 1)]) + last_pairs
            pairs[dead] = new_pairs

print(sum(pairs.values()) - pairs[1])