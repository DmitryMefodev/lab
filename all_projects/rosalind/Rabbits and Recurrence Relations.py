month = 28
offspring = 2

a = [[1], [1]]
i = 1
while i < month - 1:
    a.append([])
    for j in a[i]:
        if j == 1:
            a[i + 1].append(1)
            a[i + 1].append(offspring)
        if j == offspring:
            a[i + 1].extend([1] * offspring)
    i += 1

print(sum(a[-1]))