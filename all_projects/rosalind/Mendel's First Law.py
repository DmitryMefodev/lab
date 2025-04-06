k = 22
m = 17
n = 18
all_var = (k + m + n) * (k + n + m - 1) / 2
good_var = k * (m + n) + k * (k - 1) / 2 + (m * n) / 2 + m * (m - 1) * 3/8
print(good_var/all_var)