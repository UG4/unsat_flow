p_w = 998.23
p_s = 1195.0
w_max = 0.5448
for i in range(0, 101):
    c = i/100
    print(1/(1/p_s + (1/p_w - 1/p_s) * (c/w_max)))
