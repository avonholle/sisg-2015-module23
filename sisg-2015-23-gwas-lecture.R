
alpha = 5e-8
threshold = qchisq(1-alpha, 1)
q2 = 0.005
n = 10000
ncp = n*q2/(1-q2)
power = 1-pchisq(threshold, 1, ncp)
threshold
ncp
power

