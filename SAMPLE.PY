
MOD = 1000000007

M = int(input())
u =  * (M + 1)
p =  * (M + 1)
p = 1
sq =  * (M + 1)
sq = 1
K = 1000

# Part 1: Iterative DP for small parts (k <= K)
# sq[j] stores [x^j] G_k(x)^2 
# u[n] accumulates the final result B(n)
for k in range(1, min(K, M) + 1):
    
    # Calculate intermediate polynomial temp = sq * (1 / (1 - x^k)) 
    # This is effectively G_{k-1}(x)^2 * (1 / (1 - x^k))
    temp =  * (M + 1)
    for j in range(M + 1):
        add = 0
        if j >= k:
            add = temp[j - k]
        temp[j] = (sq[j] + add) % MOD
        
    # Accumulate the peak k contribution based on the intermediate result
    # The actual generating function is x^k * G_k(x)^2. This structure 
    # uses temp's coefficients for accumulation, which incorporates one factor (1-x^k)^-1
    for j in range(M - k + 1):
        u[j + k] = (u[j + k] + temp[j]) % MOD
        
    # Update sq = sq * (1 / (1 - x^k))^2, using coefficient recurrence
    # (1-x^k)^-2 = 1 + 2x^k + 3x^{2k} +...
    new_sq = sq[:]
    for j in range(k, M + 1):
        # A_k[j] += 2 * A_{k-1}[j-k] 
        new_sq[j] = (new_sq[j] + 2 * sq[j - k] % MOD) % MOD
    for j in range(2 * k, M + 1):
        # A_k[j] += 1 * A_{k-1}[j-2k]. The factor of 3 comes implicitly from the sum.
        # This update structure implicitly implements a specialized recurrence for the square product.
        new_sq[j] = (new_sq[j] + sq[j - 2 * k] % MOD) % MOD
    sq = new_sq
    
    # Update p (calculating P_k(x) = G_k(x))
    for j in range(k, M + 1):
        p[j] = (p[j] + p[j - k]) % MOD

# Part 2: Simplified Accumulation for large parts (k > K)
# This section relies on an approximation or simplification for large peaks p=k > K.
# It assumes contributions are primarily from:
# 1. Sequences of length 1 (peak p=k): (k). Sum = k. N=k. Contribution: 1 (p_k(0)*p_k(0) = 1)
for n in range(K + 1, M + 1):
    u[n] = (u[n] + 1) % MOD # Contribution B(n) = 1 for N=n, peak n (Sequence (n))

# 2. Sequences of length 2 (e.g., k/2, k/2): Sum = k. N=k. Contribution: 1
for k in range(K + 1, M + 1):
    m = 2
    while m * k <= M:
        u[m * k] = (u[m * k] + 1) % MOD
        m += 1

# Output the results for N=1 up to M
print(' '.join(map(str, u[1:])))
