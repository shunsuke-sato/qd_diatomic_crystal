import numpy as np
import matplotlib.pyplot as plt



def fibonacci_substitution(spins):
# A -> +1, B -> -1 
    N = len(spins)
    new_spins = []

    for j in range(N):
        if spins[j] == 1:
            new_spins.append(1)
            new_spins.append(-1)
        else:
            new_spins.append(1)

    return new_spins

def fibonacci_sequence(spins, nstep):
    for j in range(nstep):
        spins = fibonacci_substitution(spins)
#        print(spins)

    return spins

def write_spins_to_file(filename, spins):


    N = len(spins)
    n_mean = np.sum(spins)
    n_eq = (N-n_mean)//2
    n_p = n_eq + n_mean
    n_n = n_eq
    print(N,n_p,n_n)
    with open(filename, 'w') as f:
        f.write(f"{len(spins)}\n")
        for i in range(N):
            f.write(f"{spins[i]}\n")


    return

def calc_correlation(spins):
    N = len(spins)
    corr = np.zeros(N)
    dist = np.zeros(N)

    for i in range(N):
        dist[i] = i
        corr[i] = (-1)**i*np.sum(spins * np.roll(spins, i))/N

    return dist, corr


spins = np.array([-1])
nstep = 16



final_spins = fibonacci_sequence(spins, nstep)
dist, correlation = calc_correlation(final_spins)
write_spins_to_file("final_spins.dat", final_spins)

# Plot the energy over time
#plt.plot(dist, np.abs(correlation))
#plt.xlabel('Length')
#plt.ylabel('Correlation')

#plt.show()
