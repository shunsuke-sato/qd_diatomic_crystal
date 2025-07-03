import numpy as np
import matplotlib.pyplot as plt

def initialize_spins(N, n_p, n_n):
    """Initialize spins randomly to +1 or -1"""

    #spins = initialize_spins_read()
    #spins = initialize_spins_fibonacci(N, n_p, n_n)
    spins = initialize_spins_zero(N)



    #read_spins = read_spins_from_file("final_spins.dat")
    #spins = np.array(read_spins)
    #spins = np.zeros(N, dtype=int)
    #spins[0:n_p-1] = 1
    #spins[n_p:] = -1

    #for i in range(N):
    #    spins[i] = (-1)**i
    #
    return spins

def initialize_spins_read():
    """Initialize spins randomly to +1 or -1"""

    read_spins = read_spins_from_file("final_spins.dat")
    spins = np.array(read_spins)
    #spins = np.zeros(N, dtype=int)
    #spins[0:n_p-1] = 1
    #spins[n_p:] = -1

    #for i in range(N):
    #    spins[i] = (-1)**i
    #
    return spins


def initialize_spins_fibonacci(N, n_p, n_n):
    """Initialize spins randomly to +1 or -1"""

    spins = np.zeros(N, dtype=int)
    spins[0:n_p-1] = 1
    spins[n_p:] = -1

    #for i in range(N):
    #    spins[i] = (-1)**i
    #
    return spins

def initialize_spins_zero(N):
    """Initialize spins randomly to +1 or -1"""

    spins = np.zeros(N, dtype=int)
    for i in range(N):
        spins[i] = (-1)**i
    
    return spins

def energy(spins, J):
    """Compute the total energy of the 1D spin configuration with periodic boundary"""
    return -J * np.sum(spins * np.roll(spins, 1))

def metropolis_step(spins, J, beta):
    """Perform one Metropolis update step"""
    N = len(spins)
    for _ in range(N):
        i = np.random.randint(N)
        
        # Periodic boundary neighbors
        left_i = spins[i - 1]
        right_i = spins[(i + 1) % N]
        dE_i = 2 * J * spins[i] * (left_i + right_i)
        dE = dE_i

        # Metropolis criterion
        if dE < 0 or np.random.rand() < np.exp(-beta * dE):
            spins[i] *= -1

    return spins


def metropolis_step_exchange(spins, J, beta):
    """Perform one Metropolis update step"""
    N = len(spins)
    for _ in range(N):
        i = np.random.randint(N)
        j = np.random.randint(N)
        if(np.abs(i-j) <= 2):
            j = (j + N//2)%N

        if(spins[i]*spins[j] == -1):
            # Periodic boundary neighbors
            left_i = spins[i - 1]
            right_i = spins[(i + 1) % N]
            dE_i = 2 * J * spins[i] * (left_i + right_i)
            left_j = spins[j - 1]
            right_j = spins[(j + 1) % N]
            dE_j = 2 * J * spins[j] * (left_j + right_j)
            dE = dE_i + dE_j

            # Metropolis criterion
            if dE < 0 or np.random.rand() < np.exp(-beta * dE):
                spins[i] *= -1
                spins[j] *= -1
    return spins

def simulate(N, J, T, steps, n_p, n_n):
    beta = 1.0 / T
    spins = initialize_spins(N, n_p, n_n)
    energies = []

    for step in range(steps):
        spins = metropolis_step(spins, J, beta)
        if step % 10 == 0:
            E = energy(spins, J)
            energies.append(E / N)

    return spins, energies

def calc_correlation(spins):
    N = len(spins)
    corr = np.zeros(N)
    dist = np.zeros(N)

    for i in range(N):
        dist[i] = i
        corr[i] = (-1)**i*np.sum(spins * np.roll(spins, i))/N

    return dist, corr

def write_spins_to_file(filename, spins):


    N = len(spins)
    with open(filename, 'w') as f:
        f.write(f"{len(spins)}\n")
        for i in range(N):
            f.write(f"{spins[i]}\n")


    return


def read_spins_from_file(filename):
    with open(filename, 'r') as f:
        # Read the number of spins
        N = int(f.readline().strip())

        # Read the next N lines as spins
        spins = []
        for _ in range(N):
            spin = f.readline().strip()
            try:
                # Try to interpret as int, fallback to float if needed
                spin = int(spin)
            except ValueError:
                spin = float(spin)
            spins.append(spin)

    return spins

# Parameters
n_p = 1024
n_n = 1024
N = n_p + n_n           # Number of spins
J = -1.0          # Antiferromagnetic interaction (J < 0)
T = 0.4           # Temperature
steps = 4000      # Monte Carlo steps

final_spins, energy_trace = simulate(N, J, T, steps, n_p, n_n)
dist, correlation = calc_correlation(final_spins)

#for i in range(N):
#    if(i < N//2):
#        final_spins[i] = 1
#    else:
#        final_spins[i] = -1
write_spins_to_file("final_spins.dat", final_spins)


# Plot the energy over time
plt.plot(dist, correlation)
plt.xlabel('Length')
plt.ylabel('Correlation')
plt.show()
