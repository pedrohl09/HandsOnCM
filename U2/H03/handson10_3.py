# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 09:39:14 2024

@author: pedro
"""

import numpy as np
import matplotlib.pyplot as plt

# Parâmetros
n_bits = 1000             # Número de bits
T = 500                   # Tempo de símbolo OFDM
Ts = 2                    # Tempo de símbolo em portadora única
K = T // Ts               # Número de subportadoras independentes
N = 2 * K                 # DFT de N pontos
sigmas = [0, 0.1, 1]      # Vetor de variâncias do ruído

# Gerar bits aleatórios
dataIn = np.random.rand(n_bits)
dataIn = np.sign(dataIn - 0.5)  # Sequência de -1 e 1

# Conversor serial paralelo
dataInMatrix = dataIn.reshape(n_bits // 4, 4)

# Gerar constelação 16-QAM
seq16qam = (2 * dataInMatrix[:, 0] + dataInMatrix[:, 1] +
            1j * (2 * dataInMatrix[:, 2] + dataInMatrix[:, 3]))
seq16 = seq16qam.T

# Garantir propriedade da simetria
X = np.concatenate([seq16, np.conj(seq16[::-1])])

# Construindo xn
xn = np.zeros(N, dtype=complex)
for n in range(N):
    for k in range(N):
        xn[n] += (1 / np.sqrt(N)) * X[k] * np.exp(1j * 2 * np.pi * n * k / N)

# Loop de variâncias
for variance in sigmas:
    # Adição de ruído
    noise = np.sqrt(variance) * (np.random.randn(N) + 1j * np.random.randn(N))
    # Sinal recebido = xn + ruído
    rn = xn + noise

    # DFT de rn
    Y = np.zeros(K, dtype=complex)
    for k in range(K):
        for n in range(N):
            Y[k] += (1 / np.sqrt(N)) * rn[n] * np.exp(-1j * 2 * np.pi * k * n / N)

    # Plots
    plt.figure()
    plt.scatter(Y.real, Y.imag, label='Sinal Recebido')
    plt.scatter(seq16.real, seq16.imag, color='red', marker='+', label='Constelação Original')
    plt.title(f'Sinal com ruído de variância {variance}')
    plt.legend()
    plt.grid()
    plt.show()

    # Demodulação
    Z = np.zeros_like(Y, dtype=complex)
    for k in range(len(Y)):
        real_part = 3 if Y[k].real > 2 else 1 if Y[k].real > 0 else -3 if Y[k].real < -2 else -1
        imag_part = 3j if Y[k].imag > 2 else 1j if Y[k].imag > 0 else -3j if Y[k].imag < -2 else -1j
        Z[k] = real_part + imag_part

    # Contagem de erro
    error = np.sum(Z[1:K] != X[1:K])
    print(f'Para variância de {variance}, houve {error} símbolos errados.')
