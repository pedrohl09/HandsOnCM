# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 09:13:48 2024

@author: pedro
"""

import numpy as np
import matplotlib.pyplot as plt

# Parâmetros
n_bits = 100             # Número de bits
T = 50                   # Tempo de símbolo
Ts = 2                   # Tempo de símbolo em portadora única
K = T // Ts              # Número de subportadoras independentes
N = 2 * K                # N pontos da IDFT

# Gerar bits aleatórios
dataIn = np.random.rand(1, n_bits)  # Sequência de números entre 0 e 1 uniformemente distribuídos
dataIn = np.sign(dataIn - 0.5)      # Sequência de -1 e 1

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

# Construindo xt
xt = np.zeros(T + 1, dtype=complex)
for t in range(T + 1):
    for k in range(N):
        xt[t] += (1 / np.sqrt(N)) * X[k] * np.exp(1j * 2 * np.pi * k * t / T)

# Plots
plt.plot(np.abs(xt), label='x(t)')
plt.stem(range(len(xn)), np.abs(xn), linefmt='r', markerfmt='ro', basefmt=" ", label='x_n')
plt.title('Sinais OFDM')
plt.legend()
plt.xlabel('Tempo')
plt.show()