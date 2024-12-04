# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 18:30:46 2024

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

A = 1  # Amplitude do pulso
Tb = 2  # Duração do pulso que representa o bit
Eb = (A**2) * Tb  # Energia por bit
Eb2N_db = np.arange(0, 15, 1)  # Energia do sinal pelo ruído
Eb2N = np.power(10, Eb2N_db / 10)  # Converte Eb2N para linear
No = Eb / Eb2N  # Potência do ruído

# Gerar bits aleatórios
dataIn = np.random.rand(n_bits)
dataIn = np.sign(dataIn - 0.5)  # Sequência de -1 e 1

# Conversor serial paralelo
dataInMatrix = dataIn.reshape(-1, 1)

# Gerar constelação BPSK
seqbpsk = dataInMatrix[:, 0]
seqb = seqbpsk.T

# Garantir propriedade da simetria
X = np.concatenate((seqb, np.conj(seqb[::-1])))

# Construindo xn
xn = (N / np.sqrt(K)) * np.fft.ifft(X, N)

# Loop de variâncias
for i in range(len(Eb2N)):
    # Adição de ruído
    noise = np.sqrt(No[i]) * (np.random.randn(len(xn)) + 1j * np.random.randn(len(xn)))
    # Sinal recebido = xn + ruído
    rn = xn + noise

    # DFT de rn
    Y = np.fft.fft(rn, n=N) / np.sqrt(N)

    # Plot do scatterplot
    plt.figure()
    plt.scatter(np.real(Y), np.imag(Y), label='Y (scatterplot)')
    plt.scatter(np.real(seqb), np.imag(seqb), color='r', marker='+', label='seqb')
    plt.axis([-2, 2, -2, 2])
    plt.title(f"Sinal com Eb/N0 = {Eb2N_db[i]} dB")
    plt.legend()
    plt.show()

    # Demodulação
    Z = np.where(np.real(Y) > 0, 1, -1)

    # Contagem de erros
    # errors = np.sum(seqb != Z)
    # print(f"Para BPSK Eb/N0 de {Eb2N_db[i]} dB, o número de erros é {errors}")

    # Variância
    variancia = np.var(noise)
    print(f"Para BPSK Eb/N0 de {Eb2N_db[i]} dB, a variância do ruído é {variancia}\n")
