# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 18:36:11 2024

@author: pedro
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

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
error_16qam = np.zeros(len(Eb2N_db), dtype=int)
error_bpsk = np.zeros(len(Eb2N_db), dtype=int)

# Gerar bits aleatórios
dataIn = np.random.rand(n_bits)
dataIn = np.sign(dataIn - 0.5)  # Sequência de -1 e 1

print('BPSK\n')
# Conversor serial paralelo
dataInMatrix_bpsk = dataIn.reshape(-1, 1)

# Gerar constelação BPSK
seqbpsk = dataInMatrix_bpsk[:, 0]
seqb = seqbpsk.T

# Garantir propriedade da simetria
X_bpsk = np.concatenate((seqb, np.conj(seqb[::-1])))

# Construindo xn
xn_bpsk = (N / np.sqrt(K)) * np.fft.ifft(X_bpsk, N)

# Loop de variâncias
for i in range(len(Eb2N)):
    # Adição de ruído
    noise_bpsk = np.sqrt(No[i]) * (np.random.randn(len(xn_bpsk)) + 1j * np.random.randn(len(xn_bpsk)))
    # Sinal recebido = xn + ruído
    rn_bpsk = xn_bpsk + noise_bpsk

    # DFT de rn
    Y_bpsk = np.fft.fft(rn_bpsk, n=N) / np.sqrt(N)

    # Plot do scatterplot
    plt.figure()
    plt.scatter(np.real(Y_bpsk), np.imag(Y_bpsk), label='Y_bpsk (scatterplot)')
    plt.scatter(np.real(seqb), np.imag(seqb), color='r', marker='+', label='seqb')
    plt.axis([-2, 2, -2, 2])
    plt.title(f"Sinal com Eb/N0 = {Eb2N_db[i]} dB")
    plt.legend()
    plt.show()

    # Demodulação
    Z_bpsk = np.where(np.real(Y_bpsk) > 0, 1, -1)

    # Contagem de erro
    error_bpsk[i] = np.sum(Z_bpsk[1:K] != X_bpsk[1:K])
    print(f'Contagem de erro: {error_bpsk[i]}')

    # Variância
    variancia_bpsk = np.var(noise_bpsk)
    print(f"Para BPSK Eb/N0 de {Eb2N_db[i]} dB, a variância do ruído é {variancia_bpsk}\n")

print('\n', '-'*50, '\n')

print('16-QAM')    
# Conversor serial paralelo
dataInMatrix_16qam = dataIn.reshape(n_bits // 4, 4)

# Gerar constelação 16-QAM
seq16qam = (2 * dataInMatrix_16qam[:, 0] + dataInMatrix_16qam[:, 1] +
            1j * (2 * dataInMatrix_16qam[:, 2] + dataInMatrix_16qam[:, 3]))
seq16 = seq16qam.T

# Garantir propriedade da simetria
X_16qam = np.concatenate([seq16, np.conj(seq16[::-1])])

# Construindo xn
xn_16qam = np.zeros(N, dtype=complex)
for n in range(N):
    for k in range(N):
        xn_16qam[n] += (1 / np.sqrt(N)) * X_16qam[k] * np.exp(1j * 2 * np.pi * n * k / N)

# Loop de variâncias
for i in range(len(Eb2N)):
    # Adição de ruído
    noise_16qam = np.sqrt(No[i]) * (np.random.randn(len(xn_bpsk)) + 1j * np.random.randn(len(xn_bpsk)))
    # Sinal recebido = xn + ruído
    rn_16qam = xn_16qam + noise_16qam

    # DFT de rn
    Y_16qam = np.zeros(K, dtype=complex)
    for k in range(K):
        for n in range(N):
            Y_16qam[k] += (1 / np.sqrt(N)) * rn_16qam[n] * np.exp(-1j * 2 * np.pi * k * n / N)

    # Plot do scatterplot
    plt.figure()
    plt.scatter(np.real(Y_16qam), np.imag(Y_16qam), label='Y_16-QAM (scatterplot)')
    plt.scatter(np.real(seq16), np.imag(seq16), color='r', marker='+', label='seq16')
    plt.axis([-2, 2, -2, 2])
    plt.title(f"Sinal com Eb/N0 = {Eb2N_db[i]} dB")
    plt.legend()
    plt.show()

     # Demodulação
    Z_16qam = np.zeros_like(Y_16qam, dtype=complex)
    for k in range(len(Y_16qam)):
        real_part = 3 if Y_16qam[k].real > 2 else 1 if Y_16qam[k].real > 0 else -3 if Y_16qam[k].real < -2 else -1
        imag_part = 3j if Y_16qam[k].imag > 2 else 1j if Y_16qam[k].imag > 0 else -3j if Y_16qam[k].imag < -2 else -1j
        Z_16qam[k] = real_part + imag_part

    # Contagem de erro
    error_16qam[i] = np.sum(Z_16qam[1:K] != X_16qam[1:K])
    print(f'Contagem de erro: {error_16qam[i]}')

    # Variância
    variancia_16qam = np.var(noise_16qam)
    print(f"Para 16-QAM Eb/N0 de {Eb2N_db[i]} dB, a variância do ruído é {variancia_16qam}\n")

# Criando o gráfico de BER vs. Eb/N0 para BPSK
plt.figure()
plt.plot(Eb2N_db, error_bpsk, marker='o')
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('Taxa de erro de bit')
plt.title('BER vs. Eb/N0 (dB) - BPSK')
plt.grid(True)
plt.show()

# Criando o gráfico de BER vs. Eb/N0 para 16QAM
plt.figure()
plt.plot(Eb2N_db, error_16qam, marker='o')
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('Taxa de erro de bit')
plt.title('BER vs. Eb/N0 (dB) - 16-QAM')
plt.grid(True)
plt.show()

#Pe
def Q_function(x):
    return 0.5 * erfc(x / np.sqrt(2))

Pe_bpsk = Q_function(np.sqrt(2*Eb2N_db))
# Criando o gráfico de Pe vs. Eb/N0 para BPSK
plt.figure()
plt.plot(Eb2N_db, Pe_bpsk, marker='o')
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('Pe')
plt.title('Pe vs. Eb/N0 (dB) - BPSK')
plt.grid(True)
plt.show()
