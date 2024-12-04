# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 18:36:11 2024

@author: pedro
"""

import numpy as np
import cupy as cp
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
Eb2N_db = cp.arange(0, 15, 1)  # Energia do sinal pelo ruído
Eb2N = cp.power(10, Eb2N_db / 10)  # Converte Eb2N para linear
No = Eb / Eb2N  # Potência do ruído
error_16qam = cp.zeros(len(Eb2N_db), dtype=int)
error_bpsk = cp.zeros(len(Eb2N_db), dtype=int)

# Gerar bits aleatórios
dataIn = cp.random.rand(n_bits)
dataIn = cp.sign(dataIn - 0.5)  # Sequência de -1 e 1

print('BPSK\n')
# Conversor serial paralelo
dataInMatrix_bpsk = dataIn.reshape(-1, 1)

# Gerar constelação BPSK
seqbpsk = dataInMatrix_bpsk[:, 0]
seqb = seqbpsk.T

# Garantir propriedade da simetria
X_bpsk = cp.concatenate((seqb, cp.conj(seqb[::-1])))

# Construindo xn
xn_bpsk = (N / cp.sqrt(K)) * cp.fft.ifft(X_bpsk, N)

# Loop de variâncias
for i in range(len(Eb2N)):
    # Adição de ruído
    noise_bpsk = cp.sqrt(No[i]) * (cp.random.randn(len(xn_bpsk)) + 1j * cp.random.randn(len(xn_bpsk)))
    # Sinal recebido = xn + ruído
    rn_bpsk = xn_bpsk + noise_bpsk

    # DFT de rn
    Y_bpsk = cp.fft.fft(rn_bpsk, n=N) / cp.sqrt(N)

    # Converte dados CuPy para NumPy para uso com Matplotlib
    Y_bpsk_np = Y_bpsk.get()
    seqb_np = seqb.get()

    # Plot do scatterplot
    plt.figure()
    plt.scatter(np.real(Y_bpsk_np), np.imag(Y_bpsk_np), label='Y_bpsk (scatterplot)')
    plt.scatter(np.real(seqb_np), np.imag(seqb_np), color='r', marker='+', label='seqb')
    plt.axis([-2, 2, -2, 2])
    plt.title(f"Sinal com Eb/N0 = {Eb2N_db[i]} dB")
    plt.legend()
    plt.show()

    # Demodulação
    Z_bpsk = cp.where(cp.real(Y_bpsk) > 0, 1, -1)

    # Contagem de erro
    error_bpsk[i] = cp.sum(Z_bpsk[1:K] != X_bpsk[1:K])
    print(f'Contagem de erro: {error_bpsk[i]}')

    # Variância
    variancia_bpsk = cp.var(noise_bpsk)
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
X_16qam = cp.concatenate([seq16, cp.conj(seq16[::-1])])

# Construindo xn
xn_16qam = cp.zeros(N, dtype=complex)
for n in range(N):
    for k in range(N):
        xn_16qam[n] += (1 / cp.sqrt(N)) * X_16qam[k] * cp.exp(1j * 2 * cp.pi * n * k / N)

# Loop de Eb/No
for i in range(len(Eb2N)):
    # Adição de ruído
    noise_16qam = cp.sqrt(No[i]) * (cp.random.randn(len(xn_bpsk)) + 1j * cp.random.randn(len(xn_bpsk)))
    # Sinal recebido = xn + ruído
    rn_16qam = xn_16qam + noise_16qam

    # DFT de rn
    Y_16qam = cp.zeros(K, dtype=complex)
    for k in range(K):
        for n in range(N):
            Y_16qam[k] += (1 / cp.sqrt(N)) * rn_16qam[n] * cp.exp(-1j * 2 * cp.pi * k * n / N)

    # Converte dados CuPy para NumPy para uso com Matplotlib
    Y_16qam_np = Y_bpsk.get()
    seq16_np = seqb.get()    

    # Plot do scatterplot
    plt.figure()
    plt.scatter(np.real(Y_16qam_np), np.imag(Y_16qam_np), label='Y_bpsk (scatterplot)')
    plt.scatter(np.real(seq16_np), np.imag(seq16_np), color='r', marker='+', label='seqb')
    plt.axis([-4, 4, -4, 4])
    plt.title(f"Sinal com Eb/N0 = {Eb2N_db[i]} dB")
    plt.legend()
    plt.show()

     # Demodulação
    Z_16qam = cp.zeros_like(Y_16qam, dtype=complex)
    for k in range(len(Y_16qam)):
        real_part = 3 if Y_16qam[k].real > 2 else 1 if Y_16qam[k].real > 0 else -3 if Y_16qam[k].real < -2 else -1
        imag_part = 3j if Y_16qam[k].imag > 2 else 1j if Y_16qam[k].imag > 0 else -3j if Y_16qam[k].imag < -2 else -1j
        Z_16qam[k] = real_part + imag_part

    # Contagem de erro
    error_16qam[i] = cp.sum(Z_16qam[1:K] != X_16qam[1:K])
    print(f'Contagem de erro: {error_16qam[i]}')

    # Variância
    variancia_16qam = cp.var(noise_16qam)
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

#Pe BPSK
def Q_function(x):
    return 0.5 * erfc(x / cp.sqrt(2))

Pe_bpsk = Q_function(cp.sqrt(2*Eb2N))
# Criando o gráfico de Pe vs. Eb/N0 para BPSK
plt.figure()
plt.plot(Eb2N_db, Pe_bpsk, marker='o')
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('Pe')
plt.title('Pe vs. Eb/N0 (dB) - BPSK')
plt.grid(True)
plt.show()

#Pe 16-QAM
M = 16
Pe_16qam = ((4*(cp.sqrt(M)-1))/cp.sqrt(M)) * Q_function(cp.sqrt(((3 * cp.log2(M)) / (M - 1)) * Eb2N))
# Criando o gráfico de Pe vs. Eb/N0 para BPSK
plt.figure()
plt.plot(Eb2N_db, Pe_16qam, marker='o')
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('Pe')
plt.title('Pe vs. Eb/N0 (dB) - 16-QAM')
plt.grid(True)
plt.show()