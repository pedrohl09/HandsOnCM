# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 17:03:24 2024

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

A = 1 #Amplitude do pulso
Tb = 2 #Duração do pulso que representa o bit
Eb = (A**2)*Tb #Energia por bit
Eb2N_db = np.arange(0, 15, 1) #Energia do sinal pelo ruído
Eb2N = np.zeros(len(Eb2N_db))
Eb2N = np.power(10, Eb2N_db/10) #Convete Eb2N para linear
No = Eb / Eb2N #Potência do ruído

# Gerar bits aleatórios
dataIn = np.random.rand(n_bits)
dataIn = np.sign(dataIn - 0.5)  # Sequência de -1 e 1

# Conversor serial paralelo
dataInMatrix = dataIn.reshape(n_bits // 1, 1)

# # Gerar constelação 16-QAM
# seq16qam = (2 * dataInMatrix[:, 0] + dataInMatrix[:, 1] +
#             1j * (2 * dataInMatrix[:, 2] + dataInMatrix[:, 3]))
# seq16 = seq16qam.T

# Gerar constelação BPSK
seqbpsk = dataInMatrix[:, 0]
seqb = seqbpsk.T

# Garantir propriedade da simetria
X = np.concatenate((seqb, np.conj(seqb[::-1])))

# Construindo xn
xn = np.zeros(N)  # Inicializa xn com zeros
xn = (N / np.sqrt(K)) * np.fft.ifft(X, N)  # Calcula xn com a IFFT

# Loop de variâncias
for i in range(len(Eb2N)):
    # Adição de ruído
    noise = np.sqrt(No[i]) * (np.random.randn(len(xn)) + 1j * np.random.randn(len(xn)))
    # Sinal recebido = xn + ruído
    rn = xn + noise

    # DFT de rn
    Y = np.zeros(K, dtype=complex)
    Y = np.sqrt(K)/(N*np.fft.fft(rn, n=N))

   
    # Plot do scatterplot
    plt.figure()
    plt.scatter(np.real(Y), np.imag(Y), label='Y (scatterplot)')  # Scatterplot de Y
    
    # Adiciona os pontos reais e imaginários de seqb
    plt.scatter(np.real(seqb), np.imag(seqb), color='r', marker='+', label='seqb')
    
    # Define os limites do gráfico
    plt.axis([-2, 2, -2, 2])
    
    # Adiciona a legenda
    plt.legend()
    
    # Mostra o gráfico
    plt.show()

    
    # Demodulação
    Z = np.zeros(len(Y))  # Inicializa o vetor Z com zeros
    for k in range(len(Y)):
        if np.real(Y[k]) > 0:  # Verifica a parte real de Y
            Z[k] = 1
        else:
            Z[k] = -1
    
    # Variância
    variancia = np.var(noise)  # Calcula a variância do ruído

    print(f"Para BPSK Eb2N de {Eb2N_db[i]} dB, a variância é de {variancia}\n")
