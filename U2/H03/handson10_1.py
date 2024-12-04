# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 07:58:43 2024

@author: pedro
"""

import numpy as np
import random

# Passo 1: Geração de fases aleatórias
phi_k = 2 * np.pi * random.random()
phi_j = 2 * np.pi * random.random()

# Passo 2: Geração de sinais amostrados
M = 50
m = np.arange(M)
x_k = np.sin(4 * np.pi * m / 5 + phi_k)

n = 1
x_j_1 = np.sin(4 * np.pi * m / 5 + 2 * np.pi * m * n / M + phi_j)
n = 2
x_j_2 = np.sin(4 * np.pi * m / 5 + 2 * np.pi * m * n / M + phi_j)
n = 3
x_j_3 = np.sin(4 * np.pi * m / 5 + 2 * np.pi * m * n / M + phi_j)

# Passo 3: Verificação de ortogonalidade
Sum1 = np.sum(x_k * x_j_1)
print(f"O resultado para n=1 é: {Sum1}")

Sum2 = np.sum(x_k * x_j_2)
print(f"O resultado para n=2 é: {Sum2}")

Sum3 = np.sum(x_k * x_j_3)
print(f"O resultado para n=3 é: {Sum3}")