# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 14:36:34 2024

@author: Pedro Lira
"""


import numpy as np
import matplotlib.pyplot as plt
from fitter import Fitter

# Gerando dados de exemplo (substitua com seus dados reais)
data = np.loadtxt(r'C:/Users/pedro/Documents/UFRN/2024.2/Comunicações Móveis/HandsOn_02/DesvanecimentoPequenaEscala_W10.txt')

# Configurar e ajustar distribuições aos dados
f = Fitter(
    data,
    timeout=60  # Tempo limite em segundos para cada ajuste de distribuição
)
f.fit()

# Exibir o sumário com as melhores distribuições
print(f.summary())

# Obter as três melhores distribuições com base no sumário
best_distributions = f.summary().index[:3]

# Criar subplots para as três melhores distribuições
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8, 6), sharex=True)

# Plotar o histograma dos dados e a PDF de cada uma das três melhores distribuições
for ax, dist_name in zip(axes, best_distributions):
    params = f.fitted_param[dist_name]
    dist = getattr(__import__('scipy.stats', fromlist=[dist_name]), dist_name)

    # Plot do histograma
    ax.hist(data, bins=30, density=True, alpha=0.5, color='lightgray', edgecolor='black', label="Dados")

    # Plot da PDF ajustada
    x = np.linspace(min(data), max(data), 100)
    ax.plot(x, dist.pdf(x, *params), label=f"{dist_name}", lw=2)

    # Legenda e título
    ax.legend(loc="upper right")
    ax.set_ylabel("Densidade")
    ax.set_title(f"Ajuste da Distribuição: {dist_name}")

# Configurações finais
plt.xlabel("Valores")
plt.tight_layout()
plt.show()