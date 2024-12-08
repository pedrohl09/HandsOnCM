# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 00:40:30 2024

@author: pedro
"""

import numpy as np
import matplotlib.pyplot as plt

def fDrawSector(dR, dCenter):
    vtHex = np.zeros(0, dtype=complex)
    
    # Cria os pontos do hexágono
    for ie in range(1, 7):
        vtHex = np.append(vtHex, dR * (np.cos((ie-1) * np.pi / 3) + 1j * np.sin((ie-1) * np.pi / 3)))
    
    # Desloca o hexágono para o centro desejado
    vtHex += dCenter
    
    # Adiciona o primeiro ponto no final para fechar o hexágono
    vtHexp = np.append(vtHex, vtHex[0])
    
    # Plota o hexágono
    plt.plot(vtHexp.real, vtHexp.imag, 'k')

def fDrawDeploy(dR, vtBs):
    plt.gca().set_aspect('equal', adjustable='box')
    
    # Desenha setores hexagonais
    for iBsD in range(len(vtBs)):
        fDrawSector(dR, vtBs[iBsD])
    
    # Plota as posições das BSs (estações base)
    plt.plot(vtBs.real, vtBs.imag, 'sk')
    plt.axis('equal')

# Parâmetros de entrada
dR = 5e3  # Raio do hexágono
dFc = 800  # Frequência da portadora em MHz
dPasso = int(np.ceil(dR / 10))  # Resolução do grid: distância entre pontos de medição
dRMin = dPasso  # Raio de segurança
dIntersiteDistance = 2 * np.sqrt(3 / 4) * dR  # Distância entre ERBs (somente para informação)
dDimX = 5 * dR  # Dimensão X do grid
dDimY = 6 * np.sqrt(3 / 4) * dR  # Dimensão Y do grid
dPtdBm = 57  # EIRP (incluindo ganho e perdas)
dPtLinear = 10**(dPtdBm / 10) * 1e-3  # EIRP em escala linear
dHMob = 5  # Altura do receptor
dHBs = 30  # Altura do transmissor
dAhm = 3.2 * (np.log10(11.75 * dHMob))**2 - 4.97  # Modelo Okumura-Hata: Cidade grande e fc >= 400MHz

# Vetor com posições das BSs (grid Hexagonal com 7 células, uma célula central e uma camada de células ao redor)
vtBs = np.array([0], dtype=complex)
dOffset = np.pi / 6

# Cria as posições das BSs ao redor da célula central
for iBs in range(2, 8):
    vtBs = np.append(vtBs, dR * np.sqrt(3) * np.exp(1j * ((iBs-2) * np.pi / 3 + dOffset)))

# Ajusta a posição das bases (posição relativa ao canto inferior esquerdo)
vtBs += (dDimX / 2 + 1j * dDimY / 2)

# Ajuste de dimensão para medir toda a dimensão do grid
dDimY += dDimY % dPasso
dDimX += dDimX % dPasso

# Matriz de referência com a posição de cada ponto do grid
mtPosx, mtPosy = np.meshgrid(np.arange(0, dDimX + 1, dPasso), np.arange(0, dDimY + 1, dPasso))

# Iniciação da matriz de potência recebida máxima em cada ponto medido
mtPowerFinaldBm = -np.inf * np.ones(mtPosy.shape)

# Calcular O REM de cada ERB e acumular a maior potência em cada ponto de medição
mtPosEachBS = np.zeros((mtPosx.shape[0], mtPosx.shape[1], len(vtBs)), dtype=complex)
mtPowerEachBSdBm = np.zeros_like(mtPosEachBS.real)

for iBsD in range(len(vtBs)):  # Loop nas 7 ERBs
    # Matriz 3D com os pontos de medição de cada ERB
    mtPosEachBS[:, :, iBsD] = (mtPosx + 1j * mtPosy) - vtBs[iBsD]
    mtDistEachBs = np.abs(mtPosEachBS[:, :, iBsD])  # Distância entre cada ponto de medição e a sua ERB
    mtDistEachBs[mtDistEachBs < dRMin] = dRMin  # Implementação do raio de segurança
    
    # Okumura-Hata (cidade urbana) - dB
    mtPldB = 69.55 + 26.16 * np.log10(dFc) + (44.9 - 6.55 * np.log10(dHBs)) * np.log10(mtDistEachBs / 1e3) - 13.82 * np.log10(dHBs) - dAhm
    mtPowerEachBSdBm[:, :, iBsD] = dPtdBm - mtPldB  # Potências recebidas em cada ponto de medição
    
    # Plot da REM de cada ERB individualmente
    plt.figure()
    plt.pcolor(mtPosx, mtPosy, mtPowerEachBSdBm[:, :, iBsD], shading='auto')
    plt.colorbar()
    
    # Desenha setores hexagonais
    fDrawDeploy(dR, vtBs)
    plt.axis('equal')
    plt.title(f'ERB {iBsD + 1}')
    plt.show()
    
    # Cálculo da maior potência em cada ponto de medição
    mtPowerFinaldBm = np.maximum(mtPowerFinaldBm, mtPowerEachBSdBm[:, :, iBsD])

# Plot da REM de todo o grid (composição das 7 ERBs)
plt.figure()
plt.pcolor(mtPosx, mtPosy, mtPowerFinaldBm, shading='auto')
plt.colorbar()
fDrawDeploy(dR, vtBs)
plt.axis('equal')
plt.title('Todas as 7 ERB')
plt.show()