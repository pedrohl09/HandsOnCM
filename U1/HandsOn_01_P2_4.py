# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 20:49:50 2024

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

# Vetor de frequências de portadora
vtFc = [800, 1800, 2100]

for dFc in vtFc:
    # Entrada de parâmetros
    dR = 500  # Raio do hexágono
    dPasso = int(np.ceil(dR / 50))  # Resolução do grid: distância entre pontos de medição
    dRMin = dPasso  # Raio de segurança
    dIntersiteDistance = 2 * np.sqrt(3 / 4) * dR  # Distância entre ERBs (somente para informação)
    dDimX = 5 * dR  # Dimensão X do grid
    dDimY = 6 * np.sqrt(3 / 4) * dR  # Dimensão Y do grid
    dPtdBm = 21  # EIRP (incluindo ganho e perdas)
    dPtLinear = 10 ** (dPtdBm / 10) * 1e-3  # EIRP em escala linear
    dSensitivity = -90  # Sensibilidade do receptor
    dHMob = 1.5  # Altura do receptor
    dHBs = 32  # Altura do transmissor
    dAhm = 3.2 * (np.log10(11.75 * dHMob)) ** 2 - 4.97  # Modelo Okumura-Hata: Cidade grande e fc >= 400MHz

    # Vetor com posições das BSs (grid Hexagonal com 7 células, uma célula central e uma camada de células ao redor)
    vtBs = np.array([0], dtype=complex)
    dOffset = np.pi / 6

    for iBs in range(2, 8):
        vtBs = np.append(vtBs, dR * np.sqrt(3) * np.exp(1j * ((iBs - 2) * np.pi / 3 + dOffset)))   
    
    # Vetor com posições das Microcélulas (grid Hexagonal com 6 células, uma célula central e uma camada de células ao redor)    
    vtBsMicro = np.array([0], dtype=complex)    
     
    for iBs in range(1, 7):
        vtBsMicro = np.append(vtBsMicro, dR*np.exp(1j * ((iBs - 2) * np.pi/3 + np.pi)))

    # Ajuste de posição das bases (posição relativa ao canto inferior esquerdo)
    vtBs += (dDimX / 2 + 1j * dDimY / 2)
    
    # Ajuste de posição das bases das microcéululas (posição relativa ao canto inferior esquerdo)
    vtBsMicro += (dDimX / 2 + 1j * dDimY / 2)

    # Ajuste de dimensão para medir toda a dimensão do grid
    dDimY = np.ceil(dDimY + dDimY % dPasso)
    dDimX = np.ceil(dDimX + dDimX % dPasso)

    # Matriz de referência com posição de cada ponto do grid
    mtPosx, mtPosy = np.meshgrid(np.arange(0, dDimX + 1, dPasso), np.arange(0, dDimY + 1, dPasso))

    # Iniciação da matriz de potência recebida máxima em cada ponto medido
    mtPowerFinaldBm = -np.inf * np.ones(mtPosy.shape)

    # Calcular o REM de cada ERB e acumular a maior potência em cada ponto de medição
    for iBsD in range(len(vtBs)):  # Loop nas 7 ERBs
        # Matriz 3D com os pontos de medição de cada ERB
        mtPosEachBS = (mtPosx + 1j * mtPosy) - vtBs[iBsD]
        mtDistEachBs = np.abs(mtPosEachBS)  # Distância entre cada ponto de medição e a sua ERB
        mtDistEachBs[mtDistEachBs < dRMin] = dRMin  # Implementação do raio de segurança

        # Okumura-Hata (cidade urbana) - dB
        mtPldB = 69.55 + 26.16 * np.log10(dFc) + (44.9 - 6.55 * np.log10(dHBs)) * np.log10(mtDistEachBs / 1e3) - 13.82 * np.log10(dHBs) - dAhm
        mtPowerEachBSdBm = dPtdBm - mtPldB  # Potências recebidas em cada ponto de medição

        # Cálculo da maior potência em cada ponto de medição
        mtPowerFinaldBm = np.maximum(mtPowerFinaldBm, mtPowerEachBSdBm)
        
    # Calcular o REM de cada Microcélula e acumular a maior potência em cada ponto de medição
    for iBsD in range(len(vtBsMicro)):  # Loop nas 7 ERBs
        # Matriz 3D com os pontos de medição de cada ERB
        mtPosEachBSmicro = (mtPosx + 1j * mtPosy) - vtBsMicro[iBsD]
        mtDistEachBsmicro = np.abs(mtPosEachBSmicro)  # Distância entre cada ponto de medição e a sua ERB
        mtDistEachBsmicro[mtDistEachBsmicro < dRMin] = dRMin  # Implementação do raio de segurança

        # PL microcélula - dB
        mtPldBmicro = 55 + 38*np.log10(mtDistEachBsmicro/1000) + (24.5 + (1.5*(dFc)/925)) * np.log10(dFc)
        mtPowerEachBSdBmMicro = 20 - mtPldBmicro  # Potências recebidas em cada ponto de medição

        # Cálculo da maior potência em cada ponto de medição
        mtPowerFinaldBm = np.maximum(mtPowerFinaldBm, mtPowerEachBSdBmMicro)    
        
    mtOutRate = np.ones_like(mtPowerFinaldBm) 
    dSizel, dSizec = mtPowerFinaldBm.shape
    
    for il in range(dSizel):
        for ic in range(dSizec):
            if mtPowerFinaldBm[il, ic] > dSensitivity:
                mtOutRate[il, ic] = 0
            else:
                mtOutRate[il, ic] = 1

    # Plot da REM de todo o grid (composição das 7 ERBs)
    plt.figure()
    plt.pcolor(mtPosx, mtPosy, mtOutRate, shading='auto')
    plt.colorbar()
    fDrawDeploy(dR, vtBs)
    plt.axis('equal')
    plt.title(f'Microcélulas para {dFc} MHz')
    plt.show()