# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 12:03:18 2024

@author: pedro
"""

import numpy as np
import matplotlib.pyplot as plt

# Função para calcular o shadowing correlacionado
def fCorrShadowing(mtPontosMedicao, dShad, dAlphaCorr, dSigmaShad, dDimXOri, dDimYOri):
    dDimY = int(np.ceil(dDimYOri + dDimYOri % dShad))
    dDimX = int(np.ceil(dDimXOri + dDimXOri % dShad))
    mtPosxShad, mtPosyShad = np.meshgrid(np.arange(0, dDimX + 1, dShad), np.arange(0, dDimY + 1, dShad))
    mtShadowingSamples = dSigmaShad * np.random.randn(*mtPosyShad.shape)
    shadowing_corr = np.zeros((mtPontosMedicao.shape[0], mtPontosMedicao.shape[1], len(vtBs)))

    for iBs in range(len(vtBs)):
        for il in range(mtPontosMedicao.shape[0]):
            for ic in range(mtPontosMedicao.shape[1]):
                dshadPoint = mtPontosMedicao[il, ic]
                dXIndexP1 = int(np.real(dshadPoint) / dShad)
                dYIndexP1 = int(np.imag(dshadPoint) / dShad)
                shadowing_corr[il, ic, iBs] = mtShadowingSamples[dYIndexP1, dXIndexP1] if dXIndexP1 < mtShadowingSamples.shape[1] and dYIndexP1 < mtShadowingSamples.shape[0] else 0
    return shadowing_corr

# Função para desenhar as ERBs
def fDrawDeploy(dR, vtBs):
    for pos in vtBs:
        plt.plot(np.real(pos), np.imag(pos), 'ro')

for dAlphaCorr in np.arange(0, 1, 0.05):
    # Parâmetros de entrada
    dR = 200
    dFc = 800
    dShad = 50
    dSigmaShad = 8
    #dAlphaCorr = 0.5
    dPasso = 10
    dRMin = dPasso
    dIntersiteDistance = 2 * np.sqrt(3 / 4) * dR

    dDimXOri = 5 * dR
    dDimYOri = 6 * np.sqrt(3 / 4) * dR
    dPtdBm = 57
    dPtLinear = 10 ** (dPtdBm / 10) * 1e-3
    dHMob = 5
    dHBs = 30
    dAhm = 3.2 * (np.log10(11.75 * dHMob)) ** 2 - 4.97

    # Posições das BSs
    vtBs = [0]
    dOffset = np.pi / 6
    for iBs in range(2, 8):
        vtBs.append(dR * np.sqrt(3) * np.exp(1j * ((iBs - 2) * np.pi / 3 + dOffset)))
    vtBs = np.array(vtBs) + (dDimXOri / 2 + 1j * dDimYOri / 2)

    dDimY = int(np.ceil(dDimYOri + dDimYOri % dPasso))
    dDimX = int(np.ceil(dDimXOri + dDimXOri % dPasso))
    mtPosx, mtPosy = np.meshgrid(np.arange(0, dDimX + 1, dPasso), np.arange(0, dDimY + 1, dPasso))
    mtPontosMedicao = mtPosx + 1j * mtPosy

    mtPowerFinaldBm = -np.inf * np.ones(mtPosy.shape)
    mtPowerFinalShaddBm = -np.inf * np.ones(mtPosy.shape)
    mtPowerFinalShadCorrdBm = -np.inf * np.ones(mtPosy.shape)

    # Calcula o sombreamento correlacionado
    mtShadowingCorr = fCorrShadowing(mtPontosMedicao, dShad, dAlphaCorr, dSigmaShad, dDimXOri, dDimYOri)

    for iBsD in range(len(vtBs)):
        mtPosEachBS = mtPontosMedicao - vtBs[iBsD]
        mtDistEachBs = np.abs(mtPosEachBS)
        mtDistEachBs[mtDistEachBs < dRMin] = dRMin
        mtPldB = 69.55 + 26.16 * np.log10(dFc) + (44.9 - 6.55 * np.log10(dHBs)) * np.log10(mtDistEachBs / 1e3) - 13.82 * np.log10(dHBs) - dAhm
        mtShadowing = dSigmaShad * np.random.randn(*mtPosy.shape)
        mtPowerEachBSdBm = dPtdBm - mtPldB
        mtPowerEachBSShaddBm = dPtdBm - mtPldB + mtShadowing
        mtPowerEachBSShadCorrdBm = dPtdBm - mtPldB + mtShadowingCorr[:, :, iBsD]

        mtPowerFinaldBm = np.maximum(mtPowerFinaldBm, mtPowerEachBSdBm)
        mtPowerFinalShaddBm = np.maximum(mtPowerFinalShaddBm, mtPowerEachBSShaddBm)
        mtPowerFinalShadCorrdBm = np.maximum(mtPowerFinalShadCorrdBm, mtPowerEachBSShadCorrdBm)

    print(f"Valor de dAlphaCorr = {round(dAlphaCorr, 3)}")
    print(f"Desvio padrão = {round(np.std(mtShadowingCorr), 2)}")
    print("-------------------------")