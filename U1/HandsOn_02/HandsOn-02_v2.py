# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 10:22:58 2024

@author: pedro
"""

import numpy as np
import scipy.io as sio
from scipy.stats import nakagami, gamma
import matplotlib.pyplot as plt
from scipy.special import gammainc

def gera_canal(sPar):
    # Parâmetros de entrada
    nPoints = sPar['nPoints']
    totalLength = sPar['totalLength']
    P0 = sPar['P0']
    d0 = sPar['d0']
    n = sPar['n']
    sigma = sPar['sigma']
    shadowingWindow = sPar['shadowingWindow']
    m = sPar['m']
    dMed = sPar['dMed']
    txPower = sPar['txPower']
    
    # Distâncias do transmissor
    d = np.arange(d0, totalLength + dMed, dMed)
    nSamples = len(d)

    # Perda de percurso (determinística)
    vtPathLoss = P0 + 10 * n * np.log10(d / d0)

    # Sombreamento
    nShadowSamples = nSamples // shadowingWindow
    shadowing = sigma * np.random.randn(nShadowSamples)
    restShadowing = sigma * np.random.randn(1)[0] * np.ones(nSamples % shadowingWindow)
    shadowing = np.repeat(shadowing, shadowingWindow)
    shadowing = np.concatenate((shadowing, restShadowing))

    # Filtragem para evitar variação abrupta do sombreamento
    jan = shadowingWindow // 2
    vtShadCorr = np.array([np.mean(shadowing[i - jan:i + jan]) for i in range(jan, nSamples - jan)])
    
    # Ajuste do desvio padrão depois do filtro de correlação do sombreamento
    vtShadCorr = vtShadCorr * (np.std(shadowing) / np.std(vtShadCorr))
    vtShadCorr = vtShadCorr - np.mean(vtShadCorr) + np.mean(shadowing)

    # Desvanecimento de pequena escala (Nakagami)
    nakagamiSamp = 20 * np.log10(nakagami.rvs(m, size=nSamples))
    
    # Ajuste do número de amostras devido à filtragem
    txPower = np.full(nSamples, txPower)[jan:nSamples - jan]
    vtPathLoss = vtPathLoss[jan:nSamples - jan]
    vtFading = nakagamiSamp[jan:nSamples - jan]
    vtDist = d[jan:nSamples - jan]

    # Potência recebida
    vtPrxdBm = txPower - vtPathLoss + vtShadCorr + vtFading

    # Salva as variáveis do canal em um arquivo .mat
    sio.savemat(sPar['chFileName'] + '.mat', {
        'vtDist': vtDist,
        'vtPathLoss': vtPathLoss,
        'vtShadCorr': vtShadCorr,
        'vtFading': vtFading,
        'vtPrxdBm': vtPrxdBm
    })

    # Salva as variáveis do canal em um arquivo .txt
    data = np.column_stack((vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm))
    np.savetxt(sPar['chFileName'] + '.txt', data, header="vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm", delimiter=",", comments='')

    return vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm




# Parâmetros para geração do canal sintético
sPar = {
    'd0': 5,
    'P0': 0,
    'nPoints': 50000,
    'totalLength': 100,
    'n': 4,
    'sigma': 6,
    'shadowingWindow': 200,
    'm': 4,
    'txPower': 0,
    'nCDF': 40,
    'dW': 100,
    'chFileName': 'Prx_sintetico'
}

# Distância entre pontos de medição
sPar['dMed'] = sPar['totalLength'] / sPar['nPoints']

# Gera o canal sintético (usando a função `gera_canal` criada anteriormente)
vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm = gera_canal(sPar)

# Transformação de potência em mW
vtPtrxmW = 10 ** (vtPrxdBm / 10)
nSamples = len(vtPtrxmW)

# Vetores para o canal estimado
vtDesLarga = []
vtDesPequeEst = []

# Cálculo do desvanecimento lento e rápido
dMeiaJanela = round((sPar['dW'] - 1) / 2)
for ik in range(dMeiaJanela, nSamples - dMeiaJanela):
    # Desvanecimento de larga escala (dB)
    vtDesLarga.append(10 * np.log10(np.mean(vtPtrxmW[ik - dMeiaJanela:ik + dMeiaJanela])))
    # Desvanecimento de pequena escala (dB)
    vtDesPequeEst.append(vtPrxdBm[ik] - vtDesLarga[-1])

# Cálculo da envoltória normalizada para o cálculo do fading
indexes = range(dMeiaJanela, nSamples - dMeiaJanela)
vtPtrxmWNew = 10 ** (vtPrxdBm[indexes] / 10)
desLarga_Lin = 10 ** (np.array(vtDesLarga) / 10)
envNormal = np.sqrt(vtPtrxmWNew) / np.sqrt(desLarga_Lin)

# Ajuste no tamanho dos vetores devido à filtragem
vtDistEst = vtDist[dMeiaJanela:nSamples - dMeiaJanela]
vtPrxdBmEst = vtPrxdBm[dMeiaJanela:nSamples - dMeiaJanela]

# Cálculo reta de perda de percurso
vtDistLog = np.log10(vtDist)
vtDistLogEst = np.log10(vtDistEst)
dCoefReta = np.polyfit(vtDistLogEst, vtPrxdBmEst, 1)
dNEst = -dCoefReta[0] / 10
print(f"Estimação dos parâmetros de larga escala (W = {sPar['dW']}):")
print(f"   Expoente de perda de percurso estimado n = {dNEst}")

# Perda de percurso estimada
vtPathLossEst = np.polyval(dCoefReta, vtDistLogEst)
vtShadCorrEst = np.array(vtDesLarga) - vtPathLossEst
stdShad = np.std(vtShadCorrEst)
meanShad = np.mean(vtShadCorrEst)
print(f"   Desvio padrão do sombreamento estimado = {stdShad}")
print(f"   Média do sombreamento estimado = {meanShad}")

vtPathLossEst = -vtPathLossEst
vtPrxEst = sPar['txPower'] - vtPathLossEst + vtShadCorrEst + vtDesPequeEst

# Estimação da CDF do desvanecimento de pequena escala
vtn = np.arange(1, sPar['nCDF'] + 1)
xCDF = 1.2 ** (vtn - 1) * 0.01
cdffn = np.array([np.sum(envNormal <= x) for x in xCDF]) / len(envNormal)

# Monta estrutura do histograma
xccdfEst = 20 * np.log10(xCDF)
yccdfEst = cdffn

# Gráficos do canal estimado
plt.figure()
plt.plot(vtDistLogEst, vtPrxEst, label='Prx canal completo')
plt.plot(vtDistLogEst, sPar['txPower'] - vtPathLossEst, linewidth=2, label='Prx (somente perda de percurso)')
plt.plot(vtDistLogEst, sPar['txPower'] - vtPathLossEst + vtShadCorrEst, linewidth=2, label='Prx (perda de percurso + sombreamento)')
plt.xlabel('log10(d)')
plt.ylabel('Potência [dBm]')
plt.legend()
plt.title('Prx original vs estimada')
plt.show()

# Gráfico da Perda de percurso
plt.figure()
plt.plot(vtDistLog, -vtPathLoss, label='Path Loss original')
plt.plot(vtDistLogEst, -vtPathLossEst, label='Path Loss estimado')
plt.legend()
plt.title('Perda de percurso original vs estimada')
plt.show()

# Gráfico do Sombreamento
plt.figure()
plt.plot(vtDistLog, vtShadCorr, label='Shadowing original')
plt.plot(vtDistLogEst, vtShadCorrEst, label='Shadowing estimado')
plt.legend()
plt.title('Sombreamento original vs estimado')
plt.show()

# Gráfico do Fading
plt.figure()
plt.plot(vtDistLog, vtFading, label='Fading original')
plt.plot(vtDistLogEst, vtDesPequeEst, label='Fading estimado')
plt.legend()
plt.title('Fading original vs estimado')
plt.show()

# Plot das CDFs normalizadas Nakagami
plt.figure()
plt.plot(xccdfEst, yccdfEst, '--', label='CDF das amostras')
vtm = [1, 2, 4, 6]
for m_value in vtm:
    cdf_nakagami = gammainc(m_value, m_value * xCDF**2)
    plt.semilogy(20 * np.log10(xCDF), cdf_nakagami, label=f'm = {m_value}')
plt.legend()
plt.title('Estudo do fading com o conhecimento da distribuição')
plt.xlabel('x')
plt.ylabel('F(x)')
plt.axis([-30, 10, 1e-5, 1])
plt.show()
