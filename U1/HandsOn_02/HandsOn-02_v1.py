# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 21:59:50 2024

@author: pedro
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import nakagami
from scipy.special import gamma
from scipy.optimize import curve_fit
import scipy.io as sio
from sklearn.metrics import mean_squared_error
from scipy.io import savemat

def fGeraCanal(sPar):
    """
    Purpose: Generate a channel consisting of path loss, shadowing, and Nakagami fading.

    Parameters in sPar:
    - nPoints: Number of samples
    - totalLength: Maximum distance of the route
    - P0: Reference power measured at distance d0
    - d0: Reference distance d0
    - n: Path loss exponent
    - sigma: Standard deviation of the log-normal shadowing [dB]
    - shadowingWindow: Correlation window size for shadowing [samples]
    - m: Nakagami parameter
    - dMed: Distance between measurement points (totalLength/nPoints)
    - txPower: Transmission power in dBm

    Returns:
    - vtDist: Measurement points [m]
    - vtPathLoss: Path loss samples
    - vtShadCorr: Shadowing samples
    - vtFading: Small-scale fading samples
    - vtPrxdBm: Received power with complete channel
    """

    # Extract parameters
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

    # Distance from transmitter
    d = np.arange(d0, totalLength, dMed)
    nSamples = len(d)

    # Path Loss generation (deterministic)
    vtPathLoss = P0 + 10 * n * np.log10(d / d0)

    # Shadowing generation
    nShadowSamples = nSamples // shadowingWindow
    shadowing = sigma * np.random.randn(nShadowSamples)

    # Samples for the last window
    restShadowing = sigma * np.random.randn(1) * np.ones(nSamples % shadowingWindow)

    # Repeat the same shadowing value during the correlation window
    shadowing = np.tile(shadowing, (shadowingWindow, 1)).flatten()
    shadowing = np.concatenate((shadowing, restShadowing))

    # Moving average filter for smoother shadowing
    jan = shadowingWindow // 2
    vtShadCorr = np.convolve(shadowing, np.ones((2 * jan + 1,)) / (2 * jan + 1), mode='valid')
    vtShadCorr = vtShadCorr * np.std(shadowing) / np.std(vtShadCorr)
    vtShadCorr = vtShadCorr - np.mean(vtShadCorr) + np.mean(shadowing)

    # Nakagami fading generation
    nakagami_pdf = lambda x: (2. * m ** m / gamma(m)) * x ** (2 * m - 1) * np.exp(-m * x ** 2)
    nakagamiNormEnvelope = nakagami.rvs(m, size=nSamples)
    nakagamiSamp = 20 * np.log10(nakagamiNormEnvelope)

    # Received power calculation
    vtTxPower = txPower * np.ones(nSamples)
    vtTxPower = vtTxPower[jan:-jan]
    vtPathLoss = vtPathLoss[jan:-jan]
    vtFading = nakagamiSamp[jan:-jan]
    vtDist = d[jan:-jan]

    # Complete channel received power
    vtPrxdBm = vtTxPower - vtPathLoss + vtShadCorr + vtFading

    #Save channel variables to a .mat file
    savemat(f"{sPar['chFileName']}.mat", {
        'vtDist': vtDist,
        'vtPathLoss': vtPathLoss,
        'vtShadCorr': vtShadCorr,
        'vtFading': vtFading,
        'vtPrxdBm': vtPrxdBm
    })

    return vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm

def fEstimaCanal(sPar):
    """
    Estima parâmetros do canal: path loss, shadowing e desvanecimento plano.

    Entradas:
    - sPar: dicionário com parâmetros de entrada
        - P0: potência de referência medida na distância d0
        - d0: distância d0 de referência
        - dW: janela para média móvel do desvanecimento de larga escala
        - txPower: potência de transmissão em dBm
        - nCDF: Número de pontos da CDF normalizada

    Saídas:
    - sOut: dicionário com saídas estimadas
        - vtDistEst: pontos de medição [m]
        - vtPathLossEst: amostras estimadas da perda de percurso
        - vtShadCorrEst: amostras do sombreamento
        - vtDesPequeEst: amostras do desvanecimento de pequena escala
        - vtPrxEst: potência recebida estimada (canal completo)
        - dNEst: Expoente de perda de percurso estimado
        - dStdShadEst: Desvio padrão do sombreamento estimado
        - dStdMeanShadEst: Média do sombreamento estimado
        - vtXCcdfEst: Valores dos Bins da CCDF do desvanecimento de pequena escala
        - vtYCcdfEst: Quantidade de amostras em cada bin da CCDF do desvanecimento de pequena escala
    """

    # Lê canal gerado
    data = sio.loadmat(sPar['chFileName'])
    vtPrxdBm = np.array(data['vtPrxdBm'].flatten())
    vtDist = data['vtDist'].flatten()

    vtPtrxmW = 10 ** (vtPrxdBm / 10)
    nSamples = len(vtPtrxmW)

    # Vetores para canal estimado
    vtDesLarga = []
    vtDesPequeEst = []

    # Cálculo do desvanecimento lento e rápido
    dMeiaJanela = (sPar['dW'] - 1) // 2  # Meia janela
    ij = 0

    for ik in range(dMeiaJanela, nSamples - dMeiaJanela):
        # Desvanecimento de larga escala: perda de percurso + sombreamento [dB]
        vtDesLarga.append(10 * np.log10(np.mean(vtPtrxmW[ik - dMeiaJanela: ik + dMeiaJanela + 1])))
        # Desvanecimento de pequena escala [dB]
        vtDesPequeEst.append(vtPrxdBm[ik] - vtDesLarga[ij])
        ij += 1

    # Cálculo da envoltória normalizada (para efeitos de cálculo do fading)
    indexes = range(dMeiaJanela, nSamples - dMeiaJanela)
    vtPtrxmWNew = 10 ** (vtPrxdBm[list(indexes)] / 10)
    desLarga_Lin = 10 ** (np.array(vtDesLarga)[:len(indexes)] / 10)
    vtEnvNorm = np.sqrt(vtPtrxmWNew) / np.sqrt(desLarga_Lin)

    # Ajuste no tamanho dos vetores devido a filtragem
    vtDistEst = vtDist[dMeiaJanela:nSamples - dMeiaJanela]
    vtPrxdBm = vtPrxdBm[dMeiaJanela:nSamples - dMeiaJanela]

    # Cálculo reta de perda de percurso
    vtDistLog = np.log10(vtDist)
    vtDistLogEst = np.log10(vtDistEst)

    # Cálculo do coeficiente da reta que melhor se caracteriza a perda de percurso
    dCoefReta = np.polyfit(vtDistLogEst, vtPrxdBm, 1)
    # Expoente de perda de percurso estimado
    dNEst = -dCoefReta[0] / 10

    # Perda de percurso estimada para os pontos de medição
    vtPathLossEst = np.polyval(dCoefReta, vtDistLogEst)

    # Sombreamento
    vtShadCorrEst = np.array(vtDesLarga) - vtPathLossEst
    # Calcula a variância do sombreamento estimado
    dStdShadEst = np.std(vtShadCorrEst)
    dStdMeanShadEst = np.mean(vtShadCorrEst)

    vtPathLossEst = -vtPathLossEst
    vtPrxEst = sPar['txPower'] - vtPathLossEst + vtShadCorrEst + vtDesPequeEst

    # Estimação da CDF do desvanecimento de pequena escala
    vtn = np.arange(1, sPar['nCDF'] + 1)
    xCDF = 1.2 ** (vtn - 1) * 0.01

    # Cálculo da CDF
    cdffn = np.zeros(sPar['nCDF'])
    for ik in range(sPar['nCDF']):
        cdffn[ik] = np.sum(vtEnvNorm <= xCDF[ik])

    # Monta estrutura do histograma
    vtXCcdfEst = 20 * np.log10(xCDF)
    vtYCcdfEst = cdffn / cdffn[-1]  # Normaliza a CDF

    # Cria a estrutura de saída
    sOut = {
        'vtDistEst': vtDistEst,
        'vtPathLossEst': vtPathLossEst,
        'dNEst': dNEst,
        'vtShadCorrEst': vtShadCorrEst,
        'dStdShadEst': dStdShadEst,
        'dStdMeanShadEst': dStdMeanShadEst,
        'vtDesPequeEst': vtDesPequeEst,
        'vtPrxEst': vtPrxEst,
        'vtXCcdfEst': vtXCcdfEst,
        'vtYCcdfEst': vtYCcdfEst,
        'vtEnvNorm': vtEnvNorm
    }

    return sOut

# Parâmetros para geração do canal sintético
sPar = {
    'd0': 5,                         # distância de referência d0
    'P0': 0,                         # Potência medida na distância de referência d0 (em dBm)
    'nPoints': 50000,                # Número de amostras da rota de medição
    'totalLength': 100,              # Distância final da rota de medição
    'n': 4,                          # Expoente de perda de percurso
    'sigma': 6,                      # Desvio padrão do shadowing em dB
    'shadowingWindow': 200,          # Tamanho da janela de correlação do shadowing
    'm': 4,                          # Parâmetro de Nakagami
    'txPower': 0,                    # Potência de transmissão em dBm
    'nCDF': 40,                      # Número de pontos da CDF normalizada
    'dW': 100,                       # Janela de estimação do sombreamento
    'chFileName': 'Prx_sintetico'   # Nome do arquivo
}

# Distância entre pontos de medição
sPar['dMed'] = sPar['totalLength'] / sPar['nPoints']

# Chama função que gera o canal sintético
vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm = fGeraCanal(sPar)

# Mostra informações do canal sintético
print('Canal sintético:')
print(f'   Média do sombreamento: {np.mean(vtShadCorr)}')
print(f'   Std do sombreamento: {np.std(vtShadCorr)}')
print(f'   Janela de correlação do sombreamento: {sPar["shadowingWindow"]} amostras')
print(f'   Expoente de path loss: {sPar["n"]}')
print(f'   m de Nakagami: {sPar["m"]}')

# Várias janelas de filtragem para testar a estimação
vtW = [10, 50, 150, 200]
vtMSEShad = []
vtMSEFad = []

for iw in range(len(vtW)):
    # Configura valor da janela de filtragem
    sPar['dW'] = vtW[iw]
    # Chama função que estima o canal sintético
    sOut = fEstimaCanal(sPar)

    # Parser de variáveis
    vtDistEst = sOut['vtDistEst']
    vtPathLossEst = sOut['vtPathLossEst']
    dNEst = sOut['dNEst']
    vtShadCorrEst = sOut['vtShadCorrEst']
    dStdShadEst = sOut['dStdShadEst']
    dStdMeanShadEst = sOut['dStdMeanShadEst']
    vtDesPequeEst = sOut['vtDesPequeEst']
    vtPrxEst = sOut['vtPrxEst']
    vtXCcdfEst = sOut['vtXCcdfEst']
    vtYCcdfEst = sOut['vtYCcdfEst']

    vtDistLogEst = np.log10(vtDistEst)
    vtDistLog = np.log10(vtDist)

    # MSE com Shadowing conhecido
    dMeiaJanela = round((sPar['dW'] - 1) / 2)
    # Garante que ambos os vetores tenham o mesmo tamanho para cálculo do MSE
    min_len_shadowing = min(len(vtShadCorr[dMeiaJanela: -dMeiaJanela]), len(vtShadCorrEst))
    vtMSEShad.append(mean_squared_error(vtShadCorr[dMeiaJanela: dMeiaJanela + min_len_shadowing], vtShadCorrEst[:min_len_shadowing]))
    #vtMSEShad.append(mean_squared_error(vtShadCorr[dMeiaJanela: -dMeiaJanela], vtShadCorrEst))

    # MSE com Fading conhecido
    min_len_fading = min(len(vtFading[dMeiaJanela:-dMeiaJanela]), len(vtDesPequeEst))
    vtMSEFad.append(mean_squared_error(vtFading[dMeiaJanela:dMeiaJanela + min_len_fading], vtDesPequeEst[:min_len_fading]))
    #vtMSEFad.append(mean_squared_error(vtFading[dMeiaJanela: -dMeiaJanela], vtDesPequeEst))

    # Exibe resultados da estimativa
    print(f'Estimação dos parâmetros de larga escala (W = {sPar["dW"]}):')
    print(f'   Expoente de perda de percurso estimado n = {dNEst}')
    print(f'   Desvio padrão do sombreamento estimado = {dStdShadEst}')
    print(f'   Média do sombreamento estimado = {dStdMeanShadEst}')
    print(f'   MSE Shadowing = {vtMSEShad[iw]}')
    print('-------------\n')

# Display informação sobre o estudo das janelas
print('Estudo na melhor janela de filtragem')
print(f'   Janelas utilizadas = {vtW}')

# Melhor janela com Shadowing conhecido
valBestShad = min(vtMSEShad)
posBestShad = vtMSEShad.index(valBestShad)
print('   Melhor MSE relativo aos valores reais do Shadowing (melhor janela):')
print(f'      Melhor janela W = {vtW[posBestShad]}: MSE Shadowing = {valBestShad}')

# Melhor janela com Fading conhecido
valBestFad = min(vtMSEFad)
posBestFad = vtMSEFad.index(valBestFad)
print('   Melhor MSE relativo aos valores reais do Fading:')
print(f'      Melhor janela W = {vtW[posBestFad]}: MSE Shadowing = {valBestFad}')

print('----------------------------------------------------------------------------------')
print(' ')