# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 10:59:06 2024

@author: pedro
"""

import numpy as np
import scipy.io as sio
from scipy.stats import nakagami, gamma
import matplotlib.pyplot as plt
from scipy.special import gammainc
from scipy.io import loadmat
from scipy import stats
from sklearn.metrics import mean_squared_error as mse


def fGeraCanal(sPar):
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

    # Salva as variáveis do canal em um arquivo .npz
    np.savez(sPar['chFileName'] + '.npz', 
             vtDist=vtDist, 
             vtPathLoss=vtPathLoss, 
             vtShadCorr=vtShadCorr, 
             vtFading=vtFading, 
             vtPrxdBm=vtPrxdBm)

    return vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm




def fEstimaCanal(sPar):
    # Lê canal gerado
    data = np.load(f"{sPar['chFileName']}.npz")  # Supondo que o arquivo seja .npz ou similar
    vtPrxdBm = data['vtPrxdBm']
    vtDist = data['vtDist']
    
    vtPtrxmW = 10 ** (vtPrxdBm / 10)
    nSamples = len(vtPtrxmW)
    
    # Vetores para canal estimado
    vtDesLarga = []
    vtDesPequeEst = []
    
    # Cálculo do desvanecimento lento e rápido
    dMeiaJanela = round((sPar['dW'] - 1) / 2)  # Meia janela
    ij = 0
    for ik in range(dMeiaJanela, nSamples - dMeiaJanela):
        # Desvanecimento de larga escala: perda de percurso + sombreamento [dB]
        vtDesLarga.append(10 * np.log10(np.mean(vtPtrxmW[ik - dMeiaJanela:ik + dMeiaJanela])))
        # Desvanecimento de pequena escala [dB]
        vtDesPequeEst.append(vtPrxdBm[ik] - vtDesLarga[ij])
        ij += 1
    
    # Cálculo da envoltória normalizada (para efeitos de cálculo do fading)
    indexes = range(dMeiaJanela, nSamples - dMeiaJanela)
    vtPtrxmWNew = 10 ** (vtPrxdBm[indexes] / 10)
    desLarga_Lin = 10 ** (np.array(vtDesLarga)[:len(indexes)] / 10)
    vtEnvNorm = np.sqrt(vtPtrxmWNew) / np.sqrt(desLarga_Lin)
    
    # Ajuste no tamanho dos vetores devido a filtragem
    vtDistEst = vtDist[dMeiaJanela:nSamples - dMeiaJanela]
    vtPrxdBm = vtPrxdBm[dMeiaJanela:nSamples - dMeiaJanela]
    
    # Cálculo reta de perda de percurso
    vtDistLog = np.log10(vtDist)
    vtDistLogEst = np.log10(vtDistEst)
    
    # Cálculo dos coeficientes da reta que melhor se caracteriza a perda de percurso
    dCoefReta = np.polyfit(vtDistLogEst, vtPrxdBm, 1)
    
    # Expoente de perda de percurso estimado
    dNEst = -dCoefReta[0] / 10
    
    # Perda de percurso estimada para os pontos de medição
    vtPathLossEst = np.polyval(dCoefReta, vtDistLogEst)
    
    # Sombreamento
    vtShadCorrEst = np.array(vtDesLarga) - vtPathLossEst
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
        den = 0
        for ij in range(len(vtEnvNorm)):
            if vtEnvNorm[ij] <= xCDF[ik]:
                den += 1
        cdffn[ik] = den
    
    # Monta estrutura do histograma
    vtXCcdfEst = 20 * np.log10(xCDF)
    vtYCcdfEst = cdffn / cdffn[-1]
    
    # Montando a saída
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



# Parâmetros de geração do canal
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

# Chamada da função que gera o canal sintético
vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm = fGeraCanal(sPar)

# Exibe informações do canal sintético
print('Canal sintético:')
print(f'   Média do sombreamento: {np.mean(vtShadCorr)}')
print(f'   Desvio padrão do sombreamento: {np.std(vtShadCorr)}')
print(f'   Janela de correlação do sombreamento: {sPar["shadowingWindow"]} amostras')
print(f'   Expoente de path loss: {sPar["n"]}')
print(f'   m de Nakagami: {sPar["m"]}')

# Janelas de filtragem para testar a estimação
vtW = [10, 50, 150, 200]
vtMSEShad = []
vtMSEFad = []

for w in vtW:
    sPar['dW'] = w
    # Estimação do canal sintético
    sOut = fEstimaCanal(sPar)
    
    # Parsing dos resultados da estimativa
    vtDistEst = sOut['vtDistEst']
    vtPathLossEst = sOut['vtPathLossEst']
    dNEst = sOut['dNEst']
    vtShadCorrEst = sOut['vtShadCorrEst']
    dStdShadEst = sOut['dStdShadEst']
    dStdMeanShadEst = sOut['dStdMeanShadEst']
    vtDesPequeEst = sOut['vtDesPequeEst']
    vtPrxEst = sOut['vtPrxEst']
    
    # Cálculo do MSE com Shadowing conhecido
    dMeiaJanela = round((sPar['dW']-1)/2)
    vtMSEShad.append(mse(vtShadCorr[dMeiaJanela: -dMeiaJanela], vtShadCorrEst))
    
    # Cálculo do MSE com Fading conhecido
    vtMSEFad.append(mse(vtFading[dMeiaJanela: -dMeiaJanela], vtDesPequeEst))
    
    print(f'Estimação dos parâmetros de larga escala (W = {sPar["dW"]}):')
    print(f'   Expoente de perda de percurso estimado n = {dNEst}')
    print(f'   Desvio padrão do sombreamento estimado = {dStdShadEst}')
    print(f'   Média do sombreamento estimado = {dStdMeanShadEst}')
    print(f'   MSE Shadowing = {vtMSEShad[-1]}')
    print('----')

# Exibição da melhor janela de filtragem
print('Estudo na melhor janela de filtragem:')
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
print(f'      Melhor janela W = {vtW[posBestFad]}: MSE Fading = {valBestFad}')
print('----------------------------------------------------------------------------------')