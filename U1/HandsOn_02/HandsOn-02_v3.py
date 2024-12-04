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

    # # Salva as variáveis do canal em um arquivo .mat
    # sio.savemat(sPar['chFileName'] + '.mat', {
    #     'vtDist': vtDist,
    #     'vtPathLoss': vtPathLoss,
    #     'vtShadCorr': vtShadCorr,
    #     'vtFading': vtFading,
    #     'vtPrxdBm': vtPrxdBm
    # })

    # Salva as variáveis do canal em um arquivo .txt
    data = np.column_stack((vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm))
    np.savetxt(sPar['chFileName'] + '.txt', data, header="vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm", delimiter=",", comments='')

    return vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm



def f_estima_canal(sPar):
    # Verificar a extensão do arquivo para escolher a forma de leitura
    if sPar['chFileName'].endswith('.mat'):
        # Carregar os dados do canal gerado a partir do arquivo .mat
        canal_data = loadmat(sPar['chFileName'])
        vtDist = canal_data['vtDist'].flatten()
        vtPrxdBm = canal_data['vtPrxdBm'].flatten()
    elif sPar['chFileName'].endswith('.txt'):
        # Carregar dados do canal a partir de um arquivo .txt
        data = np.loadtxt(sPar['chFileName'])
        vtDist = data[:, 0]    # Coluna 0 é a distância
        vtPrxdBm = data[:, 1]  # Coluna 1 é a potência recebida em dBm
    else:
        raise ValueError("Formato de arquivo não suportado. Use '.mat' ou '.txt'.")

    # Converter potência recebida em mW
    vtPtrxmW = 10 ** (vtPrxdBm / 10)
    nSamples = len(vtPtrxmW)
    
    # Inicializar vetores de saída
    vtDesLarga = []
    vtDesPequeEst = []
    
    # Calcular o desvanecimento lento e rápido
    dMeiaJanela = round((sPar['dW'] - 1) / 2)
    for ik in range(dMeiaJanela, nSamples - dMeiaJanela):
        # Desvanecimento de larga escala (path loss + shadowing)
        mean_val = np.mean(vtPtrxmW[ik - dMeiaJanela:ik + dMeiaJanela + 1])
        vtDesLarga.append(10 * np.log10(mean_val))
        
        # Desvanecimento de pequena escala
        vtDesPequeEst.append(vtPrxdBm[ik] - vtDesLarga[-1])

    # Normalizar a envoltória para cálculo do fading
    indexes = np.arange(dMeiaJanela, nSamples - dMeiaJanela)
    vtPtrxmWNew = 10 ** (vtPrxdBm[indexes] / 10)
    desLarga_Lin = 10 ** (np.array(vtDesLarga) / 10)
    vtEnvNorm = np.sqrt(vtPtrxmWNew) / np.sqrt(desLarga_Lin)

    # Ajuste nos tamanhos dos vetores
    vtDistEst = vtDist[indexes]
    vtPrxdBm = vtPrxdBm[indexes]

    # Calcular coeficientes de reta para path loss
    vtDistLog = np.log10(vtDist)
    vtDistLogEst = np.log10(vtDistEst)
    dCoefReta = np.polyfit(vtDistLogEst, vtPrxdBm, 1)
    dNEst = -dCoefReta[0] / 10
    vtPathLossEst = np.polyval(dCoefReta, vtDistLogEst)

    # Calcular o shadowing
    vtShadCorrEst = np.array(vtDesLarga) - vtPathLossEst
    dStdShadEst = np.std(vtShadCorrEst)
    dStdMeanShadEst = np.mean(vtShadCorrEst)
    vtPathLossEst = -vtPathLossEst
    vtPrxEst = sPar['txPower'] - vtPathLossEst + vtShadCorrEst + np.array(vtDesPequeEst)

    # Cálculo dos pontos do eixo x da CDF
    vtn = np.arange(1, sPar['nCDF'] + 1)
    xCDF = 1.2 ** (vtn - 1) * 0.01

    # Calcular a CDF para fading
    cdffn = np.zeros(sPar['nCDF'])
    for ik in range(sPar['nCDF']):
        cdffn[ik] = np.sum(vtEnvNorm <= xCDF[ik])
    cdffn = cdffn / cdffn[-1]  # Normalizar a CDF

    # Converter xCDF para escala dB
    vtXCcdfEst = 20 * np.log10(xCDF)
    vtYCcdfEst = cdffn

    # Estrutura de saída
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
vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm = gera_canal(sPar)

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
    sOut = f_estima_canal(sPar)
    
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