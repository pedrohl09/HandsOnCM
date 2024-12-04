# Nome do Projeto
Hands On 02 da Primeira unidade
## Descrição

[Uma breve descrição do projeto]
Este projeto consiste na estimação do expoente de perda de percurso e o desvio padrão do sombreamento. Também é plotado os gráficos das curvas: potência recebida completa (sujeita ao desvanecimento de larga e pequena escalas) vs distância; potência recebida somente sujeita ao path loss estimado vs distância; potência recebida somente sujeita ao path loss e ao sombreamento estimados vs distância. 

Também foi feito um programa de PDF fitting para mostrar qual é a melhor distribuição (e seus parâmetros) para cada janela de filtragem simulada.

## Requisitos

* MATLAB Desktop ou MATLAB Online
* Python 3.10 ou superior
* IDE python (e.g. VS Code, pycharm, spyder, etc.)


Para executar o arquivo 'Entrega.m' basta colocar o arquivo 'Prx_Real_2021_1.mat' no mesmo diretório e executar o código através do Matlab Desktop ou Matlab Online. 

Para executar Fit Tool Python - H02 é preciso instalar a biblioteca fitter (veja o tutorial para instalar https://fitter.readthedocs.io/en/latest/), biblioteca numpy e matplotlib. Troque o caminho da varável data para o caminho do arquivo DesvanecimentoPequenaEscala_WX.

