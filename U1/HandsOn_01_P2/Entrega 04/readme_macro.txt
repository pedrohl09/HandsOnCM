Este projeto simula a cobertura de macrocélulas e ERBs (Estações Rádio Base) em uma grade hexagonal, usando o modelo de propagação Okumura-Hata para calcular a potência de sinal recebida em uma área urbana e determinar a cobertura de ERBs em frequências específicas.

O código desenha setores hexagonais para as ERBs, distribuídos em uma grade central e uma camada ao redor. Cada ERB calcula a perda de sinal com base em parâmetros como a potência de transmissão (EIRP), sensibilidade do receptor, altura das antenas e a distância entre pontos da grade e ERBs. A simulação gera um mapa de cobertura que indica, ponto a ponto, se há ou não sinal suficiente, além de um gráfico com a distribuição de pontos dentro e fora da cobertura.

Para usar o projeto, instale Python 3 com numpy e matplotlib (pip install numpy matplotlib), ajuste os parâmetros no código, como raio de cobertura e frequência, e execute o script. A simulação exibirá o mapa e o gráfico de cobertura, possibilitando a análise da eficiência para cada frequência configurada.

Com as bibliotecas instaladas, navegue até a pasta onde o arquivo "HandsOn_01_P2_5(Macro).py" está localizado utilizando o comando "cd /caminho/para/sua/pasta". Após acessar o diretório correto, execute o script com o comando "python HandsOn_01_P2_5(Macro).py" (ou "python3 HandsOn_01_P2_5(Macro).py" se estiver utilizando Python 3). Caso ocorra algum erro relacionado ao pip, utilize o comando "python -m ensurepip --upgrade" para atualizá-lo.

Se tudo estiver configurado corretamente, o código será executado sem problemas.
