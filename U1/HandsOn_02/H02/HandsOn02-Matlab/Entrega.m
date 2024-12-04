load('/MATLAB Drive/Prx_Real_2021_1.mat')

txPower = 0;
vtDist = dPath';
vtPrxdBm = dPrx; 

% Transforma potência em mWatts
vtPtrxmW = 10.^(vtPrxdBm/10);
nSamples = length(vtPtrxmW);
%
vtW = [2 5 10]; %Janela de filtragem
for iw = 1 : length(vtW)
    %Config. da janela de filtramgem
    dW = vtW(iw);
    %Vetores para Ch. estimado
    vtDesLarga = []; %Des. Larga escala
    vtDesPequeEst = []; %Des. Pequena escala
    %
    %Cálculo do des. lento e rápido
    dMeiaJanela = round((dW-1)/2);  % Meia janela
    ij = 1;
    for ik = dMeiaJanela + 1 : nSamples - dMeiaJanela
        % Desvanecimento de larga escala: perda de percurso + sombreamento [dB]
        vtDesLarga(ij) = 10*log10(mean(vtPtrxmW(ik-dMeiaJanela:ik+dMeiaJanela)));
        % Desvanecimento de pequena escala [dB]
        vtDesPequeEst(ij) = vtPrxdBm(ik)-vtDesLarga(ij);
        ij = ij + 1;
    end
    %
    %Armazenamento dos dados do vetor de desvanecimento em pequena escala
    chFileName = ['DesvanecimentoPequenaEscala', num2str(dW)];
    dlmwrite([chFileName '.txt'], vtDesPequeEst', 'delimiter', '\t');
    %
    % Cálculo da envoltória normalizada (para efeitos de cálculo do fading)
    indexes = dMeiaJanela+1 : nSamples-dMeiaJanela;
    %vtPrxW = ((10.^(vtPrxdBm(indexes)./10))/1000);
    vtPtrxmWNew = 10.^(vtPrxdBm(indexes)/10);
    desLarga_Lin = (10.^(vtDesLarga(1:length(indexes))./10));
    envNormal = sqrt(vtPtrxmWNew)./sqrt(desLarga_Lin);
    %
    % Ajuste no tamanho dos vetores devido a filtragem
    vtDistEst = vtDist( dMeiaJanela+1 : nSamples-dMeiaJanela );
    vtPrxdBm = vtPrxdBm( dMeiaJanela+1 : nSamples-dMeiaJanela );
    %
    % Cálculo reta de perda de percurso
    vtDistLog = log10(vtDist);
    vtDistLogEst = log10(vtDistEst);
    % Cálculo do coeficientes da reta que melhor se caracteriza a perda de percurso
    dCoefReta = polyfit(vtDistLogEst,vtPrxdBm,1);
    % Expoente de perda de percurso estimado
    dNEst = -dCoefReta(1)/10;
    disp(['Estimação dos parâmetros de larga escala (W = ' num2str(dW) '):'])
    disp(['   Expoente de perda de percurso estimado n = ' num2str(dNEst)]);
    % Perda de percurso estimada para os pontos de medição
    vtPathLossEst = polyval(dCoefReta,vtDistLogEst);  
    %
    % Sombreamento
    vtShadCorrEst = vtDesLarga - vtPathLossEst;
    % Calcula a variância do sombreamento estimado
    dStdShadEst = std(vtShadCorrEst);
    dStdMeanShadEst = mean(vtShadCorrEst);
    vtPathLossEst = - vtPathLossEst;
    vtPrxEst = txPower - vtPathLossEst + vtShadCorrEst + vtDesPequeEst;
    %

    %plot dos gráficos 
    figure;
    %Potência recebido do canal sem Path Loss e Shadowing
    plot(vtDistLogEst, vtPrxEst); hold all;
    %Potência recebido do canal com Path Loss
    plot(vtDistLogEst, txPower - vtPathLossEst); 
    %Potência recebido do canal com Path Loss e Shadowing
    plot(vtDistLogEst, txPower - vtPathLossEst + vtShadCorrEst)
    %
    title(['Potência recebido Rx x Log distância para W = ' num2str(dW)])
    xlabel('log_{10}(d)');
    ylabel('Potência [dBm]');
    legend('Prx canal completo', 'Prx (somente perda de percurso)', 'Prx (perda de percurso + sombreamento)');
end