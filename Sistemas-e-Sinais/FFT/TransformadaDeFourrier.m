
%   Enrico Giannobile - 19.00610-0

%   Análise de Fourrier de sinais de áudio
%   Boas práticas

% clear all;
% close all;
% clc;
% 
% %%Lendo a faixa GaitaBlues
% 
% [Y,FS] = audioread("AudioGaita.wav");   %Abertura de arquivo
% 
% Criando a variavel tempo
% 
% TS = 1/FS;                      % Tempo entre amostras
% Npontos = length(Y);            % Numero de pontos do vetor Y
% Vfinal = (Npontos - 1)*TS;      % Valor final
% 
% tempo = linspace(0,Vfinal,Npontos);
% 
% figure(1);
% 
% subplot(2,1,1);
% plot(tempo,Y(:,1));
% xlabel("Tempo em segundos");
% ylabel("Amplitude do sinal");
% title("Canal único");

% x = Y;
% 
% soma = 0;
% N = length(Y);
% for a = 1:N
%     n = a-1;
%     soma = 0;
%     for b = 1:N
%         k = b-1;
% 
%     soma = soma + Y(a)*exp(-j*2*pi*n*k/N);
% 
% 
%     end
%     x(a) = soma/N;
% end

    

[g_k,Fs] = audioread('AudioGaita.wav'); %abertura de arquivo

g_k = g_k(1:1000); %p/ transformada matricial;
N = length(g_k);                   % Numero de pontos do vetor
Ts = inv(Fs);                      % Tempo de amostragem
Tfinal = N*Ts;                     % Tempo Final
t = linspace (0,Tfinal,N);

% figure(1);

% plot(t,g_k,"LineWidth",3);
% grid;
% title("Som de uma gaita blues");
% xlabel("Tempo em segundos");
% ylabel("Amplitude do Sinal");

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 2 - Implementar série de Fourrier usando estrutura for
% 
% somatoria = 0;
% 
% %variavel computacional p --> variavel n
% %variavel computacional q --> variavel k
% 
% tic;                       %inicia o contador
% for p=1:N
% 
%     somatoria = 0;                % valor inicial da somatoria para cada n
%     n = p-1;                      % determina n dado p
%     
%     for q = 1:N
%         k = q-1;                  % determina k dado q
% 
%         somatoria = somatoria + g_k(q)*exp(-j*2*pi*n*k/N);
%     end
%     X(p) = somatoria;
% 
%     tempo_for = toc;       % para o contador e grava o valor;
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 - Transformada de Fourrier de forma matricial;



%%%%%%%%% Criando uma função para cálculo de WN

WN = @(N) exp(-j*2*pi/N);

W = WN(N);                    % determina a matriz para N pontos

MatrizFourrier = W*ones(N,N); % Matriz de Fourrier


%%%%%%%% Criandoo os vetores n e k

Vetor_n = [0:1:N-1]';
Vetor_k = [0:1:N-1]';

%%%%% Criando a matriz nk

Matriz_nk = Vetor_n*Vetor_k';

X_matriz =( MatrizFourrier.^Matriz_nk) *g_k*(1/N);

plot(t,X_matriz);

