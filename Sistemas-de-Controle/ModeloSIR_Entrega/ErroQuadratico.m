function [Eq] = ErroQuadraticoSIR(X,Rd,Id)
% Esta função determina o erro quadrático
% entre a saída real e a saída do modelo.
% u - entrada real
% y - saída real
% t - tempo real
% X = [Km pm];

beta = X(1);
r = X(2);


%%% Tamanho da população - N
N = 500;

%%% Valores iniciais do modelo
i0 = 1;
r0 = 0;
s0 = N - i0;

options = simset('solver','ode45','MaxStep',0.01,'InitialStep',0.01,'ReturnWorkspaceOutputs', 'on');
y = sim("SimulacaoSIR.slx",1000,options);

%%% Saida do modelo

I = y.logsout{1}.Values.Data; % Saída I
R = y.logsout{2}.Values.Data; % Saída R

I = resample(I,length(I),length(Id));
R = resample(R,length(R),length(Rd));

%%% Erro quadrático

Eq = sum((I - Id).^2) + sum((R - Rd).^2); %% cálculo do erro

% plot(t,y,t,ym,"linewidth", 2)

end