function [Eq] = ErroQuadraticoSIR(X,Sd,Id,Rd)
% Esta função determina o erro quadrático
% entre a saída real e a saída do modelo.
% u - entrada real
% y - saída real
% t - tempo real
% X = [Km pm];
global beta;
global r;

beta = X(1);
r = X(2);


%%% Tamanho da população - N
N = 500;

%%% Valores iniciais do modelo
i0 = 1;
r0 = 0;
s0 = N - i0;

options = simset('solver','ode45','MaxStep',0.1,'InitialStep',0.1,'ReturnWorkspaceOutputs', 'on');
simOut = sim("SimulacaoSIR.slx",1000,options);

%%% Saida do modelo

S = simOut.logsout{1}.Values.Data;  
I = simOut.logsout{2}.Values.Data;
R = simOut.logsout{3}.Values.Data;

S = resample(S,length(S),length(Sd));
I = resample(I,length(I),length(Id));
R = resample(R,length(R),length(Rd));

%%% Erro quadrático

Eq = sum((I - Id).^2) + sum((R - Rd).^2) + sum((S - Sd).^2) %% cálculo do erro

% plot(t,y,t,ym,"linewidth", 2)

end