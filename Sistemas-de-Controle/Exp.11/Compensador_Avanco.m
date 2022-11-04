%pkg load control %comando necessário no Octave 
clear all 
clc 
k = 22.5; 
tau = 0.18; 
G = tf([k], [tau 1]); 
I = tf(1, [1 0]); 
Gma = G*I; 
figure(1) 
subplot(2,1,1) 
rlocus(Gma) 
title('Lugar das Raízes Gma(s)') 
Gmf = feedback(Gma,1);    
PO = 2; 
%porcentagem de overshoot 
zeta = sqrt(log(PO/100)^2/(pi^2+log(PO/100)^2)); 
ts = 0.5; 
%tempo de acomodação em segundos 
Wn = 4/(ts*zeta); 
%polo desejado que atende aos requisitos 
polo_d = -zeta*Wn + Wn*sqrt(1-zeta^2)*i; %apenas o complexo positivo    
%pzmap(Gma)   
%polos e zeros do sistema em malha aberta 
polos_ma = pole(Gma)'; zeros_ma = zero(Gma); 
 
%escolheu-se o zero igual a parte real do polo desejado 
zc = -real(polo_d); 
zeros = [-zc zeros_ma];   
%ângulos dos polos 
for n=1:length(polos_ma)     
    if (real(polos_ma(n))>= real(polo_d))         
        ang_p(n) = pi - atan(imag(polo_d)/(abs( real(polo_d))-abs(real(polos_ma(n)))));     
    else
        ang_p(n) = atan(imag(polo_d)/(abs( real(polo_d))-abs(real(polos_ma(n)))));     
    end
end
ang_p = ang_p*180/pi; 
%ângulos dos zeros 
for n=1:length(zeros)     
    if ( real(zeros(n)) >= real(polo_d) )         
    ang_z(n) = pi - atan(imag(polo_d) /(abs(real(polo_d)) - abs(real(zeros(n)))) );     
    else
    ang_z(n) = atan( imag(polo_d) / (-abs ( real(polo_d)) + abs (real(zeros(n)))));     
    end
end
ang_z = ang_z*180/pi;
%Condição de fase 
sum_pz = sum(ang_p)- sum(ang_z); 
if sum_pz > 180     
    angulo_polo_graus = sum_pz - 180; 
else
    angulo_polo_graus = 180 - sum_pz; 
end
%convertendo para radianos
angulo_polo_rad = angulo_polo_graus*pi/180; 
%polo do compensador 
pc = (imag(polo_d)/( tan(angulo_polo_rad) ))  + zc ; 
C = tf([1 zc],[1 pc]);   
figure(1) 
subplot(2,1,2) 
rlocus(C*Gma) 
title('Lugar das Raízes C(s)*Gma(s)')   
%Condição de módulo 
kc_num = abs([polo_d+pc polo_d polo_d*tau+1]); 
kc_den = abs((polo_d+zc))*k; 
kc = prod(kc_num)/prod(kc_den);   
%imprime os resultados na tela 
fprintf('\nPolos da função em malha aberta\n') 
fprintf('%.2f\n', polos_ma) 
fprintf("\nPolo desejado\n") 
polo_d 
fprintf("\nZero do compensador\n") 
fprintf('%.2f\n', -zc) 
fprintf('\nPolo do compensador\n') 
fprintf('%.2f\n', -pc)
fprintf('\nFunção de transferência do compensador (sem o ganho)\n') 
C 
fprintf('\nGanho do compensador\n') 
fprintf('%.2f\n', kc)   
Gcmf = feedback(kc*C*Gma,1); 
figure(2) 
step(pi*Gcmf)