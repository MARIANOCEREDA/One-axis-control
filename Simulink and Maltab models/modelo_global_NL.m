%function modelo_global_NL

Jl = 0.2520; %kg/m2 (+-0.1260)
bl = 0; %Nm/rad.s (+-0.630)
Tl = 6.28; %Nm
Tl_max = 6.28;
Tl_min = -6.28;
Jm = 3.1e-6; %kg.m2
bm = 1.5e-5; %Nm/rad.s
Pp = 3; %3 pares de polos
LambdaM = 0.01546; %Wb
Lq = 0.0058; %H
Ld = 0.0066; %H
Lls = 0.0008; %H
r = 314;
Jeq = Jl/r^2 + Jm;
beq = bl/r^2 + bm;
Tleq = Tl/r;
%Datos subsistema térmico
Rs = 1.02; %ohm a 40° y 1.32ohm a 115°
Cts = 0.818; %1/°C
Rts_amb = 146.7; %°C/W
Tts_amb = Rts_amb*Cts;
alfa_cu = 3.9e-3;
Rsref = 1.02; 
Tsref = 40; 
Rts = 146.7;
Tamb = 40; %Rango de -15 a 40°C


%Simulaciones
Vqs_nom = 19.596;

%Valores de ganancias, control vectorial
Rq = 5000*Lq;
Rd = 5000*Ld;
R0 = 5000*Lls;

%Valores controlador PID
wpos = 800;
n = 2.5;
Ksa = Jeq*n*wpos^2;
Ksia = Jeq*wpos^3;
ba = Jeq*n*wpos;

%Valores observador
kw_sin_mejora = (3200)^2;
kt_sin_mejora = 6400;

%Observador con accion integral
kt = 12800;
kw = 5.12e7;
ki = 6.5536e10;

syms s

p = (s+3200)^2*(s+6400)
expand(p)

%Polos funcion de transferencia filtros pasa bajos sensores corrientes
wn_iabc = 9000;
xita_iabc = 1;
num_iabc = [wn_iabc^2];
den_iabc = [1 2*wn_iabc*xita_iabc wn_iabc^2];

ft = tf(num_iabc,den_iabc);

%Polos funcion de transferencia filtros pasa bajos sensor posición
wn_t = 7000;
xita_t = 1;
num_t = [wn_t^2];
den_t = [1 2*wn_t*xita_t wn_t^2];

%Polos funcion de transferencia filtros pasa bajos sensor temperatura
tau = 20;
num_temp = [1];
den_temp = [tau 1];
ft = tf(num_temp,den_temp);

%Polos con PID
p1_max = -613.65+0i;
p2_max = -508.87+649.33i;
p3_max = -508.87-649.33i;

p1_min = -1499.5+0i;
p2_min = -542.20+383.6i;
p3_min = -542.20-383.6i;

p1_nom = -800+0i;
p2_nom = -600+529.15i;
p3_nom = -600-529.15i;

p_modulador_torque = -5000+0i;

hold on
grid on
plot(real(p2_max),imag(p2_max),'*g');
plot(real(p3_max),imag(p3_max),'*g');
plot1 = plot(real(p1_max),imag(p1_max),'*g');

plot(real(p2_min),imag(p2_min),'*b');
plot(real(p3_min),imag(p3_min),'*b');
plot2 = plot(real(p1_min),imag(p1_min),'*b');

plot(real(p2_nom),imag(p2_nom),'r*');
plot(real(p3_nom),imag(p3_nom),'r*');
plot3 = plot(real(p1_nom),imag(p1_nom),'r*');


legend([plot1 plot2 plot3],'Parámetros máximos','Parámetros mínimos','Parámetros nominales');
title('Polos PID ante varaiación de parámetros de carga');
ylabel('Imag');
xlabel('Real');

xlim([-1600 0]);
ylim([-800 800]);

hold off

%{
%% 3 - ANÁLISIS DE ESTABILIDAD - ANÁLISIS DE AUTOVALORES DE A (CON  PRIMERA TANDA DE VALORES)
A = [0 1 0;
    0 (-beq/Jeq) (3*Pp*LambdaM)/(2*Jeq);
    0 -(LambdaM*Pp)/Lq (-Rs/Lq)];

autovals = eig(A);
  
plot(autovals,'x');
xlabel("Eje real");
ylabel("Eje imaginario");
xlim([-400 400]);
ylim([-400 400]);
grid on

num = [0,0,(Lq),(3/2*Pp*LambdaM+Rs)];
den = [(Jeq*Lq),(Rs*Jeq + Lq*beq),(3/2*Pp^2*LambdaM^2+Rs*beq),0];
tf = tf(num,den)
zeros = zero(tf)
poles = pole(tf)


%Al ser polos complejos conjugados, estamos ante un sistema SUBAMORTIGUADO
%Frecuencias naturales de oscilación, las obtenemos teniendo en cuenta que
%el denominador del denominador de la ft es de la forma:
%(s^2+2*xita*wn*s+wn^2)*s 

wn = ((3/2*Pp^2*LambdaM^2+Rs*beq)/(Jeq*Lq))^0.5
xita = ((Rs*Jeq+Lq*beq)/(Jeq*Lq))/(2*wn)

%#################################################################################################
%############################### PARA MODELO LTI EQUIVALENTE #####################################
%#################################################################################################
%% ANÁLISIS DE OBSERVABILIDAD - ANÁLISIS DE AUTOVALORES/POLOS/FRECUENCIAS (CON SEGUNDA TANDA DE VALORES)

syms beq Jeq Pp LambdaM Lq Rs Ld ids0

fprintf("#########################################################");
fprintf("############MODELO LTI EQUIVALENTE ######################");
fprintf("#########################################################");
A = [0 1 0 0;
    0 (-beq/Jeq) (3*Pp*LambdaM)/(2*Jeq) 0;
    0 -(LambdaM*Pp)/Lq (-Rs/Lq) 0
    0 0 0 -Rs/Ld];
B = [0;0;1/Lq;0];
C_t = [1 0 0 1]; %Salida desde tita
C_o = [0 1 0]; %Salida desde omega

%Observabilidad desde tita de forma simbolica
obs_t = [C_t;C_t*A;C_t*A^2;C_t*A^3]
rank(obs_t)

disp(rank(obs_t));


%% ANÁLISIS DE CONTROLABILIDAD desde entrada manipulada Vqs sin considerar perturbacion Tleq
%Agregamos el eje desacoplado q, ya que trabajamos con el sistema LTI
%equivalente aumentado
syms beq Jeq Pp LambdaM Lq Rs Ld ids0
A = [0 1 0 0;
    0 (-beq/Jeq) (3*Pp*LambdaM)/(2*Jeq) 0;
    0 -(LambdaM*Pp)/Lq (-Rs/Lq) 0
    0 0 0 -Rs/Ld];
B = [0,0;0,0;1/Lq,0;0,1/Ld];
contr_vq = [B,A*B,(A^2)*B,(A^3)*B]
rank(contr_vq)

%#################################################################################################
%######################### PARA MODELO LTI EQUIVALENTE AUMENTADO #################################
%#################################################################################################
%% ANÁLISIS DE OBSERVABILIDAD - ANÁLISIS DE AUTOVALORES/POLOS/FRECUENCIAS
fprintf("#########################################################");
fprintf("############MODELO LTI EQUIVALENTE ######################");
fprintf("#########################################################");

A = [0 1 0 0;
    0 (-beq/Jeq) (3*Pp*LambdaM)/(2*Jeq) 0;
    0 -(LambdaM*Pp)/Lq (-Rs/Lq) 0;
    0 0 0 -Rs/Ld];
C_t = [1 0 0 0]; %Salida desde tita
C_o = [0 1 0 0]; %Salida desde omega

%Observabilidad desde tita de forma simbolica
obs_t = [C_t;C_t*A;C_t*A^2];
rank(obs_t)
%Observabilidad desde omega de forma simbolica
obs_o = [C_o;C_o*2;C_o*A^2];
rank(obs_o)

fprintf("Matriz de observabilidad (tita): ");
disp(obs_t)
%}


%end