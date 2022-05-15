%function Evaluacion_polos_5_2_2

%Evaluamos los polos de la ecuación J*s^3 + ba*s^2 + Ksa*s + Ksia
Jl = 0.2520+0.1260; %kg/m2 (+-0.1260)
bl = 0+0.0630; %Nm/rad.s (+-0.0630)
Tl = 6.28; %Nm
Jm = 3.1e-6; %kg.m2
bm = 1.5e-5; %Nm/rad.s
Pp = 3; %3 pares de polos
LambdaM = 0.01546; %Wb
Lq = 0.0058; %H
Ld = 0.0066; %H
Lls = 0.0008; %H
r = 314;
Rs = 1.02; %ohm a 40° y 1.32ohm a 115°
Jeq = Jl/r^2 + Jm;
beq = bl/r^2 + bm;
wpos = 800;
n = 2.5;

%Ksa = Jeq*n*wpos^2;
%Ksia = Jeq*wpos^3;
%ba = Jeq*n*wpos;

%Valores de ganancias, control vectorial
Rq = 29;
Rd = 33;
R0 = 4;

%POLOS DEL CONTROLADOR
tf_pid = tf([1,0],[Jeq,ba,Ksa,Ksia])
pole(tf_pid)

%POLOS MODULADOR DE TORQUE
tf_mt = tf([0 1],[Lq/Rq 1])
pole(tf_mt)

%POLOS DE LA PLANTA
num = [0,0,(Lq),(3/2*Pp*LambdaM+Rs)];
den = [(Jeq*Lq),(Rs*Jeq + Lq*beq),(3/2*Pp^2*LambdaM^2+Rs*beq),0];
tf_pl = tf(num,den)

%PLOTEO DE TODOS LOS POLOS
ax = axes();
pzmap(ax,tf_pid,'r',tf_mt,'y',tf_pl,'b');
legend("Polos(X) del controlador PID","Polos(X) modulador de torque","Polos(X) de la planta");
%aux_1 = findall(ax, 'tag', 'PZ_Pole');
%aux_1.MarkerSize = 12;
%aux_1.LineWidth = 4;
grid on

a = findobj(gca,'type','line');
for i = 1:length(a)
    set(a(i),'markersize',12); %change marker size
    set(a(i), 'linewidth',2);  %change linewidth
end

%,tf_mt,'b',tf_pl,'g'


%end