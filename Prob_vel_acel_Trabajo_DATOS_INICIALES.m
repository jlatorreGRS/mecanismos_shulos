%TRABAJO TEORÍA DE MECANISMOS 
%GRUPO 5 - AMIEKE ROSELLA HARTOG Y JAVIER LATORRE RODRÍGUEZ
%PROBLEMAS DE VELOCIDAD Y ACELERACIÓN

clc
clear
close all

%-----------------------
% 1- Datos del mecanismo
%-----------------------
L2=25; %[mm]
L4=130;
L5=175;
L6=110;
x1=17.67;
y1=17.67;
x2=-16.15;
y2=74.91;
x3=156.433;
y3=45.95;
xA=0;
yA=0;
xB=50;
yB=-37;
xC=186;
yC=-60;
w2=4;  %[rad/s]
alpha=deg2rad(45);
beta=deg2rad(49.8863);
sigma=deg2rad(120.5891); %sigma=59.4109;
ro=deg2rad(115.117);
tita=deg2rad(180-74.4079);
delta=deg2rad(9.5245);


%---------------------------------------
% 2- Vector de coordenadas generalizadas
%---------------------------------------
q=[x1,y1,x2,y2,alpha,sigma,x3,y3,tita]';


%----------------------------------------------------
% 3- Evaluación de ecuaciones de restricción de "phi"
%----------------------------------------------------
phi=zeros(8,1);
phi(1)=x1^2+y1^2-L2^2;
phi(2)=(x2-xB)^2+(y2-yB)^2-L4^2;
phi(3)=(x2-xB)*(y1-yB)-(x1-xB)*(y2-yB);
if abs(sin(alpha))>0.7
    phi(4)=x1-L2*cos(alpha);            
else
    phi(4)=y1-L2*sin(alpha);            
end
phi(5)=x2-xB-L4*cos(sigma);
phi(6)=(x3-x2)^2+(y3-y2)^2-L5^2;
phi(7)=(x3-xC)^2+(y3-yC)^2-L6^2;
phi(8)=x3-xC-L6*cos(tita);      


%----------------------------
% 4- Evaluación del Jacobiano
%----------------------------
phiq=zeros(9,9);
phiq(1,:) = [ 2*x1 2*y1 0 0 0 0 0 0 0 ];
phiq(2,:) = [ 0 0 2*(x2-xB) 2*(y2-yB) 0 0 0 0 0 ];
phiq(3,:) = [-(y2-yB) (x2-xB) (y1-yB) -(x1-xB) 0 0 0 0 0 ];
if abs(sin(alpha))>0.7
    phiq(4,:) = [ 1 0 0 0  L2*sin(alpha) 0 0 0 0  ];            
else
    phiq(4,:) = [ 0 1 0 0 -L2*cos(alpha) 0 0 0 0  ];            
end
phiq(5,:) = [ 0 0 1 0 0 L4*sin(sigma) 0 0 0  ];
phiq(6,:) = [ 0 0 -2*(x3-x2) -2*(y3-y2) 0 0 2*(x3-x2) 2*(y3-y2) 0 ];
phiq(7,:) = [ 0 0 0 0 0 0 2*(x3-xC) 2*(y3-yC) 0 ];
phiq(8,:) = [ 0 0 0 0 0 0 1 0 L6*sin(tita) ];
phiq(9,:) = [ 0 0 0 0 1 0 0 0 0 ];


%---------------------------------------------------------------
% 5- Vector derecho de la ecuación de velocidades
%---------------------------------------------------------------
rhsv=[0 0 0 0 0 0 0 0 w2]';


%---------------------------------
% 6- Resolución de las velocidades
%---------------------------------
qp=phiq\rhsv;


%-----------------------------------------
% 7- Extracción de valores del vector "qp"
%-----------------------------------------
x1p=qp(1); y1p=qp(2);
x2p=qp(3);y2p=qp(4);
w2=qp(5); w4=qp(6);         
x3p=qp(7);y3p=qp(8);
w6=qp(9);              


%----------------------------------------------------
% 8- Evaluación de la derivada temporal del Jacobiano
%----------------------------------------------------
phiqp=zeros(8,9);
phiqp(1,:) = [ 2*x1p 2*y1p 0 0 0 0 0 0 0 ];
phiqp(2,:) = [ 0 0 2*x2p 2*y2p 0 0 0 0 0 ];
phiqp(3,:) = [-y2p x2p y1p -x1p 0 0 0 0 0 ];
if abs(sin(alpha))>0.7
    phiqp(4,:) = [ 0 0 0 0  L2*sin(alpha)*w2 0 0 0 0  ];            
else     
    phiqp(4,:) = [ 0 0 0 0 L2*cos(alpha)*w2 0 0 0 0  ];            
end
phiqp(5,:) = [ 0 0 0 0 0 L4*cos(sigma)*w4 0 0 0  ];
phiqp(6,:) = [ 0 0 -2*(x3p-x2p) -2*(y3p-y2p) 0 0 2*(x3p-x2p) 2*(y3p-y2p) 0 ];
phiqp(7,:) = [ 0 0 0 0 0 0 2*x3p 2*y3p 0 ];
phiqp(8,:) = [ 0 0 0 0 0 0 0 0 L6*cos(tita)*w6 ];


%------------------------
% 9- Producto de phiqp*qp
%------------------------
phiqpqp = phiqp*qp;


%---------------------------------------
% 10- Solución del sistema phiq·qpp=rhsa
%---------------------------------------
rhsa=[-phiqpqp; 0];
qpp = phiq\rhsa;


%VELOCIDADES
fprintf('\n-----------\nVelocidades\n-----------\nv1x = %f mm/s\nv1y = %f mm/s\nv2x = %f mm/s\nv2y = %f mm/s\nw2 = %f rad/s\nw4 = %f rad/s\nv3x = %f mm/s\nv3y = %f mm/s\nw6 = %f rad/s\n',...
    qp(1), qp(2), qp(3), qp(4), qp(5), qp(6),qp(7),qp(8),qp(9))  

%ACELERACIONES
fprintf('\n-------------\nAceleraciones\n-------------\na1x = %f mm/s^2\na1y = %f mm/s^2\na2x = %f mm/s^2\na2y = %f mm/s^2\naang2 = %f rad/s^2\naang4 = %f rad/s^2\na3x = %f mm/s^2\na3y = %f mm/s^2\naang6 = %f rad/s^2\n',...
    qpp(1), qpp(2), qpp(3), qpp(4), qpp(5), qpp(6),qpp(7),qpp(8),qpp(9))   
