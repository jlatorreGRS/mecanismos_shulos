%TRABAJO TEORÍA DE MECANISMOS
%GRUPO 5 - AMIEKE ROSELLA HARTOG Y JAVIER LATORRE RODRÍGUEZ
%REPRESENTACIÓN

clc
clear
close all

%-----------------------
% 1- Datos del mecanismo
%-----------------------
L2=25; %[Mm]
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
w2=2;  %[rad/s]
alpha=deg2rad(45); 
beta=deg2rad(49.8863);
sigma=deg2rad(120.5891);
ro=deg2rad(115.117);
tita=deg2rad(180-74.4079);
delta=deg2rad(9.5245);

%-------------------------
% 2-Datos de la simulación
%-------------------------
duracion=5;
At=0.01;
npasos=duracion/At;
t=zeros(npasos);
avance=zeros(npasos);

%-----------------------
% 3- Bucle principal
%-----------------------
for i=1:npasos
    
    %-----------------------
    % 4- Problema de posicion
    %-----------------------
    alpha=alpha+w2*At;
    error=10000; niter=0;
    
    %-----------------
    % 5- Paso
    %-----------------
    while error > 1e-6 && niter < 50
        
        %-----------------
        % 6- Restricciones
        %-----------------
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
        
        %-----------------
        % 7- Jacobiano
        %-----------------
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
        
        %-----------------
        % 8- Aq
        %-----------------
        Aq=-phiq\[phi;0];
        
        %-----------------
        % 9- Actualizar coordenadas
        %-----------------
        x1=x1+Aq(1);
        y1=y1+Aq(2);
        x2=x2+Aq(3);
        y2=y2+Aq(4);
        sigma=sigma+Aq(6);      
        x3=x3+Aq(7);
        y3=y3+Aq(8);
        tita=tita+Aq(9);
        
        %-----------------
        % 10- Calcular error
        %-----------------
        error=norm(Aq);
        niter=niter+1;
    end
    
    clf
    line([xA,x1],[yA,y1],'color','r','LineWidth', 4,'Marker', 'o','markers',6)
    line([x2,xB],[y2,yB],'color','k','LineWidth', 4,'Marker', 'o','markers',6)
    line([x2,x3],[y2,y3],'color','k','LineWidth', 4,'Marker', 'o','markers',6)
    line([x3,xC],[y3,yC],'color','k','LineWidth', 4,'Marker', 'o','markers',6)
    xlim([-400 400])
    ylim([-300 300])
    pause(1e-10)
    t(i)=(i-1)*0.01;
    avance(i)=y2;
    avancez2(i)=x2;
end

figure
plot(t,avance,t,avancez2,'--')
title('Avance/retroceso de herramienta')
xlabel('tiempo (s)')
ylabel('Coordenada X [- -](mm), Y [-](mm)')
