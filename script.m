clear; clc; clear variables; close all;

%%Carga de Datos

for i = 0:8
    filename = strcat("data/", num2str(i), '0.txt');
    y = dlmread(filename);
    u(:,i+1) = y(:,2);
    v(:,i+1) = y(:,3);
    if i == 0
      t=y(:,1);
    end
end

%% 1)

figure (1);

subplot(2,1,1)
plot(t,u(:,1));
hold on;
plot(t,v(:,1));
grid on;
xlabel("Tiempo [s]");ylabel("Velocidad [m/s]");
title("y00");
legend('u','v');
hold on;

subplot(2,1,2)
plot(t,u(:,6));
hold on;
plot(t,v(:,6));
grid on;
xlabel("Tiempo [s]");ylabel("Velocidad [m/s]");
title("y50");
legend('u','v');

axes( 'visible', 'off', 'title', 'Punto 1)' );






%% 2)
%%a)
um=mean(u);
vm=mean(v);

y=0:10:80;

figure (2);

subplot(2,2,1)
plot(y,um);
hold on;
plot(y,vm);
grid on;
xlabel("y [mm]");ylabel("Velocidad Media [m/s]");title("");
legend('um','vm');

%b)
uf=u-um; %%Matriz fluctuación u'(t,y)
vf=v-vm; %%Matriz fluctuación v'(t,y)

V = sqrt(um.^2 + vm.^2);

desviacionU = sqrt(var(u));
desviacionV = sqrt(var(v));

subplot(2,2,2)
plot(y,desviacionU);
hold on;
plot(y,desviacionV);
grid on;
xlabel("y [mm]");ylabel("Desviación de la fluctuación");title("");
legend('σu','σv');

% c)
ITu = desviacionU./V;
ITv = desviacionV./V;

subplot(2,2,3)
plot(y,ITu);
hold on;
plot(y,ITv);
grid on;
xlabel("y [mm]");ylabel("Intensidad de turbulencia");title("");
legend('ITu','ITv');


% d)
rho = 1.204;

txx = -rho*mean(uf.^2);
txy = -rho*mean(uf.*vf);
tyy = -rho*mean(vf.^2);

subplot(2,2,4)
hold on;
plot(y,txx);
plot(y,txy);
plot(y,tyy);
grid on;
xlabel("y [mm]");ylabel("[Pa]");title("Tensiones de Reynolds");
legend('Txx','Txy','Tyy');

axes( 'visible', 'off', 'title', 'Punto 2)' );






%3)

mu = 0.0000174;

dU_dY = diff(um) ./ diff(y);
dU_dY_extended = [dU_dY, 0];

Txy = mu*dU_dY_extended;

figure (3);
hold;
subplot(3, 1, 1);
plot(y,txy);
grid on;
xlabel("y [mm]");ylabel("[Pa]");title("Tensiones turbulentas");

subplot(3, 1, 2);
plot(y,Txy);
grid on;
xlabel("y [mm]");ylabel("[Pa]");title("Tensiones laminares");

subplot(3, 1, 3);
hold;
plot(y,dU_dY_extended);
grid on;
#TODO: PONERLE TITULO A ESTO
xlabel("y [mm]");ylabel("du/dy");title("");

axes( 'visible', 'off', 'title', 'Punto 3)' );


% se puede ver que las tensiones laminares son insignificantes frente
% a las turbulentas









%4)

function Ru = Autocorrelacion(p, u, uf)
  N=length(u);
  sumatoria = 0;
  for i = 1 : N-p
    j = uf(i,5)*uf(i+p,5);
    sumatoria = sumatoria + j;
  end
  Ru = (N/(N-1))*(1/(N-p))*sumatoria;
end

%La cantidad de pasos que grafico es arbitraria para visualizar bien los
%primeros

Ruvector = [];
for i = 1:100
    Ruvector(i) = Autocorrelacion(i, u, uf);
end

paso = 1:100;

figure(4);
subplot(2, 1, 1);
hold;
plot(paso, Ruvector);
grid on;
xlabel("pasos");ylabel("autocorrelacion");title("");

CoefAutocorrelacion = Ruvector./var(u(:,5));

subplot(2, 1, 2);
hold;
plot(paso, CoefAutocorrelacion);
grid on;
xlabel("pasos");ylabel("Coeficiente de autocorrelacion");title("");

axes( 'visible', 'off', 'title', 'Punto 4)' );
