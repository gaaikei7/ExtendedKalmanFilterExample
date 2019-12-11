clc;
clear all;
close all;


dt = 0.1;
t = 0:dt:17;
x = 1000;
y = 1000;

V = 100000/3600; % km/h a m/s 2 km/h
L = 2.92; % m
theta = 90 - atand(x/y);
Phi = 65; % grados

xpos = zeros(1, length(t));
ypos = zeros(1, length(t));
orientacion = zeros(1,length(t));

xpos(1) = x;
ypos(1) = y;

orientacion(1) = theta;


varianza = 9;
gps = zeros(2, length(t));
gps(1,1) = xpos(1) + sqrt(varianza)*randn;
gps(2,1) = ypos(1) + sqrt(varianza)*randn;
%{
gps(1,1) = awgn(xpos(1),db);
gps(2,1) = awgn(ypos(1),db);
%}


for n = 2:length(t)
    xpos(n) = xpos(n-1) + (  dt * V * cosd( 90 - orientacion(n-1) )  );
    ypos(n) = ypos(n-1) + (  dt * V * sind( 90 - orientacion(n-1) )  );
    gps(1,n) = xpos(n) + sqrt(varianza)*randn;
    gps(2,n) = ypos(n) + sqrt(varianza)*randn;
    orientacion(n) = orientacion(n-1) + (dt * V * tand(Phi) ) / L;
end

figure
plot(xpos,ypos,'b')
hold
plot(gps(1,:),gps(2,:),	'-o')


varianza = 0.01;
Vel(1:length(xpos)) = V + sqrt(varianza)*randn(1,length(xpos));

Xvec = [gps ; Vel];

Estados = [xpos; ypos; orientacion; V*ones(1,length(t))];

xhatpriori = [x; y; theta; V];
Ppriori = 10 * eye(4);
I = eye(4);

xhatArray = xhatpriori;

H = [1 0 0 0;0 1 0 0; 0 0 0 1];
R = 10*[9 0 0; 0 9 0; 0 0 0.01];
Q = 0.01*[0.01 0 0 0; 0 0.01 0 0; 0 0 0.01 0; 0 0 0 0.01];


%   Extended Kalman Filter

for i = 1:length(t)-1
F = [1 0 0 cosd(xhatpriori(3));0 1 0 sind(xhatpriori(3));0 0 1 tand(Phi)/L; 0 0 0 1];

%Prediction stage
%{
    xpos(n) = xpos(n-1) + (  dt * V * cosd( 90 - orientacion(n-1) )  );
    ypos(n) = ypos(n-1) + (  dt * V * sind( 90 - orientacion(n-1) )  );
    gps(1,n) = xpos(n) + sqrt(varianza)*randn;
    gps(2,n) = ypos(n) + sqrt(varianza)*randn;
    orientacion(n) = orientacion(n-1) + (dt * V * tand(Phi) ) / L;
%}
xhatpriori(1) = xhatpriori(1) + (  dt * xhatpriori(4) * cosd( 90 - xhatpriori(3) )  );
xhatpriori(2) = xhatpriori(2) + (  dt * xhatpriori(4) * sind( 90 - xhatpriori(3) )  );
xhatpriori(3) = xhatpriori(3) + (  dt * xhatpriori(4) * tand(Phi) ) / L;
Ppriori = F*Ppriori*F' + Q;

%Update
K = Ppriori*H'/(H*Ppriori*H'+R);
xhat = xhatpriori + K*(Xvec(:,i+1)-H*xhatpriori);
P = (I - K*H)*Ppriori;
xhatArray = [xhatArray xhat];

xhatpriori = xhat;
end

xArray = [xpos;ypos;orientacion;Vel];

xArray2 = [xpos; ypos];


Error = sqrt(sum((xhatArray - xArray).^2, 1));
ErrorGPS = sqrt(sum((gps - xArray2).^2, 1));




plot(xhatArray(1,:),xhatArray(2,:),'m')
legend('Real', 'Mediciones', 'Filtro')



figure(2)
plot(Error);
title('Error Filtro');

figure(3)
plot(ErrorGPS);
title('Error GPS');




%{
x = [0 0.2 0.4 0.8 1.2 1.6 2.0];
y = [0 0.155 0.240 0.328 0.450 0.582 0.692];
p = polyfit(x, y, 1);
v = polyval(p, x);
figure()
plot(x,y,'x','MarkerEdgeColor','black')
hold on
plot(x, v)
hold off
grid on;
xlabel('Protein standard concentration (µg/µl)');


figure()
C = (x(:).^(1:4))\y(:);
xint = linspace(0,2,100)';
plot(x,y,'o',xint,(xint.^(1:4))*C,'r-')

y = C(1)*x + C(2)*x.^2 + C(3)*x.^3 + C(4)*x.^4
spl = spline(x,y);
figure()
plot(x,y,'o',xint,ppval(spl,xint),'r-')
%}
