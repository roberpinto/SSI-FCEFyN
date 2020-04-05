%Pile in friction layered soil.
%Class Soil-Fluid-Structure Interaction in Universidad Nacional de Córdoba
%of Dr. Eng. Federico Pinto.
%This code calculate the axial force in the pile for each layer of soil
%considering SSI.
%Created date: th April 2020.
%Last review: 5th April 2020.
%Eng. Roberto Enrique Pinto Villegas.

clear all;
close all;
format short g;
clc;

%Elastic modulus of Pile.
E=input('Input the elastic modulus of Pile [N/mm^2] = ');

%Geometric of Pile.
R=input('Input the radius of Pile [m] = ');
A=pi*R^2;
I=pi*R^4/4;

%Axial load.
Ph=input('Input the axial load [kN] = ');

%Define properties for each layer of Soil.
n=input('Input the number of layer of Soil = ');

L=zeros(n,1);
kv=zeros(n,1);

for i=1:n

    disp('For the layer');disp(i);

    L(i)=input('Input the length [m] = ');

    kv(i)=input('Input the friction resistence [N/mm^2] = ');

    beta(i,1)=sqrt(kv(i)/(E*A));
    
end

kp=zeros(n,1);
kp(n)=input('Input the base resistence [Nm/mm^2] = ');

for i=n-1:-1:1
       
    W=(1+kp(i+1)/(E*A*beta(i+1)));
    
    X=exp(beta(i+1)*L(i+1));
    
    Y=(1-kp(i+1)/(E*A*beta(i+1)));
    
    Z=exp(-(beta(i+1)*L(i+1)));
    
    kp(i)=E*A*beta(i+1)*((W*X-Y*Z)/(Y*Z+W*X));
    
end

P=zeros(n+1,1);
P(1)=Ph;

for i=1:n
    
    W=(1+kp(i)/(E*A*beta(i)));
    
    X=exp(beta(i)*L(i));
    
    Y=(1-kp(i)/(E*A*beta(i)));
    
    Z=exp(-(beta(i)*L(i)));
    
    P(i+1)=P(i)*kp(i)/(E*A*beta(i))*2/(W*X-Y*Z);
    
end

disp('Axial force [kN] in each layer are: ');P

v=zeros(n+1,1);

W=(1+kp(1)/(E*A*beta(1)));
    
X=exp(beta(1)*L(1));
    
Y=(1-kp(1)/(E*A*beta(1)));
    
Z=exp(-(beta(1)*L(1)));

v(1)=P(1)/(E*A*beta(1))*(Y*Z+W*X)/(W*X-Y*Z);

for i=1:n
    
    v(i+1)=P(i+1)/kp(i);
    
end

disp('Displacement [mm] in each layer are: ');v

H=zeros(n+1,1);

for i=2:n+1;
    
    H(i)=H(i-1)-L(i-1);
    
end
    
figure;
plot(v,H)
title('Displacement vs Depth');
xlabel('Displacement [mm]');
ylabel('Depth [m]');
grid on;

figure;
plot(P,H)
title('Force vs Depth');
xlabel('Force [kN]');
ylabel('Depth [m]');
grid on;