%NxN complete Piles group.
%Class Soil-Fluid-Structure Interaction in National University of Córdoba
%of Dr. Federico Pinto.
%This code calculate the axial force in each pile considering SSI.
%Based in Randolph Theory and Elastic behaviour.
%Consider Rigid head and constant Shear modulus of Soil.
%Works only for NxN number of piles and constant spacing in each direction.
%Created date: 4th April 2020.
%Last review: 5th April 2020.
%Eng. Roberto Enrique Pinto Villegas.

clear all
format short g
clc

disp('Check the assumptions in the code!')

%Elasticity modulus of piles.
E=input('Elasticity modulus of piles [t/m^2] = ');

%Shear modulous of Sand.
Gz=input('Shear modulus of Soil [t/m^2] = ');

%Poisson modulus of soil.
nu=input('Poisson modulus of soil = ');

%Piles diameter.
R=input('Diameter of piles [m] = ')/2;

%Piles length.
L=input('Length of piles [m] = ');

%Axial Stiffness.
xi=log(2.5*50*(1-0.5));
lambda=E/Gz;
muL=sqrt(2/xi/lambda)*L/R;
k=R*Gz*((4/(1-nu))+2*pi/xi*(tanh(muL)/muL)*L/R)...
    /(1+(4/(1-nu)*(1/lambda/pi)*(tanh(muL)/muL)*L/R));

%Number of Piles.
nx=input('Number of piles in x = ');
ny=input('Number of piles in y = ');

%Piles spacing.
x=input('Separation of piles in x [m] = ');
y=input('Separation of piles in y [m] = ');

%Position of piles.
sx=-(nx*x-x)/2:x:(nx*x-x)/2;
sy=-(ny*y-y)/2:y:(ny*y-y)/2;

%Axial load of the group.
Pc=input('Input the axial load in [t] (+ load goes down) = ')

for i=1:nx;

    Sx(i,:)=sx;
    
end

for i=1:ny;

    Sy(:,i)=-sy;
    
end

%Auxiliar value.
Size=max(size(Sx,2)*nx,size(Sy,2)*ny);
D=zeros(Size);

for i=1:Size;

    for j=1:Size;

        d=sqrt((Sx(i)-Sx(j))^2+(Sy(i)-Sy(j))^2);

        D(i,j)=D(i,j)+d;

    end

end

Alpha=zeros(size(D));

for i=1:Size;

    for j=1:Size;

        if i==j;

            alpha=1;

            Alpha(i,j)=Alpha(i,j)+alpha;

        else

            alpha=sqrt(0.8/D(i,j)/(4/3));

            Alpha(i,j)=Alpha(i,j)+alpha;

        end

    end

end

%Uniform axial load for each pile. 
for i=1:Size;
    
    A(i)=1;
    
end

%Group stiffness.
disp('Group stiffness in t/m.')
Kc=A*k*inv(Alpha)*A'

%Elastic displacement of the group.
disp('Elastic displacement of the group in m.')
Uc=Pc/Kc

%Settlement relation.
disp('Settlement relation.')
Sr=Uc*Size*k/Pc

%Axial force for each pile.
disp('Axial force of each pile in t.')
P=k*inv(Alpha)*A'*Uc



