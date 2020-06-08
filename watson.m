%clear;
%clc;
%%%% Code Key_Final
syms ne ni Vp
syms Te positive
%%%%%%%% Known Parameters
rp= 7E-9;                           
S= 4*pi*rp^2;                        
me= 9.10938E-31;                
mi= 6.6335209E-26;                  
kb= 1.38064853E-23;                  
% Te= 2;                            
Ti= 300;                       
Tg= 300;                                 
eo= 8.8541878176E-12;       
e= 1.6021765E-19;                     
p= 14;                             
np= 5E12;                              
d= 10.16E-2;                          
l= d/2;                                 
RA= (18/2)*10^-2;                        
Ae= 2*pi*RA*d;                        
V= RA^2*pi*d;                           
N= p/(Ti*kb);                                        
% Interpolating Data
%%% Argon Excitaiton Linear Interpolation to calculate exitation frequency 
data= load('Ar-Excitation.txt');
E1= data(:,1);                            % Column 1- Energy- Ar [V]
Q1= data(:,2);                            % Column 2- Cross Section- Ar [m^2]
E1i= 0:1.5:60; %0:0.01:60; 
Q1i= interp1(E1,Q1,E1i);
fe1= (2/(pi^0.5))*((1/(kb*Te*11604.52500617))^1.5)*(e*E1i).^0.5.*exp((-e*E1i)./(kb*Te*11604.52500617));
fent1= trapz(E1i,fe1);
fe1n= fe1/fent1;
vth1= (e*E1i).^0.5.*((2/me)^0.5);
vexx= N*Q1i.*vth1.*fe1n;
vex= trapz(E1i,vexx);                     % Excitation Collision Frequency of Ar [1/s]
data2= load('Ar-Ionization.txt');
E2= data2(:,1);                           % Column 1- Energg- Ar [V]
Q2= data2(:,2);                           % Column 2- Cross Section- Ar [m^2]
E2i= 0:1.5:60;  %0:0.01:60;
Q2i= interp1(E2,Q2,E2i);
fe2= (2/(pi^0.5))*((1/(kb*Te*11604.52500617))^1.5)*(e*E2i).^0.5.*exp((-e*E2i)./(kb*Te*11604.52500617));
fent2= trapz(E2i,fe2);
fe2n= fe2/fent2;
vth2= (e*E2i).^0.5.*((2/me)^0.5);
viion= N*Q2i.*vth2.*fe2n;
vion= trapz(E2i,viion);                   % Ionization Collision Frequnecy of Ar [1/s]

data3= load('Ar-Elastic.txt');
E3= data3(:,1);                           % Column 1- Energg- Ar [V]
Q3= data3(:,2);                           % Column 2- Cross Section- Ar [m^2]
E3i= 0:1.5:60; %0:0.01:60;
Q3i= interp1(E3,Q3,E3i);
fe3= (2/(pi^0.5))*((1/(kb*Te*11604.52500617))^1.5)*(e*E3i).^0.5.*exp((-e*E3i)./(kb*Te*11604.52500617));
fent3= trapz(E3i,fe3);
fe3n= fe3/fent3;
vth3= (e*E3i).^0.5.*((2/me)^0.5);
vmm= N*Q3i.*vth3.*fe3n;
vm= trapz(E3i,vmm);           
ve= (((kb*Te*11604.52500617)/(2*pi*me))^0.5)*S*ne*exp((Vp*e)/(kb*Te*11604.52500617)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vi= (((kb*Ti)/(2*pi*mi))^0.5)*S*ni*(1-((Vp*e)/(kb*Ti)));
eqn1= ve-vi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Equation (2)     
k= Vp*(4*pi*eo*rp)/(e);               
eqn2= ne-ni-k*np;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Equation (3)
Eex= 11.5*e;                           
Ei= 15.8*e;                              
RHS1= ne*(vion*Ei+vex*Eex)*V;           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RHS2
phi= 40.9*10^-20;                                 
% lambi= kb*Ti/(p*phi);                  
vbar= ((8*kb*Ti)/(pi*mi))^0.5;          
mubar= N*phi*vbar;                    
Di= ((kb*Ti)/(mi*mubar));            
Da= (Di)*((Te*11604.52500617)/Ti);                    
%Vc= abs((log((abs(ni/ne))*(2*pi*(me/mi))^0.5)));     
Vsh= 1500;                   
RHS2= ni*(Da/l)*(Ae)*(Vsh/2)*e;         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RHS3- prof. LM Expression
RHS3= ((2*me)/mi)*vm*1.5*kb*((Te*11604.52500617)-Tg)*ne*V;% Third term on the RHS of the power equation [W]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RHS4
h1=(np*V)*(pi*rp^2)*(2^0.5)*(1/(me^0.5))*(2/(pi^0.5))*(1/(kb*Te*11604.52500617))^1.5;
syms E
h2=int((1+((Vp*e)/E))*(E)^2*exp(-(E)/(kb*Te*11604.52500617)), E, [-Vp*e, Inf]);
RHS4= h1*ne*(h2);                          % Fourth term on the RHS of the power equation [W]
eqn4= (ne*vion*V)-(np*vi*V)-(ni*(Da/(l))* Ae);
i = 1;
 for PP=100:20:140                         % Input Power [W]
eqn3= PP-RHS1-RHS2-RHS3-RHS4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 Equations, 4 Variables
result = vpasolve(eqn1,eqn2,eqn3,eqn4,ne,ni,Vp,Te,[10^14 10^18; 10^14 10^18;-7 1;1 5]);
LM1(1,i)= result.ne
LM2(1,i)= result.ni
LM3(1,i)= result.Vp
LM4(1,i)= (LM1(1,i)/LM2(1,i));
LM5(1,i)= result.Te
disp(i)
i = i+1;
end
 figure(1)
 plot(100:20:140,LM1,'g');
 xlabel('Power [W]', 'FontSize', 30);
 ylabel('Electron Density [m^-3]', 'FontSize', 30);
 title('Electron Density vs. Input Power', 'FontSize', 30);
 figure (2)
 plot(100:20:140,LM2,'r--*');
 xlabel('Power [W]', 'FontSize', 30);
 ylabel('Ion Density [m^-3]', 'FontSize', 30);
 title('Ion Density vs. Input Power [m^-3]', 'FontSize', 30);
 figure(3)
 plot(100:20:140,LM3,'b');
 xlabel('Power [W]', 'FontSize', 30);
 ylabel('Particle Potential [eV]', 'FontSize', 30)
 title('Particle Potential vs. Input Power', 'FontSize', 30);
 figure(4)
 plot(100:20:140,LM4,'r--*');
 xlabel('Power [W]', 'FontSize', 30);
 ylabel('Electron/Ion Density Ration [1/s]', 'FontSize', 30)
 title('Electron/Ion Density Ratio vs. Input Power', 'FontSize', 30);
 figure(5)
 plot(100:20:140,LM5,'r--*');
 xlabel('Power [W]', 'FontSize', 30);
 ylabel('Electron Temperature [K]', 'FontSize', 30)
 title('Electron Temperature [K] vs. Input Power', 'FontSize', 30);