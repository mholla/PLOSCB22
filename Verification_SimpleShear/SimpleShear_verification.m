%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this program verifies the mechanical portion of the uel 
%
% Shuolun Wang, 2020 @ND 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

load('abaqus_verfication_mech.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mechanical verification 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameter 
mu = 100;


gamma = 0:0.01:1;
T12 = mu*gamma;
Tdiff = mu*gamma.^2;

l0 = 1e-3;
gamma_aba = dis_aba(:,2)/l0;

figure(1)
h1 = plot(gamma,T12/mu,'k-','linewidth',2); 
hold on
h2 = plot(gamma,Tdiff/mu,'k--','linewidth',2);
hold on 
h3 = plot(gamma_aba,T12_aba(:,2)/mu,'ko','markersize',8);
hold on 
plot(gamma_aba,Tdiff_aba(:,2)/mu,'ko','markersize',8)
hold on 

legend([h1 h2 h3],'Shear stress','Normal stress difference',...
          'U3D8 Abaqus','interpreter','latex','location','northwest')
legend('boxoff')

xlabel('$\gamma$[-]','interpreter','latex')
ylabel('Normalized stress [-]','interpreter','latex')


set(gca,'fontsize', 23)
