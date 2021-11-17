%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this program shows the model calibration against the experiments 
%
% Shuolun Wang, 2021 @ND
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

% python script to obtain the abaqus results 
%
! abaqus cae noGUI=GetData.py


load('kawasaki_exp.mat')
load('NT11_1.txt')
load('NT12_1.txt')
load('NT13_1.txt')

load('NT11_2.txt')
load('NT12_2.txt')
load('NT13_2.txt')

load('NT11_3.txt')
load('NT12_3.txt')
load('NT13_3.txt')

load('disp_1.txt')
load('disp_2.txt')
load('disp_3.txt')




l0 = 239.5; 
%prep@E39 
figure(1)
h1 = plot(E31E39Length_micro,E31E39Density1d_micro,'ko','markersize',8);
hold on
h2 = plot(linspace(0,l0+disp_1,length(NT11_1)),NT11_1,'k-');
hold on 
h3 = plot(E33E40Length_micro,E33E40Density1d_micro,'kv','markersize',8);
hold on 
h4 = plot(linspace(0,l0+disp_1,length(NT12_1)),NT12_1,'k--');
hold on
title('Prepared at E39-40')
xlabel('Length [$\mu$m]','interpreter','latex')
ylabel('density [$\#$/$\mu$m$^3$]','interpreter','latex')
set(gca,'fontsize', 18)
legend([h1 h2 h3 h4],'IUE at 31 exp.','IUE at 31 sim.',...
    'IUE at 33-34 exp.','IUE at 33-34 sim.','location','best')
legend('boxoff')


%%prep@P6
figure(2)
h1 = plot(E31P6Length_micro,E31P6Density1d_micro,'bo','markersize',8);
hold on
h2 = plot(linspace(0,l0+disp_2,length(NT11_2)),NT11_2,'b-');
hold on 
h3 = plot(E33P6Length_micro,E33P6Density1d_micro,'bv','markersize',8);
hold on 
h4 = plot(linspace(0,l0+disp_2,length(NT12_2)),NT12_2,'b--');
hold on 
h5 = plot(E36P5Length_micro,E36P5Density1d_micro,'bs','markersize',8);
hold on 
h6 = plot(linspace(0,l0+disp_2,length(NT13_2)),NT13_2,'b:');
title('Prepared at P5-6')
xlabel('Length [$\mu$m]','interpreter','latex')
ylabel('density [$\#$/$\mu$m$^3$]','interpreter','latex')
set(gca,'fontsize', 18)
legend([h1 h2 h3 h4 h5 h6],'IUE at 31 exp.','IUE at 31 sim.',...
                           'IUE at 33-34 exp.','IUE at 33-34 sim.',...
                           'IUE at 36-37 exp.','IUE at 36-37 sim.',...
                           'location','best')
legend('boxoff')


%%prep@p16
figure(3)
h1 = plot(E31P16Length_micro,E31P16Density1d_micro,'ro','markersize',8);
hold on 
h2 = plot(linspace(0,l0+disp_3,length(NT11_3)),NT11_3,'r-');
hold on 
h3 = plot(E34P16Length_micro,E34P16Density1d_micro,'rv','markersize',8);
hold on 
h4 = plot(linspace(0,l0+disp_3,length(NT12_3)),NT12_3,'r--');
hold on 
h5 = plot(E37P16Length_micro,E37P16Density1d_micro,'rs','markersize',8);
hold on 
h6 = plot(linspace(0,l0+disp_3,length(NT13_3)),NT13_3,'r:'); 
hold on 
title('Prepared at P16')
xlabel('Length [$\mu$m]','interpreter','latex')
ylabel('density [$\#$/$\mu$m$^3$]','interpreter','latex')
set(gca,'fontsize', 18)
legend([h1 h2 h3 h4 h5 h6],'IUE at 31 exp.','IUE at 31 sim.',...
                           'IUE at 33-34 exp.','IUE at 33-34 sim.',...
                           'IUE at 36-37 exp.','IUE at 36-37 sim.',...
                           'location','best')
legend('boxoff')




%IUE31
figure(4)
a1 = plot(E31E39Length_micro,E31E39Density1d_micro,'ko','markersize',8);
hold on
a2 = plot(linspace(0,l0+disp_1,length(NT11_1)),NT11_1,'k-');
hold on 
a3 = plot(E31P6Length_micro,E31P6Density1d_micro,'bo','markersize',8);
hold on
a4 = plot(linspace(0,l0+disp_2,length(NT11_2)),NT11_2,'b-');
hold on 
a5 = plot(E31P16Length_micro,E31P16Density1d_micro,'ro','markersize',8);
hold on 
a6 = plot(linspace(0,l0+disp_3,length(NT11_3)),NT11_3,'r-');
title('IUE at E31')
xlabel('Length [$\mu$m]','interpreter','latex')
ylabel('density [$\#$/$\mu$m$^3$]','interpreter','latex')
set(gca,'fontsize', 18)
legend([a1 a2 a3 a4 a5 a6],'Prep at E39-40 exp.','Prep at E39-40 sim.',...
                           'Prep at P5-6 exp.','Prep at P5-6 sim.',...
                           'Prep at P16 exp.','Prep at P16 sim.',...
                           'location','best')
legend('boxoff')


%IUE34
figure(5)
b1 = plot(E33E40Length_micro,E33E40Density1d_micro,'kv','markersize',8);
hold on 
b2 = plot(linspace(0,l0+disp_1,length(NT12_1)),NT12_1,'k--');
hold on
b3 = plot(E33P6Length_micro,E33P6Density1d_micro,'bv','markersize',8);
hold on 
b4 = plot(linspace(0,l0+disp_2,length(NT12_2)),NT12_2,'b--');
hold on 
b5 = plot(E34P16Length_micro,E34P16Density1d_micro,'rv','markersize',8);
hold on 
b6 = plot(linspace(0,l0+disp_3,length(NT12_3)),NT12_3,'r--');
hold on 
title('IUE at E33-34')
xlabel('Length [$\mu$m]','interpreter','latex')
ylabel('density [$\#$/$\mu$m$^3$]','interpreter','latex')
set(gca,'fontsize', 18)
legend([b1 b2 b3 b4 b5 b6],'Prep at E39-40 exp.','Prep at E39-40 sim.',...
                           'Prep at P5-6 exp.','Prep at P5-6 sim.',...
                           'Prep at P16 exp.','Prep at P16 sim.',...
                           'location','best')
legend('boxoff')


 %IUE37
 
figure(6)
c1 = plot(E36P5Length_micro,E36P5Density1d_micro,'bs','markersize',8);
hold on 
c2 = plot(linspace(0,l0+disp_2,length(NT13_2)),NT13_2,'b:');
hold on
c3 = plot(E37P16Length_micro,E37P16Density1d_micro,'rs','markersize',8);
hold on 
c4 = plot(linspace(0,l0+disp_3,length(NT13_3)),NT13_3,'r:'); 
title('IUE at E36-37')
xlabel('Length [$\mu$m]','interpreter','latex')
ylabel('density [$\#$/$\mu$m$^3$]','interpreter','latex')
set(gca,'fontsize', 18)
legend([c1 c2 c3 c4],'Prep at P5-6 exp.','Prep at P5-6 sim.',...
                     'Prep at P16 exp.','Prep at P16 sim.',...
                     'location','best')
legend('boxoff')


%total 
figure(7)
h1 = plot(E31E39Length_micro,E31E39Density1d_micro,'ko','markersize',8);
hold on
h2 = plot(linspace(0,l0+disp_1,length(NT11_1)),NT11_1,'k-');
hold on 
h3 = plot(E33E40Length_micro,E33E40Density1d_micro,'kv','markersize',8);
hold on 
h4 = plot(linspace(0,l0+disp_1,length(NT12_1)),NT12_1,'k--');
hold on

h1 = plot(E31P6Length_micro,E31P6Density1d_micro,'bo','markersize',8);
hold on
h2 = plot(linspace(0,l0+disp_2,length(NT11_2)),NT11_2,'b-');
hold on 
h3 = plot(E33P6Length_micro,E33P6Density1d_micro,'bv','markersize',8);
hold on 
h4 = plot(linspace(0,l0+disp_2,length(NT12_2)),NT12_2,'b--');
hold on 
h5 = plot(E36P5Length_micro,E36P5Density1d_micro,'bs','markersize',8);
hold on 
h6 = plot(linspace(0,l0+disp_2,length(NT13_2)),NT13_2,'b:');
hold on 

h1 = plot(E31P16Length_micro,E31P16Density1d_micro,'ro','markersize',8);
hold on 
h2 = plot(linspace(0,l0+disp_3,length(NT11_3)),NT11_3,'r-');
hold on 
h3 = plot(E34P16Length_micro,E34P16Density1d_micro,'rv','markersize',8);
hold on 
h4 = plot(linspace(0,l0+disp_3,length(NT12_3)),NT12_3,'r--');
hold on 
h5 = plot(E37P16Length_micro,E37P16Density1d_micro,'rs','markersize',8);
hold on 
h6 = plot(linspace(0,l0+disp_3,length(NT13_3)),NT13_3,'r:'); 
hold on 

title('All in one')
xlabel('Length [$\mu$m]','interpreter','latex')
ylabel('density [$\#$/$\mu$m$^3$]','interpreter','latex')
set(gca,'fontsize', 18)
%legend([h1 h2 h3 h4 h5 h6],'E31_exp','E31_sim',...
%                           'E34_exp','E34_sim',...
%                           'E37_exp','E37_sim',...
%                           'interpreter','latex','location','best')
%legend('boxoff')











