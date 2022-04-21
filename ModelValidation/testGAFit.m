function [] = testGAFit()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this program calibrate FEM constitutive model via GA.
% 
% Shuolun Wang, May 2021 @ ND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

% load experimental data from kawasaki 
load('kawasaki_exp.mat')


% material parameter
%
% smooth parameter for neuron cell migration  
alpha = 0.1; 
alpha2 = 0.1;
alpha3 = 0.1; 
 
%smooth parameter for source term 
alpha_Gc = 0.05;
alpha_Gc2 = 0.05;
alpha_Gc3 = 0.05;

%smooth parameter for velocity 
alpha_v = 0.05;
alpha_v2 = 0.05; 
alpha_v3 = 0.05;


%smooth parameter for diffusivity 
alpha_d = 1.0;
alpha_d2 = 1.0; 
alpha_d3 = 1.0;


%smooth parameter for stiffness 
alpha_mu = 1.0;

%smooth parameter for growth parameter 
alpha_k = 20.0; 

%base value for source term (#/mm^3 day)
Gc = 5.0;
Gc2 = 5.0; 
Gc3 = 5.0;

%base value for velocity (mm/day)
vel = 0.1;
vel2 = 0.1;
vel3 = 0.1;


%diffusivity (mm^2/day)
diff = 0.5;
diff2 = 0.5;
diff3 = 0.5;


%cell migration threshold 
c0  = 1.0e-6;
c02 = 1.0e-6;
c03 = 1.0e-6;

%shape parameter for growth
k_s = 600000.0;
betak = 1.0;
k_c_normal = k_s*betak;
k_c_para = k_s/betak;
alpha_normal = 1.0; 
alpha_para = 1.0;


%mechanical parameter(Pa) 
lam_s = 930.0;
mu_s = 100.0;


%stiffness ratio 
betaE = 3.0; 
lam_c = lam_s*betaE;
mu_c = mu_s*betaE;


%delta function parameter for Gc as function of time
epsilon1 = 20.0;
epsilon2 = 20.0;
epsilon3 = 20.0; 



% known element parameter  
nInt = 8; nlSdv = 11; ngSdv = 11;



% plot the experimental data from kawasaki profiles
figure(1)

plot(E31E39Length_micro,E31E39Density1d_micro,'ko','linewidth',2)
hold on 
plot(E33E40Length_micro,E33E40Density1d_micro,'kv','linewidth',2)
hold on 

plot(E31P6Length_micro,E31P6Density1d_micro,'bo','linewidth',2)
hold on 
plot(E33P6Length_micro,E33P6Density1d_micro,'bv','linewidth',2)
hold on 
plot(E36P5Length_micro,E36P5Density1d_micro,'bs','linewidth',2)
hold on

plot(E31P16Length_micro,E31P16Density1d_micro,'ro','linewidth',2)
hold on 
plot(E34P16Length_micro,E34P16Density1d_micro,'rv','linewidth',2)
hold on 
plot(E37P16Length_micro,E37P16Density1d_micro,'rs','linewidth',2)
hold on

xlabel('Length ($\mu$m)','interpreter','latex')
ylabel('Cell density ($\#/\mu$m$^3$)','interpreter','latex')
set(gca,'fontsize', 18)





% set bounds of parameter
%        %Gc,      vel,      diff,        k_s       delta_v1 delta_v2  deltav3
ub = [0.00002,   1600,     30000.0,    250000.0,   200.0,     215.0    230.0]; 
lb = [0.00001,   1400,     25000.0,    150000.0,   195.0,     200.0    220.0];
% here I start the GA calibration 
% this is the fitness function 

%ftns = @(p) norm(stress_measured - FEM_constitutive(p(1),0.001*(1/p(1))));
ftns = @(p) objfunc(E31E39Length_micro,E31E39Density1d_micro,...
                    E33E40Length_micro,E33E40Density1d_micro,...
                    E31P6Length_micro,E31P6Density1d_micro,...
                    E33P6Length_micro,E33P6Density1d_micro,...
                    E36P5Length_micro,E36P5Density1d_micro,...
                    E31P16Length_micro,E31P16Density1d_micro,...
                    E34P16Length_micro,E34P16Density1d_micro,...
                    E37P16Length_micro,E37P16Density1d_micro,...
                    alpha,alpha2,alpha3,alpha_Gc,alpha_Gc2,alpha_Gc3,...
                    alpha_v,alpha_v2,alpha_v3,alpha_d,alpha_d2,alpha_d3,...
                    alpha_mu,alpha_k,p(1),p(1),p(1),p(2),p(2),p(2),...
                    p(3),p(3),p(3),c0,c02,c03,...
                    p(4),p(4),p(4),alpha_normal,alpha_para,...
                    lam_c,lam_s,mu_c,mu_s,... 
                    epsilon1,epsilon2,epsilon3,p(5),p(6),p(7),...
                    nlSdv,ngSdv,nInt); 

                
                
                
                
 





% number of parameter to be calibrated 
Parms = 7;

options = optimoptions(@ga,'MutationFcn',@mutationadaptfeasible);
options = optimoptions(options,'PlotFcn',{@gaplotbestf,@gaplotmaxconstr}, ...
    'Display','iter');
options.PopulationSize = 10;


%[B] = ga(ftns, Parms, [],[],[],[],[],[],[],[],options)
tic
[B] = ga(ftns, Parms, [],[],[],[],lb,ub,[],options)
toc

% evaluate the fitted results 
                     
[NT11_1,NT11_2,NT11_3,...
 NT12_1,NT12_2,NT12_3,...
 NT13_1,NT13_2,NT13_3,...
 disp_1,disp_2,disp_3] = FEM_constitutive(alpha,alpha2,alpha3,...
                alpha_Gc,alpha_Gc2,alpha_Gc3,...
                alpha_v,alpha_v2,alpha_v3,alpha_d,alpha_d2,alpha_d3,...
                alpha_mu,alpha_k,B(1),B(1),B(1),B(2),B(2),B(2),...
                B(3),B(3),B(3),c0,c02,c03,...
                B(4),B(4),B(4),alpha_normal,alpha_para,...
                lam_c,lam_s,mu_c,mu_s,... 
                epsilon1,epsilon2,epsilon3,B(5),B(6),B(7),...
                nlSdv,ngSdv,nInt);              
             
             
             
             
%plot the fitted results              
figure(1)             
l0 = 239.5; 
plot(linspace(0,l0+disp_1,30),NT11_1,'k-')
hold on 
plot(linspace(0,l0+disp_1,30),NT12_1,'k--')
hold on 

plot(linspace(0,l0+disp_2,30),NT11_2,'b-')
hold on 
plot(linspace(0,l0+disp_2,30),NT12_2,'b--')
hold on 
plot(linspace(0,l0+disp_2,30),NT13_2,'b:')
hold on 

plot(linspace(0,l0+disp_3,30),NT11_3,'r-')
hold on 
plot(linspace(0,l0+disp_3,30),NT12_3,'r--')
hold on
plot(linspace(0,l0+disp_3,30),NT13_3,'r:')
hold on 





end


function [NT11_1,NT11_2,NT11_3,...
          NT12_1,NT12_2,NT12_3,...
          NT13_1,NT13_2,NT13_3,...
          disp_1,disp_2,disp_3] = FEM_constitutive(alpha,alpha2,alpha3,...
                alpha_Gc,alpha_Gc2,alpha_Gc3,...
                alpha_v,alpha_v2,alpha_v3,alpha_d,alpha_d2,alpha_d3,...
                alpha_mu,alpha_k,Gc,Gc2,Gc3,vel,vel2,vel3,...
                diff,diff2,diff3,c0,c02,c03,...
                k_s,k_c_normal,k_c_para,alpha_normal,alpha_para,...
                lam_c,lam_s,mu_c,mu_s,... 
                epsilon1,epsilon2,epsilon3,delta_v1,delta_v2,delta_v3,...
                nlSdv,ngSdv,nInt)           

% %%%%%%%%function counter%%%%%%%%%%%%
% persistent  counter
% if isempty( counter )
%     counter=0; %Initializing counter
% end
% counter = counter + 1   
% 
% jobname = ['parameter' num2str(counter) '.txt'];
% fileID = fopen(jobname,'w');
% fprintf(fileID,'%f,%f,%f,%f',Gc,vel,diff,k_s);
% fclose(fileID);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% write parameter to the input file
fileID = fopen('parameter.inp','w');
fprintf(fileID,'%f,%f,%f,%f,%f,%f,%f,%f,\r\n',alpha,alpha2,alpha3,...
    alpha_Gc,alpha_Gc2,alpha_Gc3,alpha_v,alpha_v2);
fprintf(fileID,'%f,%f,%f,%f,%f,%f,%f,%f,\r\n',alpha_v3,alpha_d,...
    alpha_d2,alpha_d3,alpha_mu,alpha_k,Gc,Gc2);
fprintf(fileID,'%f,%f,%f,%f,%f,%f,%f,%f,\r\n',Gc3,vel,vel2,vel3,...
    diff,diff2,diff3,c0);
fprintf(fileID,'%f,%f,%f,%f,%f,%f,%f,%f,\r\n',c02,c03,k_s,...
    k_c_normal,k_c_para,alpha_normal,alpha_para,lam_c);
fprintf(fileID,'%f,%f,%f,%f,%f,%f,%f,%f,\r\n',lam_s,mu_c,mu_s,...
    epsilon1,epsilon2,epsilon3,delta_v1,delta_v2);
fprintf(fileID,'%f,%d,%d,%d\r\n',delta_v3,nlSdv,ngSdv,nInt);
fclose(fileID);

% run abaqus job 
! abaqus double job=j1 input=bar_exp.inp user=uel_MultiCell_deformation_v4.for ask_delete=OFF

        % wait until job is done
        a = 0;
        while a ==0
            data = fileread('j1.log');
            k = strfind(data,'COMPLETED');
            
            if k>600
                a=1;
            else
                a=0;
            end
            pause(1);
        end


% Get the odb results  
! abaqus cae noGUI=GetData.py
NT11_1 = load('NT11_1.txt');
NT11_2 = load('NT11_2.txt');
NT11_3 = load('NT11_3.txt');

NT12_1 = load('NT12_1.txt');
NT12_2 = load('NT12_2.txt');
NT12_3 = load('NT12_3.txt');

NT13_1 = load('NT13_1.txt');
NT13_2 = load('NT13_2.txt');
NT13_3 = load('NT13_3.txt');

disp_1 = load('disp_1.txt');
disp_2 = load('disp_2.txt');
disp_3 = load('disp_3.txt');


end 




function [err] = objfunc(E31E39Length_micro,E31E39Density1d_micro,...
                    E33E40Length_micro,E33E40Density1d_micro,...
                    E31P6Length_micro,E31P6Density1d_micro,...
                    E33P6Length_micro,E33P6Density1d_micro,...
                    E36P5Length_micro,E36P5Density1d_micro,...
                    E31P16Length_micro,E31P16Density1d_micro,...
                    E34P16Length_micro,E34P16Density1d_micro,...
                    E37P16Length_micro,E37P16Density1d_micro,...
                    alpha,alpha2,alpha3,alpha_Gc,alpha_Gc2,alpha_Gc3,...
                    alpha_v,alpha_v2,alpha_v3,alpha_d,alpha_d2,alpha_d3,...
                    alpha_mu,alpha_k,Gc,Gc2,Gc3,vel,vel2,vel3,...
                    diff,diff2,diff3,c0,c02,c03,...
                    k_s,k_c_normal,k_c_para,alpha_normal,alpha_para,...
                    lam_c,lam_s,mu_c,mu_s,... 
                    epsilon1,epsilon2,epsilon3,delta_v1,delta_v2,delta_v3,...
                    nlSdv,ngSdv,nInt)
             
             
            % call constutive 
          [NT11_1,NT11_2,NT11_3,...
          NT12_1,NT12_2,NT12_3,...
          NT13_1,NT13_2,NT13_3,...
          disp_1,disp_2,disp_3] = FEM_constitutive(alpha,alpha2,alpha3,...
                alpha_Gc,alpha_Gc2,alpha_Gc3,...
                alpha_v,alpha_v2,alpha_v3,alpha_d,alpha_d2,alpha_d3,...
                alpha_mu,alpha_k,Gc,Gc2,Gc3,vel,vel2,vel3,...
                diff,diff2,diff3,c0,c02,c03,...
                k_s,k_c_normal,k_c_para,alpha_normal,alpha_para,...
                lam_c,lam_s,mu_c,mu_s,... 
                epsilon1,epsilon2,epsilon3,delta_v1,delta_v2,delta_v3,...
                nlSdv,ngSdv,nInt);
    

            
             
             
             % calculate the error 
             
             err1_E39 = norm(NT11_1 - E31E39Density1d_micro)/norm(E31E39Density1d_micro);
             err2_E39 = norm(NT12_1 - E33E40Density1d_micro)/norm(E33E40Density1d_micro); 
             
             err1_P6 = norm(NT11_2 - E31P6Density1d_micro)/norm(E31P6Density1d_micro); 
             err2_P6 = norm(NT12_2 - E33P6Density1d_micro)/norm(E33P6Density1d_micro); 
             err3_P6 = norm(NT13_2 - E36P5Density1d_micro)/norm(E36P5Density1d_micro); 
             
             err1_P16 = norm(NT11_3 - E31P16Density1d_micro)/norm(E31P16Density1d_micro);
             err2_P16 = norm(NT12_3 - E34P16Density1d_micro)/norm(E34P16Density1d_micro); 
             err3_P16 = norm(NT13_3 - E37P16Density1d_micro)/norm(E37P16Density1d_micro); 
             
             l_E39_exp = 1390; % micro 
             l_P6_exp  = 1997; % micro 
             l_P16_exp = 2390.3; % micro 
             
             l_E31_sim = 239.5; %micro 
             
             err_disp_E39 = norm(l_E31_sim + disp_1 - l_E39_exp)/norm(l_E39_exp);
             err_disp_P6 = norm(l_E31_sim + disp_2 - l_P6_exp)/norm(l_P6_exp); 
             err_disp_P16 = norm(l_E31_sim + disp_3 - l_P16_exp)/norm(l_P16_exp);
             
             %weight_c = 11/8; 
             %weight_dis = 11/3; 
             % time weights
             weight_c = 1/10; 
             weight_dis = 1/4;

%              err = weight_c*(err1_E39 + err2_E39 + ...
%                    err1_P6 + err2_P6 + err3_P6 + ...
%                    err1_P16 + err2_P16 + err3_P16) + ...
%                    weight_dis*(err_disp_E39 + err_disp_P6 + err_disp_P16); 
              % error for time validation
              err = weight_c*(err1_E39 + err2_E39 + ...
                              err1_P6 + err2_P6 + err3_P6) + ...
                    weight_dis*(err_disp_E39 + err_disp_P6);              

                         
             
end 
