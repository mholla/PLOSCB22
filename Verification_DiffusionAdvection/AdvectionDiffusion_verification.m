%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% One dimensional finite element implementation for transient, nonlinear
% advection-diffusion equation for multiple neuron cell migration (3 cells)
%
% Shuolun Wang 2021 @ ND
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [T] = AdvectionDiffusion_Multicell_v3_verification()
clear % clear the workspace variables
clc % clear the command window

% load abaqus data
%
load('abaqus_verfication_cell.mat')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Problem data
%
% total time
FinalTime = 500.0;
%
% time increment
dtime = 10.0;
%
% time increment for plotting the temperature
dtOut = 50.0;
%
% length of ``bar''
%
Length = 1.0;
%
% cross sectional area (possibly a function of position)
%
Area = @(x) 1.0;
%
%
% smooth parameter 
% cell1 
alpha1 = 0.1;
alpha_Gc1 = 100;
alpha_v1 = 100;
alpha_d1 = 100;
% cell2 
alpha2 = 0.1;
alpha_Gc2 = 100;
alpha_v2 = 100;
alpha_d2 = 100;
% cell3 
alpha3 = 0.1;
alpha_Gc3 = 100;
alpha_v3 = 100;
alpha_d3 = 100;



% cell division rate

% 
epsilon1 = 20;
epsilon2 = 20;
epsilon3 = 20;


offset1 = 100; 
offset2 = 200; 
offset3 = 300; 




% cell1 
Gc1 = 25.0;
source1 = @(x,t) ((epsilon1^2)/((t-offset1)^2 + epsilon1^2))*...
                 Gc1*(exp(alpha_Gc1*(-(x - 0.2))))/(1 + exp(alpha_Gc1*(-(x - 0.2))));
% cell2
Gc2 = 25.0;
source2 = @(x,t) ((epsilon2^2)/((t-offset2)^2 + epsilon2^2))*...
                 Gc2*(exp(alpha_Gc2*(-(x - 0.2))))/(1 + exp(alpha_Gc2*(-(x - 0.2))));
% cell3
Gc3 = 25.0;
source3 = @(x,t) ((epsilon3^2)/((t-offset3)^2 + epsilon3^2))*...
                 Gc3*(exp(alpha_Gc3*(-(x - 0.2))))/(1 + exp(alpha_Gc3*(-(x - 0.2))));



% Cell migration speed 
%
% cell1
offset = -0.005;

v1 = 0.01;
velocity1 = @(x,t) offset + v1*(exp(alpha_v1*(-(x - 0.7))))/(1 + exp(alpha_v1*(-(x - 0.7))));
%velocity1 = @(x,t) 0.0;

% cell2
v2 = 0.01;
velocity2 = @(x,t) offset + v2*(exp(alpha_v2*(-(x - 0.8))))/(1 + exp(alpha_v2*(-(x - 0.8))));

% cell3
v3 = 0.01;
velocity3 = @(x,t) offset + v3*(exp(alpha_v3*(-(x - 0.9))))/(1 + exp(alpha_v3*(-(x - 0.9))));




% Diffusivity (possibly a function of position and temperature)
%
base1 = 0.0005;
base2 = 0.0005;
base3 = 0.0005;

dd1 = 0.005;
dd2 = 0.005; 
dd3 = 0.005; 

%D1 = @(x,t) base1 + dd1*(1 - (exp(alpha_d1*(-(x - 0.95))))/(1 + exp(alpha_d1*(-(x - 0.95))))); 
%D2 = @(x,t) base2 + dd2*(1 - (exp(alpha_d2*(-(x - 0.95))))/(1 + exp(alpha_d2*(-(x - 0.95))))); 
%D3 = @(x,t) base3 + dd3*(1 - (exp(alpha_d3*(-(x - 0.95))))/(1 + exp(alpha_d3*(-(x - 0.95))))); 

D1 = @(x,t) base1;
D2 = @(x,t) base2;
D3 = @(x,t) base3;



% Cell migration threshold
c01 = 100;
c02 = 100;
c03 = 100;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create the mesh you want
%
% number of elements
%
nElem = 50;
%
% nodes per element
%
nNode = 3;
%
% number of integration points per element
%
nInt = 4;
%
% the total number of nodes in this mesh
%
nNodes = nElem*nNode - (nElem - 1);
%
% generate the mesh
%
coord = zeros(nNodes,1);
for i=2:nNodes
    coord(i) = coord(i-1) + Length/(nNodes-1);
end
%
% here is the element to nodal connectivity in the form
%  node on left  = connect(elem,1)
%  next node moving right = connect(elem,2)
%  ....
%  node on the right = connect(elem,nNode)
%
connect = zeros(nElem,nNode);
for i=1:nNode
    connect(1,i) = i;
end
for i=2:nElem
    for j=1:nNode
        connect(i,j) = connect(i-1,j) + nNode - 1;
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intial condition on the temperature field
%
T = zeros(3*nNodes,1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Boundary conditions
%
% flux boundary condition
%
flux = 0.0;      % magnitude AND direction of flux BC
elemFlux = nElem; % element with the flux BC
nodeFlux = nNode; % local node number to apply the flux BC
%
% temperature boundary condition (possible function of time)
%
%Tbc = @(t) 100.0*(1.0 - cos(-t/30));
%Tbc = @(t) 300+50*sin(t/30);
Tbc = @(t) 100;

nodeT = 1; % global node number to apply the temperature BC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin the loop over time
%
timeOut = 0.0;
for time=0:dtime:FinalTime

    
% obtain the temperature field from the last converged step
%
Told = T;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin the nonlinear iteration loop
%
iter = 0;
while(iter<=100)
        
    iter = iter + 1;
    
    K = zeros(3*nNodes,3*nNodes);
    R = zeros(3*nNodes,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop over elements
    %
    for elem=1:nElem

        % get the coordinates for the nodes on this element
        %
        nodeCoords = zeros(nNode,1);
        for i=1:nNode
            nodeCoords(i,1) = coord(connect(elem,i));
        end


        % get the nodal temperatures for this element
        %
        Te1 = zeros(nNode,1);
        Te2 = zeros(nNode,1);
        Te3 = zeros(nNode,1);
  
        
        Te1Old = zeros(nNode,1);
        Te2Old = zeros(nNode,1);
        Te3Old = zeros(nNode,1);
        
        
        for i=1:nNode
            
            row = 3*(connect(elem,i)-1)+1;
            Te1(i,1) = T(row,1);
            Te2(i,1) = T(row+1,1);
            Te3(i,1) = T(row+2,1);
            Te1Old(i,1) = Told(row,1);
            Te2Old(i,1) = Told(row+1,1);
            Te3Old(i,1) = Told(row+2,1);
        end
     

        
        % compute the element residual and tangent
        %
        [Re1,Re2,Re3,Ke11,Ke22,Ke33] = element(time,nInt,nNode,nodeCoords,Te1,Te2,Te3,...
            Te1Old,Te2Old,Te3Old,Area,D1,D2,D3,source1,source2,source3,...
            dtime,velocity1,velocity2,velocity3,...
            c01,c02,c03,alpha1,alpha2,alpha3);


        % check for any flux boundary conditions on the right side
        %
        Rbc = zeros(nNode,1);
        if(elem==elemFlux)
            Rbc(nodeFlux,1) = Area(Length)*flux;
        end


        % assemble this elements tangent into the gloabl
        %
        for i=1:nNode
            for j=1:nNode
                row = 3*(connect(elem,i)-1)+1;
                col = 3*(connect(elem,j)-1)+1;
                K(row,col) = K(row,col) + Ke11(i,j);
                K(row+1,col+1) = K(row+1,col+1) + Ke22(i,j); 
                K(row+2,col+2) = K(row+2,col+2) + Ke33(i,j); 
            end
        end

        % assemble this elements residual into the global
        %
        for i=1:nNode
            row = 3*(connect(elem,i) - 1) + 1;
            R(row) = R(row) + Re1(i);
            R(row + 1) = R(row + 1) + Re2(i);
            R(row + 2) = R(row + 2) + Re3(i);
  %          R(row) = R(row) + Re(i) + Rbc(i);
        end


    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % impose temperature boundary condition
    %
    % enforce the temperature BC
    %
%     T(nodeT) = Tbc(time);
    %
    % modify the tangent and residual
    %
%     K(nodeT,nodeT) = (1.e6)*(trace(K)/nNodes);
%     R(nodeT,1) = 0.0;
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check for convergence on the residual
    %
    if(norm(R)<1.e-10)
        fprintf('Converged on the residual in %i iterations \n', iter);
        break;
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Global system solve for the temperature feild
    %
    % compute the temperature field correction
    %
    Delta = K\R;
    %
    % compute the updated temperature field
    %
    T = T + Delta
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check for convergence on the corrections
    %
    if(norm(Delta)<1.e-10)
        fprintf('Converged on the corrections in %i iterations \n', iter);
        break;
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
end
% end the nonlinear iteration loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Postprocessing to plot the temperature field
%
if(time>=timeOut)
    postprocess(Length,coord,T,nNodes,c01,c02,NT2d,NT3d); 
    timeOut = timeOut + dtOut;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
% end the loop over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



function postprocess(Length,coord,T,nNodes,c01,c02,NT2d,NT3d)

% sort the all degree of freedom into different cell species 
            
%init
c1 = zeros(nNodes,1);
c2 = zeros(nNodes,1);
c3 = zeros(nNodes,1);



for i = 1:1:nNodes
    
    c3(i) = T(3*i);
    
    c2(i) = T(3*i-1);
    
    c1(i) = T(3*i-2);
end            
            
            

% cell density profile 
figure(1); 



% c1 c2 c3 profile
%subplot(4,2,7);

%plot(datac1(:,1),datac1(:,2),'ko')% abaqus 


%h1 = plot(coord/Length,c1/c01,'k-','LineWidth',2);
%hold on
h1 = plot(coord/Length,c2/c02,'k-','LineWidth',2);
hold on 
%h2 = plot(NT2d(:,1),NT2d(:,2)/c01,'ko','markersize',8);
%hold on 
h3 = plot(NT3d(:,1),NT3d(:,2)/c02,'ko','markersize',8);
hold off
xlim([0 1]);
%ylim([0 2500]);
xlabel('$x/l$ [-]','interpreter','latex');
ylabel('$c_1/c_{0}$ [-]','interpreter','latex');
legend([h1 h3],'1D FEM Matlab','3D U3D8 Abaqus','interpreter','latex','location','northwest')
legend('boxoff')
set(gca,'FontSize',23);








drawnow




end



function [Re1,Re2,Re3,Ke11,Ke22,Ke33] = element(time,nInt,nNode,nodeCoords,Te1,Te2,Te3,...
            Te1Old,Te2Old,Te3Old,Area,D1,D2,D3,source1,source2,source3,...
            dtime,velocity1,velocity2,velocity3,...
            c01,c02,c03,alpha1,alpha2,alpha3)

% obtain gauss points and weights
%
if(nInt==1);
    [xi,w] = GaussInt1Pt();
elseif(nInt==2);
    [xi,w] = GaussInt2Pt();
elseif(nInt==3);
    [xi,w] = GaussInt3Pt();
elseif(nInt==4);
    [xi,w] = GaussInt4Pt();
elseif(nInt==5);
    [xi,w] = GaussInt5Pt();
else
    error('nInt is not programmed');
end


% obtain nodal coordinates
%
if(nNode==2);
    x1 = nodeCoords(1,1);
    x2 = nodeCoords(2,1);
elseif(nNode==3);
    x1 = nodeCoords(1,1);
    x2 = nodeCoords(2,1);
    x3 = nodeCoords(3,1);
elseif(nNode==4);
    x1 = nodeCoords(1,1);
    x2 = nodeCoords(2,1);
    x3 = nodeCoords(3,1);
    x4 = nodeCoords(4,1);
else
    error('nNode is not programmed');
end


% init
%
Re1 = zeros(nNode,1);
Re2 = zeros(nNode,1);
Re3 = zeros(nNode,1);



Ke11 = zeros(nNode,nNode);
Ke22 = zeros(nNode,nNode); 
Ke33 = zeros(nNode,nNode); 


%
% loop over integration points
%
for intPt=1:nInt
    
    % compute shape functions and derivatives
    %
    if(nNode==2)
        [N,B,Jac] = shapeLinear(x1,x2,xi(intPt));
    elseif(nNode==3)
        [N,B,Jac] = shapeQuadratic(x1,x2,x3,xi(intPt));
    elseif(nNode==4)
        [N,B,Jac] = shapeCubic(x1,x2,x3,x4,xi(intPt));
    else
        error('nNode is not programmed');
    end
    
    % current location of this integ point
    %
    x = N*nodeCoords;
    
    
    % compute the temperature and dTdX at this integ point based on nodal
    % values of the temperature
    %
    Told1 = N*Te1Old;
    T1 = N*Te1;
    dTdX1 = B*Te1;
    dTdtime1 = (T1 - Told1)/dtime;
    
    Told2 = N*Te2Old;
    T2 = N*Te2;
    dTdX2 = B*Te2;
    dTdtime2 = (T2 - Told2)/dtime;
    
    Told3 = N*Te3Old;
    T3 = N*Te3;
    dTdX3 = B*Te3;
    dTdtime3 = (T3 - Told3)/dtime;    
    
    
    % update the element residual
    %


    Hfunc1 = (exp(alpha1*(T1 - c01)))/(1.0 + exp(alpha1*(T1 - c01)));
    Hfunc2 = (exp(alpha2*(T2 - c02)))/(1.0 + exp(alpha2*(T2 - c02)));
    Hfunc3 = (exp(alpha3*(T3 - c03)))/(1.0 + exp(alpha3*(T3 - c03)));
    
    
    flux1 = - T1 * Hfunc1 * velocity1(x,T1) + D1(x,T1) * dTdX1;
    flux2 = - T2 * Hfunc2 * velocity2(x,T2) + D2(x,T2) * dTdX2;
    flux3 = - T3 * Hfunc3 * velocity3(x,T3) + D3(x,T3) * dTdX3;


    Re1 = Re1 + Jac*w(intPt)*...
        (...
        transpose(N)*dTdtime1...
       +transpose(B)*flux1...
       -transpose(N)*source1(x,time)...
        );
    
    
    Re2 = Re2 + Jac*w(intPt)*...
        (...
        transpose(N)*dTdtime2...
       +transpose(B)*flux2...
       -transpose(N)*source2(x,time)...
        );   
    
    Re3 = Re3 + Jac*w(intPt)*...
        (...
        transpose(N)*dTdtime3...
       +transpose(B)*flux3...
       -transpose(N)*source3(x,time)...
        );    
    
    
    
    % this is just a hack for now
    %
%    dKdT1 = (D1(x,T1+.1)-D1(x,T1))/0.1;  
%    dKdT2 = (D2(x,T2+.1)-D2(x,T2))/0.1;
%    dKdT3 = (D3(x,T3+.1)-D3(x,T3))/0.1;

    
%    dSdT1 = (source1(x,T1+.1)-source1(x,T1))/0.1;
%    dSdT2 = (source2(x,T2+.1)-source2(x,T2))/0.1;
%    dSdT3 = (source3(x,T3+.1)-source3(x,T3))/0.1;


    % update the element tangent
    %
    

    fac11 = alpha1*exp(alpha1*(T1-c01));
    fac12 = alpha2*exp(alpha2*(T2-c02));
    fac13 = alpha3*exp(alpha3*(T3-c03));
    
    
    fac21 = alpha1*exp(2.0*alpha1*(T1-c01));
    fac22 = alpha2*exp(2.0*alpha2*(T2-c02));
    fac23 = alpha3*exp(2.0*alpha3*(T3-c03));
    
    
    fac31 = exp(alpha1*(T1-c01)) + 1.0;
    fac32 = exp(alpha2*(T2-c02)) + 1.0;
    fac33 = exp(alpha3*(T3-c03)) + 1.0;
    
    
    dHdc1 = fac11/fac31 - fac21/fac31^2;
    dHdc2 = fac12/fac32 - fac22/fac32^2;
    dHdc3 = fac13/fac33 - fac23/fac33^2;

    
    
    dqdc1 = -(Hfunc1 + T1*dHdc1)*velocity1(x,T1);
    dqdc2 = -(Hfunc2 + T2*dHdc2)*velocity2(x,T2);
    dqdc3 = -(Hfunc3 + T3*dHdc3)*velocity3(x,T3);

    
    
    Ke11 = Ke11 + Jac*w(intPt)*...
        (...
        -(1.0/dtime)*transpose(N)*N...
        - D1(x,T1)*transpose(B)*B...
        - dqdc1*transpose(B)*N...
        );    
    
    Ke22 = Ke22 + Jac*w(intPt)*...
        (...
        -(1.0/dtime)*transpose(N)*N...
        - D2(x,T2)*transpose(B)*B...
        - dqdc2*transpose(B)*N...
        );

    Ke33 = Ke33 + Jac*w(intPt)*...
        (...
        -(1.0/dtime)*transpose(N)*N...
        - D3(x,T3)*transpose(B)*B...
        - dqdc3*transpose(B)*N...
        );     

    
        
end % loop over integration points

return;
end

function [xi,w] = GaussInt1Pt()
% Gauss integration locations and weights for 1pt integration
xi = 0.0;
%
w = 2.0;
return;
end

function [xi,w] = GaussInt2Pt()
% Gauss integration locations and weights for 2pt integration
xi(1) = -sqrt(1.0/3.0);
xi(2) =  sqrt(1.0/3.0);
%
w(1) = 1.0;
w(2) = 1.0;
return;
end

function [xi,w] = GaussInt3Pt()
% Gauss integration locations and weights for 3pt integration
xi(1) = -0.7745966692;
xi(2) =  0.0;
xi(3) =  0.7745966692;
%
w(1) = 0.5555555556;
w(2) = 0.8888888889;
w(3) = 0.5555555556;
return;
end

function [xi,w] = GaussInt4Pt()
% Gauss integration locations and weights for 4pt integration
xi(1) = -0.8611363116;
xi(2) = -0.3399810436;
xi(3) =  0.3399810436;
xi(4) =  0.8611363116;
%
w(1) = 0.3478548451;
w(2) = 0.6521451549;
w(3) = 0.6521451549;
w(4) = 0.3478548451;
return;
end

function [xi,w] = GaussInt5Pt()
% Gauss integration locations and weights for 5pt integration
xi(1) = -0.9061798459;
xi(2) = -0.5384693101;
xi(3) =  0.0;
xi(4) =  0.5384693101;
xi(5) =  0.9061798459;
%
w(1) = 0.2369268851;
w(2) = 0.4786286705;
w(3) = 0.5688888889;
w(4) = 0.4786286705;
w(5) = 0.2369268851;
return;
end

function [N,B,Jac] = shapeLinear(x1,x2,xi)
% shape functions and derivatives for a 2-node 1D element

% element length
%
Le = x2 - x1;

% the shape function matrix
%
N = (1.0/2.0)*[1.0-xi 1.0+xi];

% derivatives of shape functions
%
B = (1.0/Le)*[-1.0 1.0];

% the mapping jacobian
%
Jac = Le/2.0;

return;
end

function [N,B,Jac] = shapeQuadratic(x1,x2,x3,xi)
% shape functions and derivatives for a 3-node 1D element

% the shape function matrix
%
N1 = (1/2)*xi*(xi - 1);
N2 = (1 + xi)*(1 - xi);
N3 = (1/2)*xi*(xi + 1);
N = [N1 N2 N3];


% derivatives of shape functions in the xi coordinate
%
dN1dXi = xi - 1/2;
dN2dXi = -2*xi;
dN3dXi = 1/2 + xi;
%
% the mapping jacobian
%
Jac = dN1dXi*x1 + dN2dXi*x2 + dN3dXi*x3;
%
% derivatives of shape functions in the x coordinate
%
B = [dN1dXi dN2dXi dN3dXi]/Jac;

return;
end

function [N,B,Jac] = shapeCubic(x1,x2,x3,x4,xi)
% shape functions and derivatives for a 4-node 1D element

% the shape function matrix
%
N1 = (-9/16)*(xi+1/3)*(xi-1/3)*(xi-1);
N2 = (27/16)*(xi+1)*(xi-1/3)*(xi-1);
N3 = (-27/16)*(xi+1)*(xi+1/3)*(xi-1);
N4 = (9/16)*(xi+1)*(xi+1/3)*(xi-1/3);
N = [N1 N2 N3 N4];

% derivatives of shape functions in the xi coordinate
%
dN1dXi = (9*xi)/8 - (27*xi^2)/16 + 1/16;
dN2dXi = (81*xi^2)/16 - (9*xi)/8 - 27/16;
dN3dXi = 27/16 - (81*xi^2)/16 - (9*xi)/8;
dN4dXi = (27*xi^2)/16 + (9*xi)/8 - 1/16;
%
% the mapping jacobian
%
Jac = dN1dXi*x1 + dN2dXi*x2 + dN3dXi*x3 + dN4dXi*x4;
%
% derivatives of the shape function in the x coordinate
%
B = [dN1dXi dN2dXi dN3dXi dN4dXi]/Jac;

return;
end

