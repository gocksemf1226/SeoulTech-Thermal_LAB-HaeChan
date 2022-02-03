%%% 0. MAIN-LOOP %%%
function top(nelx,nely,volfrac,penal,rmin);     
% INITIALIZE
x(1:nely,1:nelx) = volfrac; % distributing the material evenly in the design domain
loop = 0; 
change = 1.;
iteration_vector = []; %HC
change_vector = []; %HC
objective_vector = []; %HC
design_variable_vector = []; %HC
% START ITERATION
while change > 1e-2  % the main-loop is terminated if the change in design variables is less than 1 percent
  loop = loop + 1;
  xold = x;
% FE-ANALYSIS
  [U]=FE(nelx,nely,x,penal);     % a call to the Finite Element subroutine %%%%% 2. FE-ANALYSIS %%%%%
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [KE] = lk; % element stiffness matrix subroutine is called only once
  c = 0.;
  for ely = 1:nely % 19-28_a loop over all EL determines obj fnc & sensitivities
    for elx = 1:nelx
      n1 = (nely+1)*(elx-1)+ely; % upper left EL N# in global node numbers
      n2 = (nely+1)* elx   +ely; % upper right EL N# in global node numbers 
      % --> 'n1','n2' are used to extract the EL disp vector Ue from the global disp vector U
      Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1); % EL disp vector
      c = c + x(ely,elx)^penal*Ue'*KE*Ue;
      dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;      % Eq(4)
    end
  end
% FILTERING OF SENSITIVITIES
  [dc]   = check(nelx,nely,rmin,x,dc);  % sensitivity analysis is followed by a call to the [MESH INDEPENDENCY FILTER]  
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
  [x]    = OC(nelx,nely,x,volfrac,dc);  % sensitivity analysis is followed by a call to the [OPTIMALITY CRITERIA optimizer]. 'x' is the vector of design variables.
% PRINT RESULTS                     % 34-37_the current compliance as well as other parameters are printed
  change = max(max(abs(x-xold)));   % 기존 design variable 'xold'와 updated design variable 'x' 행렬의 차이의 절댓값의 최댓값
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
  iteration_vector(loop) = loop; %HC
  change_vector(loop) = change; %HC
  objective_vector(loop) = c; %HC
% PLOT DENSITIES  
  colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);    % the resulting density distribution is plotted
  fig1 = figure(1); %HC
  set(fig1, 'OuterPosition', [0, 0, 420, 420]) %HC
  assignin('base','design_variable_vector',x); %HC
end
figure(2); %HC
fig2 = figure(2); %HC
set(fig2, 'OuterPosition', [1, 1, 1000, 800]); %HC
subplot(211); %HC
plot(iteration_vector, change_vector, 'b-'); %HC
grid on; %HC
title('change'); %HC
xlabel('iteration'); %HC
ylabel('change'); %HC
subplot(212); %HC
plot(iteration_vector, objective_vector, 'y-'); %HC
grid on; %HC
title('Objective Function'); %HC
xlabel('iteration'); %HC
ylabel('Obj values'); %HC
% final design variable을 xlsx에 저장 ※속도 느려짐 - 사용시 line48 - ctrl+t
writematrix(x, 'design_variable_vector.xlsx','Sheet', 'Sheet1', 'Range', 'A1:BH20'); %HC
%%%%%%%%%% 4. OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dc)  % the updated design variables are found by the optimizer
l1 = 0; l2 = 100000; move = 0.2;    % lower & upper bound for the Lagrangian multiplier
while (l2-l1 > 1e-4)    % [bi-sectioning algorithm] : 64-73_the value of the Lagrangian multiplier that satisfies the volume constraint can be found. 
                        % 64_The interval which bounds the Lagrangian multiplier is repeatedly halved until its size is less than the convergence criteria 
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
  if sum(sum(xnew)) - volfrac*nelx*nely > 0;    % material volume(monotonously decreasing fnc of the Lagrangian multiplier-lag)
    l1 = lmid;
  else
    l2 = lmid;
  end
end
%%%%%%%%%% 3. MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc)
dcn=zeros(nely,nelx);
for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
        fac = rmin-sqrt((i-k)^2+(j-l)^2);
        sum = sum+max(0,fac);
        dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
      end
    end
    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
  end
end
%%%%%%%%%% 2. FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,x,penal)
[KE] = lk; 
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));       % 'sparse' : 메모리 낭비 방지를 위해 행렬값 대부분이 0으로 채워진 행렬
F = sparse(2*(nely+1)*(nelx+1),1); 
U = zeros(2*(nely+1)*(nelx+1),1);
for elx = 1:nelx    %96-104_ The global stiffness matrix is formed by a loop over all EL. 
                    % 'n1','n2' are used to insert the EL stiffness matrix at the right places in the global stiffness matrix.
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
  end
end
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM) % N & EL are numbered column wise  % Boundary Conditions
F(2,1) = -1;    % column wise from left to right, 2 dof --> This command applies a vertical unit force force in the upper left corner.
fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);    % it is easier to define the dof that are fixed
alldofs     = [1:2*(nely+1)*(nelx+1)];     
freedofs    = setdiff(alldofs,fixeddofs);   %'freedofs' indicate the dof which are unconstrained.
                                            %84_'freedofs' are found automatically using the MATLAB operator 'setdiff' which finds the dof as the difference between all dof and the fixed dof
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);     % Supports are implemented by eliminating fixed dof from the linear EQs. 
U(fixeddofs,:)= 0;
%%%%%%%%%% 1. ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [KE]=lk    % The 8 by 8 matrix for a square bi-linear 4-node EL was determined analytically 
%                     using a symbolic manipulation software.
E = 1.;     % Young's modulus
nu = 0.3;   % Poisson's ratio
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 -1/4+nu/12 -1/8-nu/8  nu/6  1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Ole Sigmund, Department of Solid         %
% Mechanics, Technical University of Denmark, DK-2800 Lyngby, Denmark.     %
% Please sent your comments to the author: sigmund@fam.dtu.dk              %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "A 99 line topology optimization code written in Matlab"                 %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but does not guaranty that the code is    %
% free from errors. Furthermore, he shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
