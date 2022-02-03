%%%2D Conduction Topology Optimization%%%
%%% Made by HaeChan Jeong - HC &&&
function Conduction_2D_top(nelx,nely,volfrac,penal,rmin);
% INITIALIZE
x(1:nely,1:nelx) = volfrac; 
loop = 0; 
change = 1.;
iteration_vector = []; %HC
change_vector = []; %HC
objective_vector = []; %HC
mean_temp_vector = []; %HC
Ue_vector = zeros(nely,nelx); %HC

% START ITERATION
while change > 1e-2
  loop = loop + 1;
  xold = x;
% FE-ANALYSIS
  [U]=FE(nelx,nely,x,penal);         
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [KE] = lk;
  c = 0.;
  for ely = 1:nely
    for elx = 1:nelx
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      Ue = U([n1; n2; n2+1; n1+1],1);
      Ue_mean(ely,elx) = mean(Ue); %HC
      c = c + (0.001+0.999*x(ely,elx)^penal)*Ue'*KE*Ue;
      dc(ely,elx) = -0.999*penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
    end
  end
% FILTERING OF SENSITIVITIES
  [dc]   = check(nelx,nely,rmin,x,dc);    
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
  [x]    = OC(nelx,nely,x,volfrac,dc); 


% PRINT RESULTS
  change = max(max(abs(x-xold)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change ) ...
        ' mean Temp.: ' sprintf('%6.3f', mean(U))])
  iteration_vector(loop) = loop; %HC
  change_vector(loop) = change; %HC
  objective_vector(loop) = c; %HC
  mean_temp_vector(loop) = mean(U); %HC
% PLOT DENSITIES
%   top_plot = subplot(211); colormap(top_plot, gray); %HC
  colormap(gray); imagesc(-x); axis equal; axis tight; axis off; 
  fig1 = figure(1); %HC
  set(fig1, 'OuterPosition', [0, 0, 400,400]) %HC
%   mean_temp_plot = subplot(212); colormap(mean_temp_plot, jet); %HC
  colormap(jet); imagesc(Ue_mean); colorbar; axis equal; axis tight; axis off; pause(1e-1); %HC
  fig2 = figure(2); %HC
  set(fig2, 'OuterPosition', [0, 400, 400,400]) %HC
end

% Plot
figure(3); %HC
fig3 = figure(3); %HC
set(fig3, 'OuterPosition', [1000, 0, 1000, 800]); %HC
subplot(311); plot(iteration_vector, change_vector, 'b-'); %HC
grid on; title('change'); xlabel('iteration'); ylabel('change'); %HC
subplot(312); plot(iteration_vector, objective_vector, 'y-'); %HC
grid on; title('Objective Function'); xlabel('iteration'); ylabel('Obj values'); %HC
subplot(313); plot(iteration_vector, mean_temp_vector, 'r-');
grid on; title('Mean Temperature'); xlabel('interation'); ylabel('Mean T'); %HC



%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dc)  
l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
  if sum(sum(xnew)) - volfrac*nelx*nely > 0;
    l1 = lmid;
  else
    l2 = lmid;
  end
end



%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,x,penal)
[KE] = lk; 
K = sparse((nelx+1)*(nely+1), (nelx+1)*(nely+1));
F = sparse((nely+1)*(nelx+1),1); U = zeros((nely+1)*(nelx+1),1);
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [n1;n2;n2+1;n1+1];
    K(edof,edof) = K(edof,edof) + (0.001+0.999*x(ely,elx)^penal)*KE;
  end    
end
% DEFINE LOADS AND SUPPORTS (SQURE PLATE WITH HEAT SINK)
F(:,1)=0.01;
fixeddofs = [nely/2+1-(nely/20):2:nely/2+1+(nely/20)];
alldofs = [1:(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;



%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
KE = [2/3 -1/6 -1/3 -1/6
    -1/6 2/3 -1/6 -1/3
    -1/3 -1/6 2/3 -1/6
    -1/6 -1/3 -1/6 2/3];



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
