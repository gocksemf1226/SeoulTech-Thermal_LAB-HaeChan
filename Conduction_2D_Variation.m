%%%2D Conduction Topology Optimization%%%
%%% Made by HaeChan Jeong - HC &&&
% nely = 40; nelx = 40; volfrac = 0.4; penal = 3.0; rmin=1.2;
% 40,40,0.4,3.0,1.2

function Conduction_2D_Variation(nelx,nely,volfrac,penal,rmin);
% INITIALIZE
x(1:nely,1:nelx) = volfrac; 
loop = 0; 
change = 1.;
iteration_vector = []; %HC
change_vector = []; %HC
objective_vector = []; %HC
mean_temp_vector = []; %HC
Ue_vector = zeros(nely,nelx); %HC
% Force = []; %HC
% START ITERATION

% while loop > 4
while change > 1e-3

  loop = loop + 1;
  xold = x;
% FE-ANALYSIS
  [K, U, U_avg, lambda_av, lambda_vr]=FE(nelx,nely,x,penal);  

  U_var = var(U, 1);

% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [KE] = lk;
  c = 0.;
  w1 = 1e-1; w2 = 3.097e-4;
  N = 2*2;
  unit_vec = ones(N,1); 
  for ely = 1:nely
    for elx = 1:nelx
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      Ue = U([n1; n2; n2+1; n1+1],1);
      Ue_mean(ely,elx) = mean(Ue); %HC
      dKE = 0.999*penal*x(ely,elx)^(penal-1)*KE;
      lambda_av_e = lambda_av([n1; n2; n2+1; n1+1],1);
      lambda_vr_e = lambda_vr([n1; n2; n2+1; n1+1],1);

      %%%% Average temperature %%%%
      av = 1 / N * (unit_vec' * Ue);
      dav(ely,elx) = lambda_av_e' * dKE * Ue;
      phi_av = w1 * dav;

      %%%% Variance temperature %%%%
      vr = 1 / N * (Ue - unit_vec * U_avg)' * (Ue - unit_vec * U_avg);
      dvr(ely,elx) = 2/ N * dav(ely,elx) * unit_vec' * (Ue-unit_vec * U_avg) + lambda_vr_e' * dKE * Ue;
      phi_vr = w2 * dvr;


      dphi(ely,elx) = phi_av(ely,elx) + phi_vr(ely,elx);
    end
  end

      assignin('base', 'dphi', dphi);


% FILTERING OF SENSITIVITIES
  [dphi]   = check(nelx,nely,rmin,x,dphi);    

% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
  [x]    = OC(nelx,nely,x,volfrac,dphi); 

% PRINT RESULTS
  change = max(max(abs(x-xold)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj1.: ' sprintf('%10.4f', av) ' Obj2.: ' sprintf('%10.4f', vr) ...
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
  set(fig1, 'OuterPosition', [0, 100, 500,500]) %HC
%   mean_temp_plot = subplot(212); colormap(mean_temp_plot, jet); %HC
  colormap(jet); imagesc(Ue_mean); colorbar; axis equal; axis tight; axis off; pause(1e-8); %HC
  fig2 = figure(2); %HC
  set(fig2, 'OuterPosition', [500, 100, 500,500]) %HC

end


figure(3); %HC
fig3 = figure(3); %HC
set(fig3, 'OuterPosition', [400, 0, 800, 600]); %HC
subplot(311); plot(iteration_vector, change_vector, 'b-'); %HC
grid on; title('change'); xlabel('iteration'); ylabel('change'); %HC
subplot(312); plot(iteration_vector, objective_vector, 'y-'); %HC
grid on; title('Objective Function'); xlabel('iteration'); ylabel('Obj values'); %HC
subplot(313); plot(iteration_vector, mean_temp_vector, 'r-');
grid on; title('Mean Temperature'); xlabel('interation'); ylabel('Mean T'); %HC

U_var = var(U, 1)
mean(dav)
mean(dvr)

% % Prescribed Heat Plot
% fig4 = figure(4); %HC
% set(fig4, 'OuterPosition', [0, 400, 400,800]) %HC
% rshForce = reshape(Force,[nely+1 nelx+1]);
% colormap(jet); imagesc(rshForce); colorbar; axis equal; axis tight; axis off; pause(1e-2); %HC

% end
%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dphi)  
l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dphi./lmid)))));
  if sum(sum(xnew)) - volfrac*nelx*nely > 0;
    l1 = lmid;
  else
    l2 = lmid;
  end
end
% end
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dphin]=check(nelx,nely,rmin,x,dphi)
dphin=zeros(nely,nelx);
for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
        fac = rmin-sqrt((i-k)^2+(j-l)^2);
        sum = sum+max(0,fac);
        dphin(j,i) = dphin(j,i) + max(0,fac)*x(l,k)*dphi(l,k);
      end
    end
    dphin(j,i) = dphin(j,i)/(x(j,i)*sum);
  end
end
% end
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K, U, U_avg, lambda_av, lambda_vr]=FE(nelx,nely,x,penal)
[KE] = lk; 
K = sparse((nelx+1)*(nely+1), (nelx+1)*(nely+1));
F = sparse((nely+1)*(nelx+1),1); U = zeros((nely+1)*(nelx+1),1);

N = (nely+1)*(nelx+1); unit_vec = ones((nely+1)*(nelx+1),1);
lambda_av = zeros((nely+1)*(nelx+1),1); lambda_vr = zeros((nely+1)*(nelx+1),1);
% U_avg = zeros((nely+1)*(nelx+1),1);

for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [n1;n2;n2+1;n1+1];
    K(edof,edof) = K(edof,edof) + (0.001+0.999*x(ely,elx)^penal)*KE;
  end    
end
% DEFINE LOADS AND SUPPORTS (SQURE PLATE WITH HEAT SINK)

assignin('base','Stiffness_vector',K); %HC

%%% Hot_Spot Zone Setting %%%
% Distributed Heat Zone
F(:,1)=0.01;
% Hot_Spot Zone Center Position
xi = 11; xf = 31; yi = 11; yf = 31;

% % Hot_Spot Zone - SQUARE
% for i = xi:xf
%     for j = yi:yf
%         hotdofs = [1+(i-1)*(nely+1)+(j-1)];
%         F(hotdofs,1) = 0.1;
%     end
% end

% % Hot_Spot Zone - CONTOUR
% for m = 0:(xf-xi)/2
%     n = m;
%     for i = [(xi+m) (xf-m)]
%         for j = (yi+n):(yf-n)
%             hotdofs = [1+(i-1)*(nely+1)+(j-1)];
%             F(hotdofs,1) = 0.01*(m+1);
%         end
%     end
%     for j = [(yi+n) (yf-n)]
%         for i = (xi+n):(xf-n)
%             hotdofs = [1+(i-1)*(nely+1)+(j-1)];
%             F(hotdofs,1) = 0.01*(m+1);
%         end
%     end
% end

% % Force term to 'base workspace'
% assignin('base','Force',F); %HC

% Boundary Condition Position
Left = 1+nely/2-(nely/20):1+nely/2+(nely/20);
Right = ((nely+1)*nelx)+1+nely/2-(nely/20):((nely+1)*nelx)+1+nely/2+(nely/20);
Top = 1+(nely+1)*((nelx/2)-(nelx/20)):nely+1:1+(nely+1)*((nelx/2)+(nelx/20));
Bottom = (nely+1)*(1+(nelx/2)-(nelx/20)):nely+1:(nely+1)*(1+(nelx/2)+(nelx/20));

fixeddofs = [Left];
alldofs = [1:(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;
assignin('base', 'average_temp', U);

U_avg = 1/N * unit_vec' * U;

lambda_av(freedofs,:) = -1 / N * K(freedofs,freedofs) \ unit_vec(freedofs,:);      
lambda_av(fixeddofs,:)= 0;

lambda_vr(freedofs,:) = -2 / N * K(freedofs,freedofs) \ (U(freedofs,:) - unit_vec(freedofs,:) * U_avg);
lambda_vr(fixeddofs,:)= 0;

%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
KE = [2/3 -1/6 -1/3 -1/6
    -1/6 2/3 -1/6 -1/3
    -1/3 -1/6 2/3 -1/6
    -1/6 -1/3 -1/6 2/3];
% end
