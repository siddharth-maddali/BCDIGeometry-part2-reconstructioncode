%% Demo script for the phase-retrieval with non-orthogonal TO orthogonal geometry
%% conversion AFTER the phase retrieval is performed. The orthogonal frame chosen 
%% for the final display is the one that is defined by the detection plane and the 
%% vector orthogonal to this plane.

clear
close all

addpath('../../Display/')

% Loading of the RC dataset, the support in the orthogonal / non orthogonal
% real-space representation, and the spatial and reciprocal space axis
load('dp.mat')

% We use the support in the non-orthogonal real-space representation
supp                = supp_non_ortho;                           

%% Generation of the initial guess in the non-orthogonal real space
rng('default')
psi_non_ortho       = supp .*exp(1i*2*pi*(rand(size(dp))-.5));

%% ER: perform the phase retrieval using iterative algorithm (via 3DFFT)
%% NOTE: the result is in the NON-ORTHO real-space frame (\bar{e}_1,\bar{e}_2, \bar{e}_3}. 
iter_num    = 10;           % Number of ER update
alpha       = 1;            % Updating stepsize for ER

[psi_non_ortho,dp_error] = ER(sqrt(dp), psi_non_ortho, supp, [r2; r3], alpha, iter_num);

%% DISPLAY: 
%% We use the appropriate transformation to produce a reconstruction in the ORTHO frame (e_1,e_2,e_3}.

[R1,R2,R3]          = meshgrid(r1,r2,r3);
[Q1,Q2,Q3]          = meshgrid(q1,q2,q3);

psi_non_ortho       = psi_non_ortho .* supp;
param               = [r1(2)-r2(1), r2(2)-r2(1), r3(2)-r3(1), theta_B];
phase_ramp          = exp(1i*2*pi*R2.*Q3*tan(theta_B));
PSI_non_ortho    	= fftn(psi_non_ortho);
psi_ortho           = NonOrthoFourier_TO_Real(PSI_non_ortho,phase_ramp,param);

figure(10)
h = displayisosurf(abs(supp), 0.1*max(abs(supp(:))), 'g',r1,r2,r3);
axis(.5*[min(r1) max(r1) min(r2) max(r2) min(r3) max(r3)]),
grid
xlabel('r1 [um]'), ylabel('r2 [um]'), zlabel('r3 [um]')
title('Support in the NON ORTHO FRAME')

figure(11)
h = displayisosurf(abs(psi_non_ortho), 0.1*max(abs(psi_non_ortho(:))), 'g',r1,r2,r3);
axis([min(r1) max(r1) min(r2) max(r2) min(r3) max(r3)]),
grid
xlabel('r1 [um]'), ylabel('r2 [um]'), zlabel('r3 [um]')
title('|psi| in the NON ORTHO FRAME')

figure(12)
subplot(121)
imagesc(r2, r3, abs(squeeze(psi_non_ortho(:,length(r1)/2,:))))
axis image, axis xy
xlabel('r2 [um]'), ylabel('r3 [um]')
title('|psi| in the NON ORTHO FRAME')
colorbar, 
grid
subplot(122)
imagesc(r2, r3, angle(squeeze(psi_non_ortho(:,length(r1)/2,:))))
axis image, axis xy
xlabel('r2 [um]'), ylabel('r3 [um]')
title('Angle(psi) in the NON ORTHO FRAME')
colorbar, 
grid

figure(13)
h = displayisosurf(abs(psi_ortho), 0.1*max(abs(psi_ortho(:))), 'g',r1,r2,r3);
axis([min(r1) max(r1) min(r2) max(r2) min(r3) max(r3)]),
grid
xlabel('r1 [um]'), ylabel('r2 [um]'), zlabel('r3 [um]')
title('|psi| in the ORTHO FRAME')

figure(14)
subplot(121)
imagesc( r2, r3, abs(squeeze(psi_ortho(:,length(r1)/2,:))))
axis image, axis xy
xlabel('r2 [um]'), ylabel('r3 [um]')
title('|psi| in the ORTHO FRAME')
colorbar, 
grid
subplot(122)
imagesc( r2, r3, angle(squeeze(psi_ortho(:,length(r1)/2,:))) )
axis image, axis xy
axis([min(r2) max(r2) min(r3) max(r3)]),
title('Angle(psi) in the ORTHO FRAME')
xlabel('r2 [\mu m]'), ylabel('r3 [\mu m]')
colorbar, %colormap(1-pink)
grid