%% Demo script for the phase-retrieval with non-orthogonal TO orthogonal geometry
%% conversion WHILE the phase retrieval is performed. The orthogonal frame chosen 
%% for the final display is the one that defined by the detection plane and the 
%% vector orthogonal to this plane. The reconsruction strategy uses the RC slices
%% one by one to retrieve the sample via projection/backprojections operations...
clear
close all

addpath('../../Display/')

%% Loading of the RC dataset, the support in the orthogonal / non orthogonal
% real-space representation, and the spatial and reciprocal space axis
load('dp.mat')

% We use the support in the non-orthogonal real-space representation
supp                        = supp_ortho;                           

%% Generation of the initial guess in the non-orthogonal real space
rng('default')
psi_ortho                   = supp .*exp(1i*2*pi*(rand(size(dp))-.5));

% . The mesh in the ORTHO real/reciprocal space...
[R1,R2,R3]                  = meshgrid(r1,r2,r3);
[Q1,Q2,Q3]                  = meshgrid(q1,q2,q3);

figure(10)
h = displayisosurf(abs(supp), 0.1*max(abs(supp(:))), 'g',r1,r2,r3);
axis([min(r1) max(r1) min(r2) max(r2) min(r3) max(r3)]),
grid
xlabel('r1 [um]'), ylabel('r2 [um]'), zlabel('r3 [um]')
title('Support in the ORTHO FRAME')

figure(11)
h = displayisosurf(psi_ortho, -.01, 'g',r1,r2,r3);
axis([min(r1) max(r1) min(r2) max(r2) min(r3) max(r3)]),
xlabel('r1 [um]'), ylabel('r2 [um]'), zlabel('r3 [um]')
title('|psi| in the  ORTHO FRAME')
grid

figure(12)
subplot(121)
imagesc(r2, r3, abs(squeeze(psi_ortho(:,length(r1)/2,:))))
axis image, axis xy
xlabel('r2 [um]'), ylabel('r3 [um]')
title('|psi| in the  ORTHO FRAME')
colorbar, 
grid
subplot(122)
imagesc(r2, r3, angle(squeeze(psi_ortho(:,length(r1)/2,:))))
axis image, axis xy
xlabel('r2 [um]'), ylabel('r3 [um]')
title('Angle(psi) in the  ORTHO FRAME')
colorbar, 
grid


%% ER: perform the phase retrieval using iterative algorithm (via 2DFFT + projections/backprojections)
%% NOTE: the result is in the ORTHO real-space frame (e_1,e_2, e_3}. 
iter_num                    = 20;           % Number of ER update
alpha                       = 1;            % Updating stepsize for ER

[psi_ortho,supp,dp_error] = ER_ortho_BP(sqrt(dp),...
    psi_ortho, supp, ...
    [r1; r2; r3; q1; q2; q3], R3, R2, theta_B,...    
    alpha, iter_num);
