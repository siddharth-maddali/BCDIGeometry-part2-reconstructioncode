clear
close all

addpath('../../Display/')

% The numerical data can be plagued by shot-noise noise
NOISE       = 0;
max_counts  = 1e4;

%% Defnition of the Bragg angle...
theta_B     = 30*pi/180;
detB     	= cos(theta_B);

%% Definition of the square shape in the orthogonal real space
s1          = 2;
s2          = 2;
s3          = 2;
phi         = 15*pi/180;

%% Definition of the \ortho spatial-vectors 
r1_min      = -25; 
r1_max      = 25; 
d_r1        = s1/10;
r2_min      = -25; 
r2_max      = 25; 
d_r2        = s2/10;
r3_min      = -25; 
r3_max      = 25; 
d_r3        = s3/10;

r1          = r1_min:d_r1:(r1_max-d_r1);
r2          = r2_min:d_r2:(r2_max-d_r2);
r3          = r3_min:d_r3:(r3_max-d_r3);

%% Definition of the probe in the orthogonal real space
Probe_ortho = ones(length(r2), length(r1), length(r3));

%% Generation of the square shape in the orthogonal real space 
% (we allow a rotation of the square in the (e2,e3) plan)
[R1,R2,R3]  = meshgrid(r1,r2,r3);
R3_phi      = (cos(phi)*R2 -sin(phi)*R3);
R2_phi      = (sin(phi)*R2 +cos(phi)*R3);
R1_phi      = R1;
Sq_ortho    = (abs(R1_phi) < s1 & abs(R2_phi) < s2 & abs(R3_phi) < s3);

% smooth
G_kernel    = exp(-.5*(R1.^2/(2*(.5*d_r1)^2) + R2.^2/(2*(.5*d_r2)^2) + R3.^2/(2*(.5*d_r3)^2) ));
Sq_ortho    = ifftn(fftn(ifftshift(Sq_ortho)) .* fftn(G_kernel));
Sq_ortho    = Sq_ortho/max(abs(Sq_ortho(:)));

% Generation of the exit-field in the orthogonal real-space
psi_ortho   = Sq_ortho .* Probe_ortho;

%% definition of the orthogonal recirpocal-space

% The _frequency_ sampling rate along (q1,q2,q3) is chosen arbitrarly (on the contrar2 to any real-life experimental condition)
d_q1        = (1/d_r1)/length(r1);
q1          = 0:d_q1:((1/d_r1)-d_q1) ;
q1          = q1 - .5/d_r1;
d_q2        = (1/d_r2)/length(r2);
q2          = 0:d_q2:((1/d_r2)-d_q2) ;
q2          = q2 - .5/d_r2;
d_q3        = (1/d_r3)/length(r3);
q3          = 0:d_q3:((1/d_r3)-d_q3) ;
q3          = q3 - .5/d_r3;

% . The _frequency_ sampling rate along q3_hat is such that q3_hat*cos(\theta) = q3 (for the sake of simlicity) 
alpha               = 1/cos(theta_B);                               
d_q3_hat            = alpha*d_q1;
q3_hat              = 0:d_q3_hat:((1/d_r1)-d_q3_hat);   
q3_hat              = q3_hat - .5/d_r1;

%% calculate the diffraction wavefield (via a 3DFT)
% we apply the coordinate transform via the Fourier Transform and a phase ramp

% The RC measurements are the intensities of the FT of the exit field wrt non-orthogonal reciprocal representation
[Q1,Q2,Q3]              = meshgrid(q1,q2,q3);
phase_ramp              = exp(sqrt(-1)*2*pi*R2.*(Q3*tan(theta_B)));
PSI_non_ortho           = Real_TO_NonOrthoFourier(psi_ortho,phase_ramp,[d_r1, d_r2, d_r3, theta_B]);

% Sample in the non-orthogonal frame
psi_non_ortho           = ifftn(ifftshift(PSI_non_ortho));


%% perform the measurement and add statistic noise
% intensity measurement
dp          = abs(PSI_non_ortho).^2;

% add poisson noise
if NOISE  
    dp          = dp/max(dp(:))*max_counts;
    
    for n = 1 : size(dp,3)
        temp        = dp(:,:,n)/1e12;
        temp        = imnoise(temp,'poisson');
        dp(:,:,n)   = temp*1e12;
    end
end


% Extract the orth/non_orthogonal support
supp_ortho      = double(abs(psi_ortho) > 1e-2 * max(abs(psi_ortho(:))));
supp_non_ortho  = double(abs(psi_non_ortho) > 1e-2 * max(abs(psi_non_ortho(:))));

% We save the data and supports...
save('dp.mat','dp', 'supp_ortho', 'psi_ortho', 'supp_non_ortho', 'psi_non_ortho', 'theta_B', 'r1', 'r2', 'r3', 'q1', 'q2', 'q3', 'q3_hat')


%% For display purpose only: we want a display that fits with the view of the 
%% paper => we need to permute dim [r2 r1 r3] to [r2 r3 r1]

psi_ortho_disp  = permute(psi_ortho, [2 3 1]);
PSI_ortho_disp  = fftshift(fftn(psi_ortho_disp));
dp_disp         = permute(dp, [2 3 1]);
psi_non_ortho_disp = permute(psi_non_ortho, [2 3 1]);


figure(1)
h = displayisosurf(psi_ortho, -.45, 'g',r1,r2,r3);
axis(.5*[min(r1) max(r1) min(r2) max(r2) min(r3) max(r3)]),
grid
xlabel('r1 [\mu m]'), ylabel('r2 [\mu m]'), zlabel('r3 [\mu m]')

figure(2),
imagesc(r2,r3,squeeze(psi_ortho(:,length(r1)/2,:)))
xlabel('r_2 [\mu m]')
ylabel('r_3 [\mu m]')
title('psi in orth frame')
axis xy
grid, colorbar, colormap(1-gray)

figure(3)
h = displayisosurf(abs(fftshift(fftn(psi_ortho))), -.05, 'g',q1,q2,q3);
axis(.5*[min(q1) max(q1) min(q2) max(q2) min(q3_hat) max(q3_hat)]),
grid
xlabel('q1 [1/um]'), ylabel('q2 [1/um]'), zlabel('q3 [1/um]')
title('|PSI| in the ORTH FRAME')

figure(4)
h = displayisosurf(dp, -.002, 'g',q1,q2,q3);
axis(.5*[min(q1) max(q1) min(q2) max(q2) min(q3_hat) max(q3_hat)]),
grid
xlabel('q1 [1/um]'), ylabel('q2 [1/um]'), zlabel('q3 [1/um]')
title('|PSI| in the DETECTION FRAME (FT METH)')

figure(11)
h = displayisosurf(dp_disp, -.002, 'g',q3,q1,q2);
axis(.5*[min(q3_hat) max(q3_hat) min(q1) max(q1) min(q2) max(q2)]),
grid
xlabel('q3 [1/um]'), ylabel('q1 [1/um]'), zlabel('q2 [1/um]')
title('|PSI| in the DETECTION FRAME')

figure(5)
imagesc( q2, q3, log10(abs(squeeze(dp(:,length(q2)/2,:)))) )
axis image, axis xy
axis([min(q3_hat) max(q3_hat) min(q3) max(q3)]),
title('|PSI| in the DETECTION FRAME')
xlabel('q_2 [1/\mu m]'), ylabel('q_3 [1/\mu m]')
colorbar, 
grid

figure(6)
h = displayisosurf(abs(psi_non_ortho), -.01, 'g',r1,r2,r3);
axis(.5*[min(r1) max(r2) min(r3) max(r2) min(r3) max(r3)]),
grid
xlabel('r1 [um]'), ylabel('r2 [um]'), zlabel('r3 [um]')
title('|psi| in the NON ORTHO REAL SPACE')

figure(7),
imagesc(r2,r3,squeeze(imag(psi_non_ortho(:,length(r2)/2,:))))
xlabel('r_3 [\mu m]')
ylabel('r_2 [\mu m]')
title('psi in the non orth frame')
axis xy
grid, colorbar, colormap(1-gray)

%% Les commandes d'affichage ci-dessous permettent les affichages présentés dans le papier 

AFFICHAGE_PAPIER = 0
if(AFFICHAGE_PAPIER)
figure(8)
h = displayisosurf(psi_ortho_disp, -.45, 'g',r3,r1,r2);
axis(.5*[min(r3) max(r3) min(r1) max(r1) min(r2) max(r2)]),
grid
xlabel('r3 [\mu m]'), ylabel('r1 [\mu m]'), zlabel('r2 [\mu m]')

figure(9),
imagesc(r3,r2,squeeze(psi_ortho_disp(length(r1)/2,:,:)).')
xlabel('r_3 [\mu m]')
ylabel('r_2 [\mu m]')
title('psi in orth frame')
axis xy
grid, colorbar, colormap(1-gray)

figure(10)
h = displayisosurf(abs(PSI_ortho_disp), -.05, 'g',q3,q1,q2);
axis(.5*[min(q3_hat) max(q3_hat) min(q1) max(q1) min(q2) max(q2)]),
grid
xlabel('q3 [1/um]'), ylabel('q1 [1/um]'), zlabel('q2 [1/um]')
title('|PSI| in the ORTH FRAME')

figure(12)
imagesc(q3, q2, log10(abs(squeeze(dp_disp(length(q2)/2,:,:)).')))
axis image, axis xy
axis([min(q3_hat) max(q3_hat) min(q2) max(q2)]),
title('|PSI| in the DETECTION FRAME')
xlabel('q_3 [1/\mu m]'), ylabel('q_2 [1/\mu m]')
colorbar, 
grid

figure(13)
h = displayisosurf(psi_non_ortho_disp, -.45, 'g',r3,r1,r2);
axis(.5*[min(r3) max(r3) min(r3) max(r3) min(r2) max(r2)]),
grid
xlabel('r3 [\mu m]'), ylabel('r1 [\mu m]'), zlabel('r2 [\mu m]')
title('|psi| in the NON ORTHO REAL SPACE')

figure(14),
imagesc(r3,r2,squeeze(abs(psi_non_ortho_disp(length(r1)/2,:,:))).')
xlabel('r_3 [\mu m]')
ylabel('r_2 [\mu m]')
title('psi in non orth frame')
axis xy
grid, colorbar, colormap(1-gray)
end

