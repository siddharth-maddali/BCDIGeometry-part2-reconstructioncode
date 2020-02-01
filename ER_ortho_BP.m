function [psi_ortho,supp,dp_error] = ER_ortho_BP(dp_sqrt,... 
    psi_ortho, supp, ...
    geo_param, R3, R2, theta_B,... 
    alpha, iter_num)


%
% [psi_nonortho,dp_error] = ER_ortho(dp_sqrt,psi_nonortho, supp, geo_param, alpha, iter_num)
%
% Description: 
%   Support-based, 3D phase-retrieval routine via the standard error-reduction
%   routine.  
%   
%   If the 3D Fourier components measurement are obtained in a non-orthogonal 
%   geometry,  the sample will be retrieved in an ORTHOGHONAL real-space frame. 
%   This frame is conjugate with the orthogonal reciprocal-space frame (k_1, k_2, k_3) 
%   defined by the detection plane and the unitary vector perpendicular to the detection plane.
%
% Inputs:
%
%   dp_sqrt         \in IN^{N2 x N1 x N3}   :   square root of the 3D intensity measurement stack (i.e., the sqrt of the intensity extracted along the RC)
%   psi_ortho       \in IC^{N2 x N1 x N3}   :   initial-guess of the 3D retrieved field 
%   supp            \in IN^{N2 x N1 x N3}   :   binary map defining the support of the real-space exit-field IN 
%                                               THE FRAME CONJUGATE WITH THE ORTHOGONAL MEASUREMENT GEOMETRY (k_1, k_2, k_3)
%   geo_param       \in IR^{6 X N1}         :   the axis along e1, e2, e3 and k_1, k_2, k_3
%   R2              \in IR^{N2 x N1 x N3}   :   Coordinate matrix as provided by [R1,R2,R3] = meshgrid(r1,r2,r3);
%   Q3              \in IR^{N2 x N1 x N3}   :   Coordinate matrix as provided by [Q1,Q2,Q3] = meshgrid(q1,q2,q3);
%   theta_B         \in IR                  :   Bragg angle [rad]
%   alpha           \in IR_+                :   updating step-size for the ER routine
%   iter_num        \in IN                  :   total number of ER updates
%
% Outputs: 
%
%   psi_ortho       \in IC^{N2 x N1 x N3}   :   retrieved 3D exit-field in the real-space ORTHOGONAL frame
%   dp_error        \IR^{iter_num}          :   error metric values

DISPLAY = 1
DEBUG   = 1;
PAPIER  = 1;

[~,N1,N3]  = size(dp_sqrt);

r1          = geo_param(1,:);
r2          = geo_param(2,:);
r3          = geo_param(3,:);
q1          = geo_param(4,:);
q2          = geo_param(5,:);
q3          = geo_param(6,:);

d_r1        = r1(2)-r1(1);
d_r2        = r2(2)-r2(1);
d_r3        = r3(2)-r3(1);

if DISPLAY
    % initiate the reconstruction display
    figure('Position', [500,500,700,300], 'Name', 'Reconstruction', 'NumberTitle', 'off', 'Color', [1,1,1]);
    h1	= imagesc( r2, r3, squeeze(abs(psi_ortho(:,fix(N1/2)+1,:))),'Parent',subplot(121)); title('abs(\psi)'); xlabel('r2'), ylabel('r3'), axis image; axis xy; colorbar
    h2	= imagesc( r2, r3, squeeze(angle(psi_ortho(:,fix(N1/2)+1,:))),'Parent',subplot(122)); title('angle(\psi)'); xlabel('r2'), ylabel('r3'), axis image; axis xy; colorbar
    
    % initiate the error plot
    figure('Position',[1400,500,700,300], 'Name', 'Error Plot', 'NumberTitle', 'off', 'Color', [1,1,1]);
    err_axes = axes('tag', 'err_ax', 'XLim', [1,max(2,iter_num)],'YLim',[0,1]);
    title('log10(Error metric) plot'); xlabel('Iterations'); ylabel('log10(Error)'); grid on; box on
    err  = line(0,0,'Parent',gca,'Color','b','linewidth',2);
end

%% We sort the angle indices so that we update from lower to higher
%% frequency components
[~,index]   = sort(abs(q3));
telapsed    = zeros(1,N3);

%% Main ER/BP iteration...
for iter = 1 : iter_num
    disp(iter)    
    dp_error_n              = ones(1,length(q3));    
    h                       = waitbar(0,'Loop over the RC angles...');
    
    %% loop through the rocking angles  
    counter                 = 0;
    psi_ortho_new           = 0;    
    
    for n_q3 = index 
        
        waitbar(counter/N3,h,sprintf('angle : %d over %d',counter,N3))
        
        %% Current diffraction pattern
        sqrtI               = squeeze(dp_sqrt(:,:,n_q3));
        
        %% -- Computation of the current slice in the 3D expected far-field...
        tstart              = tic;        
        phase_ramp          = exp(sqrt(-1)*2*pi*q3(n_q3)*(R3 - tan(theta_B)*R2));
        psi_mod_integ       = cos(theta_B)*squeeze(sum(psi_ortho.*conj(phase_ramp),3))*d_r3;        
        PSI_non_ortho       = fftshift(fft2(fftshift(psi_mod_integ)))*d_r1*d_r3;
        telapsed(n_q3)      = toc(tstart);
        %%--
        
        %% fitting error for the 2D diffracted field in qx(n_qx)
        dp_error_n(1,n_q3)	= sum(sum(abs(sqrtI - abs(PSI_non_ortho)).^2));
                
        %% -- Contribution of the current slice in the update of the 3D expected far-field
        tstart              = tic;        
        PSI_non_ortho       = sqrtI.*exp(1i*angle(PSI_non_ortho));
        psi_mod_integ       = ifftshift(ifft2(ifftshift(PSI_non_ortho)))/(cos(theta_B)*d_r1*d_r2*d_r3);
        psi_ortho_new       = psi_ortho_new + repmat(psi_mod_integ, [1,1,length(q3)]).*phase_ramp/N3;        
        telapsed(n_q3)      = telapsed(n_q3) + toc(tstart);
        %%--                
        counter = counter +1;
    end
    
    close(h)
    
    %% apply the support constraint using ER algorithm
     psi_ortho           = psi_ortho - alpha*supp.*(psi_ortho - psi_ortho_new);
    
    dp_error(1,iter) = sum(sum(dp_error_n));
        
    if(DISPLAY)
        % update the display
        set(h1,'CData',squeeze(abs(psi_ortho(:,fix(N1/2)+1,:))));
        set(h2,'CData',squeeze(angle(psi_ortho(:,fix(N1/2)+1,:))));
        set(err_axes,'YLim',[min(log10(dp_error)) - 1,max(log10(dp_error))]);
        set(err,'XData',1:iter);
        set(err,'YData',log10(dp_error));
        drawnow
    end
    
    averageTime = sum(telapsed)
    
end

