function [psi_nonortho,dp_error] = ER(dp_sqrt,psi_nonortho, supp, geo_param, alpha, iter_num)

%
% [psi_nonortho,dp_error] = ER(dp_sqrt,psi_nonortho, supp, geo_param, alpha, iter_num)
%
% Description: 
%   Support-based, 3D phase-retrieval routine via the standard error-reduction
%   routine.  If the 3D Fourier components measurement are obtained in a
%   non-orthogonal geometry, then the sample will be retrieved in the
%   frame, conjugate with the non-orthogonal reciprocal-space frame, and it is 
%   not non-orthogonal.
%
% Inputs:
%
%   dp_sqrt         \in IN^{N2 x N1 x N3}   :   square root of the 3D intensity measurement stack (i.e., the sqrt of the intensity extracted along the RC)
%   psi_nonortho    \in IC^{N2 x N1 x N3}   :   initial-guess of the 3D retrieved field 
%   supp            \in IN^{N2 x N1 x N3}   :   binary map defining the support of the real-space exit-field IN 
%                                               THE FRAME CONJUGATE WITH THE (non-orthogonal) MEASUREMENT GEOMETRY (\bar{k}_1, \bar{k}_2, \bar{k}_3)
%   geo_param       \in IR^{2 X N1}         :   -FOR DISPLAY PURPOSE ONLY- the axis vector along ex and ez
%   alpha           \in IR_+                :   updating step-size for the ER routine
%   iter_num        \in IN                  :   total number of ER updates
%
% Outputs: 
%
%   psi_nonortho    \in IC^{N2 x N1 x N3}   :   retrieved 3D exit-field in the real-space NON-ORTHOGONAL frame
%   dp_error        \IR^{iter_num}          :   error metric values


DISPLAY = 1;

[~,N1,~]  = size(dp_sqrt);
r2          = geo_param(1,:);
r3          = geo_param(2,:);

if DISPLAY
    % initiate the reconstruction display
    figure('Position', [500,500,700,300], 'Name', 'Reconstruction', 'NumberTitle', 'off', 'Color', [1,1,1]);
    h1	= imagesc( r2, r3, squeeze(abs(psi_nonortho(:,fix(N1/2)+1,:))),'Parent',subplot(121)); title('abs(\psi)'); xlabel('r2'), ylabel('r3'), axis image; axis xy; colorbar
    h2	= imagesc( r2, r3, squeeze(angle(psi_nonortho(:,fix(N1/2)+1,:))),'Parent',subplot(122)); title('angle(\psi)'); xlabel('r2'), ylabel('r3'), axis image; axis xy; colorbar    
    
    % initiate the error plot
    figure('Position',[1400,500,700,300], 'Name', 'Error Plot', 'NumberTitle', 'off', 'Color', [1,1,1]);
    err_axes = axes('tag', 'err_ax', 'XLim', [1,iter_num],'YLim',[0,1]);
    title('log10(Error metric) plot'); xlabel('Iterations'); ylabel('log10(Error)'); grid on; box on
    err  = line(0,0,'Parent',gca,'Color','b','linewidth',2);
end

telapsed    = zeros(1,iter_num);

for iter = 1 : iter_num
       
    disp(iter)
    
    % forward calculation
    tstart                  = tic;    
    PSI_non_ortho           =  fftshift(fftn(fftshift(psi_nonortho)));
    telapsed(iter)          = toc(tstart);
    
    % Computation of the Error metric one aims at minimizing with ER
    dp_error(1,iter)        = sum((dp_sqrt(:)-abs(PSI_non_ortho(:))).^2);

    % replace the modulus and keep the phase
    PSI_non_ortho           = dp_sqrt.*exp(1i*angle(PSI_non_ortho));
    
    % backward calculation
    tstart                  = tic;   
    psi_nonortho_new        = ifftshift(ifftn(ifftshift(PSI_non_ortho))); 
    telapsed(iter)          = telapsed(iter) + toc(tstart);
    
    % apply the support constraint using ER algorithm
    psi_nonortho       = psi_nonortho - alpha*supp.*(psi_nonortho - psi_nonortho_new) ;
              
    if DISPLAY
        % update the display
        set(h1,'CData',squeeze(abs(psi_nonortho_new(:,fix(N1/2)+1,:))));
        set(h2,'CData',squeeze(angle(psi_nonortho_new(:,fix(N1/2)+1,:))));
        set(err_axes,'YLim',[min(log10(dp_error)) - 1,max(log10(dp_error))]);
        set(err,'XData',1:iter);
        set(err,'YData',log10(dp_error));
        drawnow                
    end
    
end

averageTime     = mean(telapsed)


