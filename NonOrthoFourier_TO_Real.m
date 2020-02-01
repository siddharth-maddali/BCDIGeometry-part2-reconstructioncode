function psi_ortho           = NonOrthoFourier_TO_Real(PSI_non_ortho,phase_ramp,param)

% psi_ortho             = transform_NonOrthoRecip_TO_OrthoReal(PSI_non_ortho,phase_ramp,param)
%
% Description: 
%   Provides the ORTHOGONAL representation of the 3D exit field (noted \psi) FROM the
%   NON-ORTHOGONAL representation of the 3D far-field (noted \bar{\Psi}).
%
% Inputs:
%
%   psi_ortho       \in IC^{N2 x N1 x N3}   :   3D exit-field in the orthogonal frame (noted \psi)
%   phase_ramp      \in IN^{N2 x N1 x N3}   :   3D matrix involved in the  mapping \bar{\Psi} -> \psi 
%   param           \in IR^{4}              :   geometrical parameters involved in the transformation 
%
% Output:
%
%   PSI_non_ortho   \in IC^{N2 x N1 x N3}   :   3D far-field in the NON-orthogonal frame (noted \Psi)

d_r1                = param(1); 
d_r2                = param(2);
d_r3                = param(3);
theta_B             = param(4);
detB                = cos(theta_B);
zeta_ortho_rx       = fftshift(ifft(ifft(PSI_non_ortho,[],1),[],2),3) .* conj(phase_ramp); % Phase ramp
psi_ortho           = ifft(ifftshift(zeta_ortho_rx,3),[],3)/(detB *d_r1*d_r2*d_r3); % ifft along the rows...
