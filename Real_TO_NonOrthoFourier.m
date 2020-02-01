function PSI_non_ortho           = Real_TO_NonOrthoFourier(psi_ortho,phase_ramp,param)

% PSI_non_ortho           = transform_OrthoReal_TO_NonOrthoRecip(psi_ortho,phase_ramp,param)
%
% Description: 
%   Provides the NON-ORTHOGONAL representation of the 3D far field \bar{\Psi} FROM the
%   ORTHOGONAL representation of the real-space exit-field \psi.
%
% Inputs:
%
%   psi_ortho       \in IC^{N2 x N1 x N3}   :   3D exit-field in the orthogonal frame (noted \psi)
%   phase_ramp      \in IN^{N2 x N1 x N3}   :   3D matrix involved in the  mapping \psi -> \bar{\Psi}
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
XI_ortho_rx         = fftshift(fft(psi_ortho,[],3),3) .* phase_ramp; % Phase ramp
PSI_non_ortho       = fftshift(fft(fftshift(fft(XI_ortho_rx,[],1),1),[],2),2)*detB*d_r1*d_r2*d_r3;