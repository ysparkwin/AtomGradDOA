function [f,g] = mCoFit(theta,Y,noisepower)
%% CoFit
if size(theta,1)~=1, theta = theta.'; end
Nsensor = size(Y,1);
% SCM
Syy        = Y * Y' / size(Y,2);
noisepower = max(1e-7*trace(Syy)/Nsensor,noisepower);
% steering vector
avec    = exp(-1j*pi*(0:(Nsensor-1))'*sind(theta));
% steering vector d/d theta
avecD   = (-1j*pi*(0:(Nsensor-1))'*cosd(theta)) .* avec;
% steering vector pseudo-inverse
avecP   = pinv(avec);
% A * A^+
M       = avec*avecP;

% S_y^RN
SyR       = ( M*(Syy - noisepower*eye(Nsensor))*M' + (noisepower)*eye(Nsensor) );

% Objective function
f       = norm( (Syy - SyR),'fro')^2;

if nargout > 1 % gradient required
    g = zeros(numel(theta),1);
    for iVar = 1:numel(theta)
        dMdthdAAp   = avecD(:,iVar) * avecP(iVar,:);
        dMdth       = (eye(Nsensor) - M) * dMdthdAAp ...
                      + dMdthdAAp' * (eye(Nsensor) - M);
        g(iVar)     = 4*real( trace( (SyR-Syy)'*dMdth*Syy*M' ) );
    end

    f = f * 1e8; % 1000 for step-sizing. You can comment this.
end
end