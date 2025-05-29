function [f,g] = mWCoFit(theta,Y,noisepower)
%% WCoFit
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
SigmayInv = pinv(SyR);

% Objective function
f       = norm( sqrtm(SigmayInv) * (Syy - SyR),'fro')^2;

if nargout > 1 % gradient required
    g = zeros(numel(theta),1);
    for iVar = 1:numel(theta)
        dMdthdAAp   = avecD(:,iVar) * avecP(iVar,:);
        dMdth       = (eye(Nsensor) - M) * dMdthdAAp ...
                      + dMdthdAAp' * (eye(Nsensor) - M);
        % dSyRdth     = dMdth * Syy * M' + M * Syy * dMdth';
        % g(iVar)     = trace( (eye(Nsensor)-SigmayInv*Syy*Syy*SigmayInv)*dSyRdth );
        g(iVar)     = 2*real( trace( (eye(Nsensor)-SigmayInv*Syy*Syy*SigmayInv)*dMdth*Syy*M' ) );
    end

    f = f * 1000; % 1000 for step-sizing. You can comment this.
end
end