function [f,g] = CBFwithGrad(theta,Y,Nsensor)
%% CBF
if size(theta,1)~=1, theta = theta.'; end
% SCM
Ryy = Y * Y' / size(Y,2);
% steering vector
avec   = exp(-1j*pi*(0:(Nsensor-1))'*sind(theta));
% steering vector d/d theta
avecD  = (-1j*pi*(0:(Nsensor-1))'*cosd(theta)) .* avec;
% Calculate objective f [ - to change finding PEAKS to min. problem]
% f      = -real(avec(theta)'*Ryy*avec(theta));
f      = -trace(real(avec'*Ryy*avec)) / Nsensor / Nsensor;

if nargout > 1 % gradient required
    g = zeros(numel(theta),1);
    for iVar = 1:numel(theta)
        dAdth         = zeros(Nsensor,numel(theta));
        dAdth(:,iVar) = avecD(:,iVar);
        % g(iVar)       = - real( trace( dAdth'*Ryy*avec(theta) + avec(theta)'*Ryy*dAdth ) );
        g(iVar)       = - 2*real( trace( dAdth'*Ryy*avec ) ) / Nsensor / Nsensor;
    end
end
end