function [sp] = mCBF(sensingMatrix,receivedSignal)
% Function for estimating the spatial spectrum P = a^H R a / a^H a
% This is called the conventional beamformer.
%
% sp - Estimated spatial spectrum
% a - Matrix of steering vectors

Ryy = receivedSignal * receivedSignal' / size(receivedSignal,2);

na = size(sensingMatrix, 2);
sp = zeros(na, 1);

for i = 1:na
    aa = sensingMatrix(:, i);
    sp(i) = aa'*Ryy*aa/(aa'*aa);
end

end