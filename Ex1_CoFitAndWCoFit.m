% Version 1.0: (11/11/2024)
% written by Yongsung Park

% Yongsung Park & Peter Gerstoft & Christoph F. Mecklenbräuker
% AOPE/WHOI
% MPL/SIO/UCSD
% yongsung.park@whoi.edu / gerstoft@ucsd.edu / christoph.mecklenbraeuker@tuwien.ac.at
% noiselab.ucsd.edu

% Citation
% Y. Park, P. Gerstoft, and C. F. Mecklenbräuker, “Atom-Constrained Gridless DOA Refinement With Wirtinger Gradients,” IEEE Open J. Signal Process. 5, 1134–1146 (2024).
% https://doi.org/10.1109/OJSP.2024.3496815

% P. Gerstoft and Y. Park, “Atom-Constrained Maximum Likelihood Gridless DOA with Wirtinger Gradients,” Proc. IEEE ICASSP (2025).
% https://doi.org/10.1109/ICASSP49660.2025.10889232

%%
clear; clc;
close all;

dbstop if error;

addpath([cd,'/_common'])
% addpath(['../_common'])

%% Noise model
model = 'Gaussian',                      nu_model_string = '';

%% Environment parameters
freq   =  2.0E+03;    % frequency (Hz)
c0     =  343;        % speed of sound (m/s) in dry air at 20°C
lambda =  c0/freq;    % wavelength (m)
wavenum=  2*pi/lambda;% wave number (rad/m)

%% Array configuration
Mlist=[181];
for mlist = 1:length(Mlist)
M = Mlist(mlist);
antenna_array.type = 'ULA';  % or 'UCA'
switch(antenna_array.type)
    case 'ULA' % definition of uniform linear array geometry
        theta  = 0.0;          % irrelevant elevation angle of incoming wave [degrees]
        antenna_array.N = 20;  % no. of sensors in ULA
        antenna_array.d = 0.5; % sensor spacing of ULA measured in wavelengths
        for n=0:(antenna_array.N-1)
            antenna_array.x(n+1) = n * antenna_array.d * lambda;
            antenna_array.y(n+1) = 0.0;
            antenna_array.z(n+1) = 0.0;
        end
        % array steering matrix of size N x M for all azimuth, the elevation is irrelevant.
        % M = 181;   % standard dictionary, 1 deg azimuth resolution
        % M = 18001; % high resolution dictionary, 0.01 deg azimuth resolution
        dphi=180/(M-1);
        phi_vec = [-90:dphi:90];
end

% Design/steering matrix (Sensing matrix)
sensingMatrix = zeros(antenna_array.N,M);
sensingMatrixD = zeros(antenna_array.N,M); %-- CRB-YP
for m=1:M
    kvec = wavenum * [sind(phi_vec(m))*cosd(theta);cosd(phi_vec(m))*cosd(theta);sind(theta)];
    kvecD = wavenum * [cosd(phi_vec(m))*cosd(theta);-sind(phi_vec(m))*cosd(theta);sind(theta)];
    sensingMatrix(:,m)   = exp(-1j * kvec.' * [antenna_array.x;antenna_array.y;antenna_array.z]); % normalization to |a_n|=1 or ||a||_2^2 = N.
    sensingMatrixD(:,m)  = (-1j * kvecD.' * [antenna_array.x;antenna_array.y;antenna_array.z])...
        .* exp(-1j * kvec.' * [antenna_array.x;antenna_array.y;antenna_array.z]); % normalization to |a_n|=1 or ||a||_2^2 = N.
end

%% Number of sensors / grid-points / snapshots
Nsensor     = antenna_array.N;  % number of sensors
Ntheta      = M;                % number of angular-search grid
Nsnapshot   = 25;               % number of snapshots

%% Simulation parameters
% noise standard deviation sigma

SNRs = 36:-3:-6;

% number of sources
Number_of_DOAs = 1;

NmonteCarlo = 250;
LSnapshot = Nsnapshot;
% LSnapshot = Nsnapshot * NmonteCarlo; % Number of array data vector observations "Large"


%% loop over various SNR levels

% for isnr = 1 %1:length(SNRs)
for isnr = 1:length(SNRs)

%     rng('default') % YP: We need to be careful where to set 'rng'
%     rng(1,'twister')

    SNR  = SNRs(isnr);

% evaluate SBL
    options = SBLSet;
    options.Nsource = Number_of_DOAs + 0;
%     options.gamma_range=10^-3;

    errorDOAseparation = 1; % [deg.]
    errorDOAsepP = floor(errorDOAseparation/dphi) - 1;
    errorDOApeak = Number_of_DOAs + 2;
    errCut = 10; % Maximum RMSE cut-off.

% obtain active indices --------------------------%
    options.activeIndices = 1;
    options.activeIndRepN = 10;
    options.convergence.min_iteration = options.activeIndRepN;
%-------------------------------------------------%

%     for n_monteCarlo = 1
    for n_monteCarlo = 1:NmonteCarlo  % parfor loop over snapshot realizations

rng(n_monteCarlo,'twister')

        disp(' ')
        disp(['SNR',num2str(SNR),'#Sim : ',num2str(n_monteCarlo)])

        % number of sources
        x_src   = ones(Number_of_DOAs,1);
        switch(Number_of_DOAs)
            case 1
                %         DOA_src = (asind(gen_bs(sind(-75), sind(75), Number_of_DOAs, asin(2/Nsensor))));
                DOA_src = -10 + dphi*rand - dphi/2;
            case 2
                DOA_src = [-10 + dphi*rand - dphi/2; 10 + dphi*rand - dphi/2];
            case 3
                DOA_src = [-3 + dphi*rand - dphi/2; 2 + dphi*rand - dphi/2; 75 + dphi*rand - dphi/2];
            otherwise
                error('this Number_of_DOAs is not implemented')
        end

        DOA_MC(:,n_monteCarlo) = DOA_src;

%         for k=1:Number_of_DOAs
%             [~,m_src(k)] = min(abs(phi_vec - DOA_src(k)));
%         end
%         a_src = sensingMatrix(:,m_src);
%         SNRmaxTmp = -10*log10((1-diag(abs(a_src'*sensingMatrix(:,m_src+1)))/Nsensor)*2);
%         if exist('SNRmax','var')==0, SNRmax = []; end
%         SNRmax = [SNRmax;SNRmaxTmp];

        % Steering vectors for true sources
%         for k=1:Number_of_DOAs
%             m_src(k) = find(phi_vec == DOA_src(k));
%         end
%         a_src = sensingMatrix(:,m_src);
        a_src  = zeros(antenna_array.N,Number_of_DOAs);
        a_srcD = zeros(antenna_array.N,Number_of_DOAs);
        for k=1:Number_of_DOAs
            kvec = wavenum * [sind(DOA_src(k))*cosd(theta);cosd(DOA_src(k))*cosd(theta);sind(theta)];
            kvecD = wavenum * [cosd(DOA_src(k))*cosd(theta);-sind(DOA_src(k))*cosd(theta);sind(theta)];
            a_src(:,k)   = exp(-1j * kvec.' * [antenna_array.x;antenna_array.y;antenna_array.z]); % normalization to |a_n|=1 or ||a||_2^2 = N.
            a_srcD(:,k)  = (-1j * kvecD.' * [antenna_array.x;antenna_array.y;antenna_array.z])...
                        .* exp(-1j * kvec.' * [antenna_array.x;antenna_array.y;antenna_array.z]); % normalization to |a_n|=1 or ||a||_2^2 = N.
        end

        % Noise modeling
        sigma = 1 * norm(a_src*x_src,'fro') / (10^(SNR/20));
        %     SNR_gen = 10*log10(norm(a_src*x_src,'fro')^2 ./ (sigma.^2));
        %     % check the generated SNR

        switch(model)
            case 'Laplace-like'
                [y,xAmp] = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                    sigma,model);
            case 'Gaussian'
                [y,xAmp] = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                    sigma,model);
            case 'epscont'
                [y,xAmp] = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                    sigma,model,epsilon_model,lambda_model);
            case 'Complex-Student' % method in Ollila & Koivunen, PIMRC 2003
                [y,xAmp] = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                    sigma,model,nu_model);
            case 'Heteroscedastic'
                [y,xAmp] = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                    sigma,model);
            otherwise
                error(['unknown model ', model]);
        end

        Y = y;
%         Y = y(:,(n_monteCarlo-1)*Nsnapshot+(1:Nsnapshot));

%% CRB-YP Van Trees Book Eq.(8.106) & (8.110)
        XAMP  = xAmp;
%         XAMP  = xAmp(:,(n_monteCarlo-1)*Nsnapshot+(1:Nsnapshot));
        vanTreeV  = a_src;
        vanTreeD  = a_srcD;          % D Eq.(8.100)

        vanTreeSf = diag(diag(XAMP*XAMP'/Nsnapshot)); % S_f
        Pn   = sigma^2;
        
        % H Eq.(8.101) where P_V Eq.(8.96)
        H = vanTreeD'...
            *(eye(Nsensor) - vanTreeV/(vanTreeV'*vanTreeV)*vanTreeV')...
            *vanTreeD;

        % Eq.(8.110)
        CRBa = real(H .* (vanTreeSf.'));
        CRBa = eye(size(XAMP,1)) / CRBa * (Pn / Nsnapshot / 2);

        if exist('outputsCRBa','var')==0, outputsCRBa = []; end
        outputsCRBa = [outputsCRBa;mean(diag(CRBa))];
        
        % Eq.(8.106)
        CRBaux1 = vanTreeV' * vanTreeV * (vanTreeSf / Pn);
        CRBaux2 = eye(size(XAMP,1)) / ( eye(size(XAMP,1)) + CRBaux1 );
        CRB = real( vanTreeSf * (CRBaux2 * CRBaux1) .* (H.') );
        CRB = eye(size(XAMP,1)) / CRB * (Pn / Nsnapshot / 2);

        if exist('outputsCRBd','var')==0, outputsCRBd = []; end
        outputsCRBd = [outputsCRBd;mean(diag(CRB))];

%% CBF w/ fminunc w/ trust-region (Analytic gradient)
        [Pcbf] = abs(mCBF(sensingMatrix,Y));
        [~, Ilocs] = findpeaks(Pcbf,'SORTSTR','descend','Npeaks', Number_of_DOAs);
% CBF w/ fminunc
        % options = optimoptions('fminunc','Algorithm','quasi-newton','Display','iter');
        optionsCBF = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'Display','off');

        theta0 = phi_vec(Ilocs);
        thetaEst = fminunc(@(theta)CBFwithGrad(theta,Y,Nsensor),theta0,optionsCBF);

        DoA_error_cbf = errorDOAcutoff(thetaEst,DOA_src,10);
        disp(['RMSE CBF fminunc      : ',num2str(sqrt(mean(power(DoA_error_cbf,2))))])

        [f,g] = CBFwithGrad(thetaEst,Y,Nsensor);
        disp(['obj.     :                           ',num2str(-f)])
        disp(['gradient :                           ',num2str(-g.')])

        if exist('outputscbf','var')==0, outputscbf = []; end
        outputcbf = struct('theta',thetaEst,'error',DoA_error_cbf);
        outputscbf = [outputscbf; outputcbf];

%% SBL Initialization
        loss_param = inf;
        [gammaInd81,report81] = SBL_v5p12(sensingMatrix, Y, 'SBL-G', loss_param, options, errorDOApeak, errorDOAsepP);

%% Results
        optionsOPT = optimset('Display','off');
%         optionsOPT = optimset('TolX',1e-15, 'TolFun',1e-6, 'Display','off');

        % SBL model to post-process, gamma 81: Gauss / 82: Huber / 83: MVT / 84: Tyler
        for modelN = 81%[81 82 83 84]
        evalText = ['gammaModel = gammaInd',num2str(modelN),';']; eval(evalText);
        evalText = ['reportModel = report',num2str(modelN),';']; eval(evalText);

        % DOA estimates separation: 2 deg.
        [~,gammaModel]    = findpeaks(reportModel.results.final_iteration.gamma,'SortStr','descend','NPeaks',errorDOApeak,'MinPeakDistance',2/dphi);
        [~,gammaModelInd] = sort(SBL_v4(sensingMatrix(:,gammaModel),Y,options),'descend');
        gammaModel        = gammaModel(gammaModelInd(1:Number_of_DOAs));
        % gammaModel = gammaModel(1:Number_of_DOAs);


% starting point
theta0     = [phi_vec(gammaModel(1:Number_of_DOAs))];
noisepower = reportModel.results.final_iteration.noisepower;

disp([' '])

disp(['Covariance fitting'])
disp([' '])
%%        % Method : CoFit w/ fminunc w/ quasi-Newton (Numerical gradient)
        % noisepower = reportModel.results.final_iteration.noisepower;
        % noisepower = 1e-1;

        optionsOPTunc = optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
        % optionsOPTunc = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'Display','iter');

        DOAsbl7 = fminunc(@(theta)mCoFit(theta,Y,noisepower),theta0,optionsOPTunc);
        DoA_error_sbl7 = errorDOAcutoff(DOAsbl7,DOA_src,errCut);
        disp(['RMSE CoFit  grad-Num : ',num2str(sqrt(mean(power(DoA_error_sbl7,2))))])

        [f,g] = mCoFit(DOAsbl7,Y,noisepower);
        disp(['obj.     :                           ',num2str(f)])
        disp(['gradient :                           ',num2str(g.')])

%%        % Method : COVfit w/ fminunc w/ trust-region (Analytic gradient)
        % noisepower = reportModel.results.final_iteration.noisepower;
        % noisepower = 1e-1;

        % optionsOPTunc = optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
        optionsOPTunc = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'Display','off');

        DOAsbl8 = fminunc(@(theta)mCoFit(theta,Y,noisepower),theta0,optionsOPTunc);
        DoA_error_sbl8 = errorDOAcutoff(DOAsbl8,DOA_src,errCut);
        disp(['RMSE CoFit  grad-A   : ',num2str(sqrt(mean(power(DoA_error_sbl8,2))))])

        [f,g] = mCoFit(DOAsbl8,Y,noisepower);
        disp(['obj.     :                           ',num2str(f)])
        disp(['gradient :                           ',num2str(g.')])

        [f,g] = mCoFit(DOA_src.',Y,noisepower);
        disp(['If we compute obj. and gradient at the TRUE DOA'])
        disp(['obj.     :                           ',num2str(f)])
        disp(['gradient :                           ',num2str(g.')])

%%        % Method : WCoFit w/ fminunc w/ quasi-Newton (Numerical gradient)
        % noisepower = reportModel.results.final_iteration.noisepower;
        % noisepower = 1e-1;

        optionsOPTunc = optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
        % optionsOPTunc = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'Display','iter');

        DOAsbl5 = fminunc(@(theta)mWCoFit(theta,Y,noisepower),theta0,optionsOPTunc);
        DoA_error_sbl5 = errorDOAcutoff(DOAsbl5,DOA_src,errCut);
        disp(['RMSE WCoFit grad-Num : ',num2str(sqrt(mean(power(DoA_error_sbl5,2))))])

        [f,g] = mWCoFit(DOAsbl5,Y,noisepower);
        disp(['obj.     :                           ',num2str(f)])
        disp(['gradient :                           ',num2str(g.')])

%%        % Method : WCOVfit w/ fminunc w/ trust-region (Analytic gradient)
        % noisepower = reportModel.results.final_iteration.noisepower;
        % noisepower = 1e-1;

        % optionsOPTunc = optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
        optionsOPTunc = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'Display','off');

        DOAsbl6 = fminunc(@(theta)mWCoFit(theta,Y,noisepower),theta0,optionsOPTunc);
        DoA_error_sbl6 = errorDOAcutoff(DOAsbl6,DOA_src,errCut);
        disp(['RMSE WCoFit grad-A   : ',num2str(sqrt(mean(power(DoA_error_sbl6,2))))])

        [f,g] = mWCoFit(DOAsbl6,Y,noisepower);
        disp(['obj.     :                           ',num2str(f)])
        disp(['gradient :                           ',num2str(g.')])

        [f,g] = mWCoFit(DOA_src.',Y,noisepower);
        disp(['If we compute obj. and gradient at the TRUE DOA'])
        disp(['obj.     :                           ',num2str(f)])
        disp(['gradient :                           ',num2str(g.')])

%%      save
        text1 = ['if exist(','''outputsSBL',num2str(modelN),'''',',','''var''',')==0, outputsSBL',num2str(modelN),' = []; end'];
        eval(text1)
        saveCharVar = who('DoA_error_sbl*');
        text2 = ['outputSBL',num2str(modelN),' = struct(',''''];
        for iVar = 1:numel(saveCharVar)
            text2 = [text2,'error',num2str(iVar),''''...
                ,',',saveCharVar{iVar},',',''''];
        end
        text2(end-1:end) = []; text2 = [text2,');'];
        eval(text2)
        text3 = ['outputsSBL',num2str(modelN),' = [outputsSBL',num2str(modelN),'; outputSBL',num2str(modelN),'];'];
        eval(text3)

        end
    end % end of the for-loop
    saveCharVar = who('outputs*');
    saveChar = ['save([ ''p12Prand_'', model(1), ''mode_'', ''s'', num2str(Number_of_DOAs), ''MC'' , num2str(NmonteCarlo) , ''SNRn'' , num2str(isnr), ''g'', num2str(M)], ''SNRs'' , ''NmonteCarlo'' '];
    for ichar = 1:numel(saveCharVar)
        saveChar = [saveChar,',''',char(saveCharVar{ichar}),''''];
    end
    saveChar = [saveChar,');'];
    eval(saveChar)

    if isnr > 1
        delete( [ 'p12Prand_', model(1), 'mode_', 's', num2str(Number_of_DOAs), 'MC' , num2str(NmonteCarlo) , 'SNRn' , num2str(isnr-1), 'g', num2str(M), '.mat' ] )
    end

end % end of for isnr=1:length(sigma_vec) loop

if mlist > 1
    delete( [ 'p12Prand_', model(1), 'mode_', 's', num2str(Number_of_DOAs), 'MC' , num2str(NmonteCarlo) , 'SNRn' , num2str(isnr), 'g', num2str(Mlist(mlist-1)), '.mat' ] )
end

end

%% Plot
ss = Number_of_DOAs;

plotColor = lines(6);

saveCharVar = who('outputs*');
saveCharVar(2) = [];

figure(1);

for n_output = [1 2 3]

    dataLoadchar = ['dataLoad = ',char(saveCharVar(n_output)),';'];
    eval(dataLoadchar);

    Nsnrs= numel(SNRs);
    MInd = numel(outputsSBL81) / numel(SNRs) / NmonteCarlo;

    for mInd=1:MInd
        for ind=1:Nsnrs
            if      n_output == 1
                rmseSNR(ind) = sqrt( mean( dataLoad( ...
                    1+ (ind-1)*NmonteCarlo+ (mInd-1)*NmonteCarlo*Nsnrs:NmonteCarlo+ (ind-1)*NmonteCarlo+ (mInd-1)*NmonteCarlo*Nsnrs ) ...
                    )*180/pi*180/pi );
            elseif  n_output == 2
                totET1 = [];
                totET2 = [];
                for index=1:NmonteCarlo
                    totET1 = [totET1;dataLoad((ind-1)*NmonteCarlo+index+ (mInd-1)*NmonteCarlo*Nsnrs).error2];
                    totET2 = [totET2;dataLoad((ind-1)*NmonteCarlo+index+ (mInd-1)*NmonteCarlo*Nsnrs).error4];
                end
                Nout = 0.0; % Portion of Outliers, (ignore)
                totET1 = sort(abs(totET1));
                totET2 = sort(abs(totET2));
                rmseSNR1(ind) = sqrt(mean(power(totET1(1:length(totET1)-floor(length(totET1)*Nout)),2)))+10^-10;
                rmseSNR2(ind) = sqrt(mean(power(totET2(1:length(totET2)-floor(length(totET2)*Nout)),2)))+10^-10;
            else
                totET1 = [];
                for index=1:NmonteCarlo
                    totET1 = [totET1;dataLoad((ind-1)*NmonteCarlo+index+ (mInd-1)*NmonteCarlo*Nsnrs).error];
                end
                Nout = 0.0; % Portion of Outliers, (ignore)
                totET1 = sort(abs(totET1)); totET2 = sort(abs(totET2));
                rmseSNR1(ind) = sqrt(mean(power(totET1(1:length(totET1)-floor(length(totET1)*Nout)),2)))+10^-10;
            end
        end

        pcolor = lines;
        if      n_output == 1 && mInd == 1
            hold on; plot(SNRs,rmseSNR,'k','linewidth',1);
        elseif  n_output == 2 && mInd == 1
            hold on; plot(SNRs,rmseSNR1,'+-','Color',pcolor(2,:),'linewidth',1.5,'MarkerSize',8);
            hold on; plot(SNRs,rmseSNR2,'x','Color',pcolor(3,:),'linewidth',1.5,'LineStyle','--','MarkerSize',8);
            legend('CRB','CoFit','WCoFit','Interpreter','latex','FontSize',14);

        elseif  n_output == 3 && mInd == 1 && ss == 1
            hold on; plot(SNRs,rmseSNR1,'bs:','linewidth',1.5);
            legend('CRB','CoFit','WCoFit','CBF','Interpreter','latex','FontSize',14);
        end

    end
end

box on; grid on;
axis([-6 36 7e-3 15])
set(gca,'fontsize',16,'TickLabelInterpreter','latex','yscale','log')
set(gca,'YTick',[1e-2 1e-1 1e0 1e1])

xlabel('ASNR~[dB]','interpreter','latex');
ylabel('RMSE~[$^\circ$]','interpreter','latex');

%%
rmpath([cd,'/_common'])
% rmpath(['../_common'])
%% End------------------------------------------------------------------------------------------------------------------------ %%


%% Signal generation
function [receivedSignal,s_src] = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
    sigma,model,model_param1,model_param2)
% function to generate sensor observations
if (strcmpi(model,'Gaussian') || isempty(model) )
    noise_realization = sigma * complex(randn(Nsensor,LSnapshot),randn(Nsensor,LSnapshot))/sqrt(2);
    if 0
        % deterministric source
        receivedSignal = (a_src * x_src * ones(1,LSnapshot)) + noise_realization;
    else
        % stochastic source
        s_src = x_src .* complex(randn(Number_of_DOAs,LSnapshot),randn(Number_of_DOAs,LSnapshot))/sqrt(2 * Number_of_DOAs);
        receivedSignal = ( a_src * s_src ) + noise_realization;
    end
else
    error(['please specify noise model as a string equal to Laplace-like, ...' ...
        'Gaussian, epscont, Complex-Student or Heteroscedastic\n']);
end
end


