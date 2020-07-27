%% ------------------------------------------------------------------------
% Estimation.m
% 
% -------------------------------------------------------------------------

%% Reset Matlab

clear all; close all; clc;


%% Select the manager

imgr        = 90457;  % 23000: DFA. 90457: Vanguard.


%% Load data

mdata       = xlsread('DataEstimation.xlsx');
mdata       = mdata(mdata(:,1)==imgr,:);
vrweight    = mdata(:,3);


%% Linear IV

% Structure the data

vLNrweight  = log(vrweight);

vloc        = find(vrweight>0); 
mdata0      = mdata(vloc,:);
mchars      = mdata0(:,4:8);
vLNme       = mdata0(:,9);
vLNmeIV     = mdata0(:,10);
vLNrweight  = vLNrweight(vloc);

dT          = length(vLNme);
vones       = ones(dT,1);

mX          = [vLNme mchars vones];
mZ          = [vLNmeIV mchars vones];

% Unconstrained estimation

vb_linearIV     = inv(mZ' * mX) * mZ' * vLNrweight;

% Constrained estimation

if vb_linearIV(1) > .99
    
    mX          = [mchars vones];
    vb_linearIV = inv(mX' * mX) * mX' * (vLNrweight - .99*vLNme);
    vb_linearIV = [.99; vb_linearIV];
    
end


%% Non-linear IV

% Structure the data

mchars          = mdata(:,4:8);
vLNme           = mdata(:,9);
vLNmeIV         = mdata(:,10);

[dT,dN]         = size(mchars);
dN              = dN+2;
vones           = ones(dT,1);

mX              = [vLNme mchars vones];
mZ              = [vLNmeIV mchars vones];

% Unconstrained estimation

vb_nonlinearIV  = vb_linearIV;
dchange         = 1;

while dchange > 1E-4
    
    vlatent                 = vrweight .* exp(-mX * vb_nonlinearIV);
    mZ_tilde                = (vlatent * ones(1,dN)) .* mZ;

    vb_nonlinearIV_new      = vb_nonlinearIV + inv(mZ_tilde' * mX) * mZ' * (vlatent - 1);
    
    dchange                 = max(abs(vb_nonlinearIV - vb_nonlinearIV_new));    
    
    vb_nonlinearIV          = vb_nonlinearIV_new;
    
end

% Constrained estimation

if vb_nonlinearIV(1) > .99
    
    vb_nonlinearIV  = vb_linearIV(2:end);
    dchange         = 1;
    mX              = [mchars vones];
    dN              = dN-1;
    
    while dchange > 1E-4 
    
        vlatent                 = vrweight .* exp(-.99 * vLNme - mX * vb_nonlinearIV);
        mX_tilde                = (vlatent * ones(1,dN)) .* mX;

        vb_nonlinearIV_new      = vb_nonlinearIV + inv(mX_tilde' * mX) * mX' * (vlatent - 1);

        dchange                 = max(abs(vb_nonlinearIV - vb_nonlinearIV_new));    

        vb_nonlinearIV          = vb_nonlinearIV_new; 
        
    end
    
    vb_nonlinearIV  = [.99;vb_nonlinearIV];
    
end



%% Display estimates

disp([vb_linearIV vb_nonlinearIV])



























































