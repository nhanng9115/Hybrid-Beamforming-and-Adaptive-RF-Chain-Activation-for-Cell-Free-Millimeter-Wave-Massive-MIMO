function power = get_power(L,N,Nr,n_AP,nbar,RFA_scheme,K,Pt,sigma2,~)

Nt = 4;
P_LNA = 20e-3; P_M = 0.6e-3; P_PS = 30e-3; P_ADC = 200e-3; P_RF = 40e-3; P_SW = 5e-3;
% eta = 0.4; 
% Pt = 1000*K*rho*sigma2/eta; Pbt = 250; P0 = 825; 
% Pmin = P_RF + P_ADC + P_LNA + Nr*P_PS + P_M; Psleep = 0.667*Pmin;
BW = 100e6;
Pue = 0.1; Pfix = 0.825; eta = 0.3; alpha = 2; Tc = 2e-3; tau_p = K-1; tau_d = 180;  Pfhtotal = 50; Cfh = 100e6; 
P0 = K*Pt*sigma2/eta + K*Pue + L*Pfix;
Rfh = 2*Nt*K*tau_d*alpha/Tc;
Pfh = Pfhtotal/Cfh * Rfh;
P1 = P_LNA + 2*P_M; P2 = Nr*P_PS + P_RF + P_ADC;
switch RFA_scheme
    case 'ARFA'
        power = P0 + n_AP*Pfh + n_AP*Nr*P1 + L*nbar*P2;
    case 'fix_N'
        power = P0 + L*Pfh + L*Nr*P1 + L*N*P2;
    case 'fix_nbar'
        power = P0 + L*(Pfh + Nr*P1 + nbar*P2);
    case 'APS'
        power = P0 + L*nbar/N*(Pfh + Nr*P1 + N*P2);
    case 'AS'
        N_AS = N;
        power = P0 + L*(Pfh + Nr*P_SW + N_AS*(P_RF + P_ADC + P1));
    otherwise
        power = 0;
end

end %eof