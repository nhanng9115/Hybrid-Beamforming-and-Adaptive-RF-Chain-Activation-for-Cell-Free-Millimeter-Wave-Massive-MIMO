clc;
clear;
close all;
warning off
%% simulation setup
Nr = 64; Nt = 4; K = 8; N = 8; L0 = 32; P0 = 50; BW = 100; nbar0 = 2;
% sim_mode = 1;
for sim_mode = 5
    switch sim_mode
        case 1
            disp('SE and EE vs Pt......')
            file_name = strcat(num2str(Nr),'x',num2str(Nt),'x',num2str(L0),'x',num2str(K));
            load(file_name,'H_reli','H_all','LS_all','Ph_all','rx_array_all');
            x_vec = 10:5:40; Pt_dBm_vec = x_vec; L_vec = L0*ones(length(x_vec),1); nbar_vec = nbar0*ones(length(x_vec),1);
            xaxis = 'Pt';
        case 2
            disp('SE and EE vs L......')
            x_vec = [8,16,32,48,64]; L_vec = x_vec; Pt_dBm_vec = P0*ones(length(x_vec),1); nbar_vec = nbar0*ones(length(x_vec),1);
            xaxis = 'L';
        case 3
            disp('SE and EE vs nbar......')
            file_name = strcat(num2str(Nr),'x',num2str(Nt),'x',num2str(L0),'x',num2str(K));
            load(file_name,'H_reli','H_all','LS_all','Ph_all','rx_array_all');
            x_vec = [1:8]; nbar_vec = x_vec; L_vec = L0*ones(length(x_vec),1); Pt_dBm_vec = P0*ones(length(x_vec),1);
            xaxis = 'n_bar';
        case 4
            disp('SE and EE vs K......')
            x_vec = [2,4,8,12,16]; L_vec = L0*ones(length(x_vec),1); Pt_dBm_vec = P0*ones(length(x_vec),1); nbar_vec = nbar0*ones(length(x_vec),1);
            xaxis = 'K';
        case 5
            disp('SE and EE vs n_path......')
            x_vec = [1,2,5,10,15,20]; L_vec = L0*ones(length(x_vec),1); Pt_dBm_vec = P0*ones(length(x_vec),1); nbar_vec = nbar0*ones(length(x_vec),1);
            xaxis = 'n_path';
        otherwise
            error('Choose simulation mode......')
    end
    
    %% Simulation initialization
    BF_scheme = [1,...%ERFA_UB
        1,...%ERFA_C
        1,...%GS-RFS
        1,...%SC-RFS
        1,...%LS-RFS
        1,... %APS
        1];% AS
    
    ITER = 20;
    
    R_Fixed_N = zeros(ITER,length(x_vec));
    R_Fixed_nbar = zeros(ITER,length(x_vec));
    R_ARFA_FS = zeros(ITER,length(x_vec));
    R_ARFA_SC = zeros(ITER,length(x_vec));
    R_ARFA_LS = zeros(ITER,length(x_vec));
    R_APS = zeros(ITER,length(x_vec));
    R_AS = zeros(ITER,length(x_vec));
    
    EE_Fixed_N = zeros(ITER,length(x_vec));
    EE_Fixed_nbar = zeros(ITER,length(x_vec));
    EE_ARFA_FS = zeros(ITER,length(x_vec));
    EE_ARFA_SC = zeros(ITER,length(x_vec));
    EE_ARFA_LS = zeros(ITER,length(x_vec));
    EE_APS = zeros(ITER,length(x_vec));
    EE_AS = zeros(ITER,length(x_vec));
    
    P_Fixed_N = zeros(ITER,length(x_vec));
    P_Fixed_nbar = zeros(ITER,length(x_vec));
    P_ARFA_FS = zeros(ITER,length(x_vec));
    P_ARFA_SC = zeros(ITER,length(x_vec));
    P_ARFA_LS = zeros(ITER,length(x_vec));
    P_APS = zeros(ITER,length(x_vec));
    P_AS = zeros(ITER,length(x_vec));
    
    %% Start simulation
    var_e = 0.01;
    for ss = 1:length(x_vec)
        % Load simulation data
        if sim_mode == 2
            file_name = strcat(num2str(Nr),'x',num2str(Nt),'x',num2str(x_vec(ss)),'x',num2str(K));
            load(file_name,'H_all','LS_all','Ph_all','rx_array_all');
        elseif sim_mode == 4
            K = x_vec(ss);
            file_name = strcat(num2str(Nr),'x',num2str(Nt),'x',num2str(L0),'x',num2str(K));
            load(file_name,'H_all','LS_all','Ph_all','rx_array_all');
        elseif sim_mode == 5
            n_path = x_vec(ss);
            file_name = strcat(num2str(Nr),'x',num2str(Nt),'x',num2str(L0),'x',num2str(K),'x',num2str(n_path));
            load(file_name,'H_all','LS_all','Ph_all','rx_array_all');
        end
        disp(strcat('x = ', num2str(x_vec(ss))))
        
        Pt_dBm = Pt_dBm_vec(ss); rho = db2pow(Pt_dBm - 30); sigma2 = db2pow(-87 - 30); %rho = db2pow(Pt_dBm)/db2pow(-87);
        L = L_vec(ss);
        nbar = nbar_vec(ss);
        
        parfor ii = 1:ITER
            % Load channel relizations
            LS = LS_all(:,ii); Ph = Ph_all(:,ii); a_Rx = rx_array_all(:,:,:,ii);
            Htmp = H_all(:,:,:,:,ii);
            evar_scaled = mean(LS) * var_e;
            
            H1 = []; H = zeros(Nr,K*Nt,L);
            for ll = 1:L
                H_l = [];
                for kk = 1:K
                    E_kl = sqrt(evar_scaled/(K*Nt)) * sqrt(1/2) * (randn(Nr,Nt) + 1i*randn(Nr,Nt));
                    H_kl = Htmp(:,:,kk,ll) - E_kl;
                    H_l = cat(2, H_l, H_kl);
                end
                H1 = cat(1,H1,H_l);
                H(:,:,ll) = H_l;
            end
            gamma = rho/(sigma2 + rho*Nt*evar_scaled);
            
            %% Start stimulation
            if BF_scheme(1)
                tic
                [F,W] = HBF_D(H,Nr,Nt,N*ones(L,1),N,L,K,gamma); X = F*W;
                R_Fixed_N(ii,ss) = log2(det(eye(L*K*Nt) + gamma*pinv(X)*(H1*H1')*(X)));
                power = get_power(L,N,Nr,0,0,'fix_N',K,rho,sigma2,R_Fixed_N(ii,ss));
                EE_Fixed_N(ii,ss) = R_Fixed_N(ii,ss)/power;
                P_Fixed_N(ii,ss) = power;
                t_Fixed_N(ii,ss) = toc;
            end
            
            if BF_scheme(2)
                tic
                %nbar = N;
                [F,W] = HBF_D(H,Nr,Nt,nbar*ones(L,1),nbar,L,K,gamma); X = F*W;
                R_Fixed_nbar(ii,ss) = log2(det(eye(L*K*Nt) + gamma*pinv(X)*(H1*H1')*(X)));
                power = get_power(L,N,Nr,0,nbar,'fix_nbar',K,rho,sigma2,R_Fixed_nbar(ii,ss));
                EE_Fixed_nbar(ii,ss) = R_Fixed_nbar(ii,ss)/power;
                P_Fixed_nbar(ii,ss) = power;
                t_Fixed_nbar(ii,ss) = toc;
            end
            
            if BF_scheme(3)
                tic
                [R_ARFA_FS(ii,ss),EE_ARFA_FS(ii,ss),P_ARFA_FS(ii,ss)] = ARFA_FS(H1,H,Nr,Nt,nbar,L,K,gamma,N,rho,sigma2);
                t_ARFA_FS(ii,ss) = toc;
            end
            
            if BF_scheme(4)
                tic
                [R_ARFA_SC(ii,ss),EE_ARFA_SC(ii,ss),P_ARFA_SC(ii,ss)] = ARFA_SV(H1,H,Nr,Nt,nbar,L,K,gamma,N,rho,sigma2);
                t_ARFA_SC(ii,ss) = toc;
            end
            
            if BF_scheme(5)
                tic
                [R_ARFA_LS(ii,ss),EE_ARFA_LS(ii,ss),P_ARFA_LS(ii,ss)] = ARFA_PL(H1,H,Nr,Nt,nbar,L,K,gamma,LS,N,rho,sigma2);
                t_ARFA_LS(ii,ss) = toc;
            end
            
            if BF_scheme(6)
                tic
                [R_APS(ii,ss),EE_APS(ii,ss),P_APS(ii,ss)] = APS(H1,H,Nr,Nt,nbar,L,K,gamma,Ph,N,rho,sigma2);
                t_APS(ii,ss) = toc;
            end
            
            if BF_scheme(7)
                tic
                N_AS = 32;
                [R_AS(ii,ss),EE_AS(ii,ss),P_AS(ii,ss)] = AS(H1,H,Nr,Nt,nbar,L,K,gamma,N_AS,rho,sigma2);
                t_AS(ii,ss) = toc;
            end
        end
    end
    
    %% rate
    R_Fixed_N = mean(R_Fixed_N,1).';
    R_Fixed_nbar = mean(R_Fixed_nbar,1).';
    R_ARFA_FS = mean(R_ARFA_FS,1).';
    R_ARFA_LS = mean(R_ARFA_LS,1).';
    R_ARFA_SC = mean(R_ARFA_SC,1).';
    R_APS = mean(R_APS,1).';
    R_AS = mean(R_AS,1).';
    Rate = [R_Fixed_N,R_Fixed_nbar,R_ARFA_FS,R_ARFA_LS,R_ARFA_SC,R_APS,R_AS];
    
    %% EE
    EE_Fixed_N = mean(EE_Fixed_N,1).';
    EE_Fixed_nbar = mean(EE_Fixed_nbar,1).';
    EE_ARFA_FS = mean(EE_ARFA_FS,1).';
    EE_ARFA_SC = mean(EE_ARFA_SC,1).';
    EE_ARFA_LS = mean(EE_ARFA_LS,1).';
    EE_APS = mean(EE_APS,1).';
    EE_AS = mean(EE_AS,1).';
    EE = [EE_Fixed_N,EE_Fixed_nbar,EE_ARFA_FS,EE_ARFA_SC,EE_ARFA_LS,EE_APS,EE_AS];
    
    %% Power
    P_Fixed_N = mean(P_Fixed_N,1).';
    P_Fixed_nbar = mean(P_Fixed_nbar,1).';
    P_ARFA_FS = mean(P_ARFA_FS,1).';
    P_ARFA_SC = mean(P_ARFA_SC,1).';
    P_ARFA_LS = mean(P_ARFA_LS,1).';
    P_APS = mean(P_APS,1).';
    P_AS = mean(P_AS,1).';
    Power = [P_Fixed_N,P_Fixed_nbar,P_ARFA_FS,P_ARFA_SC,P_ARFA_LS,P_APS,P_AS];
    
    %% run time
    t_Fixed_N = mean(t_Fixed_N,1).';
    t_Fixed_nbar = mean(t_Fixed_nbar,1).';
    t_ARFA_FS = mean(t_ARFA_FS,1).';
    t_ARFA_LS = mean(t_ARFA_LS,1).';
    t_ARFA_SC = mean(t_ARFA_SC,1).';
    t_APS = mean(t_APS,1).';
    t_AS = mean(t_AS,1).';
    time = [t_Fixed_N,t_Fixed_nbar,t_ARFA_FS,t_ARFA_LS,t_ARFA_SC,t_APS,t_AS]
    
    %% save results and plot figures
    save_file = strcat('result_mode',num2str(sim_mode),'_',num2str(nbar));
    save(save_file,'Rate','EE','Power');
    %load(save_file,'Rate','EE','Power');
    leg_str = {'Fixed-activation HBF, $n_l=N, \forall l$', 'Fixed-activation HBF, $n_l=\bar{n}, \forall l$',...
        'Proposed SC-ARFA', 'Proposed PL-based D-ARFA',...
        'Proposed SV-based D-ARFA', 'APS', 'AS'};
    color = {'--k', ':ko','-b*','-rs','-bd', '--gp', '--ks'};
    
    
    if sim_mode == 1
        plot_fig(x_vec,BF_scheme,Rate,leg_str,color,'rate',xaxis);
        plot_fig(x_vec,BF_scheme,EE,leg_str,color,'EE',xaxis);
    elseif sim_mode > 1
        plot_fig(x_vec,BF_scheme,Rate,leg_str,color,'rate',xaxis);
        plot_fig(x_vec,BF_scheme,EE,leg_str,color,'EE',xaxis);
        plot_fig(x_vec,BF_scheme,Power,leg_str,color,'power',xaxis);
    else
        plot_fig(x_vec,BF_scheme,Rate,leg_str,color,'rate',xaxis);
        plot_fig(x_vec,BF_scheme,Power,leg_str,color,'power',xaxis);
    end
end
