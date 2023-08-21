clear all;
clc;
D = 1000; f = 28e9; Ga = db2pow(15); n_path = 3;
% Nr = 64; Nt = 4; K = 2; L = 10; N = 8; nbar = 2;
Nr = 64; Nt = 4; L = 32; N = 8; nbar = 2; K = 8;
n_chan = 10;

for ii = 1:n_chan
    H_tmp = [];
    
    for ll = 1:L
        x_l = D*rand; y_l = D*rand;
        AoD_min = -pi/6; AoD_max = pi/6;
        
        H_K = zeros(Nr, Nt, K);  % One user channel
        a_Tx = zeros(Nt, n_path, K); % TX steering vector
        a_Rx = zeros(Nr, n_path, K); % RX steering vector
        rx_array = zeros(Nr,K);
        LS_ll = zeros(K,1);
        
        
        for kk = 1:K
            Hu = zeros(Nr, Nt);
            max_gain = 0;
            % large scale fading
            x_kk = D*rand;
            y_KK = D*rand;
            dis = sqrt((x_kk-x_l)^2 + (y_KK-y_l)^2);
            LS_KK = zeros(n_path,1);
            
            % small scale fading
            for pp = 1:n_path
                LS_KK(pp,1) = pathloss(dis, f);  
                
                AoD = AoD_min + (AoD_max-AoD_min).* rand;
                a_Tx(:,pp,kk) = array_response(AoD, Nt);
                a_TX_u = a_Tx(:,pp,kk);
                
                AoA = AoD_min + (AoD_max-AoD_min).* rand;
                a_Rx(:,pp,kk) = array_response(AoA, Nr);
                a_RX_u = a_Rx(:,pp,kk);
                
                alpha = sqrt(1/2) * (randn + 1j*randn);
                Hpp = alpha * a_RX_u * a_TX_u';
                Hu = Hu + sqrt(Nr*Nt/n_path)*Hpp;
                
                gain_tmp = norm(Hpp)^2;
                if gain_tmp > max_gain
                    max_gain = gain_tmp;
                    rx_array(:,kk) = a_RX_u;
                end
            end
            
            H_all(:,:,kk,ll,ii) = sqrt(Ga*mean(LS_KK)) *  Hu;
            LS_ll(kk,1) = Ga*mean(LS_KK);
            Ph_ll(kk,1) = norm(H_all(:,:,kk,ll,ii),'fro')^2;
        end
        
        %[H_l,LS_l,rx_array_l] = ULA_channel(K, Nr, Nt,n_path,D,f);
        %H_all(:,:,ll,ii) = H_l;
        LS_all(ll,ii) = mean(LS_ll);
        Ph_all(ll,ii) = sum(Ph_ll);
        rx_array_all(:,:,ll,ii) = rx_array;
        %H_tmp = cat(1,H_tmp,H_l);
    end
    %H_reli(:,:,ii) = H_tmp;
end
file_name = strcat(num2str(Nr),'x',num2str(Nt),'x',num2str(L),'x',num2str(K));
save(file_name,'H_all','LS_all','Ph_all','rx_array_all');