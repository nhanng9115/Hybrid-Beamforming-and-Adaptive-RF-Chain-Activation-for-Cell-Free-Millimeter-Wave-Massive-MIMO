function [H,LSmean,rx_array] = ULA_channel(K,Nr,Nt,n_path,D,f)
x_l = D*rand; y_l = D*rand;
AoD_min = -pi/6; AoD_max = pi/6;

H_K = zeros(Nr, Nt, K);  % One user channel
a_Tx = zeros(Nt, n_path, K); % TX steering vector
a_Rx = zeros(Nr, n_path, K); % RX steering vector
rx_array = zeros(Nr,K);
LSmean = 0;
Ga = db2pow(24.5);
for u = 1:K
    Hu = zeros(Nr, Nt);
    max_gain = 0;
    % large scale fading
    x_u = D*rand;
    y_u = D*rand;
    dis = sqrt((x_u-x_l)^2 + (y_u-y_l)^2);
    LS_u = 0;
    
    % small scale fading
    for pp = 1:n_path
        LS_pp = pathloss(dis, f);  LS_u = LS_u + LS_pp/n_path;
        AoD = AoD_min + (AoD_max-AoD_min).* rand;
        a_Tx(:,pp,u) = array_response(AoD, Nt);
        a_TX_u = a_Tx(:,pp,u);
        
        AoA = AoD_min + (AoD_max-AoD_min).* rand;
        a_Rx(:,pp,u) = array_response(AoA, Nr);
        a_RX_u = a_Rx(:,pp,u);
        
        alpha = sqrt(1/2) * (randn + 1j*randn);
        Hp = alpha * a_RX_u * a_TX_u';
        Hu = Hu + Hp;
        
        gain_tmp = norm(Hp)^2;
        if gain_tmp > max_gain
            max_gain = gain_tmp;
            rx_array(:,u) = a_RX_u;
        end
    end
    LSmean = LSmean + Ga*LS_u/K;
    H_K(:,:,u) = sqrt(Ga*LS_u)*sqrt(1/n_path) * Hu;
    
end

H = []; % (Nr x K*Nt)- matrix for all users
for u = K:-1:1
    Hk = zeros(Nr,Nt);
    Hk(:,:) = H_K(:,:,u);
    H = cat(2, Hk, H);
end

end