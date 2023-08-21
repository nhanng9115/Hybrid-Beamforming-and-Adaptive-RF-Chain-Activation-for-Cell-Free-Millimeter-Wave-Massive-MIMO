function R = AS_UESA(H, q, Nr, N, K, rho, collect_data, count_max, doES)
% collect_data = 0;
% count_max = 1;
I = eye(K);
if doES == 1
    if Nr == 32
        if collect_data == 1
            comb_all = combnk(1:Nr,N);
            save('.\data\comb_AS_32','comb_all')
        else
            dat = load('.\data\comb_AS_32');
            comb_all = dat.comb_all;
        end
    else
        comb_all = combnk(1:Nr,N);
    end
    n_comb = size(comb_all,1);
    R = 0;
    count = 0;
    best_comb = comb_all(1,:);
    for i = 1:n_comb
        comb = comb_all(i,:);
        Htmp = H(comb,:);
        Rtmp = log2(det(I + rho * (Htmp'*Htmp)));
        if Rtmp > R
            best_comb = comb;
            R = Rtmp;
            count = 0;
        else
            count = count + 1;
        end
        if count == count_max
            break;
        end
    end
else
    comb = 1:Nr;
    for s = 1:Nr-N
        Rvec = [];
        ant_vec = [];
        %comb_tmp = comb(comb~=1);
        
        for i = 1:length(comb)
            ant_de = comb(i);
            ant_vec = cat(2,ant_vec,ant_de);
            comb_tmp = comb(comb~=ant_de);
            Htmp = H(comb_tmp,:);
            Rtmp = log2(det(I + rho * (Htmp'*Htmp)));
            Rvec = cat(2,Rvec,Rtmp);
        end
        [~, i_min] = max(Rvec);
        ant_de_min = ant_vec(i_min);
        comb = comb(comb~=ant_de_min);
        1;
    end
    %comb
    Htmp = H(comb,:);
    R = log2(det(I + rho * (Htmp'*Htmp)));
end
end % eof