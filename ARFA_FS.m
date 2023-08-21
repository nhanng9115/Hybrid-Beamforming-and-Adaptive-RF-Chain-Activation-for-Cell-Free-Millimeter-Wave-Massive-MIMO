function [R, EE, power] = ARFA_FS(H1,H,Nr,Nt,nbar,L,K,gamma,N,Pt,sigma2)

N_UB = N;
N_vec = nbar*ones(L,1);
N_vec_best = N_vec;

[F,W, R_vec] = HBF_SC(H,Nr,Nt,N_vec,nbar,L,K,gamma);
R = log2(det(eye(L*nbar) + gamma*F'*(H1*H1')*F));
[~, ll_order] = sort(R_vec, 'descend');

comp = 1; 
while N_vec(ll_order(1)) < N_UB
    tt = 1;
    while (tt < round(L/2))
        kk = L;
        while (N_vec(ll_order(kk)) == 0) && (kk > tt)
            kk = kk - 1;
        end
        if (tt < kk) && (N_vec(ll_order(tt)) < N_UB)
            N_vec(ll_order(tt)) = N_vec(ll_order(tt)) + 1;
            N_vec(ll_order(kk)) = N_vec(ll_order(kk)) - 1;
            if sum(N_vec) == L*nbar
                comp = comp + 1;
                [F,W] = HBF_SC(H,Nr,Nt,N_vec,nbar,L,K,gamma);
                Rtmp = log2(det(eye(L*nbar) + gamma*F'*(H1*H1')*F));
                if Rtmp > R
                    N_vec_best = N_vec;
                end
            end
        else
            break
        end
        tt = tt + 1;
    end
end
if (sum(N_vec_best) ~= L*nbar) || (sum(N_vec_best > N_UB) > 0)
    error('***********');
end
[F,W] = HBF_SC(H,Nr,Nt,N_vec_best,nbar,L,K,gamma);
X = F*W;
R = log2(det(eye(L*K*Nt) + gamma*pinv(X)*(H1*H1')*(X)));

% compute EE
n_AP = L-sum(N_vec_best == 0);
power = get_power(L,N,Nr,n_AP,nbar,'ARFA',K,Pt,sigma2,R);
EE = R/power;
end % eof