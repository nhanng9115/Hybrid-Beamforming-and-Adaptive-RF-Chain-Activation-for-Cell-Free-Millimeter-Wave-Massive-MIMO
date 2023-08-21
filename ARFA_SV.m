function [R, EE, power] = ARFA_SV(H1,H,Nr,Nt,nbar,L,K,gamma,N,Pt,sigma2)

comp = 1;

N_vec_conv = N*ones(L,1);
[~,~,eig_tab] = HBF_D(H,Nr,Nt,N_vec_conv,N,L,K,gamma);
all_eig = eig_tab(:);

all_eig_sort = sort(all_eig, 'descend');
all_eig_max = all_eig_sort(1:L*nbar);
eig_min = all_eig_max(L*nbar);

for l = 1:L
    eig_l = eig_tab(:,l);
    N_vec(l) = sum(eig_l >= eig_min);
end

% sum(N_vec)
if (sum(N_vec) ~= L*nbar)
    error('***********');
end
[F,W] = HBF_D(H,Nr,Nt,N_vec,nbar,L,K,gamma); X = F*W;
R = log2(det(eye(L*K*Nt) + gamma*pinv(X)*(H1*H1')*(X)));

% compute EE
n_AP = L-sum(N_vec == 0);
power = get_power(L,N,Nr,n_AP,nbar,'ARFA',K,Pt,sigma2,R);
EE = R/power;
end %eof