function [R, EE, power] = APS(H1,H,Nr,Nt,nbar,L,K,gamma,pRx,N,Pt,sigma2)
comp = 0;
n_AP = floor(L*nbar/N);

[~, l_order] = sort(pRx,'descend');

F = zeros(n_AP*Nr,n_AP*nbar);
I = eye(K*Nt);
Q = I;
R_vec = zeros(L,1);
col_end = 0;

for ii = 1:n_AP % choose N_AP best APs
    l = l_order(ii);
    col_start = col_end + 1;
    col_end = col_start + N - 1;
    
    %H_l = H((l-1)*Nr+1:(l-1)*Nr+Nr, :);
    H_l = H(:,:,l);
    %T_l = H_l*(Q^(-1))*H_l';
    [U,D,V] = svd(H_l);
    F_l = U(:,1:N);
    F_l = quant_sub(Nr,n_AP,N,F_l); % quantize U
    %R_vec(l) = log2(det(eye(N) + gamma*(F_l' * T_l * F_l)));

    %G_l = H_l'*F_l*F_l'*H_l;
    %E_l = I + gamma*Q^(-1)*G_l;
    %Q = Q * E_l;
    F((l-1)*Nr+1:(l-1)*Nr+Nr, col_start:col_end) = F_l;
    
    % digital precoding
    W_l = F_l'*H_l;
    W(col_start:col_end,(l-1)*K*Nt+1:(l-1)*K*Nt+K*Nt) = W_l;
    
    H_APS((l-1)*Nr+1:(l-1)*Nr+Nr, :) = H1((l-1)*Nr+1:(l-1)*Nr+Nr, :);
end
% R = log2(det(eye(size(H_APS*H_APS',2)) + gamma* (H_APS*H_APS')));
X = F*W;
R = log2(det(eye(size(X,2)) + gamma*pinv(X)*(H_APS*H_APS')*(X)));

% R = sum(R_vec);
% R = log2(det(eye(size(F'*(H_APS*H_APS')*F,2)) + gamma* F'*(H_APS*H_APS')*F));
% compute EE
power = get_power(L,N,Nr,n_AP,nbar,'APS',K,Pt,sigma2,R);
EE = R/power;
end % eof