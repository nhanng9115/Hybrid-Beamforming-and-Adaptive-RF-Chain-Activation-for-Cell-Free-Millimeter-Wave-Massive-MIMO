function [F,W, R_vec, sigma_max] = HBF_SC(H,Nr,Nt,N_vec,nbar,L,K,gamma)
q = 2^4;
sigma_max = 0;
F = zeros(L*Nr,L*nbar);
I = eye(K*Nt);
Q = I;
R_vec = zeros(L,1);
col_end = 0;
for l = 1:L
    N_l = N_vec(l);
    col_start = col_end + 1;
    col_end = col_start + N_l - 1;
    
    %H_l = H((l-1)*Nr+1:(l-1)*Nr+Nr, :);
    H_l = H(:,:,l);
    T_l = H_l*(Q^(-1))*H_l';
    [U,D,V] = svd(T_l);
    F_l = U(:,1:N_l);
    %F_l*F_l'
    F_l = quant_sub(Nr,L,N_l,F_l); % quantize U
    R_vec(l) = log2(det(eye(N_l) + gamma*(F_l' * T_l * F_l)));

    G_l = H_l'*F_l*F_l'*H_l;
    E_l = I + gamma*Q^(-1)*G_l;
    Q = Q * E_l;
    F((l-1)*Nr+1:(l-1)*Nr+Nr, col_start:col_end) = F_l;
    
    % digital precoding
    %W_l = F_l'*H_l;
    J = F_l'*H_l*H_l'*F_l + 1/gamma*F_l'*F_l;
    W_l = J^-1*F_l'*H_l;
    W(col_start:col_end,(l-1)*K*Nt+1:(l-1)*K*Nt+K*Nt) = W_l;
end
% X = F*W;
if (size(F,2) ~= L*nbar)
    error('***********');
end
end