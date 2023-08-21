function [F,W, eig_tab] = HBF_D(H,Nr,Nt,N_vec,nbar,L,K,gamma)
F = zeros(L*Nr,L*nbar);
I = eye(K);
Q = I;
eig_tab = zeros(Nr,L);
col_end = 0;

for l = 1:L
    N_l = N_vec(l);
    col_start = col_end + 1;
    col_end = col_start + N_l - 1;
    
    %H_l = H((l-1)*Nr+1:(l-1)*Nr+Nr, :);
    H_l = H(:,:,l);
    [U,D,V] = svd(H_l);
    
    eigs = diag(D);
    [eig_sort, col_sort] = sort(eigs,'descend');
    eig_tab(1:N_l,l) = eig_sort(1:N_l);
    U = U(:,col_sort(1:N_l));
    F_l = U(:,1:N_l);
    F_l = quant_sub(Nr,L,N_l,F_l);

    F((l-1)*Nr+1:(l-1)*Nr+Nr, col_start:col_end) = F_l;
    % digital precoding
    %W_l = F_l'*H_l;
    J = F_l'*H_l*H_l'*F_l + 1/gamma*F_l'*F_l;
    W_l = J^-1*F_l'*H_l;
    W(col_start:col_end,(l-1)*K*Nt+1:(l-1)*K*Nt+K*Nt) = W_l;
end

if (size(F,2) ~= L*nbar)
    error('***********');
end
end