function [R, EE, power] = AS(H1,H,Nr,Nt,nbar,L,K,gamma,N_AS,Pt,sigma2)
comp = 0;
H_selected = [];

for l = 1:L
    %H_l = H(Nr*(l-1)+1:Nr*l, :);
    H_l = H(:,:,l);
    I = 1:Nr;
    W = [];
    G = H_l*H_l';
    for n = 1:N_AS
        if n == 1
            for i = 1:length(I)
                g_n(i) = sqrt(G(i,i));
            end
            [i_max, w_n] = max(g_n);
            I(w_n) = [];
            W = [W, w_n];
        else
            for i = 1:length(I)
                phi_i = 0;
                for m = 1:length(W)
                    phi_i = phi_i + sqrt(1 - G(i,W(m))*G(W(m),i)/(G(i,i)*G(W(m),W(m))));
                end
                g_n(i) = sqrt(G(i,i)/Nr) + phi_i;
            end
            [i_max, w_n] = max(g_n);
            I(w_n) = [];
            W = [W, w_n];
        end
    end
    H1_l = H1(Nr*(l-1)+1:Nr*l, :);
    H_selected = cat(1,H_selected,H1_l(W,:));
end

R = log2(det(eye(size(H_selected,2)) + gamma*(H_selected'*H_selected)));
% compute EE
power = get_power(L,N_AS,Nr,L,nbar,'AS',K,Pt,sigma2,R);
EE = R/power;
end % eof