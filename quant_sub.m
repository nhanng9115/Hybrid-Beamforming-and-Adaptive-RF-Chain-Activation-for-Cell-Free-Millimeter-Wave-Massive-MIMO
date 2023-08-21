function F_quan = quant_sub(Nr,L,N,U)
q = 2^4;

Dict = 1/sqrt(Nr) * exp(1i * 2*pi/q * [0:q-1]).';

% Conver U to andlge and argument
Phi_u = angle(U);
Phi_dict = angle(Dict);

F_quan = zeros(Nr,N);
for nn = 1:N
    for rr = 1:Nr
        delta_phi = wrapToPi(Phi_dict - Phi_u(rr,nn));
        m = abs(delta_phi);
        [~, iMin] = min(m);
        F_quan(rr,nn) = Dict(iMin);
    end
end
end