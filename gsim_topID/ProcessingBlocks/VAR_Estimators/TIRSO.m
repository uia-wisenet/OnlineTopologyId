function [t_A, v_strongcvxPar]= TIRSO (P, lambda, gamma, mu, sigma, v_alpha, m_y, m_A_initial,b_trace,b_diminishing)
% This function computes an estimate via TIRSO
% INPUT: P, lambda,  gamma, mu, sigma, v_alpha, m_y, m_A_initial
% OUTPUT: A (a_n, n=1,...,N)
% Initialization: a_n[P]=0, n=1,...,N; Phi[P-1]= sigma^2*eye(N*P); m_r= zeros(N*P,N);
% y has size N X T 
% A has a size of N X NP
% alpha has a size of T X 1
[N,T] = size(m_y);
assert(length(v_alpha)==T,'v_alpha should be a vector of size T');
assert(size(m_A_initial,1)==N && size(m_A_initial,2)==N*P,'A_initial should have of size N X NP');
Phi_prev= sigma^2*eye(N*P); % initializing Phi
r_prev= zeros(N*P,N); % r has NP X N size to avoid transpose 
a_prev= m_A_initial; 
m_r= zeros(N*P,N); 
t_A= zeros(N,N*P,T); 
v_strongcvxPar = zeros(T,1);
for t= P+1: T
    % receive data y[t]
    % form g[t] via g[t]= vec([y[t-1],...,y[t-P]]^T)
    y_prev= m_y(:,t-P:t-1);
    aux = transpose(fliplr(y_prev));
    g = aux(:);
    % update Phi
    Phi_t= gamma*Phi_prev+ mu*(g*g');
    v_strongcvxPar(t)= min(eig(Phi_t));
    if b_trace==1
        if b_diminishing==1
            v_alpha(t)= 1/(trace(Phi_t)*sqrt(t));
        else 
            v_alpha(t)= 1/(trace(Phi_t));
        end
    end
    for n=1:N
        %update r_n
        m_r(:,n)= gamma*r_prev(:,n) + mu*m_y(n,t)*g;
        grad_n= Phi_t*(a_prev(n,:))' - m_r(:,n); % v_n in the paper
        for nprime=1:N
            groupindices=(nprime-1)*P+1:nprime*P; % n,n' group indices
            af_nnprime= (a_prev(n,groupindices))' - v_alpha(t)*grad_n(groupindices);
            if n ~= nprime
                t_A(n,groupindices,t)= max(0,(1- (v_alpha(t)*lambda)/norm(af_nnprime)))*af_nnprime; %indicator rem
            else
                t_A(n,groupindices,t)=af_nnprime;
            end
        end
    end
    % to store A, Phi, and r_n
    a_prev=t_A(:,:,t); 
    Phi_prev= Phi_t; 
    r_prev= m_r; 
end
end