function t_A= OSGDTISO (P, lambda, v_alpha, m_y, m_A_initial)
% This function computes an estimate via Online Subgradient Descent method
% INPUT: P, lambda, alpha,
% OUTPUT: A (a_n[t], n=1,...,N)
% Initialization: a_n[P]=0, n=1,...,N
% y is of size N X T
% A has a size of N X NP
% alpha has a size of T X 1

a_prev= m_A_initial; % size N X NP
[N,T] = size(m_y);
assert(numel(v_alpha)==T,'v_alpha should be a vector');
assert(size(m_A_initial,1)==N && size(m_A_initial,2)==N*P,'A_initial should have of size N X NP');
t_A= zeros(N,N*P,T);
regsubgrad_n= zeros(N*P,1);
for t=P+1:T % in paper, t=P,...
    % receive data y[t]
    % form g[t] via g[t]= vec([y[t-1],...,y[t-P]]^T)
    y_prev= m_y(:,t-P:t-1);
    m_aux = transpose(fliplr(y_prev));
    g = m_aux(:);
    for n=1:N
        lossgrad_n= (g'*(a_prev(n,:))' - m_y(n,t))*g; % this is v_n in the paper
        for nprime=1:N
            groupindices=(nprime-1)*P+1:nprime*P;
            if (nprime==n) 
                regsubgrad_n(groupindices)=zeros(P,1);
            elseif (norm(a_prev(n,groupindices))==0)
                randVector= randn(P,1);
                regsubgrad_n(groupindices)=randVector./(norm(randVector));
            else
                regsubgrad_n(groupindices)= a_prev(n,groupindices)/norm(a_prev(n,groupindices));
            end
        end
        subgradient_n= lossgrad_n + lambda*regsubgrad_n;
        t_A(n,:,t)=  a_prev(n,:)- v_alpha(t)* subgradient_n';
    end
    a_prev= t_A(:,:,t); % to store A for the next time instant
end

end