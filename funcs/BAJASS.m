function scores = BAJASS(par,y,true_L)
%%% BAJASS: Baseline jammer-mitigating synchronization method
% Returns the corresponding scores for y. 
% The resulting l_hat as a function of the threshold has to be determined
% externally for reasons of computational efficiency (evaluating different
% scores while only evaluating JASS once). 
% In contrast to JASS, the score is computed by simply removing the
% strongest dimension of the windowed receive signal, and then computing a
% normalized correlation with the synchronization sequence
  scores = zeros(1,true_L+1);
  for l=0:true_L
    Y_l = y(:,l+1:l+par.seq_length);
    if par.num_power_its==0
      [V_full,D] = eig(Y_l*Y_l');
      [~,ind] = sort(diag(D),'descend');
      V = V_full(:,ind(1:par.I_est));
    else
      V = power_method(Y_l*Y_l', par);
    end
    P = eye(par.B)-V*pinv(V);
    proj_Y = P*Y_l;
    c_l = proj_Y*par.s';
    score = norm(c_l)^2/norm(proj_Y,'fro')^2;
    scores(l+1) = score;
  end
end

function V = power_method(A, par)
% returns (an approximate version of) the par.I_est eigenvectors which belong to the
% largest eigenvalues of a positive semidefinite symmetric matrix A.
  V = zeros(par.B, par.I_est);
  A_ii = A;
  for ii=1:par.I_est
    v = 1/sqrt(2)*(randn(par.B,1)+1j*randn(par.B,1));
    for kk=1:par.num_power_its
      v = A_ii*v;
      v = v/norm(v);
    end
    lambda = v'*A_ii*v;
    V(:,ii) = v;
    A_ii = A_ii - lambda*(v*v');
  end
end