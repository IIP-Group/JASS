function scores = correlation_sync(par,y,true_L)
%%% Baseline non-jammer-mitigating synchronization method
% Returns the corresponding scores for y. 
% The resulting l_hat as a function of the threshold has to be determined
% externally for reasons of computational efficiency (evaluating different
% scores while only evaluating JASS once). 
% In contrast to JASS, the score is computed computing a normalized correlation 
% with the synchronization sequence, without trying to mitigate the jammer
  scores = zeros(1,true_L+1);
  for l=0:true_L
    Y_l = y(:,l+1:l+par.seq_length);
    c_l = Y_l*par.s';
    score = norm(c_l)^2/norm(Y_l,'fro')^2;
    scores(l+1) = score;
  end
end

