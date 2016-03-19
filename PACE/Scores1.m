function [xi_est, xi_var,y_predOrig]=Scores1(y, t, mu, phi, lambda, sigma, sigma_new, noeig, error, out1, regular)

  [muSub, phiSub] = convertMuPhi(t, out1, mu, phi, regular);
  ncohort = length(y);
  LAMBDA = diag(lambda);
  xi_var = cell(1,ncohort);
  %fprintf(1,'Start calculating the PC scores \n');
  %update \xi 

  if error==1 
      sigma1 = sigma_new;
      if regular == 2
         yy= reshape(cell2mat(y), length(y{1}), ncohort)';
         error0 = diag(sigma1);
         A = LAMBDA*phiSub'*pinv(phiSub*LAMBDA*phiSub'+error0);
         MU = repmat(muSub, ncohort,1);
         B = yy-MU;
         xi_est = (A*B')';
         y_predOrig = MU+xi_est*phiSub';
         y_predOrig = num2cell(y_predOrig,2);
         C = LAMBDA-A*(LAMBDA*phiSub')';
         for i = 1:ncohort
            xi_var{i} = C;
         end     
      else
        sigmaSub = cell(1,length(t));%
        for i = 1:length(t)%
          sigmaSub{i} = mapX1d(out1,sigma,t{i});%
        end%
	     y_predOrig = cell(1,ncohort);
	     xi_est = zeros(ncohort, noeig);
	     zeta_est = xi_est;
	     phii= phiSub;
	     mu_i = muSub;
	     for i = 1:ncohort
            if regular ~= 2
	           phii = phiSub{i};
	           mu_i = muSub{i};
            end
	        yi= y{i};
            error0=diag(sigmaSub{i});
            A = LAMBDA*phii'*pinv(phii*LAMBDA*phii'+error0);
            xi_est(i,:)=(A*(yi-mu_i)')';
            xi_var{i}=LAMBDA-A*(LAMBDA*phii')';
            y_predOrig{i} = mu_i+xi_est(i,:)*phii';
         end
      end
  elseif error==0 
      if regular == 2
          yy= reshape(cell2mat(y), length(y{1}), ncohort)';
          A = LAMBDA*phiSub'*pinv(phiSub*LAMBDA*phiSub');
          MU = repmat(muSub, ncohort,1);
          B = yy-MU;
          xi_est = (A*B')';
          y_predOrig = MU+xi_est*phiSub';
          y_predOrig = num2cell(y_predOrig,2);
          C = LAMBDA-A*(LAMBDA*phiSub')';
          for i = 1:ncohort
    	     xi_var{i} = C;
          end
      else
          y_predOrig = cell(1,ncohort);
       	  xi_est = zeros(ncohort, noeig);
	      phii= phiSub;
	      mu_i = muSub;
	      for i = 1:ncohort
             if regular ~= 2
		        phii = phiSub{i};
		        mu_i = muSub{i};
             end
	         yi= y{i};
			 A = LAMBDA*phii'*pinv(phii*LAMBDA*phii');
             xi_est(i,:)=(A*(yi-mu_i)')';
             xi_var{i}=LAMBDA-A*(LAMBDA*phii')';
             y_predOrig{i} = mu_i+xi_est(i,:)*phii';
          end
      end
  end
end
