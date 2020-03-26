function [distmat,errors] = get_distance_matrix(data,timestep,alg)
    
% This function takes a dataset (cell) of size (n_vs,n_tm,2)
% where n_vs is the number of curves and n_tm is the number
% of timesteps.

    n_vs = size(data);
    n_tm = n_vs(2);
    % Number of timesteps
    n_vs = n_vs(1);
    % Number of vases

    timestep = min(timestep,n_tm);
    % In the case that a number greater than the total number
    % of timesteps is passed, function returns distance matrix
    % at last timestep.

    distmat = zeros(n_vs,n_vs);
    % this is where we store all of the distances.

    errors = zeros(2,n_vs);
    % this is to catch any errors that may occur during
    % the registration.

    k1 = 1;
    k2 = 1;
    tot = n_vs*(n_vs-1)/2;

    lq = round(tot/4);
    hlf = 2*lq;
    uq = lq+hlf;
    % calculate quartiles to display progress

    for i = 1:n_vs

        T1 = data{i,timestep,2};
        F1 = data{i,timestep,1};
 
        for j = i+1:n_vs
        % we can optimise by not calculating distances
        % between all possible pairs/orders, since
        % d(i,j) = d(j,i) and d(i,i) = 0. If we have 
        % already found the distance between vase i and 
        % vase j, we won't formally compute the distance
        % between vase j and vase i.
            T2 = data{j,timestep,2};
            F2 = data{j,timestep,1};

            try
                d = plf_reparam_distance(F1,T1,F2,T2,alg);
                % computes L2 distances after pw registration
            catch
                errors(:,k2) = [i j];
                d = 1000;
                k2 = k2+1;
                % we note down the vase pairs that caused
                % errors during registration.
            end

            distmat(i,j) = d;
            distmat(j,i) = d;
            % since for our distance metric d(x,y) = d(y,x).

%             if (k1 == lq) || (k1 == hlf) || (k1 == uq)
%                 disp(string(round((k1/tot)*100))+"% done")
%             end

            k1 = k1+1;

        end
    end
    
    errors = errors(:,1:k2-1);
end