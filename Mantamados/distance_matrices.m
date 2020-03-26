% We wish to compute pairwise-distances between 100 vases at 1000
% different timesteps.

% We will be computing the L2 distance after pairwise-registration.
% We shall be using functions in the shape analysis library, libsrvf.
% Here, the distances between two vases c1,c2, are computed as follows:
%       1) Evaluate c1,c2 to make sure they are strictly-increasing.
%       2) Compute square root vlocity representations of the
%          evaluated curves, q1,q2.
%       3) Find optimal warping function, gamma for the registration
%          of c1 to c2.
%       4) Compute L2 distance between q1 and the composition (q2,gamma)

%% 1) Load Data
% We load the data from the csv into a cell of size
% (100(no. vases),1000(no.timesteps),2(x,y)).
% This makes it much easier to use.

% 1.1. Load data from csv
pth = "C:\Users\arian\Downloads\PC_RW_large_Mantamados_open.csv\PC_RW_large_Mantamados_open.csv";
data = readtable(pth);  

n_vs = 100;
% no. vases
n_tm = 1000;
% no. timesteps

% 1.2. Create cell
data_cell = cell(n_vs,n_tm,2);

for i=1:n_vs
    % Find all table entries for vase i
    values = data(data{:,1}==i,1:5);
    for j=1:n_tm
        % Find coordinates of vase i and timestep j
        x = values{values{:,2}==j,5};
        y = values{values{:,2}==j,4};
        data_cell(i,j,1) = {x};
        data_cell(i,j,2) = {y};
    end
end


%% 2) Scale Data
% libsrvf registration algorithm works with strictly
% increasing curves. So we swap y and x around
% and ensure curves are *stricly* increasing.
% The algorithm also requires that the start and
% end points of two curves are close to each other
% therefore we scale the curves so they are all
% between 0 and 1.

% 2.1 Create cell
scaled_data = cell(n_vs,n_tm,2);

for i=1:n_vs
    for j=1:n_tm
        % Load y coordinate as T and x coordinate
        % as F for vase i and timestep j
        T = data_cell{i,j,2};
        F = data_cell{i,j,1};
        % Scale F and T so that values are
        % between 0 and 1.
        Tn = T./(max(T) - min(T));
        Tn = Tn - min(Tn);
        Fn = F./(max(F) - min(F));
        Fn = Fn - min(Fn);
        % Removing any horizontal and vertical
        % jumps. Horizontal jumps violate
        % the *strictly* increasing condition whilst
        % vertical jumps will cause issues further
        % along the line.
        [Fn_,Tn_] = remove_jumps(Fn,Tn);

        scaled_data(i,j,1) = {Fn_};
        scaled_data(i,j,2) = {Tn_};
    end
end

%% 3) Find Percentiles of Distance Matrices
% Using the function get_distance_matrix, which
% finds the distance matrix between all vases
% at a given timestep, we compute the 5th,50th
% 95th percentile of the distances at each timestep

% 3.1. Create array that will contain all percentiles.
vase_prct = zeros(n_tm,4);
% 3.2. Assign the timesteps to the first column
vase_prct(:,1) = 1:n_tm;
% 3.3. Add cell to catch any errors from the registration.
tot_errors = cell(1,n_tm);

for i=1:n_tm
    % 3.4. Get distance matrix for timestep i.
    [dm,errors] = get_distance_matrix(scaled_data,i,'dp');
    if ~isempty(errors)
        % 3.5. If there are any errors, add them to the
        % total error cell, tot_errors.
        tot_errors(1,i) = {errors};
    end
    % 3.6. Find percentiles of distances in distance matrix
    prct = prctile(dm,[5 50 95],'all');
    vase_prct(i,2:4) = prct;
    % Progress notifications:
    if (i==100) || (i==250) || (i==500) || (i==750)
        disp(string(i/10)+"% done")
    end
end

% 3.7. Save to CSV.
writetable(table(vase_prct),'rw_percentiles.csv')

        
        
        
        