% Clear and clc
clear
clc

% Load data file
load 'data'
% Number of companies is 21 in paper
num_company = 21;

% Init required data struct
X1_2000 = zeros(21,1);
X1_2001 = zeros(21,1);
X1_2002 = zeros(21,1);

X2_2000 = zeros(21,1);
X2_2001 = zeros(21,1);
X2_2002 = zeros(21,1);

Y1_2000 = zeros(21,1);
Y1_2001 = zeros(21,1);
Y1_2002 = zeros(21,1);
Y2_2000 = zeros(21,1);
Y2_2001 = zeros(21,1);
Y2_2002 = zeros(21,1);

Z1_2000 = zeros(21,1);
Z1_2001 = zeros(21,1);
Z1_2002 = zeros(21,1);
Z2_2000 = zeros(21,1);
Z2_2001 = zeros(21,1);
Z2_2002 = zeros(21,1);

% Instantiate data struct in Table4 in paper
for i=1:21
    X1_2000(i) = X1_sep(1+3*(i-1));
    X1_2001(i) = X1_sep(2+3*(i-1));
    X1_2002(i) = X1_sep(3*i);
    X2_2000(i) = X2_sep(1+3*(i-1));
    X2_2001(i) = X2_sep(2+3*(i-1));
    X2_2002(i) = X2_sep(3*i);
    Y1_2000(i) = Y1_sep(1+3*(i-1));
    Y1_2001(i) = Y1_sep(2+3*(i-1));
    Y1_2002(i) = Y1_sep(3*i);
    Y2_2000(i) = Y2_sep(1+3*(i-1));
    Y2_2001(i) = Y2_sep(2+3*(i-1));
    Y2_2002(i) = Y2_sep(3*i);
    Z1_2000(i) = Z1_sep(1+3*(i-1));
    Z1_2001(i) = Z1_sep(2+3*(i-1));
    Z1_2002(i) = Z1_sep(3*i);
    Z2_2000(i) = Z2_sep(1+3*(i-1));
    Z2_2001(i) = Z2_sep(2+3*(i-1));
    Z2_2002(i) = Z2_sep(3*i);
end

% Init required index
E_s_ori = zeros(21,1);
E_s = zeros(21,1);
E_s_2000 = zeros(21,1);
E_s_2001 = zeros(21,1);
E_s_2002 = zeros(21,1);
W_2000 = zeros(21,1);
W_2001 = zeros(21,1);
W_2002 = zeros(21,1);
MPI_s_00_to_01 = zeros(21,1);
MPI_s_01_to_02 = zeros(21,1);
MPI_s_00_to_02 = zeros(21,1);

% Main body to implement model(1) in paper page 2.
for k=1:num_company
    % 6 optimization parameters for linear programming
    u1 = optimvar('u1','LowerBound',1e-6, 'Type', 'continuous');
    u2 = optimvar('u2','LowerBound',1e-6, 'Type', 'continuous');
    v1 = optimvar('v1','LowerBound',1e-6, 'Type', 'continuous');
    v2 = optimvar('v2','LowerBound',1e-6, 'Type', 'continuous');
    % template for this ooptimization problem
    prob = optimproblem('ObjectiveSense','maximize');
    % constraints
    prob.Objective = u1 * Y1(k) + u2 * Y2(k);
    prob.Constraints.cons0 = v1 * X1(k) + v2 * X2(k) == 1;
    for j = 1:num_company
        eval(['prob.Constraints.cons',num2str(j),'=','u1 * Y1(j) + u2 * Y2(j) - v1 * X1(j) - v2 * X2(j) <= 0']);
    end
    % Solve
    sol = solve(prob);
%     % display 4 parameters
%     disp ([sol.u1, sol.u2, sol.v1, sol.v2]);
    E_s_ori(k) = (sol.u1 * Y1(k) + sol.u2 * Y2(k))/(sol.v1 * X1(k)+sol.v2 * X2(k));
end

for k=1:num_company
    % 4 optimization parameters for linear programming
    u1 = optimvar('u1','LowerBound',1e-6, 'Type', 'continuous');
    u2 = optimvar('u2','LowerBound',1e-6, 'Type', 'continuous');
    v1 = optimvar('v1','LowerBound',1e-6, 'Type', 'continuous');
    v2 = optimvar('v2','LowerBound',1e-6, 'Type', 'continuous');
    % template for this ooptimization problem
    prob = optimproblem('ObjectiveSense','maximize');
    % constraints
    prob.Objective = u1 * Y1(k) + u2 * Y2(k);
    prob.Constraints.cons0 = v1 * X1(k) + v2 * X2(k) == 1;
    for j = 1:21
        eval(['prob.Constraints.cons_2000_',num2str(j),'=','u1 * Y1_2000(j) + u2 * Y2_2000(j) - v1 * X1_2000(j) - v2 * X2_2000(j) <= 0']);
        eval(['prob.Constraints.cons_2001_',num2str(j),'=','u1 * Y1_2001(j) + u2 * Y2_2001(j) - v1 * X1_2001(j) - v2 * X2_2001(j) <= 0']);
        eval(['prob.Constraints.cons_2001_',num2str(j),'=','u1 * Y1_2002(j) + u2 * Y2_2002(j) - v1 * X1_2002(j) - v2 * X2_2002(j) <= 0']);
    end
    % Solve调用求解器进行求解
    sol = solve(prob);
    % display 4 parameters
%     fprintf('\t\tu1\t\t\t\tu2\t\t\t\tv1\t\t\t\tv2\n');
%     disp ([double(sol.u1), double(sol.u2), double(sol.v1), double(sol.v2)]);
    E_s(k) = (sol.u1 * Y1(k) + sol.u2 * Y2(k))/(sol.v1 * X1(k)+sol.v2 * X2(k));
    
    % compute index in model(5)
    E_s_2000(k) = (sol.u1 * Y1_2000(k) + sol.u2 * Y2_2000(k))/(sol.v1 * X1_2000(k)+sol.v2 * X2_2000(k));
    E_s_2001(k) = (sol.u1 * Y1_2001(k) + sol.u2 * Y2_2001(k))/(sol.v1 * X1_2001(k)+sol.v2 * X2_2001(k));
    E_s_2002(k) = (sol.u1 * Y1_2002(k) + sol.u2 * Y2_2002(k))/(sol.v1 * X1_2002(k)+sol.v2 * X2_2002(k));
    MPI_s_00_to_01(k) = E_s_2001(k)/E_s_2000(k);
    MPI_s_01_to_02(k) = E_s_2002(k)/E_s_2001(k);
    MPI_s_00_to_02(k) = E_s_2002(k)/E_s_2000(k);
end

% compute average in Table2
sum_s_ori = 0;
sum_s = 0;
fprintf('\t\tThe following are results of Multi-period Model in Paper.Table2\n\n');
fprintf('\t\t\t    Aggregate Model\t\t\t     Multi-period Model\n\n');
fprintf('\t       E_s_k     E_s_k\n\n');
for k=1:num_company
    fprintf('\t%d ',k);
    disp([E_s_ori(k), E_s(k)]);
    sum_s_ori = sum_s_ori + E_s_ori(k);
    sum_s = sum_s + E_s(k);
end
fprintf('Average');
disp([sum_s_ori/num_company,sum_s/num_company]);

% compute average in Table3
sum_E_s_00_to_01 = 0;
sum_E_s_01_to_02 = 0;
sum_E_s_00_to_02 = 0;

sum_MPI_s_00_to_01 = 0;
sum_MPI_s_01_to_02 = 0;
sum_MPI_s_00_to_02 = 0;

sum_W_00_to_01 = 0;
sum_W_01_to_02 = 0;
sum_W_00_to_02 = 0;

% Display results in Table3
fprintf('\t\t The following are results of Multi-period Model in Paper.Table3\n\n');
fprintf('  \tYear  Period\t   E_s_k     MPI_s_k\n\n');
for k=1:num_company
    fprintf('  \t2000  [00-01]  ');
    disp([E_s_2000(k), MPI_s_00_to_01(k)]);
    fprintf('  \t2001  [01-02]  ');
    disp([E_s_2001(k), MPI_s_01_to_02(k)]);
    fprintf('  \t2002  [00-02]  ');
    disp([E_s_2002(k), MPI_s_00_to_02(k)]);
    fprintf('\n');
    
    sum_E_s_00_to_01 = sum_E_s_00_to_01 + E_s_2000(k);
    sum_E_s_01_to_02 = sum_E_s_01_to_02 + E_s_2001(k);
    sum_E_s_00_to_02 = sum_E_s_00_to_02 + E_s_2002(k);
    
    sum_MPI_s_00_to_01 = sum_MPI_s_00_to_01 + MPI_s_00_to_01(k);
    sum_MPI_s_01_to_02 = sum_MPI_s_01_to_02 + MPI_s_01_to_02(k);
    sum_MPI_s_00_to_02 = sum_MPI_s_00_to_02 + MPI_s_00_to_02(k);
end

fprintf(' Average  [00-01]  ');
disp([sum_E_s_00_to_01/num_company, sum_MPI_s_00_to_01/num_company]);
fprintf(' Average  [01-02]  ');
disp([sum_E_s_01_to_02/num_company, sum_MPI_s_01_to_02/num_company]);
fprintf(' Average  [00-02]  ');
disp([sum_E_s_00_to_02/num_company, sum_MPI_s_00_to_02/num_company]);
    
