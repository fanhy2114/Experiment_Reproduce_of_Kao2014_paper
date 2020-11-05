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
E_1_ori = zeros(21,1);
E_2_ori = zeros(21,1);
E_s = zeros(21,1);
E_1 = zeros(21,1);
E_2 = zeros(21,1);
E_s_2000 = zeros(21,1);
E_s_2001 = zeros(21,1);
E_s_2002 = zeros(21,1);
E_1_2000 = zeros(21,1);
E_1_2001 = zeros(21,1);
E_1_2002 = zeros(21,1);
E_2_2000 = zeros(21,1);
E_2_2001 = zeros(21,1);
E_2_2002 = zeros(21,1);
W_2000 = zeros(21,1);
W_2001 = zeros(21,1);
W_2002 = zeros(21,1);
W_prime_2000 = zeros(21,1);
W_prime_2001 = zeros(21,1);
W_prime_2002 = zeros(21,1);
MPI_s_00_to_01 = zeros(21,1);
MPI_s_01_to_02 = zeros(21,1);
MPI_s_00_to_02 = zeros(21,1);
MPI_1_00_to_01 = zeros(21,1);
MPI_1_01_to_02 = zeros(21,1);
MPI_1_00_to_02 = zeros(21,1);
MPI_2_00_to_01 = zeros(21,1);
MPI_2_01_to_02 = zeros(21,1);
MPI_2_00_to_02 = zeros(21,1);

% Main body to implement model(2) in paper page 2.
for k=1:num_company
    % 6 optimization parameters for linear programming
    u1 = optimvar('u1','LowerBound',1e-6, 'Type', 'continuous');
    u2 = optimvar('u2','LowerBound',1e-6, 'Type', 'continuous');
    v1 = optimvar('v1','LowerBound',1e-6, 'Type', 'continuous');
    v2 = optimvar('v2','LowerBound',1e-6, 'Type', 'continuous');
    w1 = optimvar('w1','LowerBound',1e-6, 'Type', 'continuous');
    w2 = optimvar('w2','LowerBound',1e-6, 'Type', 'continuous');
    % template for this ooptimization problem
    prob = optimproblem('ObjectiveSense','maximize');
    % constraints
    prob.Objective = u1 * Y1(k) + u2 * Y2(k);
    prob.Constraints.cons0 = v1 * X1(k) + v2 * X2(k) == 1;
    for j = 1:num_company
        eval(['prob.Constraints.cons',num2str(j),'=','w1 * Z1(j) + w2 * Z2(j) - v1 * X1(j) - v2 * X2(j) <= 0']);
        eval(['prob.Constraints.consx',num2str(j),'=','u1 * Y1(j) + u2 * Y2(j) - w1 * Z1(j) - w2 * Z2(j) <= 0']);
    end
    % Solve
    sol = solve(prob);
%     % display 6 parameters
%     disp ([sol.u1, sol.u2, sol.v1, sol.v2, sol.w1, sol.w2]);
    % The left part 3 column of Table2, Aggregate Model results
    E_s_ori(k) = (sol.u1 * Y1(k) + sol.u2 * Y2(k))/(sol.v1 * X1(k)+sol.v2 * X2(k));
    E_1_ori(k) = (sol.w1 * Z1(k) + sol.w2 * Z2(k))/(sol.v1 * X1(k)+sol.v2 * X2(k));
    E_2_ori(k) = (sol.u1 * Y1(k) + sol.u2 * Y2(k))/(sol.w1 * Z1(k)+sol.w2 * Z2(k));
end

% Main boody of model(4), model(5), model(7) and model(8)
for k=1:num_company
    % 6 optimization parameters for linear programming
    u1 = optimvar('u1','LowerBound',1e-6, 'Type', 'continuous');
    u2 = optimvar('u2','LowerBound',1e-6, 'Type', 'continuous');
    v1 = optimvar('v1','LowerBound',1e-6, 'Type', 'continuous');
    v2 = optimvar('v2','LowerBound',1e-6, 'Type', 'continuous');
    w1 = optimvar('w1','LowerBound',1e-6, 'Type', 'continuous');
    w2 = optimvar('w2','LowerBound',1e-6, 'Type', 'continuous');
    % template for this ooptimization problem
    prob = optimproblem('ObjectiveSense','maximize');
    % constraints
    prob.Objective = u1 * Y1(k) + u2 * Y2(k);
    prob.Constraints.cons0 = v1 * X1(k) + v2 * X2(k) == 1;
    for j = 1:21
        eval(['prob.Constraints.cons_2000_',num2str(j),'=','w1 * Z1_2000(j) + w2 * Z2_2000(j) - v1 * X1_2000(j) - v2 * X2_2000(j) <= 0']);
        eval(['prob.Constraints.cons_2001_',num2str(j),'=','w1 * Z1_2001(j) + w2 * Z2_2001(j) - v1 * X1_2001(j) - v2 * X2_2001(j) <= 0']);
        eval(['prob.Constraints.cons_2001_',num2str(j),'=','w1 * Z1_2002(j) + w2 * Z2_2002(j) - v1 * X1_2002(j) - v2 * X2_2002(j) <= 0']);
        eval(['prob.Constraints.consx_2000_',num2str(j),'=','u1 * Y1_2000(j) + u2 * Y2_2000(j) - w1 * Z1_2000(j) - w2 * Z2_2000(j) <= 0']);
        eval(['prob.Constraints.consx_2001_',num2str(j),'=','u1 * Y1_2001(j) + u2 * Y2_2001(j) - w1 * Z1_2001(j) - w2 * Z2_2001(j) <= 0']);
        eval(['prob.Constraints.consx_2001_',num2str(j),'=','u1 * Y1_2002(j) + u2 * Y2_2002(j) - w1 * Z1_2002(j) - w2 * Z2_2002(j) <= 0']);
    end
    % Solve调用求解器进行求解
    sol = solve(prob);
    % display 6 parameters
%     fprintf('\t\tu1\t\t\t\tu2\t\t\t\tv1\t\t\t\tv2\t\t\t\tw1\t\t\t\tw2\n');
%     disp ([double(sol.u1), double(sol.u2), double(sol.v1), double(sol.v2), double(sol.w1), double(sol.w2)]);
    % The left part 3 column of Table2, Multi-period Model results
    E_s(k) = (sol.u1 * Y1(k) + sol.u2 * Y2(k))/(sol.v1 * X1(k)+sol.v2 * X2(k));
    E_1(k) = (sol.w1 * Z1(k) + sol.w2 * Z2(k))/(sol.v1 * X1(k)+sol.v2 * X2(k));
    E_2(k) = (sol.u1 * Y1(k) + sol.u2 * Y2(k))/(sol.w1 * Z1(k)+sol.w2 * Z2(k));
    
    % compute index in model(5)
    E_s_2000(k) = (sol.u1 * Y1_2000(k) + sol.u2 * Y2_2000(k))/(sol.v1 * X1_2000(k)+sol.v2 * X2_2000(k));
    E_s_2001(k) = (sol.u1 * Y1_2001(k) + sol.u2 * Y2_2001(k))/(sol.v1 * X1_2001(k)+sol.v2 * X2_2001(k));
    E_s_2002(k) = (sol.u1 * Y1_2002(k) + sol.u2 * Y2_2002(k))/(sol.v1 * X1_2002(k)+sol.v2 * X2_2002(k));
    E_1_2000(k) = (sol.w1 * Z1_2000(k) + sol.w2 * Z2_2000(k))/(sol.v1 * X1_2000(k)+sol.v2 * X2_2000(k));
    E_1_2001(k) = (sol.w1 * Z1_2001(k) + sol.w2 * Z2_2001(k))/(sol.v1 * X1_2001(k)+sol.v2 * X2_2001(k));
    E_1_2002(k) = (sol.w1 * Z1_2002(k) + sol.w2 * Z2_2002(k))/(sol.v1 * X1_2002(k)+sol.v2 * X2_2002(k));
    E_2_2000(k) = (sol.u1 * Y1_2000(k) + sol.u2 * Y2_2000(k))/(sol.w1 * Z1_2000(k)+sol.w2 * Z2_2000(k));
    E_2_2001(k) = (sol.u1 * Y1_2001(k) + sol.u2 * Y2_2001(k))/(sol.w1 * Z1_2001(k)+sol.w2 * Z2_2001(k));
    E_2_2002(k) = (sol.u1 * Y1_2002(k) + sol.u2 * Y2_2002(k))/(sol.w1 * Z1_2002(k)+sol.w2 * Z2_2002(k));
    W_2000(k) = (sol.v1 * X1_2000(k)+sol.v2 * X2_2000(k))/(sol.v1 * X1(k) + sol.v2 * X2(k));
    W_2001(k) = (sol.v1 * X1_2001(k)+sol.v2 * X2_2001(k))/(sol.v1 * X1(k) + sol.v2 * X2(k));
    W_2002(k) = (sol.v1 * X1_2002(k)+sol.v2 * X2_2002(k))/(sol.v1 * X1(k) + sol.v2 * X2(k));
    W_prime_2000(k) = (sol.w1 * Z1_2000(k)+sol.w2 * Z2_2000(k))/(sol.w1 * Z1(k) + sol.w2 * Z2(k));
    W_prime_2001(k) = (sol.w1 * Z1_2001(k)+sol.w2 * Z2_2001(k))/(sol.w1 * Z1(k) + sol.w2 * Z2(k));
    W_prime_2002(k) = (sol.w1 * Z1_2002(k)+sol.w2 * Z2_2002(k))/(sol.w1 * Z1(k) + sol.w2 * Z2(k));
    MPI_s_00_to_01(k) = E_s_2001(k)/E_s_2000(k);
    MPI_s_01_to_02(k) = E_s_2002(k)/E_s_2001(k);
    MPI_s_00_to_02(k) = E_s_2002(k)/E_s_2000(k);
    MPI_1_00_to_01(k) = E_1_2001(k)/E_1_2000(k);
    MPI_1_01_to_02(k) = E_1_2002(k)/E_1_2001(k);
    MPI_1_00_to_02(k) = E_1_2002(k)/E_1_2000(k);
    MPI_2_00_to_01(k) = E_2_2001(k)/E_2_2000(k);
    MPI_2_01_to_02(k) = E_2_2002(k)/E_2_2001(k);
    MPI_2_00_to_02(k) = E_2_2002(k)/E_2_2000(k);
end

% compute average in Table2
sum_s_ori = 0;
sum_1_ori = 0;
sum_2_ori = 0;
sum_s = 0;
sum_1 = 0;
sum_2 = 0;
fprintf('\t\tThe following are results of Multi-period Model in Paper.Table2\n\n');
fprintf('\t\t\t    Aggregate Model\t\t\t     Multi-period Model\n\n');
fprintf('\t       E_s_k     E_1_k     E_2_k     E_s_k    E_1_k     E_2_k\n\n');
for k=1:num_company
    fprintf('\t%d ',k);
    disp([E_s_ori(k), E_1_ori(k), E_2_ori(k), E_s(k), E_1(k), E_2(k)]);
    sum_s_ori = sum_s_ori + E_s_ori(k);
    sum_1_ori = sum_1_ori + E_1_ori(k);
    sum_2_ori = sum_2_ori + E_2_ori(k);
    sum_s = sum_s + E_s(k);
    sum_1 = sum_1 + E_1(k);
    sum_2 = sum_2 + E_2(k);
end
fprintf('Average');
disp([sum_s_ori/num_company,sum_1_ori/num_company,sum_2_ori/num_company, sum_s/num_company, sum_1/num_company, sum_2/num_company]);

% compute average in Table3
sum_E_s_00_to_01 = 0;
sum_E_s_01_to_02 = 0;
sum_E_s_00_to_02 = 0;
sum_E_1_00_to_01 = 0;
sum_E_1_01_to_02 = 0;
sum_E_1_00_to_02 = 0;
sum_E_2_00_to_01 = 0;
sum_E_2_01_to_02 = 0;
sum_E_2_00_to_02 = 0;

sum_MPI_s_00_to_01 = 0;
sum_MPI_s_01_to_02 = 0;
sum_MPI_s_00_to_02 = 0;
sum_MPI_1_00_to_01 = 0;
sum_MPI_1_01_to_02 = 0;
sum_MPI_1_00_to_02 = 0;
sum_MPI_2_00_to_01 = 0;
sum_MPI_2_01_to_02 = 0;
sum_MPI_2_00_to_02 = 0;

sum_W_00_to_01 = 0;
sum_W_01_to_02 = 0;
sum_W_00_to_02 = 0;
sum_W_prime_00_to_01 = 0;
sum_W_prime_01_to_02 = 0;
sum_W_prime_00_to_02 = 0;

% Display results in Table3
fprintf('\t\t The following are results of Multi-period Model in Paper.Table3\n\n');
fprintf('  \tYear  Period\t   E_s_k     MPI_s_k    E_1_k    W_1_k    MPI_1_k    E_2_k     W_2_k     MPI_2_k\n\n');
for k=1:num_company
    fprintf('  \t2000  [00-01]  ');
    disp([E_s_2000(k), MPI_s_00_to_01(k), E_1_2000(k), W_2000(k), MPI_1_00_to_01(k), E_2_2000(k), W_prime_2000(k), MPI_2_00_to_01(k)]);
    fprintf('  \t2001  [01-02]  ');
    disp([E_s_2001(k), MPI_s_01_to_02(k), E_1_2001(k), W_2001(k), MPI_1_01_to_02(k), E_2_2001(k), W_prime_2001(k), MPI_2_01_to_02(k)]);
    fprintf('  \t2002  [00-02]  ');
    disp([E_s_2002(k), MPI_s_00_to_02(k), E_1_2002(k), W_2002(k), MPI_1_00_to_02(k), E_2_2002(k), W_prime_2002(k), MPI_2_00_to_02(k)]);
    fprintf('\n');
    
    sum_E_s_00_to_01 = sum_E_s_00_to_01 + E_s_2000(k);
    sum_E_s_01_to_02 = sum_E_s_01_to_02 + E_s_2001(k);
    sum_E_s_00_to_02 = sum_E_s_00_to_02 + E_s_2002(k);
    sum_E_1_00_to_01 = sum_E_1_00_to_01 + E_1_2000(k);
    sum_E_1_01_to_02 = sum_E_1_01_to_02 + E_1_2001(k);
    sum_E_1_00_to_02 = sum_E_1_00_to_02 + E_1_2002(k);
    sum_E_2_00_to_01 = sum_E_2_00_to_01 + E_2_2000(k);
    sum_E_2_01_to_02 = sum_E_2_01_to_02 + E_2_2001(k);
    sum_E_2_00_to_02 = sum_E_2_00_to_02 + E_2_2002(k);
    
    sum_MPI_s_00_to_01 = sum_MPI_s_00_to_01 + MPI_s_00_to_01(k);
    sum_MPI_s_01_to_02 = sum_MPI_s_01_to_02 + MPI_s_01_to_02(k);
    sum_MPI_s_00_to_02 = sum_MPI_s_00_to_02 + MPI_s_00_to_02(k);
    sum_MPI_1_00_to_01 = sum_MPI_1_00_to_01 + MPI_1_00_to_01(k);
    sum_MPI_1_01_to_02 = sum_MPI_1_01_to_02 + MPI_1_01_to_02(k);
    sum_MPI_1_00_to_02 = sum_MPI_1_00_to_02 + MPI_1_00_to_02(k);
    sum_MPI_2_00_to_01 = sum_MPI_2_00_to_01 + MPI_2_00_to_01(k);
    sum_MPI_2_01_to_02 = sum_MPI_2_01_to_02 + MPI_2_01_to_02(k);
    sum_MPI_2_00_to_02 = sum_MPI_2_00_to_02 + MPI_2_00_to_02(k);
    
    sum_W_00_to_01 = sum_W_00_to_01 + W_2000(k);
    sum_W_01_to_02 = sum_W_01_to_02 + W_2001(k);
    sum_W_00_to_02 = sum_W_00_to_02 + W_2002(k);
    sum_W_prime_00_to_01 = sum_W_prime_00_to_01 + W_prime_2000(k);
    sum_W_prime_01_to_02 = sum_W_prime_01_to_02 + W_prime_2001(k);
    sum_W_prime_00_to_02 = sum_W_prime_00_to_02 + W_prime_2002(k);
end

fprintf(' Average  [00-01]  ');
disp([sum_E_s_00_to_01/num_company, sum_MPI_s_00_to_01/num_company, sum_E_1_00_to_01/num_company, sum_W_00_to_01/num_company,...
        sum_MPI_1_00_to_01/num_company, sum_E_2_00_to_01/num_company, sum_W_prime_00_to_01/num_company, sum_MPI_2_00_to_01/num_company]);
fprintf(' Average  [01-02]  ');
disp([sum_E_s_01_to_02/num_company, sum_MPI_s_01_to_02/num_company, sum_E_1_01_to_02/num_company, sum_W_01_to_02/num_company,...
        sum_MPI_1_01_to_02/num_company, sum_E_2_01_to_02/num_company, sum_W_prime_01_to_02/num_company, sum_MPI_2_01_to_02/num_company]);
fprintf(' Average  [00-02]  ');
disp([sum_E_s_00_to_02/num_company, sum_MPI_s_00_to_02/num_company, sum_E_1_00_to_02/num_company, sum_W_00_to_02/num_company,...
        sum_MPI_1_00_to_02/num_company, sum_E_2_00_to_02/num_company, sum_W_prime_00_to_02/num_company, sum_MPI_2_00_to_02/num_company]);
    
