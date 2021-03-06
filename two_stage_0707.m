% Clear and clc
clear
clc

% Load data file
load '0706_data'

% Number of cities
num_cities = 30;
years = [2004,2005,2006,2007,2008,2009,2010,2011,2013];
num_years = length(years);

X1 = zeros(num_cities, 1);
X2 = zeros(num_cities, 1);
X3 = zeros(num_cities, 1);
X4 = zeros(num_cities, 1);
Z1 = zeros(num_cities, 1);
Z2 = zeros(num_cities, 1);
Z3 = zeros(num_cities, 1);
Y1 = zeros(num_cities, 1);
Y2 = zeros(num_cities, 1);
Y3 = zeros(num_cities, 1);
Y4 = zeros(num_cities, 1);

for i=1:num_cities
    for j=1:num_years
        X1(i) = X1(i) + eval(['X1_',num2str(years(j)),'(i)']);
        X2(i) = X2(i) + eval(['X2_',num2str(years(j)),'(i)']);
        X3(i) = X3(i) + eval(['X3_',num2str(years(j)),'(i)']);
        X4(i) = X4(i) + eval(['X4_',num2str(years(j)),'(i)']);
        Z1(i) = Z1(i) + eval(['Z1_',num2str(years(j)),'(i)']);
        Z2(i) = Z2(i) + eval(['Z2_',num2str(years(j)),'(i)']);
        Z3(i) = Z3(i) + eval(['Z3_',num2str(years(j)),'(i)']);
        Y1(i) = Y1(i) + eval(['Y1_',num2str(years(j)),'(i)']);
        Y2(i) = Y2(i) + eval(['Y2_',num2str(years(j)),'(i)']);
        Y3(i) = Y3(i) + eval(['Y3_',num2str(years(j)),'(i)']);
        Y4(i) = Y4(i) + eval(['Y4_',num2str(years(j)),'(i)']);
    end
end

E_s_ori = zeros(num_cities, 1);
E_1_ori = zeros(num_cities, 1);
E_2_ori = zeros(num_cities, 1);
E_s_total = zeros(num_cities, 1);
E_1_total = zeros(num_cities, 1);
E_2_total = zeros(num_cities, 1);
E_s = zeros(num_cities, num_years);
E_1 = zeros(num_cities, num_years);
E_2 = zeros(num_cities, num_years);
W = zeros(num_cities, num_years);
W_prime = zeros(num_cities, num_years);
MPI_s = zeros(num_cities, num_years);
MPI_1 = zeros(num_cities, num_years);
MPI_2 = zeros(num_cities, num_years);

% Main boody of model(2)
for k=1:num_cities
    % 11 optimization parameters for linear programming
    u1 = optimvar('u1','LowerBound',1e-11, 'Type', 'continuous');
    u2 = optimvar('u2','LowerBound',1e-11, 'Type', 'continuous');
    u3 = optimvar('u3','LowerBound',1e-11, 'Type', 'continuous');
    u4 = optimvar('u4','LowerBound',1e-11, 'Type', 'continuous');
    v1 = optimvar('v1','LowerBound',1e-11, 'Type', 'continuous');
    v2 = optimvar('v2','LowerBound',1e-11, 'Type', 'continuous');
    v3 = optimvar('v3','LowerBound',1e-11, 'Type', 'continuous');
    w1 = optimvar('w1','LowerBound',1e-11, 'Type', 'continuous');
    w2 = optimvar('w2','LowerBound',1e-11, 'Type', 'continuous');
    w3 = optimvar('w3','LowerBound',1e-11, 'Type', 'continuous');
    w4 = optimvar('w4','LowerBound',1e-11, 'Type', 'continuous');
    % template for this ooptimization problem
    prob = optimproblem('ObjectiveSense','maximize');
    % constraints
    prob.Objective = u1 * Y1(k) + u2 * Y2(k) + u3 * Y3(k) + u4 * Y4(k);
    prob.Constraints.cons0 = v1 * X1(k) + v2 * X2(k) + v3 * X3(k) == 1;
    
    for j = 1:num_cities
        eval(['prob.Constraints.cons',num2str(j),'=','u1 * Y1(j) + u2 * Y2(j) + u3 * Y3(j) + u4 * Y4(j) - v1 * X1(j) - v2 * X2(j) - v3 * X3(j) <= 0']);
        eval(['prob.Constraints.consx',num2str(j),'=','w1 * Z1(j) + w2 * Z2(j) + w3 * Z3(j) - v1 * X1(j) - v2 * X2(j) - v3 * X3(j) <= 0']);
        eval(['prob.Constraints.consy',num2str(j),'=','u1 * Y1(j) + u2 * Y2(j) + u3 * Y3(j) + u4 * Y4(j) - w1 * Z1(j) - w2 * Z2(j) - w3 * Z3(j) - w4* X4(j) <= 0']);
    end
    % Solve调用求解器进行求解
    sol = solve(prob);
    % The left part 3 column of Table2, Aggregate Model results
    E_s_ori(k) = (sol.u1 * Y1(k) + sol.u2 * Y2(k) + sol.u3 * Y3(k) + sol.u4 * Y4(k))/(sol.v1 * X1(k) + sol.v2 * X2(k) + sol.v3 * X3(k));
    E_1_ori(k) = (sol.w1 * Z1(k) + sol.w2 * Z2(k) + sol.w3 * Z3(k))/(sol.v1 * X1(k) + sol.v2 * X2(k) + sol.v3 * X3(k));
    E_2_ori(k) = (sol.u1 * Y1(k) + sol.u2 * Y2(k) + sol.u3 * Y3(k) + sol.u4 * Y4(k))/(sol.w1 * Z1(k) + sol.w2 * Z2(k) + sol.w3 * Z3(k) + sol.w4 * X4(k));
end


% Main boody of model(4), model(5), model(7) and model(8)
for k=1:num_cities
    % 11 optimization parameters for linear programming
    u1 = optimvar('u1','LowerBound',1e-11, 'Type', 'continuous');
    u2 = optimvar('u2','LowerBound',1e-11, 'Type', 'continuous');
    u3 = optimvar('u3','LowerBound',1e-11, 'Type', 'continuous');
    u4 = optimvar('u4','LowerBound',1e-11, 'Type', 'continuous');
    v1 = optimvar('v1','LowerBound',1e-11, 'Type', 'continuous');
    v2 = optimvar('v2','LowerBound',1e-11, 'Type', 'continuous');
    v3 = optimvar('v3','LowerBound',1e-11, 'Type', 'continuous');
    w1 = optimvar('w1','LowerBound',1e-11, 'Type', 'continuous');
    w2 = optimvar('w2','LowerBound',1e-11, 'Type', 'continuous');
    w3 = optimvar('w3','LowerBound',1e-11, 'Type', 'continuous');
    w4 = optimvar('w4','LowerBound',1e-11, 'Type', 'continuous');
    % template for this ooptimization problem
    prob = optimproblem('ObjectiveSense','maximize');
    % constraints
    prob.Objective = u1 * Y1(k) + u2 * Y2(k) + u3 * Y3(k) + u4 * Y4(k);
    prob.Constraints.cons0 = v1 * X1(k) + v2 * X2(k) + v3 * X3(k) == 1;
    for i = 1:num_years
        for j = 1:num_cities
            year = num2str(years(i));
            eval(['prob.Constraints.cons_',year,'_',num2str(j),'=','u1 * Y1_',year,'(j)','+ u2 * Y2_',year, '(j)', '+ u3 * Y3_', year, '(j)', '+ u4 * Y4_', year, '(j)', '- v1 * X1_', year, '(j)', '- v2 * X2_', year, '(j)', '- v3 * X3_', year, '(j) <= 0']);
            eval(['prob.Constraints.consx_',year,'_',num2str(j),'=','w1 * Z1_',year,'(j)','+ w2 * Z2_',year, '(j)', '+ w3 * Z3_', year, '(j)', '- v1 * X1_', year, '(j)', '- v2 * X2_', year, '(j)', '- v3 * X3_', year, '(j) <= 0']);
            eval(['prob.Constraints.consy_',year,'_',num2str(j),'=','u1 * Y1_',year,'(j)','+ u2 * Y2_',year, '(j)', '+ u3 * Y3_', year, '(j)', '+ u4 * Y4_', year, '(j)', '- w1 * Z1_', year, '(j)', '- w2 * Z2_', year, '(j)', '- w3 * Z3_', year, '(j)', '- w4 * X4_', year, '(j) <= 0']);
        end
    end
    % Solve调用求解器进行求解
    sol = solve(prob);
%     disp([sol.u1, sol.u2, sol.u3, sol.u4, sol.v1, sol.v2, sol.v3, sol.w1, sol.w2, sol.w3, sol.w4]);
    E_s_total(k) = (sol.u1 * Y1(k) + sol.u2 * Y2(k) + sol.u3 * Y3(k) + sol.u4 * Y4(k))/(sol.v1 * X1(k) + sol.v2 * X2(k) + sol.v3 * X3(k));
    E_1_total(k) = (sol.w1 * Z1(k) + sol.w2 * Z2(k) + sol.w3 * Z3(k))/(sol.v1 * X1(k) + sol.v2 * X2(k) + sol.v3 * X3(k));
    E_2_total(k) = (sol.u1 * Y1(k) + sol.u2 * Y2(k) + sol.u3 * Y3(k) + sol.u4 * Y4(k))/(sol.w1 * Z1(k) + sol.w2 * Z2(k) + sol.w3 * Z3(k) + sol.w4 * X4(k));
    
    for j = 1:num_years
        year = num2str(years(j));
        E_s(k,j) = (sol.u1 * eval(['Y1_',year,'(k)']) + sol.u2 * eval(['Y2_',year,'(k)'])+ sol.u3 * eval(['Y3_',year,'(k)'])+ sol.u4 * eval(['Y4_',year,'(k)']))/(sol.v1 * eval(['X1_',year,'(k)']) + sol.v2 * eval(['X2_',year,'(k)'])+ sol.v3 * eval(['X3_',year,'(k)']));
        E_1(k,j) = (sol.w1 * eval(['Z1_',year,'(k)']) + sol.w2 * eval(['Z2_',year,'(k)'])+ sol.w3 * eval(['Z3_',year,'(k)']))/(sol.v1 * eval(['X1_',year,'(k)']) + sol.v2 * eval(['X2_',year,'(k)'])+ sol.v3 * eval(['X3_',year,'(k)']));
        E_2(k,j) = (sol.u1 * eval(['Y1_',year,'(k)']) + sol.u2 * eval(['Y2_',year,'(k)'])+ sol.u3 * eval(['Y3_',year,'(k)'])+ sol.u4 * eval(['Y4_',year,'(k)']))/(sol.w1 * eval(['Z1_',year,'(k)']) + sol.w2 * eval(['Z2_',year,'(k)'])+ sol.w3 * eval(['Z3_',year,'(k)'])+ sol.w4 * eval(['X4_',year,'(k)']));
        W(k,j) = (sol.v1 * eval(['X1_',year,'(k)']) + sol.v2 * eval(['X2_',year,'(k)'])+ sol.v3 * eval(['X3_',year,'(k)']))/(sol.v1 * X1(k) + sol.v2 * X2(k) + sol.v3 * X3(k));
        W_prime(k,j) = (sol.w1 * eval(['Z1_',year,'(k)']) + sol.w2 * eval(['Z2_',year,'(k)'])+ sol.w3 * eval(['Z3_',year,'(k)']))/(sol.w1 * Z1(k) + sol.w2 * Z2(k) + sol.w3 * Z3(k));
    end
    
    for j = 1 : num_years - 1
        year = num2str(years(j));
        MPI_s(k,j) = E_s(k,j+1)/E_s(k,j);
        MPI_1(k,j) = E_1(k,j+1)/E_1(k,j);
        MPI_2(k,j) = E_2(k,j+1)/E_2(k,j);
    end
    
    MPI_s(k,num_years) = E_s(k,num_years)/E_s(k,1);
    MPI_1(k,num_years) = E_1(k,num_years)/E_1(k,1);
    MPI_2(k,num_years) = E_2(k,num_years)/E_2(k,1); 
end

% compute average in Table2
sum_s_ori = 0;
sum_1_ori = 0;
sum_2_ori = 0;
sum_s_total = 0;
sum_1_total = 0;
sum_2_total = 0;
fprintf('\t\tThe following are results of Multi-period Model in Paper.Table2\n\n');
fprintf('\t\t\t    Aggregate Model\t\t\t     Multi-period Model\n\n');
fprintf('\t       E_s_k     E_1_k     E_2_k     E_s_k    E_1_k     E_2_k\n\n');
for k=1:num_cities
    fprintf('\t%d ',k);
    disp([E_s_ori(k), E_1_ori(k), E_2_ori(k), E_s_total(k), E_1_total(k), E_2_total(k)]);
    sum_s_ori = sum_s_ori + E_s_ori(k);
    sum_1_ori = sum_1_ori + E_1_ori(k);
    sum_2_ori = sum_2_ori + E_2_ori(k);
    sum_s_total = sum_s_total + E_s_total(k);
    sum_1_total = sum_1_total + E_1_total(k);
    sum_2_total = sum_2_total + E_2_total(k);
end
fprintf('Average');
disp([sum_s_ori/num_cities,sum_1_ori/num_cities,sum_2_ori/num_cities, sum_s_total/num_cities, sum_1_total/num_cities, sum_2_total/num_cities]);

sum_E_s = zeros(num_years, 1);
sum_E_1 = zeros(num_years, 1);
sum_E_2 = zeros(num_years, 1);
sum_MPI_s = zeros(num_years, 1);
sum_MPI_1 = zeros(num_years, 1);
sum_MPI_2 = zeros(num_years, 1);
sum_W = zeros(num_years, 1);
sum_W_prime = zeros(num_years, 1);

% Display results in Table3
fprintf('\t\t The following are results of Multi-period Model\n\n');
fprintf('  \tYear  Period\t   E_s_k     MPI_s_k    E_1_k    W_1_k    MPI_1_k    E_2_k     W_2_k     MPI_2_k\n\n');
for k=1:num_cities
    for j = 1 : num_years
        year = years(j);
        if (j == num_years)
            fprintf('  \t%d  [0%d-%d]  ',years(num_years), years(1)-2000, years(num_years)-2000);
        elseif(year -1999 < 10)
            fprintf('  \t%d  [0%d-0%d]  ',year, year-2000, year-1999);
        elseif (year -2000 < 10)
            fprintf('  \t%d  [0%d-%d]  ',year, year-2000, year-1999);
        else
            fprintf('  \t%d  [%d-%d]  ',year, year-2000, year-1999);
        end
        disp([E_s(k,j), MPI_s(k,j), E_1(k,j), W(k,j), MPI_1(k,j), E_2(k,j), W_prime(k,j), MPI_2(k,j)]);
        
        sum_E_s(j) = sum_E_s(j) + E_s(k, j);
        sum_E_1(j) = sum_E_1(j) + E_1(k, j);
        sum_E_2(j) = sum_E_2(j) + E_2(k, j);
        
        sum_MPI_s(j) = sum_MPI_s(j) + MPI_s(k, j);
        sum_MPI_1(j) = sum_MPI_1(j) + MPI_1(k, j);
        sum_MPI_2(j) = sum_MPI_2(j) + MPI_2(k, j);
        
        sum_W(j) = sum_W(j) + W(k,j);
        sum_W_prime(j) = sum_W_prime(j) + W_prime(k,j);
    end
    fprintf('\n');
end

for j = 1:num_years
    year = years(j);
    if(j == num_years)
        fprintf(' Average  [0%d-%d]  ', years(1)-2000, years(num_years)-2000);
    elseif(year -1999 < 10)
        fprintf(' Average  [0%d-0%d]  ', year-2000, year-1999);
    elseif (year -2000 < 10)
        fprintf(' Average  [0%d-%d]  ', year-2000, year-1999);
    else
        fprintf(' Average  [%d-%d]  ', year-2000, year-1999);
    end
    disp([sum_E_s(j)/num_cities, sum_MPI_s(j)/num_cities, sum_E_1(j)/num_cities, sum_W(j)/num_cities,...
        sum_MPI_1(j)/num_cities, sum_E_2(j)/num_cities, sum_W_prime(j)/num_cities, sum_MPI_2(j)/num_cities]);
end

        