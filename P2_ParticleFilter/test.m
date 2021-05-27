estConst = EstimatorConst;
N = 100000;
mu_0 = [[estConst.pA(1); estConst.pB(1)]; [estConst.pA(2); estConst.pB(2)]; 0; 0];
var_0 = [estConst.d; estConst.d; estConst.phi_0; estConst.l];
x_0 = sampleInitial(mu_0, var_0, N);
scatter(x_0(:, 1), x_0(:, 2)); 
%%  
function x_bar = sampleInitial(mu, var, N) %Sample uniform around (p0-x_bar) <= d 
    %Stupid implementation for now, fix later
    stateSize = numel(var);
    u_bar = rand(N, stateSize); %Mutual independence of the initial conditions
    x_bar = zeros(N, stateSize);
    for i = 1:N
        decision = rand;
        for j = 1:stateSize
            if j==1
                if decision < 0.5
                    x_bar(i, j) = var(j)*(2*u_bar(i, j) - 1) + mu(j);
                else
                    x_bar(i, j) = var(j)*(2*u_bar(i, j) - 1) + mu(j+1);
                end
            elseif j==2
                bound = 2*u_bar(i, j-1) - 1;
                if decision < 0.5
                    x_bar(i, j) = sqrt(var(j)*var(j)*(1 - bound*bound))*(2*u_bar(i, j) - 1) + mu(j+1);
                else
                    x_bar(i, j) = sqrt(var(j)*var(j)*(1 - bound*bound))*(2*u_bar(i, j) - 1) + mu(j+2);
                end
            else
                x_bar(i, j) = var(j)*(2*u_bar(i, j) - 1) + mu(j+2);
            end
        end
    end
end