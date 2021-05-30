estConst = EstimatorConst;
N = 10000;
mu_0 = [[estConst.pA(1); estConst.pB(1)]; [estConst.pA(2); estConst.pB(2)]; 0; 0];
var_0 = [estConst.d; estConst.d; estConst.phi_0; estConst.l];
x_0 = sampleInitial(mu_0, var_0, N);
scatter(x_0(:, 1), x_0(:, 2)); 

%% 
estConst = EstimatorConst;
x = -4*estConst.epsilon:estConst.epsilon/100:4*estConst.epsilon;
y = zeros(1, numel(x));
for i = 1:numel(x)
    y(i) = pdf_wk(x(i), estConst.epsilon);
end
plot(x, y);
%%

out = g(@(x)(f(x, 2)), 3);
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


function pwk = pdf_wk(x, epsilon)
    if x<0
        x = -x;
    end
    if x<=2*epsilon
        pwk = (2/(5*epsilon)) * (1 - (1/(2*epsilon)) * x);
    elseif (x>2*epsilon) && (x<= 2.5*epsilon)
        pwk = (2/(5*epsilon*epsilon)) * (x - 2*epsilon);
    elseif (x>2.5*epsilon) && (x<=3*epsilon)
        pwk = 1/(5*epsilon) - (2/(5*epsilon*epsilon)) * (x - 2.5*epsilon);
    else
        pwk = 0;
    end
end


function y = f(x, a)
    y = a*x^2;
end

function z = g(f, b);
    z = f(1) + b;
end