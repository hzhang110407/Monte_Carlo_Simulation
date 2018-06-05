%set parameters
s0 = 100; %today's price
k = 100; %strike
t = 1; %time to expiry
vol = 0.2; %volatility
r = 0.05; % risk-free rate

steps = 252*10; %num of time steps
n = 100000; %num of times of simulations

[call, put] = bs(100,100,0.05,1,0.2); %the closed form solution of european options
fprintf('The closed form solution of the above European oprions\n Call : %0.4f, Put : %0.4f\n', call,put);
fprintf('\n');
%simulate the price path using euler shceme
X = normrnd(0, 1, [n, steps-1]);
st = zeros(n, steps);
st(:,1) = s0;
dt = t/(steps-1);

for i = 2:steps
    st(:, i) = st(:, i-1) .* (1 + (r - 0.5*vol*vol)*dt + vol*power(dt, 0.5)*X(:,i-1)...
        + 0.5*vol*vol*X(:,i-1).^2*dt );  
    %st(:,i) = st(:,i-1).*exp((r - 0.5*vol*vol)*dt + vol*power(dt,0.5)*X(:,i-1));
end

euro = zeros(1,2);
[euro(1),euro(2)] = european(k, st, r, t);

fprintf('The solution given by monte carlo simulation of the above European options\n Call : %0.4f, Put : %0.4f\n',...
    euro(1), euro(2));
fprintf('\n');

asian_g_clo = asian_g_c(s0, k, r, t, vol, steps, dt);

Asian_a = zeros(1,4);
[Asian_a(1), Asian_a(2), Asian_a(3), Asian_a(4)] = asian_a(k, st, r, t);

Asian_g = zeros(1,4);
[Asian_g(1), Asian_g(2), Asian_g(3), Asian_g(4)] = asian_g(k, st, r, t);

fprintf('The solution given by monte carlo simulation of the above Asian options with arithmetic average\n');
fprintf('Fixed_call : %0.4f, Fixed_put : %0.4f, Floating_call : %0.4f, Floating_put : %0.4f\n',...
    Asian_a(1), Asian_a(2), Asian_a(3), Asian_a(4));
fprintf('\n');
fprintf('the solution given by monte carlo simulation of the above Asian options with  geometric average\n');
fprintf('Fixed_call : %0.4f, Fixed_put : %0.4f, Floating_call : %0.4f, Floating_put : %0.4f\n',...
    Asian_g(1), Asian_g(2), Asian_g(3), Asian_g(4));
fprintf('\n');
fprintf('The closed form solution of the asian fixed call option with geometric average is %0.4f\n', asian_g_clo);
fprintf('\n');

fprintf('Following is some simple error analysis\n Firstly, check the error between BSM and Simulation');
fprintf('The error for European call : %0.4f, the error for European put : %.4f\n\n',abs(call-euro(1)), abs(put-euro(2)));
fprintf('Then for the fixed call with geometric average : %0.4f\n\n', abs(asian_g_clo - Asian_g(1)));

fprintf('check the relationships among different contracts\n');
relation = Asian_a - Asian_g;
disp(relation);
if (relation(1) >= 0 && relation(2) <= 0 && relation(3) <= 0 && relation(4) >=0)
    disp('all satisfied');
else
    disp('not satisfied');
end

fprintf('Following is the possible optimization\n');

error = zeros(4,3);

for idx = 1:4
    n = 10*power(10, idx);
    e = zeros(10,3);
    
    for j = 1:10
        %simulate the price path using euler shceme
        X = normrnd(0, 1, [n, steps-1]);
        st = zeros(n, steps);
        st(:,1) = s0;
        dt = t/(steps-1);

        for i = 2:steps
            st(:, i) = st(:, i-1) .* (1 + (r - 0.5*vol*vol)*dt + vol*power(dt, 0.5)*X(:,i-1)...
                + 0.5*vol*vol*X(:,i-1).^2*dt );  
            %st(:,i) = st(:,i-1).*exp((r - 0.5*vol*vol)*dt + vol*power(dt,0.5)*X(:,i-1));
        end

        euro = zeros(1,2);
        [euro(1),euro(2)] = european(k, st, r, t);
        Asian_g = zeros(1,4);
        [Asian_g(1), Asian_g(2), Asian_g(3), Asian_g(4)] = asian_g(k, st, r, t);
        e(j,1) = abs(call - euro(1));
        e(j,2) = abs(put - euro(2));
        e(j,3) = abs(asian_g_clo - Asian_g(1));
    end
    error(idx,:) = mean(e);
    %disp(error(idx,:));
end

disp(error);

%appling variance reduction tech.

for idx = 1:4
    n = 10*power(10, idx);
    e = zeros(10,3);
    
    for j = 1:10
        %simulate the price path using euler shceme
        X = normrnd(0, 1, [n, steps-1]);
        st1 = zeros(n, steps);
        st1(:,1) = s0;
        st2 = st1;
        dt = t/(steps-1);

        for i = 2:steps
            st1(:, i) = st1(:, i-1) .* (1 + (r - 0.5*vol*vol)*dt + vol*power(dt, 0.5)*X(:,i-1)...
                + 0.5*vol*vol*X(:,i-1).^2*dt );  
            st2(:, i) = st2(:, i-1) .* (1 + (r - 0.5*vol*vol)*dt + vol*power(dt, 0.5)*(-X(:,i-1))...
                + 0.5*vol*vol*X(:,i-1).^2*dt );  
        end

        euro = zeros(1,2);
        price = zeros(2,2);
        [price(1,1), price(1,2)] = european(k, st1, r, t);
        [price(2,1), price(2,2)] = european(k, st2, r, t);
        %[euro(1),euro(2)] = mean(price);
        euro = mean(price);
        price = zeros(2,4);
        [price(1,1),price(1,2),price(1,3),price(1,4)] = asian_g(k, st1, r, t);
        [price(2,1),price(2,2),price(2,3),price(2,4)] = asian_g(k, st2, r, t);
        Asian_g = zeros(1,4);
        %[Asian_g(1), Asian_g(2), Asian_g(3), Asian_g(4)] = mean(price);
        Asian_g = mean(price);
        e(j,1) = abs(call - euro(1));
        e(j,2) = abs(put - euro(2));
        e(j,3) = abs(asian_g_clo - Asian_g(1));
    end
    error(idx,:) = mean(e);
    %disp(error(idx,:));
end

disp(error);

function [call,put] = european(k, st, r, t)
    call = mean(max(st(:,size(st,2)) - k, 0))*exp(-r*t);
    put = mean(max(k - st(:,size(st,2)), 0))*exp(-r*t);
end

function [call_fixed, put_fixed, call_floating, put_floating] = asian_a(k, st, r, t)
    call_fixed = mean(max(mean(st,2) - k,0))*exp(-r*t);
    put_fixed = mean(max(k - mean(st,2), 0))*exp(-r*t);
    call_floating = mean(max(st(:,size(st,2)) - mean(st,2), 0))*exp(-r*t);
    put_floating = mean(max(mean(st,2) - st(:,size(st,2)),0))*exp(-r*t);
end

function [call_fixed, put_fixed, call_floating, put_floating] = asian_g(k, st, r, t)
    call_fixed = mean(max(geomean(st,2) - k,0))*exp(-r*t);
    put_fixed = mean(max(k - geomean(st,2), 0))*exp(-r*t);
    call_floating = mean(max(st(:,size(st,2)) - geomean(st,2), 0))*exp(-r*t);
    put_floating = mean(max(geomean(st,2) - st(:,size(st,2)),0))*exp(-r*t);
end

%a closed form solution for comparasion
function call = asian_g_c(s0, k, r, t, vol, steps, dt)
vol2 = vol*vol;
v1 = 0;
v2 = 0;
steps = steps-1;
n = steps;
for i = 1:steps
   v1 = v1 + (n-i+1)*dt;
   v2 = v2 + power(n-i+1, 2)*dt;
end

v1 = v1/steps;
v2 = v2/power(steps, 2);
%sigma = vol*vol*v2;
E = s0*exp(r*v1 + 0.5*vol2*(v2-v1));
rho = (log(E/k)+0.5*vol2*v2)/(vol*power(v2,0.5));
call = E*exp(-r*t)*normcdf(rho) - k*exp(-r*t)*normcdf(rho - vol*power(v2,0.5));

end

function [call, put] = bs(stock, strike, rate, time, volatility)

d1 = (1/(volatility*power(time, 0.5))) * ( -log(strike/stock) + (rate + power(volatility,2)/2) * time);

d2 = d1 - volatility*power(time, 0.5);

call = normcdf(d1)*stock - normcdf(d2)*exp(-rate*time).*strike;

put = strike*exp(-rate*time) - stock + call;

end

