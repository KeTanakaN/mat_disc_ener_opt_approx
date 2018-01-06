clear all;
close all;
dir_smp = './DataSmp/'; % Directory to save the sampling points
dir_err = './DataErr/'; % Directory to save the errors 

%% This program uses ADVANPIX Multiprecision Computing Toolbox
% http://www.advanpix.com/documentation/users-manual/
% addpath('C:\Users\ketanaka\Documents\Multiprecision Computing Toolbox');
addpath('C:\Users\KenichiroTanaka\Documents\Multiprecision Computing Toolbox');
mp.Digits(75);

m_pi = mp('pi'); % MULTI-precision 'pi'

%% Functions etc.
% -----------------------------------
% (1) (The example in the previous paper (Single Exponential))
% d = m_pi/4;
% Q = @(x)log(cosh(2*x));
% w = @(x) exp(-Q(x));
% f = @(x) w(x);
% x = mp([-25:0.05:25]); % MULTI
% prefix = 'w1_se'; % for filenames
% -----------------------------------
% (2) (The example in the previous paper (Gaussian))
% d = m_pi/4;
% Q = @(x) (x.^2);
% w = @(x) exp(-Q(x));
% f = @(x) w(x) .* (x.^2)./((m_pi^2)/16+x.^2);
% x = mp([-10:0.02:10]); % MULTI
% prefix = 'w2_ga'; % for filenames
% -----------------------------------
% (3) (The example in the previous paper (Double Exponential))
% d = m_pi/4;
% Q = @(x) log(cosh((m_pi/2)*sinh(2*x)));
% w = @(x) exp(-Q(x));
% f = @(x) w(x);
% x = mp([-3:0.006:3]); % MULTI
% prefix = 'w3_de'; % for filenames
% -----------------------------------
% (4) (The SE (TANH) transform of (1-t^2)^(1/2)*(1+t^2))
d = m_pi;
Q = @(x) log(cosh(x/2));
w = @(x) exp(-Q(x));
f = @(x) w(x) .* (1+tanh(x/2).^2);
x = mp([-100:0.2:100]); % MULTI
prefix = 'w4_se'; % for filenames
% -----------------------------------
% (5) (The DE transform of (1-t^2)^(1/2)*(1+t^2))
% d = m_pi/2;
% Q = @(x) log(cosh((m_pi/2)*sinh(x)));
% w = @(x) exp(-Q(x));
% f = @(x) w(x) .* (1+tanh((m_pi/2)*sinh(x)).^2);
% x = mp([-6:0.012:6]); % MULTI
% prefix = 'w5_de'; % for filenames
% -----------------------------------
% (6) (The SE (TANH) transform of (1-t)^(1/2)*(1+t)^(3/2)*(1+t^2))
% d = m_pi;
% Q = @(x) ((1/2)*log(1+exp(x)) + (3/2)*log(1+exp(-x)));
% w = @(x) exp(-Q(x));
% f = @(x) w(x) .* 4 .* (1+tanh(x/2).^2); 
% x = mp([-40:0.14:100]); % MULTI
% prefix = 'w6_se'; % for filenames
% -----------------------------------
% (7) (The DE transform of (1-t)^(1/2)*(1+t)^(3/2)*(1+t^2))
% d = m_pi/2;
% Q = @(x) ((1/2)*log(1+exp(m_pi*sinh(x))) + (3/2)*log(1+exp(-m_pi*sinh(x))));
% w = @(x) exp(-Q(x));
% f = @(x) w(x) .* 4 .* (1+tanh((m_pi/2)*sinh(x)).^2);
% x = mp([-4.5:0.01:5.5]); % MULTI
% prefix = 'w7_de'; % for filenames
% -----------------------------------

%% Computation of approximants
d  = d - mp(10^(-10)); % MULTI
Cd = m_pi/(4*d);

err = zeros(10,1);
i = 1;
for n = 21:20:201 % The number of the sampling points
    filename = strcat(dir_smp, prefix, '_n_', num2str(n), '.txt');
    a = mp(dlmread(filename)); % MULTI (The samling points)

    tic;

    % Computation of the coefficients (lambda) for the barycentric formulas
    mu     = mp(ones(n,1)); % MULTI
    lambda = mp(ones(n,1)); % MULTI
    for k=1:n
        for j=1:n
            if j ~= k
                mu(k) = mu(k) * tanh(Cd*(a(k) - a(j)));
            end
        end
        lambda(k) = mu(k)^(-1);
    end

    % Computation of the approximation by the barycentric formula (1)
    % appr = APPR_tB_baryc_one(x, lambda, a, f, w, Cd);
    % postfix = 'baryc_one';

    % Computation of the approximation by the barycentric formula (2)
    appr = APPR_tB_baryc_two(x, lambda, a, f, w, Cd);
    postfix = 'baryc_two';

    toc;

    % Computation of the error
    err(i) = max(abs(f(x) - appr));
    i = i+1;
end

% Output of the errors
filename = strcat(dir_err, prefix, '_err_', postfix, '.txt');
% dlmwrite(filename, err);
% type(filename);
