clear all;
close all;
dir_smp = './DataSmp/'; % The directory to save the sampling points

%% Settings
% The external field Q (the minus logarithm of a weight function) and its derivatives
% -----------------------------------
% (1) (The example in the previous paper (Single Exponential))
% Q   = @(x)log(cosh(2*x));
% d1Q = @(x)(2*tanh(2*x));
% d2Q = @(x)(4*sech(2*x)^2);
% d = pi/4;
% prefix = 'w1_se'; % for filenames
% Num = '1'; % for a picture
% -----------------------------------
% (2) (The example in the previous paper (Gaussian))
% Q   = @(x) (x.^2);
% d1Q = @(x) (2*x);
% d2Q = @(x) 2;
% d = pi/4;
% prefix = 'w2_ga'; % for filenames
% Num = '2'; % for a picture
% -----------------------------------
% (3) (The example in the previous paper (Double Exponential))
% Q   = @(x)log(cosh((pi/2)*sinh(2*x)));
% d1Q = @(x)(pi*cosh(2*x).*tanh((pi/2)*sinh(2*x)));
% d2Q = @(x)(pi^2*cosh(2*x).^2.*sech((pi/2)*sinh(2*x)).^2 + 2*pi*sinh(2*x).*tanh((pi/2)*sinh(2*x)));
% d = pi/4;
% prefix = 'w3_de'; % for filenames
% Num = '3'; % for a picture
% -----------------------------------
% (4) (The SE (TANH) transform of (1-t^2)^(1/2))
% Q   = @(x) log(cosh(x/2));
% d1Q = @(x) (tanh(x/2)/2);
% d2Q = @(x) (sech(x/2)^2)/4;
% d = pi;
% prefix = 'w4_se'; % for filenames
% Num = '4'; % for a picture
% -----------------------------------
% (5) (The DE transform of (1-t^2)^(1/2))
% Q   = @(x) log(cosh((pi/2)*sinh(x)));
% d1Q = @(x) ((pi/2)*cosh(x).*tanh((pi/2)*sinh(x)));
% d2Q = @(x) ((pi/2)^2*cosh(x).^2.*sech((pi/2)*sinh(x)).^2 + (pi/2)*sinh(x).*tanh((pi/2)*sinh(x)));
% d = pi/2;
% prefix = 'w5_de'; % for filenames
% Num = '5'; % for a picture
% -----------------------------------
% (6) (The SE (TANH) transform of (1-t)^(1/2)*(1+t)^(3/2))
% Q   = @(x) ((1/2)*log(1+exp(x)) + (3/2)*log(1+exp(-x)));
% d1Q = @(x) (1/2 - 2./(1+exp(x)));
% d2Q = @(x) 1./(1+cosh(x));
% d = pi;
% prefix = 'w6_se'; % for filenames
% Num = '6'; % for a picture
% -----------------------------------
% (7) (The DE transform of (1-t)^(1/2)*(1+t)^(3/2))
Q   = @(x) ((1/2)*log(1+exp(pi*sinh(x))) + (3/2)*log(1+exp(-pi*sinh(x))));
d1Q = @(x) (pi/2)*cosh(x).*(1 - 4./(1+exp(pi*sinh(x))));
d2Q = @(x) (pi/2)*(pi*(cosh(x).*sech((pi/2)*sinh(x))).^2 + sinh(x) - 4*sinh(x)./(1+exp(pi*sinh(x))));
d = pi/2;
prefix = 'w7_de'; % for filenames
Num = '7'; % for a picture
% -----------------------------------

% The function K for the kernel and its derivatives
d   = d - 10^(-10);
Cd  = pi/(4*d);
K   = @(x)(-log(tanh(Cd * abs(x))));
d1K = @(x)(-2*Cd*csch(2*Cd*x));
d2K = @(x)(Cd^2*(csch(Cd*x).^2 + sech(Cd*x).^2));

%% Genaration of the sampling points
for n = 21:20:201 % The number of sampling points

    % Newton's method
    h = 0.1;
    a = ([-floor(n/2):(-floor(n/2)+n-1)]*h)'; % Initial sequence

    eps_tol = 1e-14;
    
    tic;
    for newton = 1:1000
        Grad = d1Q(a)*(n-1)/n;
        for i=1:n
            for j=1:n
                if j~=i
                    Grad(i) = Grad(i) + d1K(a(i) - a(j));
                end
            end
        end
        
        Hess = zeros(n, n);
        for i=1:n
            for j=1:n
                if j~=i
                    Hess(i,j) = - d2K(a(i) - a(j));
                end
            end
            Hess(i,i) = - sum(Hess(i,:)) + d2Q(a(i))*(n-1)/n;
        end
        
        b = Hess\Grad;
        a = a - b;
        eps = max(abs(b)) % Display the value for confirmation. 
        if eps < eps_tol
            break;
        end
    end
    toc;
    
    % The potential with the external field
    dx = 1.5*(a(n) - a(1))/(5*n);
    x = [1.5*a(1):dx:1.5*a(n)];
    one_col = ones(length(a),1);
    one_row = ones(1,length(x));
    A = kron(a,one_row);
    X = kron(x,one_col);
    
    tmp_arr = K(X - A);
    for j=1:length(x);
        dp(j) = (sum(tmp_arr(:,j)) + Q(x(j)));
    end
    
    % The optimal energy and modified Robin constant
    E = 0;
    for i=1:n
        for j=1:n
            if j~=i
                E = E + K(a(i) - a(j));
            end
        end
    end
    Qsum = sum(Q(a))*(n-1)/n;
    E = E + 2*Qsum;
    F = E - Qsum;
    
    % Output of the potential
    f = ones(1,length(x))*(F/n);
    disp('F/n =');disp(F/n);
    plot(x,dp,'Linewidth',2);
    hold on;
    plot(x,f,'--r','Linewidth', 2);
    ylim([F/n-5, F/n+15]);
    set(gca,'FontName','Times','FontSize',12,'FontWeight','bold');
    xlabel('x');
    ylabel('Function values');
    title(['Function U_{n}^{D}(x) + Q(x) for weight ', Num ,'  (n = ',num2str(n),')']);
    legend('Function', 'F_{K,Q}^{D}(n)/n','Location', 'SouthEast');
    grid on;
    hold off;
    
    pause;
    
    % Output of the sampling points
    filename = strcat(dir_smp, prefix, '_n_', num2str(n), '.txt');
    % dlmwrite(filename,a,'precision',15);

end
