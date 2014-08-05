N = 200;
mu = 30;
x = (1:1:2e2)';
a =4;

for  lambda = N:-1:1
    P(:, lambda) = exp(-a*lambda)* (a*lambda).^x ./ factorial(x) ...
        * exp(-mu)* mu.^lambda ./ factorial(lambda);
end

p = sum(P, 2);

figure
plot(x, p)