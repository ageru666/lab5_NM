x = linspace(-1, 3, 50);
y = function1(x); 
xi = linspace(min(x), max(x), 50);

lagrangePoly = polynomialLagrange(x, y);
newtonPoly = polynomialNewton(x, y);
splinePoly = evaluateCubicSpline(x, y, xi);

lagrangeStr = polyToString(lagrangePoly);
disp(['Рівняння полінома Лагранжа: ', lagrangeStr]);

newtonStr = polyToString(newtonPoly);
disp(['Рівняння поліномаа Ньютона: ', newtonStr]);

cubicStr = polyToString(splinePoly);
disp(['Рівняння полінома кубика: ', cubicStr]);

% Plotting
xi = linspace(min(x), max(x), 1000);
yi_lagrange = polyval(lagrangePoly, xi);
yi_newton = polyval(newtonPoly, xi);
yi = evaluateCubicSpline(x, y, xi);

subplot(1, 3, 1);
plot(xi, yi_lagrange);
title('Lagrange');

subplot(1, 3, 2);
plot(xi, yi_newton);
title('Newton ');

subplot(1, 3, 3);
plot(xi, yi);
title('Cubic Spline');

function y = function1(x)
    y = 2 * x.^2 - x + 1;
end

function y = function2(x)
    y = 2 * sin(3*x);
end

function y = function3(x)
    y = x .^ 2 - 3 * sin(3 * x);
end

function str = polyToString(p)
  
    n = length(p);
    str = '';
    for i = 1:n
        coef = p(i);
        if coef == 0
            continue;
        end

        if i ~= 1 && coef > 0
            str = [str, ' + '];
        elseif coef < 0
            str = [str, ' - '];
            coef = -coef;
        end

        degree = n - i;
        if degree == 0
            str = [str, num2str(coef)];
        elseif degree == 1
            str = [str, num2str(coef), 'x'];
        else
            str = [str, num2str(coef), 'x^', num2str(degree)];
        end
    end
end

function p = polynomialLagrange(x, y)
    n = length(x);
    p = 0;
    for i = 1:n
        L = 1;
        for j = 1:n
            if j ~= i
                L = conv(L, poly([x(j)])) / (x(i) - x(j));
            end
        end
        p = p + y(i) * L;
    end
end

function p = polynomialNewton(x, y)
    n = length(x);
    diffTable = dividedDifferences(x, y);
    p = diffTable(1,1);
    for i = 2:n
        term = diffTable(1,i);
        polyTerm = 1;
        for j = 1:(i-1)
            polyTerm = conv(polyTerm, poly([x(j)]));
        end
        term = term * polyTerm;
        
       
        if length(term) > length(p)
            p = [zeros(1, length(term) - length(p)), p];
        elseif length(term) < length(p)
            term = [zeros(1, length(p) - length(term)), term]; 
        end
        
        p = p + term;
    end
end

function diffTable = dividedDifferences(x, y)
    n = length(x);
    diffTable = zeros(n,n);
    diffTable(:,1) = y';
    for j = 2:n
        for i = 1:(n-j+1)
            diffTable(i,j) = (diffTable(i+1,j-1) - diffTable(i,j-1)) / (x(i+j-1) - x(i));
        end
    end
end

function z = cubicSplineCoeffs(x, y)
    n = length(x);
    h = diff(x);

    dy = diff(y) ./ h;
    
    w = zeros(n-1, 1);
    z = zeros(n, 1);

    for i = 2:n-1
        m = h(i-1) * (2 - w(i-1)) + 2 * h(i);
        w(i) = h(i) / m;
        z(i) = (6 * (dy(i) - dy(i-1)) - h(i-1) * z(i-1)) / m;
    end
    
    for i = n-1:-1:1
        z(i) = z(i) - w(i) * z(i+1);
    end
end

function result = evaluateCubicSpline(x, y, xi)
    z = cubicSplineCoeffs(x, y);
    n = length(x);
    result = zeros(size(xi));

    for k = 1:length(xi)
        i = find(x <= xi(k), 1, 'last');
        if i >= n, i = n-1; end
        if i < 1, i = 1; end

        h = x(i+1) - x(i);
        dx = xi(k) - x(i);

        result(k) = z(i)/6/h * (x(i+1) - xi(k))^3 ...
                  + z(i+1)/6/h * (xi(k) - x(i))^3 ...
                  + (y(i+1)/h - z(i+1)*h/6) * dx ...
                  + (y(i)/h - z(i)*h/6) * (x(i+1) - xi(k));
    end
end


