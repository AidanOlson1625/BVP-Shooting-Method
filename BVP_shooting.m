% Run this code to generate a graph of the solution to the spring mass system
% using both Newton's Method and the Bisection Method

    graphs(1, -100, 100)


% We also reccomend running these functions to find the correct initial 
% velocity, p, for a system with initial conditions. We encourage you to 
% experiment with different boundary values and damping coefficients.

    % bisection(0,20)
    % newtons(3)


function p = error(p0)
errTol = 0.005;
err = inf;
max_iter = 200;
it = 1;
while err > errTol && it <= max_iter
    p1 = newtons(p0);
    err = abs(p1-p0);
    p0 = p1;
    it = it + 1;
end
it = it - 1;

function p1 = newtons(p0)
    g = 9.81; L = 2; c = 0.1; alpha = 9*pi/10; beta = -2;
    a = 0; b= 2; n = 20000; h = (b-a)/n;
        f1 = @(w,th) w;
        f2 = @(w,th) -c*w - g/L*sin(th);

        th = zeros(n+1,1);
        w = zeros(n+1,1);
        th(1) = alpha; w(1) = p0;
        
            for i = 1 : n
                k1 = h*f1(w(i),th(i));
                l1 = h*f2(w(i),th(i));
                
                k2 = h*f1(w(i)+k1/2,th(i)+l1/2);
                l2 = h*f2(w(i)+k1/2,th(i)+l1/2);

                k3 = h*f1(w(i)+k2/2,th(i)+l2/2);
                l3 = h*f2(w(i)+k2/2,th(i)+l2/2);

                k4 = h*f1(w(i)+k3,th(i)+l3);
                l4 = h*f2(w(i)+k3,th(i)+l3);
                
           
                th(i+1) = th(i) + (1/6)*(l1 + 2*l2 + 2*l3 + l4);
                w(i+1) = w(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
            end
        f3 = @(z,d) d;
        f4 = @(z,d) -g/L*z - c*d;
        d = zeros(n+1,1);
        z = zeros(n+1,1);
        z(1) = 0; d(1) = 1;
            for i = 1 : n
                k1 = h*f3(z(i),d(i));
                l1 = h*f4(z(i),d(i));

                k2 = h*f3(z(i)+k1/2,d(i)+l1/2);
                l2 = h*f4(z(i)+k1/2,d(i)+l1/2);

                k3 = h*f3(z(i)+k2/2,d(i)+l2/2);
                l3 = h*f4(z(i)+k2/2,d(i)+l2/2);

                k4 = h*f3(z(i)+k3,d(i)+l3);
                l4 = h*f4(z(i)+k3,d(i)+l3);
                
                
                d(i+1) = d(i) + (1/6)*(l1 + 2*l2 + 2*l3 + l4);
                z(i+1) = z(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
            end

         p1 = p0 - (th(end) - beta)/z(end);
    
end
disp('-------------------------------------------')
disp(' ')
disp("Newton's Method:")
disp(['i = ',num2str(it),', p = ',num2str(p1)]);

p = p1;
end

function p_final = bisection(p1, p2)
tol = 0.005;
FA = BVP(p1); FB = BVP(p2);

    if sign(FA)*sign(FB) > 0
        disp(' ')
        disp(['BVP(',num2str(p1),') and BVP(',num2str(p2),') have the same sign.'])
        fprintf('BVP(%d) equals: %d\n',p1, FA)
        fprintf('BVP(%d) equals: %d\n',p2, FB)
        disp('Try a different interval.')
        disp(' ')
        return;
    end
 
 
    for j=1:250
        p = p1 + (p2-p1)/2;
        FP = BVP(p);
        p1_history(j) = p1; p2_history(j) = p2; p_history(j) = p;
        if abs(FP) < tol || FP == 0
            disp('Bisection Method:')
            disp(['i = ',num2str(j),', ','p = ',num2str(p)])
            p_final = p;
            break
        end
   
    
        if sign(FA)*sign(FP) > 0
            p1 = p;
            FA = FP;
        else
            p2 = p;
        end
    end
 
        function F = BVP(p)
        g = 9.81; L = 2; c = 0.1; alpha = 9*pi/10; beta = -2;

        th(1) = alpha; w(1) = p;
        h = 0.01; n = 200;
        
        f1 = @(w,th) w;
        f2 = @(w,th) -c*w - (g/L)*sin(th);
            for i = 1:n
                k1 = h*f1(w(i),th(i));
                l1 = h*f2(w(i),th(i));

                k2 = h*f1(w(i)+k1/2,th(i)+l1/2);
                l2 = h*f2(w(i)+k1/2,th(i)+l1/2);

                k3 = h*f1(w(i)+k2/2,th(i)+l2/2);
                l3 = h*f2(w(i)+k2/2,th(i)+l2/2);

                k4 = h*f1(w(i)+k3,th(i)+l3);
                l4 = h*f2(w(i)+k3,th(i)+l3);
                
                
                th(i+1) = th(i) + (1/6)*(l1 + 2*l2 + 2*l3 + l4);
                w(i+1) = w(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
            end
        
        F = th(n+1) - beta;       
        end
end


function graphs(p_initial_guess, bisection_p1_initial, bisection_p2_initial)
c = 0.1;

function A = final(p)
g = 9.81; L = 2; c = 0.1; b = 2; alpha = 9*pi/10; beta = -2;
        
        f1 = @(w,th) w;
        f2 = @(w,th) -c*w - (g/L)*sin(th);

        th(1) = alpha; w(1) = p;
        h = 0.01; n = 200;
        
            for i = 1:n
                k1 = h*f1(w(i),th(i));
                l1 = h*f2(w(i),th(i));

                k2 = h*f1(w(i)+k1/2,th(i)+l1/2);
                l2 = h*f2(w(i)+k1/2,th(i)+l1/2);

                k3 = h*f1(w(i)+k2/2,th(i)+l2/2);
                l3 = h*f2(w(i)+k2/2,th(i)+l2/2);

                k4 = h*f1(w(i)+k3,th(i)+l3);
                l4 = h*f2(w(i)+k3,th(i)+l3);
                
                
                th(i+1) = th(i) + (1/6)*(l1 + 2*l2 + 2*l3 + l4);
                w(i+1) = w(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
            end
        A = th;
end

n = 200;
t = linspace(0,2,n+1);
plot(t, final(error(p_initial_guess)), 'r', t, final(bisection(bisection_p1_initial, bisection_p2_initial)), 'b')
legend(["Newton's Method", "Bisection Method"]);
xlabel('t');
ylabel('θ(t)');
disp(['c = ', num2str(c)])
end