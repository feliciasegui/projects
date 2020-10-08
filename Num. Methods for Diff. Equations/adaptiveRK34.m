% Task 2.4 

function [t, y] = adaptiveRK34(f, y0, t0, tf, tol)
y = y0;  % vektor för y-värden
t = t0;  % vektor för t-värden
h = 0;
err = 0;
i = 1;      % startvärde för i
h(i) = (abs(tf-t0).*tol.^(1/4))./(100.*(1+norm(f(t0,y0)))); % startvärde på h 
k = 3;      % konstant k, order
[~, err0] = RK34step(f, y(:, i), t(i), h(i)); % tar fram första error, 
% se upp för h(0)

while (tf-t(i))>h(i)                                % while loop sådan att vi inte hamnar utanför tf.
   i=i+1;                                           % ökar med 1 varje gång, 
    if i==2
        h(i) = newstep(tol, err(i-1), err0, h(i-1), k);   %nya h skapas med hjälp av det nya felet, det gamla felet och det gamla h
        [y(:,i), err(i)] = RK34step(f, y(:,i-1), t(i-1), h(i-1)); %tar fram nya y(2) mha redan givna värden
        t(i) = t(i-1) + h(i);                             %nytt t skapas 
    else                                            % blir problem med index om den inte går vidare nedåt på h(i)
        h(i) = newstep(tol, err(i-1), err(i-2), h(i-1), k);
        [y(:,i), err(i)] = RK34step(f, y(:,i-1), t(i-1), h(i));
        t(i) = t(i-1) + h(i);
    end
end

[y(:,i), ~] = RK34step(f, y(:,i-1), t(i-1), tf-t(i));
t(i) = t(i-1) + h(i);
plot(t, y);
end
