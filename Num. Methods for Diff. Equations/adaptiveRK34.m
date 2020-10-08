% Task 2.4 

function [t, y] = adaptiveRK34(f, y0, t0, tf, tol)
y = y0;  % vektor f�r y-v�rden
t = t0;  % vektor f�r t-v�rden
h = 0;
err = 0;
i = 1;      % startv�rde f�r i
h(i) = (abs(tf-t0).*tol.^(1/4))./(100.*(1+norm(f(t0,y0)))); % startv�rde p� h 
k = 3;      % konstant k, order
[~, err0] = RK34step(f, y(:, i), t(i), h(i)); % tar fram f�rsta error, 
% se upp f�r h(0)

while (tf-t(i))>h(i)                                % while loop s�dan att vi inte hamnar utanf�r tf.
   i=i+1;                                           % �kar med 1 varje g�ng, 
    if i==2
        h(i) = newstep(tol, err(i-1), err0, h(i-1), k);   %nya h skapas med hj�lp av det nya felet, det gamla felet och det gamla h
        [y(:,i), err(i)] = RK34step(f, y(:,i-1), t(i-1), h(i-1)); %tar fram nya y(2) mha redan givna v�rden
        t(i) = t(i-1) + h(i);                             %nytt t skapas 
    else                                            % blir problem med index om den inte g�r vidare ned�t p� h(i)
        h(i) = newstep(tol, err(i-1), err(i-2), h(i-1), k);
        [y(:,i), err(i)] = RK34step(f, y(:,i-1), t(i-1), h(i));
        t(i) = t(i-1) + h(i);
    end
end

[y(:,i), ~] = RK34step(f, y(:,i-1), t(i-1), tf-t(i));
t(i) = t(i-1) + h(i);
plot(t, y);
end
