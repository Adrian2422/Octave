clc();
clear();
close all;

x0 = [-3 -2 -1 0 1 2 3];
y0 = [-15 1 5 3 1 5 21]; 
n = length(x0);

function y = lagrange(n, x, x0, y0)
## x0 - wektor z wêz³ami x
## y0 - wektor z wartoœciami dla x0
## y - obliczona wartosc
y = 0;
  for i=1:n
      p = 1;
      for j=1:n
          if j == i   # zapobiega dzieleniu przez 0
              continue;
          endif;
          p *= (x-x0(j)) / (x0(i)-x0(j));
      endfor;
      y += y0(i) * p;   
  endfor;
endfunction

## pomiêdzy -3,3 ma byæ 100 punktów, w których x0(1-7) ma byæ wyrysowane
s = linspace(-5, 5); 
m = length(s);  ## m = 100

## pêtla obliczaj¹ca y dla 100 x z linespace
yi = zeros(1,m);
for i=1:m
  y = lagrange(n, s(i), x0, y0);
  yi(i) += y;
endfor


## interpolacja punktu x
x = 2.25;
yx = lagrange(n, x, x0, y0);


plot(s, yi, 'b-;funkcja;',x0, y0, 'ro;wêz³y;', x, yx, 'g*;zadany x;');
axis([-4, 4, -17, 23]);
  title('funkcja x0');
  xlabel('x');
  ylabel('y');

yi



