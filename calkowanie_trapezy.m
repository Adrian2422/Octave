clc();
clear();
close all;

## Zadanie 6 a): 
## Oblicz wartoœæ przybli¿on¹ ca³ki 
##  ca³ka  oznaczona w przedziale <0,3>
##  funkcji -x^4 + 2.01x^2 + 1dx
## dla liczby wêz³ów równej 7

## funkcja z regu³¹ trapezów
function t = trapezy(f, a, b, n)
  h = (b-a)/n;
  s =(0.5*f(a) + f(b));
  for i=1:n-1
    s += f(a+i*h); 
  endfor
  t = h*s;
  
  printf("Wynik ca³kowania dla %d wêz³ów wynosi: %d \n", n, t); 
endfunction

##przedzia³ ca³ki
a = 0; 
b = 3;
## liczba wêz³ów
n = 7; 
## funkcja
f = @(x) (-(x.^4)+ (2.01*x.^2) + 1); 

trapezy(f, a, b, n);
