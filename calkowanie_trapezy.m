clc();
clear();
close all;

## Zadanie 6 a): 
## Oblicz warto�� przybli�on� ca�ki 
##  ca�ka  oznaczona w przedziale <0,3>
##  funkcji -x^4 + 2.01x^2 + 1dx
## dla liczby w�z��w r�wnej 7

## funkcja z regu�� trapez�w
function t = trapezy(f, a, b, n)
  h = (b-a)/n;
  s =(0.5*f(a) + f(b));
  for i=1:n-1
    s += f(a+i*h); 
  endfor
  t = h*s;
  
  printf("Wynik ca�kowania dla %d w�z��w wynosi: %d \n", n, t); 
endfunction

##przedzia� ca�ki
a = 0; 
b = 3;
## liczba w�z��w
n = 7; 
## funkcja
f = @(x) (-(x.^4)+ (2.01*x.^2) + 1); 

trapezy(f, a, b, n);
