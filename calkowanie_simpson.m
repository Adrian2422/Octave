clc();
clear();
close all;

## Zadanie 6 a): 
## Oblicz warto�� przybli�on� ca�ki 
##  ca�ka  oznaczona w przedziale <0,3>
##  funkcji -x^4 + 2.01x^2 + 1dx
## dla liczby w�z��w r�wnej 7

## funkcja z regu�� simpsona
function s = simpson(f, a, b,  n)
  h = (b-a)/(2*n);
  spa = 0;
  snp = 0;
  ## parzysty sk�adnik sumy
  for i=1:n
    x = a+h*(2*i-1);
    spa += f(x);
  endfor
  ## nieparzysty sk�adnik sumy
  for i=1:(n-1)
    x  = a+h*2*i;
    snp += f(x);
  endfor
  
## ostateczna posta� wzoru z�o�onego Simpsona
  s = h*(f(a)+f(b)+4*spa+2*snp)/3;
  printf("Wynik ca�kowania dla %d w�z��w wynosi: %d \n", n, s); 
  
endfunction;
##przedzia� ca�ki
a = 0; 
b = 3;
## liczba w�z��w
n = 7; 
## funkcja
f = @(x) (-(x.^4)+ (2.01*x.^2) + 1); 

simpson(f, a, b, n);

