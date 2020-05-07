clc();
clear();
close all;

## Zadanie 6 a): 
## Oblicz wartoœæ przybli¿on¹ ca³ki 
##  ca³ka  oznaczona w przedziale <0,3>
##  funkcji -x^4 + 2.01x^2 + 1dx
## dla liczby wêz³ów równej 7

## funkcja z regu³¹ simpsona
function s = simpson(f, a, b,  n)
  h = (b-a)/(2*n);
  spa = 0;
  snp = 0;
  ## parzysty sk³adnik sumy
  for i=1:n
    x = a+h*(2*i-1);
    spa += f(x);
  endfor
  ## nieparzysty sk³adnik sumy
  for i=1:(n-1)
    x  = a+h*2*i;
    snp += f(x);
  endfor
  
## ostateczna postaæ wzoru z³o¿onego Simpsona
  s = h*(f(a)+f(b)+4*spa+2*snp)/3;
  printf("Wynik ca³kowania dla %d wêz³ów wynosi: %d \n", n, s); 
  
endfunction;
##przedzia³ ca³ki
a = 0; 
b = 3;
## liczba wêz³ów
n = 7; 
## funkcja
f = @(x) (-(x.^4)+ (2.01*x.^2) + 1); 

simpson(f, a, b, n);

