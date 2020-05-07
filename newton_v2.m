clc();
clear();
close all;

## Poprawiona wersja programu do interpolacji Newtona, teraz mo�na zmienia� dowolnie ilo�� w�z��w.
## Ilorazy r�nicowe s� przechowywane w macierzy, co u�atwia manipulowanie nimi w p�tlach.
## Program zosta� podzielony na funkcje dla czytelno�ci i lepszej organizacji kodu.
## 
## Adrian Kloc


x0 = [-2 -1 0 2 4];
y0 = [-1 0 5 99 -55];
n = length(x0);

## funckja obliczaj�ca ilorazy - zwraca macierz z ilorazami
function [y1] = ilorazy_roz(n, x0, y0)
m = 1;
y1 = [1;1];
  
## obliczanie pierwszego rz�du iloczyn�w
  for i=1:(n-1),
      y = (y0(i+1)-y0(i))/(x0(i+1)-x0(i));
      y1(i, 1) = y;
      y = 0;
    endfor;
##obliczanie reszt rz�d�w na podstawie pierwszego
  for i=2:(n-1),
    for j=1:(n-i),
      y = ((y1(j+1, m)-y1(j, m))/(x0(j+i)-x0(j)));
      y1(j, i) = y;
      y = 0;
    endfor;
    m += 1;
  endfor;
endfunction;

[y1] = ilorazy_roz(n, x0, y0);

## funkcja obliczaj�ca warto�ci dla punkt�w w przestrzeni linspace - zwraca przestrze� i jej warto�ci
function [s, yi] = inter_linspace(n, y1, x0, y0)
  s = [-10:0.1:55]; ## w�asna przestrze�
  l = length(s);    ## d�ugo�� wektora s
  yi = zeros(1,l);  ## wektor dla obliczonych warto�ci
  
  w1 = 0;   ##
  w2 = 1;   ## to s� zmienne tymczasowe, potrzebne dla oblicze� kolejnych p�tli
  w3 = 0;   ##
  w4 = 0;   ##

  for i=1:l
    for j=1:n-1
      for k=1:j
        w1 = s(i) - x0(k);
        w2 *= w1;
      endfor
      w3 += y1(1,j) * w2;
      w1 = 0;
      w2 = 1;
    endfor
    w4 = y0(1) + w3;
    yi(i) = w4;
    w3 = 0;
    w4 = 0;
  endfor
endfunction
[s, yi] = inter_linspace(n, y1, x0, y0);
## funkcja obliczaj�ca warto�� dla zadanego x - zwraca x i y
function [x, y] = inter_x(n, y1, x0, y0)
  x = 1.45; ## zadany x przez u�ytkownika
  
  w1 = 0;
  w2 = 1;
  w3 = 0;
  w4 = 0;

    for j=1:n-1
      for k=1:j
        w1 = x - x0(k);
        w2 *= w1;
      endfor
      w3 += y1(1,j) * w2;
      w1 = 0;
      w2 = 1;
    endfor
    w4 = y0(1) + w3;
    y = w4;
endfunction

[x, y] = inter_x(n, y1, x0, y0);

## rysowanie wykresu ze wszystkimi potrzebnymi danymi
plot(x0, y0, 'ro;w�z�y;', s, yi, 'b-;f(y);', x, y, 'go;x;');
  axis([-3, 5, -150, 150]);
  title('funkcja x0');
  xlabel('x');
  ylabel('y');