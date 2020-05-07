clc();
clear();
close all;

A = [-2   -2   -3   0;
     -2    1   -3  -3;
     -2   -2   -4  -2;
     -1    2   -3  -2];
b = [5; 2; 6; 1];

function [x, Ab] = elim_Gauss(A,b)
%% A - Macierz wsp�czynnik�w
%% b - Wektor wyraz�w wolnych
%% x - Wektor z rozwi�zaniami uk�adu r�wna� liniowych
%%
%% Ab - macierz wsp�czynnik�w powi�kszona o 
%%      wektor wyraz�w wolnych

  Ab = [A,b]; 
  [w, k] = size(Ab);
  x = zeros(w,1);
    
  for j = 1:w-1
      for i = j+1:w
          Ab(i,j:k) = Ab(i,j:k)-Ab(i,j)/Ab(j,j)*Ab(j,j:k);
      endfor;
  endfor;

  x(w) = Ab(w,k)/Ab(w,w);
  for i = w-1:-1:1
      x(i)=(Ab(i,k)-Ab(i,i+1:w)*x(i+1:w))/Ab(i,i);
  endfor;
  
%% czasami w obliczonych x pojawia si� minus przy zerze
%% p�tla poprawia ten b��d
  for i=1:w
      if(x(i, 1) == -0)
        x(i, 1) = 0;
      endif;
  endfor;
  
endfunction;

[x, Ab] = elim_Gauss(A,b)