clc();
clear();
close all;

x = [0.000  0.020 0.036 0.060 0.094 0.133 0.164 0.196 0.234 0.264 0.285 0.300]; ## napiêcie
y = [0.0    4.7   7.7   10.5  11.5  10.0  7.0   6.0   8.0   12.0  16.0  19.0];  ## natê¿enie
yh = [0, y, 0];

## jako ¿e wbudowana funkcja csape nie dzia³a³a na moim komputerze
## postanowi³em wkleiæ ca³¹ funkcjê do programu.
function pp = csape (x, y, cond, valc); 

  x = x(:);
  n = length(x);
  if (n < 3) 
    error ("csape requires at least 3 points");
  endif

  ## Check the size and shape of y
  ndy = ndims (y);
  szy = size (y);
  if (ndy == 2 && (szy(1) == n || szy(2) == n))
    if (szy(2) == n)
      a = y.';
    else
      a = y;
      szy = fliplr (szy);
    endif
  else
    a = shiftdim (reshape (y, [prod(szy(1:end-1)), szy(end)]), 1);
  endif
  m = size (a, 2);

  if exist("valc", "var") && size(valc) == [1 2]
    valc = valc';
  endif

  b = c = zeros (n, m);
  h = diff (x);
  idx = ones (columns(a),1);

  if (nargin < 3 || strcmp(cond,"complete"))

    # set first derivative at end points
    if (nargin < 4)
      valc = zeros(2, m);
      n_use = min(n, 4);
      for i = 1:m
        valc(1, i) = polyval(polyder(polyfit(x(1:n_use), a(1:n_use, i), n_use-1)), x(1));
        valc(2, i) = polyval(polyder(polyfit(x(end-n_use+1:end), a(end-n_use+1:end, i), n_use-1)), x(end));
      endfor
    endif

    if (n == 3)
      warning ("off", "Octave:broadcast", "local");
      e = 2 * [h(1) h(1:(end-1))+h(2:end) h(end)];
      A = spdiags([[h; 0], e(:), [0; h]], [-1,0,1], n, n);
      d = diff(a) ./ h; #uses broadcasting if columns(a) > 1
      g = 3 * diff(d);
      A(1, 1) = 2 * h(1);
      A(1, 2) = h(1);
      A(n, n) = 2 * h(end);
      A(end, end-1) = h(end);
      g = [3*(d(1, :) - valc(1, :)); g; 3*(valc(2, :) - d(end, :))];
      c = A \ g;
    else
      dg = 2 * (h(1:n - 2) .+ h(2:n - 1));
      dg(1) = dg(1) - 0.5 * h(1);
      dg(n - 2) = dg(n-2) - 0.5 * h(n - 1);

      e = h(2:n - 2);

      g = 3 * diff (a(2:n,:)) ./ h(2:n - 1,idx)...
        - 3 * diff (a(1:n - 1,:)) ./ h(1:n - 2,idx);
      g(1,:) = 3 * (a(3,:) - a(2,:)) / h(2) ...
          - 3 / 2 * (3 * (a(2,:) - a(1,:)) / h(1) - valc(1,:));
      g(n - 2,:) = 3 / 2 * (3 * (a(n,:) - a(n - 1,:)) / h(n - 1) - valc(2,:))...
          - 3 * (a(n - 1,:) - a(n - 2,:)) / h(n - 2);

      c(2:n - 1,:) = spdiags([[e(:);0],dg,[0;e(:)]],[-1,0,1],n-2,n-2) \ g;
      
      c(1,:) = (3 / h(1) * (a(2,:) - a(1,:)) - 3 * valc(1, :) - c(2,:) * h(1)) / (2 * h(1));
      c(n,:) = - (3 / h(n - 1) * (a(n,:) - a(n - 1,:)) - 3 * valc(2, :) + c(n - 1,:) * h(n - 1)) / (2 * h(n - 1));
    end
    b(1:n - 1,:) = diff (a) ./ h(1:n - 1, idx)...
      - h(1:n - 1,idx) / 3 .* (c(2:n,:) + 2 * c(1:n - 1,:));
    d = diff (c) ./ (3 * h(1:n - 1, idx));

  elseif (strcmp(cond,"variational") || strcmp(cond,"second"))

    if ((nargin < 4) || strcmp(cond,"variational"))
      ## set second derivatives at end points to zero
      valc = zeros (2, 1);
    endif

    c(1,:) = valc(1, :) / 2;
    c(n,:) = valc(2, :) / 2;

    g = 3 * diff (a(2:n,:)) ./ h(2:n - 1, idx)...
      - 3 * diff (a(1:n - 1,:)) ./ h(1:n - 2, idx);

    g(1,:) = g(1,:) - h(1) * c(1,:);
    g(n - 2,:) = g(n-2,:) - h(n - 1) * c(n,:);

    if( n == 3)
      dg = 2 * h(1);
      c(2:n - 1,:) = g / dg;
    else
      dg = 2 * (h(1:n - 2) .+ h(2:n - 1));
      e = h(2:n - 2);
      c(2:n - 1,:) = spdiags([[e(:);0],dg,[0;e(:)]],[-1,0,1],n-2,n-2) \ g;
    end
        
    b(1:n - 1,:) = diff (a) ./ h(1:n - 1,idx)...
      - h(1:n - 1,idx) / 3 .* (c(2:n,:) + 2 * c(1:n - 1,:));
    d = diff (c) ./ (3 * h(1:n - 1, idx));
  
  elseif (strcmp(cond,"periodic"))

    h = [h; h(1)];

    ## XXX FIXME XXX --- the following gives a smoother periodic transition:
    ##    a(n,:) = a(1,:) = ( a(n,:) + a(1,:) ) / 2;
    a(n,:) = a(1,:);

    tmp = diff (shift ([a; a(2,:)], -1));
    g = 3 * tmp(1:n - 1,:) ./ h(2:n,idx)...
      - 3 * diff (a) ./ h(1:n - 1,idx);

    if (n > 3)
      dg = 2 * (h(1:n - 1) .+ h(2:n));
      e = h(2:n - 1);

      ## Use Sherman-Morrison formula to extend the solution
      ## to the cyclic system. See Numerical Recipes in C, pp 73-75
      gamma = - dg(1);
      dg(1) -=  gamma;
      dg(end) -= h(1) * h(1) / gamma; 
      z = spdiags([[e(:);0],dg,[0;e(:)]],[-1,0,1],n-1,n-1) \ ...
          [[gamma; zeros(n-3,1); h(1)],g];
      fact = (z(1,2:end) + h(1) * z(end,2:end) / gamma) / ...
          (1.0 + z(1,1) + h(1) * z(end,1) / gamma);

      c(2:n,:) = z(:,2:end) - z(:,1) * fact;
    endif

    c(1,:) = c(n,:);
    b = diff (a) ./ h(1:n - 1,idx)...
      - h(1:n - 1,idx) / 3 .* (c(2:n,:) + 2 * c(1:n - 1,:));
    b(n,:) = b(1,:);
    d = diff (c) ./ (3 * h(1:n - 1, idx));
    d(n,:) = d(1,:);

  elseif (strcmp(cond,"not-a-knot"))

    g = zeros(n - 2,columns(a));
    g(1,:) = 3 / (h(1) + h(2)) * (a(3,:) - a(2,:)...
          - h(2) / h(1) * (a(2,:) - a(1,:)));
    g(n - 2,:) = 3 / (h(n - 1) + h(n - 2)) *...
        (h(n - 2) / h(n - 1) * (a(n,:) - a(n - 1,:)) -...
         (a(n - 1,:) - a(n - 2,:)));

    if (n > 4)

      g(2:n - 3,:) = 3 * diff (a(3:n - 1,:)) ./ h(3:n - 2,idx)...
        - 3 * diff (a(2:n - 2,:)) ./ h(2:n - 3,idx);

      dg = 2 * (h(1:n - 2) .+ h(2:n - 1));
      dg(1) = dg(1) - h(1);
      dg(n - 2) = dg(n-2) - h(n - 1);

      ldg = udg = h(2:n - 2);
      udg(1) = udg(1) - h(1);
      ldg(n - 3) = ldg(n-3) - h(n - 1);
      c(2:n - 1,:) = spdiags([[ldg(:);0],dg,[0;udg(:)]],[-1,0,1],n-2,n-2) \ g;

    elseif (n == 4)

      dg = [h(1) + 2 * h(2); 2 * h(2) + h(3)];
      ldg = h(2) - h(3);
      udg = h(2) - h(1);
      c(2:n - 1,:) = spdiags([[ldg(:);0],dg,[0;udg(:)]],[-1,0,1],n-2,n-2) \ g;
      
    else # n == 3
            
      dg= [h(1) + 2 * h(2)];
      c(2:n - 1,:) = g/dg(1);

    endif

    c(1,:) = c(2,:) + h(1) / h(2) * (c(2,:) - c(3,:));
    c(n,:) = c(n - 1,:) + h(n - 1) / h(n - 2) * (c(n - 1,:) - c(n - 2,:));
    b = diff (a) ./ h(1:n - 1, idx)...
      - h(1:n - 1, idx) / 3 .* (c(2:n,:) + 2 * c(1:n - 1,:));
    d = diff (c) ./ (3 * h(1:n - 1, idx));

  else
    msg = sprintf("unknown end condition: %s",cond);
    error (msg);
  endif

  d = d(1:n-1,:); c=c(1:n-1,:); b=b(1:n-1,:); a=a(1:n-1,:);
  pp = mkpp (x, cat (2, d'(:), c'(:), b'(:), a'(:)), szy(1:end-1));

endfunction;




p = spline(x, y);  ## splajn not-a-knot
ph = spline(x, yh); ## splajn hermitowski
pn = csape(x, y, 'variational'); ## splajn naturalny

l = linspace(-0.4, 0.6, 200);

plot( l,ppval(ph,l),"g-;spl. hermitowski;",...
      l,ppval(pn,l),"c-;spl. naturalny;",...
      l, ppval(p, l), "b-;spl. not-a-knot;", x, y, "r+;wêz³y;");
axis([-0.5, 0.8, -30, 30]);
  title('Interpolacja charakterystyki pr¹dowo-napiêciowej diody tunelowej');
  xlabel('U [V]');
  ylabel('I [mA]');