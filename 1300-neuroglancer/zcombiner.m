function k = zcombiner(a, b, z, y, root)
% ZCOMBINER - Combines several z tiles, all over the place
%   ZCOMBINER(a, b, z, y, root) constructs all the combined tiles
%   at a/b/z/y using ZCOMBINE.
%   ZCOMBINER(a, b, z, [], root) finds all the y columns at a/z and
%   calls itself recursively. z may be cell containing {zhi}, in which
%   case ZCOMBINER iterates over all the pertaining zlo values.
%   ZCOMBINER(a, b, [], [], root) finds all the z columns at a and
%   calls itself recursively.
%   ZCOMBINER leaves breadcrumbs to signify that it has completed an 
%   entire level and returns the number of output images that it has
%   made.

if nargin<5 || isempty(root)
  root = '/lsi2/dw/170428/q1pyramid';
end

Zmax = 9604; % max Z at B=0
Zmax = ceil(Zmax / 2^b); % max Z at given B
Ymax = ceil(125000 / 512); % max Y tiles at A=0
Ymax = ceil(Ymax / 2^a); % max Y at given A
Xmax = ceil(45000 / 512); % max X tiles at A=0
Xmax = ceil(Xmax / 2^a); % max X at given A

k = 0;
if nargin<3 || isempty(z)
  % got a/b; iterate zmaj
  guard = sprintf('%s/A%iB%i/complete', root, a, b);
  if ~exist(guard)
    Zmaxmaj = ceil(Zmax/100);
    fprintf(2, 'Working on A%i B%i [Z=0..%i99]\n', a, b, Zmaxmaj);
    for zmaj = 0:Zmaxmaj
      k = k + zcombiner(a, b, {zmaj}, [], root);
    end
    if k>0
      unix(sprintf('touch %s', guard));
    end
  end
elseif iscell(z)
  % got a/b/zmaj; iterate zmin
  zmaj = z{1};
  guard = sprintf('%s/A%iB%i/Z%i/complete', root, a, b, zmaj);
  if ~exist(guard)
    fprintf(2, 'Working on A%i B%i Z%ixx\n', a, b, zmaj);
    for zmin = 0:99
      k = k + zcombiner(a, b, zmaj*100 + zmin, [], root);
    end
    if k>0
      unix(sprintf('touch %s', guard));
    end
  end
elseif nargin<4 || isempty(y)
  % got a/b/zmaj/zmin; iterate y
  zmin = mod(z, 100);
  zmaj = floor(z/100);
  guard = sprintf('%s/A%iB%i/Z%i/%i/complete', root, a, b, zmaj, zmin);
  if ~exist(guard)
    fprintf(2, 'Working on A%i B%i Z%i [Y=0:%i]\n', a, b, z, Ymax);
    for y = 0:Ymax
      k = k + zcombiner(a, b, zmaj*100 + zmin, y, root);
    end
    if k>0
      unix(sprintf('touch %s', guard));
    end
  end
else
  % got a/b/zmaj/zmin/y; iterate x
  zmin = mod(z, 100);
  zmaj = floor(z/100);
  guard = sprintf('%s/A%iB%i/Z%i/%i/Y%i/complete', root, a, b, ...
                  zmaj, zmin, y);
  if ~exist(guard)
    fprintf(2, 'Working on A%i B%i Z%i Y%i\n', a, b, z, y);
    for x = 0:Xmax
      if ~exist(sprintf('%s/A%iB%i/Z%i/%i/Y%i/X%i.jpg', root, a, b, ...
                  zmaj, zmin, y, x))
        if zcombine(a, b, x, y, z, root)
          k = k + 1;
          fprintf(2, 'Constructed A%i B%i Z%i Y%i X%i\n', a, b, z, y, x);
        else
          fprintf(2, 'Nothing for A%i B%i Z%i Y%i X%i\r', a, b, z, y, x);
        end
      end
    end
    if k>0
      unix(sprintf('touch %s', guard));
    end
  end
end
return

