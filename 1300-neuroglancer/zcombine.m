function n = zcombine(a, b, x, y, z, root)
% ZCOMBINE - Combines several z tiles
%   ZCOMBINE(a, b, x, y, z, b) loads 2^b tiles to construct the
%   averaged tile a/b/z/y/x from a/z0..z1/y/x.
%   ZCOMBINE is smart enough to divide by less than 2^b if some source
%   tiles are missing and doesn't save anything when there are zero source
%   tiles.
%   n = ZCOMBINE(...) returns the number of tiles.

if nargin<6
  root = '/lsi2/dw/170428/q1pyramid';
end

Z = 2^b;
n = 0;
for dz=1:Z
  z1 = z*Z + dz-1;
  zlo = mod(z1, 100);
  zhi = floor(z1/100);
  fn = sprintf('%s/Z%i/%i/A%i/Y%i/X%i.jpg', ...
               root, zhi, zlo, a, y, x);
  if exist(fn)
    img = imread(fn);
    if n==0
      ima = zeros(size(img));
    end
    ima = ima + double(img);
    n = n + 1;
  end
end

if n>0
  ima = ima / n;
  ima = uint8(ima);
  zlo = mod(z, 100);
  zhi = floor(z/100);
  pth = sprintf('%s/A%iB%i/Z%i/%i/Y%i', root, a, b, zhi, zlo, y);
  mkdir(pth);
  printf('Writing %s\n', sprintf('%s/X%i.jpg', pth, x));
  imwrite(ima, sprintf('%s/X%i.jpg', pth, x));
end
