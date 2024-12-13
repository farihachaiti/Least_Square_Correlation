%  Example call of the GEOTRANS-function
%  --------------------------------------
%  function test                    
%  %        ====
%  f = double(imread('image.tif'));                     % read grayscale image                      
%  H = [1 0.1 2; 0.2 0.7 3; 0 0 1];    % affine example transformation
%  g = geotrans (H, f);                   % transformation of iamge f using H
%  figure; imshow (g, []);                                 % show result

function g = geotrans(H, f)
%        ==================
[h w] = size(f);       
[Xp, Yp] = ndgrid(1 : w, 1 : h);
X = norm2(H \ [Xp(:) Yp(:) ones(w*h,1)]');          % indirect transformation
clear g;        % Resampling of intensity values using bilinear interpolation
xI = reshape(X(1,:), w, h)';
yI = reshape(X(2,:), w, h)';
g(:,:,1) = interp2(f(:,:,1), xI, yI, 'linear');                                        
g2 = ~isfinite(g);
g(g2) = 0;

function n = norm2(x)
%        ============
for i = 1:3
    n(i,:) = x(i,:) ./ x(3,:);
end