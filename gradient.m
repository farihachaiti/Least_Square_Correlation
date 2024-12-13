%  example call of GRADIENT-function
%  ---------------------------------------------
%  function test                    
%  %        ====
%  f = double(imread('image.tif'));                     % read greyvalue image                      
%  [fx, fy] = gradient(f, 2);                 % compute gradients in x and y
%  figure(1); imshow (fx, []);                % horizontal greyvalue changes 
%  figure(2); imshow (fy, []);                  % vertical greyvalue changes
  
function [fx, fy] = gradient(f, sigma)
%        =============================                          
r = sqrt(2)*pi*sigma;            % radius of feature mask based on standard deviation
[x, y] = meshgrid(-r:r, -r:r);                       % centering of image wirdow
                                % 1. derivative of Gauss function in x-direction
mask = -x/(2*pi*sigma.^4).*exp(-0.5*(x.^2+y.^2)/(2*sigma^2));   
fx = imfilter(f, mask , 'conv');                     % horizontal convolution and 
fy = imfilter(f, mask', 'conv');             % vertical convolution with the mask