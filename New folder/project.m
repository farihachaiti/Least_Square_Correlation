function project
  clc
  pkg load image
  
  r = 2;
  thres = 0.5; 
  dmin = 10; 
  dmax = 50; 
  %reference = double(imread('left. png')); 
  %search = double(imread('right.png'));
  %[h, w] = size(reference);
  %h = 15;
  %w = 100;

  
  #Affine Transformation
%  printf("Affine Transformation\n")
  
 % reference = double(imread('input.jpg'));
 % [h, w] = size(reference);  
  
 % H = [1 0.1 2; 0.2 0.7 3; 0 0 1];  
 % search = geotrans(H, reference);   
 % figure; imshow (search, []);   
  #Least Squares Matching
  printf("Least Squares Matching \n")
    h = 100;
    w = 100;
    D = zeros(h, w);
  reference = double(imread('target_img.png'));
  search = double(imread('distorted_img_A.png'));
  %for iterator = 1:10   
    for i = 1+r:h-r % For each row i and column j of the reference
      for j = 1+r:w-r % image in a distance r from the border
        cmax = thres;
        start = max(j-dmax,r+1);% Crop the search space
        stop = max(j-dmin,r+1);
        f = search(i-r:i+r,j-r:j+r);
        [fx,fy] = gradient(f, 2);      
        for k = start:stop
          g = reference(i-r:i+r,k-r: k+r);
          z = least_squares_correlation(fx,fy,i,j,f,g);
          H = [1+z(1,1) z(2,1) z(5,1); z(3,1) 1+z(4,1) z(6,1); 0 0 1]
          g_prime = geotrans(H,g);
          %figure; imshow (g_prime, []); 
          disparity = sum(sum([f-g_prime]))
          if disparity>thres
            min_disparity = thres
          else
            min_disparity = disparity
            D(i,j) = j-k;
          endif
        endfor
      endfor
    endfor
%endfor
figure(2); imshow(D,[]);
end

function z = least_squares_correlation(fx,fy,x,y,f,g)

  A1 = [sum(sum(fx^2*x^2)) sum(sum(fx^2*x*y)) sum(sum(fx*fy*x^2)) sum(sum(fx*fy*x*y)) sum(sum(fx^2*x)) sum(sum(fx*fy*x))];
  A2 = [sum(sum(fx^2*x*y)) sum(sum(fx^2*y^2)) sum(sum(fx*fy*x*y)) sum(sum(fx*fy*y^2)) sum(sum(fx^2*y)) sum(sum(fx*fy*y))];
  A3 = [sum(sum(fx*fy*x^2)) sum(sum(fx*fy*x*y)) sum(sum(fy^2*x^2)) sum(sum(fy^2*x*y)) sum(sum(fx*fy*x)) sum(sum(fy^2*x))];
  A4 = [sum(sum(fx*fy*x*y)) sum(sum(fx*fy*y^2)) sum(sum(fy^2*x*y)) sum(sum(fy^2*y^2)) sum(sum(fx*fy*y)) sum(sum(fy^2*y))];
  A5 = [sum(sum(fx^2*x)) sum(sum(fx^2*y)) sum(sum(fx*fy*x)) sum(sum(fx*fy*y)) sum(sum(fx^2)) sum(sum(fx*fy))];
  A6 = [sum(sum(fx*fy*x)) sum(sum(fx*fy*y)) sum(sum(fy^2*x)) sum(sum(fy^2*y)) sum(sum(fx*fy)) sum(sum(fy^2))];
  
  A = [A1;A2;A3;A4;A5;A6]
  
  b = transpose([sum(sum(f.-g*fx*x)) sum(sum(f.-g*fx*y)) sum(sum(f.-g*fy*x)) sum(sum(f.-g*fy*y)) sum(sum(f.-g*fx)) sum(sum(f.-g*fy))]);
  z = pinv(A)*b;

end
