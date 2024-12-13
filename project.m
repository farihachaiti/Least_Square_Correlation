function project
  clc
  pkg load image
   
  #Least Squares Matching
  printf("Least Squares Matching \n")
  
  H = [1 0.1 2; 0.2 0.7 3; 0 0 1];      
  f = double(imread('input.jpg'));
  g = geotrans(H,f); 
   
 % figure; imshow(g,[]);

  
  %f = double(imread('target_img.png'));
  %g = double(imread('distorted_img_A.png'));
  %g = double(imread('distorted_img_B.png'));
  %g = double(imread('distorted_img_C.png'));
  
  [w,h] = size(f);
  x = [-round(w/2):round(w/2)];
  y = [-round(h/2):round(h/2)];

  X = [];
  for i = 1:w
      X = [X; x(i) * ones(1, 100)];
  end
  Y = [];
  for j = 1:h
      Y = [Y; y(j) * ones(1, 100)];
  end

  g_prime = g;
 for iterator = 1:200
   [fx,fy] = gradient(g_prime, 2);  
    z = least_squares_correlation(fx,fy,X,Y,g_prime,f);
    save('-append','parameters.mat','z');
    H = [1+z(1) z(2) z(5); z(3) 1+z(4) z(6); 0 0 1];
    g_prime = geotrans(H,g_prime);
  endfor
  g_prime = geotrans(H,g);
  figure; imshow(g_prime,[]);
end

function z = least_squares_correlation(fx,fy,x,y,f,g)

  A1 = [sum(sum((fx.^2)*(x.^2))) sum(sum(((fx.^2)*x).*y))  sum(sum((fx*fy)*(x.^2))) sum(sum(((fx*fy)*x).*y)) sum(sum((fx.^2)*x)) sum(sum((fx*fy)*x))];
  A2 = [sum(sum(((fx^2)*x).*y)) sum(sum((fx^2)*(y.^2))) sum(sum(((fx*fy)*x).*y)) sum(sum((fx*fy)*y.^2)) sum(sum((fx^2)*y)) sum(sum((fx*fy)*y))];
  A3 = [sum(sum((fx*fy)*(x.^2))) sum(sum(((fx*fy)*x).*y)) sum(sum((fy^2)*(x.^2))) sum(sum(((fy^2)*x).*y)) sum(sum((fx*fy)*x)) sum(sum((fy^2)*x))];
  A4 = [sum(sum(((fx*fy)*x).*y)) sum(sum((fx*fy)*(y.^2))) sum(sum(((fy^2)*x).*y)) sum(sum((fy^2)*(y.^2))) sum(sum((fx*fy)*y)) sum(sum((fy^2)*y))];
  A5 = [sum(sum((fx^2)*x)) sum(sum((fx^2)*y)) sum(sum((fx*fy)*x)) sum(sum((fx*fy)*y)) sum(sum(fx^2)) sum(sum(fx*fy))];
  A6 = [sum(sum((fx*fy)*x)) sum(sum((fx*fy)*y)) sum(sum((fy^2)*x)) sum(sum((fy^2)*y)) sum(sum(fx*fy)) sum(sum(fy^2))];
  
  
  A = [A1;A2;A3;A4;A5;A6];
  
  b = [sum(sum((f-g).*sum(sum(fx)).*sum(sum(x)))) sum(sum((f-g).*sum(sum(fx)).*sum(sum(y)))) sum(sum((f-g).*sum(sum(fy)).*sum(sum(x)))) sum(sum((f-g).*sum(sum(fy)).*sum(sum(y)))) sum(sum((f-g).*sum(fx))) sum(sum((f-g).*sum(fy)))];
  z = pinv(A)*transpose(b);

end