function zi = interp2dd(x, y, z, xi, yi, varargin)
% function zi = INTERP2DD(x, y, z, xi, yi, varargin)
% interpolates the grid z created by the x and y vector by using the
% interpolation points xi and yi. 
% Eg: 
% x = 0:0.5:10; y = 0:0.5:10; temp = 1:(length(x)*length(y)); 
% z = reshape(temp, length(x), length(y));
% xi = 0:0.1:10; yi = 0:0.1:10; zi = interp2dd(x, y, z, xi, yi, 'cubic'); 
% figure; subplot(2, 1, 1); imagesc(x, y, z); 
% subplot(2, 1, 2); imagesc(xi, yi, zi); 

    [x1, y1] = meshgrid(xi, yi); 
    zi =  interp2(x, y, z, x1, y1, varargin{:}); 
    
end
