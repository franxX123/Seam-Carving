
image = imread('York.jpg');
n_img1 = mySeamCarving(image,1200);
imshow(n_img1, [])

% transposed_nimg1 = permute(n_img1, [2 1 3]);
% n_img2 = mySeamCarving(transposed_nimg1,720);
% 
% final_img = permute(n_img2, [2 1 3]);
% imshow(final_img, [])
% size(final_img)

function new_image = mySeamCarving(init_img, new_col)
    
    img = init_img;
    [init_row, init_col, init_channel] = size(init_img);
    num_iterations = init_col - new_col;
    
    for iteration = 1:num_iterations
        % Part a
        filter = fspecial('sobel');

        im_dy_red = imfilter(double(img(:,:,1)), filter, 'conv');
        im_dx_red = imfilter(double(img(:,:,1)), filter', 'conv');

        im_dy_green = imfilter(double(img(:,:,2)), filter, 'conv');
        im_dx_green = imfilter(double(img(:,:,2)), filter', 'conv');

        im_dy_blue = imfilter(double(img(:,:,3)), filter, 'conv');
        im_dx_blue = imfilter(double(img(:,:,3)), filter', 'conv');

        im_gradient_red = sqrt(im_dy_red.^2 + im_dx_red.^2);
        im_gradient_green = sqrt(im_dy_green.^2 + im_dx_green.^2);
        im_gradient_blue = sqrt(im_dy_blue.^2 + im_dx_blue.^2);

        E = im_gradient_red + im_gradient_green + im_gradient_blue;
        % 	imshow(E, []);

        % Part b
        [img_row, img_col, channels] = size(img);
        M = zeros(img_row, img_col);

        % Part c
        M(1,:) = E(1,:);

        % Part d
        M = createScoring(E, M);

        % Test 1
        % E = [0,5,1;6,8,4;9,1,7];
        % M = zeros(size(E));
        % 
        % M(1,:) = E(1,:);
        % M = createScoring(E, M)

        % Part e
        % Get the lowest values at the bottom row of M

        bottom_row = img_row;

        min = 1;
        last_row = M(bottom_row,:);

        for i = 1:img_col
            if last_row(i) < last_row(min)
                min = i;
            end
        end

        bottom_col = min;
        M(bottom_row, bottom_col); % The value at the bottom

        % Part f
        % Extract the seam from M
        seam = BFSGetSeam(bottom_row, bottom_col, M);

        % Part g
        % Remove the seam from the image

        new_img = uint8(zeros(img_row, img_col-1, channels));
        for c = 1:channels
            for i = 1:img_row 
                delta = 0;
                for j = 1:img_col  

                   if j == seam(i) 
                        delta = -1;
                   else
                       new_img(i,j + delta,c) = uint8(img(i,j,c));
                   end

                end
            end
        end
        
        img = new_img;
        % size(new_img);
        % imshow(new_img, []);
    end
    
    new_image = img;
end

% Part e Helper
% This returns a 1D array where the indexes represents the row number
% the values in the indexes represent which column along the row is the
% smallest.
function seam = BFSGetSeam(bot_x, bot_y, M)
    [row, col] = size(M);
    seam_path = zeros(1, row);
    seam_path(row,1) = bot_y;
    x = bot_x;
    y = bot_y;
    
    while x > 1
        m = 0; % m for minimum
        x = x - 1;
        if (y-1) >= 1 && M(x, y-1) < M(x,y)
            m = y - 1;
        else    
            m = y;
        end
        
        if (y+1) <= col && M(x,m) > M(x, y + 1)
            m = y + 1;
        end
        
        seam_path(x) = m;
        y = m; % We change the place along the columns as we go up the rows
    end
    
    seam = seam_path;
end


% Part d Helper
function result = createScoring(energyM,optM)
    [energy_row, energy_col] = size(energyM);
    for x = 2:energy_row
        for y = 1:energy_col
            
            top_left = Inf;
            top_mid =  optM(x-1, y);
            top_right = Inf;
            
            if y-1 > 0
                top_left = optM(x-1, y-1);
            end         
           
            if y+1 <= energy_col
                top_right = optM(x-1, y+1);
            end
                
            optM(x,y) = energyM(x,y) + min([top_left, top_mid, top_right]);
        end
    end
    
    result = optM;
end






