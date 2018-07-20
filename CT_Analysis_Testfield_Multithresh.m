clear, clc, format short g

fontsize = 20;

% Make some overall designations first

Number_Slices = input('Enter the number of slices you want to analyze per time point: ');

% Define the primary folder through directory designation

A = dir; 
Primary_Folder = A.folder;

disp('  ')

% Primary loop to obtain average grayscale values

for aa = 1:(length(dir)-4)
    
    % Selecting the folder for a given time
    
    Folder = A(aa+2).name;
    disp(strcat('You are analyzing the current folder:',{' '},Folder)) 
    cd(Folder)
    
        % Selecting the hypotonic image and isotonic image you want to analyze
        
        for bb = 1:uint16(Number_Slices)
            
            % Hypotonic image
            
            filename_designation = char(strcat('Enter hypotonic image',{' '},num2str(bb),{' '},'of',{' '},num2str(Number_Slices),':',{' '}));
            filename = input(filename_designation,'s');
            
            disp('  ')
            
            I = uint8(imread(filename));
            I_size = size(I);
            I_rows = uint16(I_size(1));
            I_columns = uint16(I_size(2));
           figure
            imshow(I) 
            title('Original')
            hold on 
            
            
            
            % Boost the image contrast to aide in boundry recognition
            
            I2 = adapthisteq(I);
            
            % Contrasted image
            
            IC = uint8(I2);
            IC_size = size (IC);
            IC_rows = uint16(IC_size(1));
            IC_columns = uint16(IC_size(2));

            % Convert contrasted image to binary, 
            
            BW = im2bw(IC);
            
            % Sharpen image borders by clearing out unneccesary noise
            
            BW_2 = imclearborder(BW);
            
            % Remove any objects of size < 600 pixels
            
            BW_3 = bwareaopen(BW_2,600);
            
            % Display the cleaned figure
            
            figure
            imshow(BW_3)
            title('Cleaned Figure')
            % Display the contrasted image to plot boundaries
            
            imshow(IC);
            hold on;
            title('Contrasted Image')
                        
            % Find the boundaries of the binary image 
            
            boundaries = bwboundaries(BW_3);	
            
            % Determine how many seperate boundaries are within the image
            
            numberOfBoundaries = size(boundaries);
            
            % Plot the boundaries upon the image based upon their reference
            % value i; in this code, the first boundary will be given to
            % the hypotonic kidney, while the second will be given to the
            % isotonic kidney
            
            for i = 1 : numberOfBoundaries
                    thisBoundary = boundaries{i};
                    plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 1);
            end
            hold off;
            
            % The boundary plotted in the previous for loop will now be
            % used to isolate both kidneys within the original contrasted
            % image
            
            [rows columns numberofcolorbands] = size(IC);
            
            maskedImage = zeros(rows, columns, 'uint8');
            maskedImage(BW_3) = IC(BW_3); 
            
            % Display the new masked image 
            
            imshow(maskedImage);
            title('New Masked Image')
            
            % Most 'holes' presented by the tubing will now be filled with
            % an averaged grayscale value
            
            maskedImage_2 = imfill(maskedImage);
            
            % This filled image is now converted to binary and displayed
            
            BW_4 = im2bw(maskedImage_2); %switch im2bw to imbinarize
            imshow(maskedImage_2)         
            %hold on;
            title('Green Boundary of Original Photo')
            
            % New boundaries are now calculated and plotted for the filled
            % image
            
%             boundaries_2 = bwboundaries(maskedImage_2);	
%             numberOfBoundaries_2 = size(boundaries_2);
%             for j = 1 : numberOfBoundaries_2
%                     thisBoundary_2 = boundaries_2{j};
%                     plot(thisBoundary_2(:,2), thisBoundary_2(:,1), 'g', 'LineWidth', 1);
%             end
%             imfill(maskedImage_2);
%             hold off;
            %note: threshold -> 165                              %Solomon's
                                                            %original code 
%             Cortex = maskedImage_2 > 165;
%             Medulla = (maskedImage_2 < 165)&(maskedImage_2 > 40);
%             figure
%             imshow(Cortex)
%             figure
%             imshow(Medulla)

%             the following 4 lines: testing adaptthresh, creates a binary image 
%             T = adaptthresh(maskedImage, 0.4);
%             BW_5 = imbinarize(maskedImage,T); 
%             figure
%             imshowpair(maskedImage, BW_5, 'montage')

%the following 6 lines: testing multithresh - for some reason, didn't
%work... 
% Z = imread('maskedImage');
% imshow(Z)
% level = multithresh(Z); 
% seg_I = imquantize(Z, level); 
% figure 
% imshow(seg_I,[])

%the following . lines: testing segmenting image into three levels with 2
%thresholds 
% Z = imread('maskedImage_2'); 
% figure
% imshow(Z)
% axis off
% title('Original Image')

% The following 6 lines: takes maskedImage_2 and makes the cortex white,
% medulla grey, and everything else black 
thresh = multithresh(maskedImage_2, 2); 
seg_I = imquantize(maskedImage_2, thresh); 
figure
imshow(seg_I, []) 
imshowpair(seg_I, maskedImage_2, 'montage') %places the seg_I image next to the maskedImage_2
title(['Multithresh Image for ', num2str(filename)]) %gives the resulting image a title with the corresponding photo number


P_value = impixel(I)
%c = [
%r = [

P = impixel(I,seg_I,c,r)

%The following 5 lines: takes seg_I figure and makes it colored
% RGB = label2rgb(seg_I);
% figure; 
% imshow(RGB)
% axis off
% title('RGB Segmented Image')


           
            % The boundaries calculated in the previous step are now
            % referenced as {1} and {2}, which corresponds to the hypotonic
            % boundary and the isotonic boundary 
            
            [B,~] = bwboundaries(maskedImage_2,'noholes');
            [r c] = size(maskedImage_2);
            [B numblobs] = bwboundaries(maskedImage_2);
            
            % The coordinates of the hypotonic kidney are now calculated,
            % and plotted in a binary polymask 
            
            firstBoundary = B{1};
            x = firstBoundary(:,2);
            y = firstBoundary(:,1);
            firstblob = poly2mask(x,y,rows,columns);
            
            % The coordinates of the isotonic kidney are then calculated
            % and plotted in the same fashion
            
            secondBoundary = B{2};
            x_2 = secondBoundary(:,2);
            y_2 = secondBoundary(:,1);
            secondblob = poly2mask(x_2,y_2,rows,columns);
            
            % Preallocate both regions
            
            Hypotonic_region_pixels = zeros(I_rows,I_columns,'uint16');
            Isotonic_region_pixels = zeros(I_rows,I_columns,'uint16');
            
            % Set the if statement count variable for both regions to zero 
            
            Loop_count_HR = uint16(0);
            Loop_count_IR = uint16(0);
            
            for dd = 1:I_rows
                for ee = 1:I_columns
                    
                    % Both the hypotonic and isotonic regions are set equal
                    % to their corresponding polymasks, with reference to
                    % the for loop indices; this ensures that ever value
                    % not within said polymask is registered as a zero
                    
                    Hypotonic_region_pixels(dd,ee) = firstblob(dd,ee);
                    Isotonic_region_pixels(dd,ee) = secondblob(dd,ee);
                    
                    % These two if statements will then scan the image for
                    % any values within the polymask, assign it a logical
                    % "yes", and store it within the loop count
                    % (specifically assigned to the original image)
                    
                    if Hypotonic_region_pixels(dd,ee) == 1
                        Loop_count_HR = Loop_count_HR+1;
                        I_hold_HR(Loop_count_HR) = I(dd,ee);
                    end
                    if Isotonic_region_pixels(dd,ee) == 1
                        Loop_count_IR = Loop_count_IR+1;
                        I_hold_IR(Loop_count_IR) = I(dd,ee);
                    end
                end
            end
            
            % The means of both loop counts are then averaged, which takes
            % the pixels registered as "yes" within the original image, and
            % averages their grayscale values
            
            Hypotonic_region_I(aa,bb) = mean(I_hold_HR);
            Isotonic_region_I(aa,bb) = mean(I_hold_IR)
            
            % NOTE: Because the imclearborder and bwareaopen functions were
            % used to clean the image, a large portion of the renal pelvis
            % was erased, therefore leading to smaller average grayscale
            % values
            
        end
        cd(Primary_Folder)
end
  