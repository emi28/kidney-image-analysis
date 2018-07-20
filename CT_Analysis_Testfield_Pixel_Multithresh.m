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
            
            % Display the contrasted image to plot boundaries
            
            imshow(IC);
            hold on;
            
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
            
            % Most 'holes' presented by the tubing will now be filled with
            % an averaged grayscale value
            
            maskedImage_2 = imfill(maskedImage);
            
            % This filled image is now converted to binary and displayed
            
            BW_4 = im2bw(maskedImage_2);
            imshow(maskedImage_2)         
            hold on;
            
            % New boundaries are now calculated and plotted for the filled
            % image
            
            boundaries_2 = bwboundaries(maskedImage_2);	
            numberOfBoundaries_2 = size(boundaries_2);
            for j = 1 : numberOfBoundaries_2
                    thisBoundary_2 = boundaries_2{j};
                    plot(thisBoundary_2(:,2), thisBoundary_2(:,1), 'g', 'LineWidth', 1);
            end
            imfill(maskedImage_2);
            hold off;
            
% Solomon's Original Code for differentiating Cortex from Medulla
%             Cortex = maskedImage_2 > 165;
%             Medulla = (maskedImage_2 < 165)&(maskedImage_2 > 40);
%             figure
%             imshow(Cortex)
%             figure
%             imshow(Medulla)



              %The following code creates two thresholds and shows an image
              %which I named, seg_I, of maskedImage_2 where the cortex is white, the medulla is
              %grey, and everything else is black. 
              thresh = multithresh(maskedImage_2, 2); %This line finds the thresholds using Otsu's method
              seg_I = imquantize(maskedImage_2, thresh); %This line creates the new image with white, black and gray
              figure
              imshow(seg_I, [])
              figure
              imhist(seg_I)
              hold on
             % P_value = impixel(seg_I)
              
              imshowpair(seg_I, maskedImage_2, 'montage') %this line shows the new image next to maskedImage_2
            
              title(['Multithresh Image for ', num2str(filename)]) %titles the image with the corresponding photo name
              
%              for 1:length(filename)
%                 if pixelvalue == 1
%                     pixelvalue = Cortex
%                 else if pixelvalue == 0.5
%                         pixelvalue = Medulla
%                     else 
%                         pixelvalue = Outside
            

              %The following code formats the seg_I image in RGB 
               RGB = label2rgb(seg_I);
               figure; 
%                 cmap = colormap(seg_I)          %my attempts to get the
%                                                   colored values 
%                 cmap_2 = colormap(RGB)
%                 mymap =                         %Here, I tried to change
%                                                   the three colors 
%                     [1 0 0 
%                     0 1 0 
%                     0 0 1]
%                 colormap(mymap)
                imshow(RGB)
                axis off
                title('RGB Segmented Image')
               
                
%                 low = 
%                 high = 
%                 Binary_Medulla_Area = roicolor(RGB, low, high)   %Tried
%                                                                   to find the color value of the regions (gray) 
%    
              
              
              %P_value = impixel(seg_I)                            %another
                                                                   %attempt to find the color values of the regions in the seg_I
                                                                   %image
            
              
              %c = [
              %r = [

              % P = impixel(I,seg_I,c,r)
              

            
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
            Isotonic_region_I(aa,bb) = mean(I_hold_IR);
            
            % NOTE: Because the imclearborder and bwareaopen functions were
            % used to clean the image, a large portion of the renal pelvis
            % was erased, therefore leading to smaller average grayscale
            % values
            
        end
        cd(Primary_Folder)
end
  