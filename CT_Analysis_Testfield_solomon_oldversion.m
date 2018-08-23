clear, clc, format short g

fontsize = 20;

% Make some overall designations first

Number_Slices = input('Enter the number of slices you want to analyze per time point: ');

% Define the primary folder through directory designation

A = dir; 
Primary_Folder = A.folder;

disp('  ')

% Primary loop to obtain average grayscale values

for aa = 1:(length(dir)-14)
    
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
            
            [rows columns numberofcolorbands] = size(I);
            
            maskedImage = zeros(rows, columns, 'uint8');
            maskedImage(BW_3) = I(BW_3); 
            
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
            title('maskedImage_2')
            hold on;

            
            [B,~] = bwboundaries(maskedImage_2,'noholes');
            [r e] = size(maskedImage_2);
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
            
            maskedImage_4 = zeros(rows, columns, 'uint8');
            maskedImage_4(firstblob) = maskedImage_2(firstblob);
            
            simul_hypo_1 = im2bw(maskedImage_4);
            simul_hypo = bwareaopen(simul_hypo_1,1000);
            measurements = regionprops(simul_hypo,'Area');
            
            boundaries_SH = bwboundaries(simul_hypo);
            numberOfBoundaries_SH = size(boundaries_SH, 1);
            figure
            imshow(I)
            title('boundaries SH')
            hold on;
% NOTE: I did not write the code from this point on, only adapted it to my
% variables
            for blobIndex = 1 : numberOfBoundaries_SH
                thisBoundary_SH = boundaries_SH{blobIndex};
                x_SH = thisBoundary_SH(:, 2); % x = columns. <-- actually gives all rows of values in column 2
                y_SH = thisBoundary_SH(:, 1); % y = rows. <-- actually gives all rows of values in column 3
                plot(x_SH,y_SH, 'y');
                hold on
                
            % Find which two boundary points are farthest from each other.
                maxDistance_SH = -inf;
                for k_SH = 1 : length(x_SH)
                    distances_SH = sqrt( (x_SH(k_SH) - x_SH) .^ 2 + (y_SH(k_SH) - y_SH) .^ 2 );
                    [thisMaxDistance, indexOfMaxDistance] = max(distances_SH);
                    if thisMaxDistance > maxDistance_SH
                        maxDistance_SH = thisMaxDistance;
                        index1 = k_SH;
                        index2 = indexOfMaxDistance;
                    end
                end
            end

            xMidPoint = mean([x_SH(index1), x_SH(index2)]);
            yMidPoint = mean([y_SH(index1), y_SH(index2)]);
            longSlope = (y_SH(index1) - y_SH(index2)) / (x_SH(index1) - x_SH(index2))
                   
            perpendicularSlope = -1/longSlope
            % Use point slope formula (y-ym) = slope * (x - xm) to get points
            y1 = perpendicularSlope * (1 - xMidPoint) + yMidPoint;
            y2 = perpendicularSlope * (columns - xMidPoint) + yMidPoint;
            %theta = atand((longSlope-perpendicularSlope)/(1-longSlope*perpendicularSlope))
            %theta finds the angle between two line segments 
            
            %gets angle of perpendicular green line
            h = imdistline(gca, [1,columns], [y1,y2]); 
            angle = getAngleFromHorizontal(h) 
                     
	
            % Get the profile perpendicular to the midpoint so we can find out when if first enters and last leaves the object.
            [cx,cy,c] = improfile(simul_hypo,[1, columns], [y1, y2], 1000);
            % Get rid of NAN's that occur when the line's endpoints go above or below the image.
            c(isnan(c)) = 0;
            firstIndex = find(c, 1, 'first');
            lastIndex = find(c, 1, 'last');
            % Compute the distance of that perpendicular width.
            perpendicularWidth = sqrt( (cx(firstIndex) - cx(lastIndex)) .^ 2 + (cy(firstIndex) - cy(lastIndex)) .^ 2 );
            % Get the average perpendicular width.  This will approximately be the area divided by the longest length.
            averageWidth = measurements(blobIndex).Area / maxDistance_SH;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Start of the code for finding the perpendicular line of the
            %boundary between two consecutive boundary coordinates 
           
            x_first = 0; 
            y_first = 0;
            x_next = 0;
            y_next = 0;
            for idx = 1:length(x_SH)-1 %NOTE: -1 gets rid of the index error
                x_first = x_SH(idx);   %x_first is the x-coordinate of a point on the boundary
                x_next = x_SH(idx+1);  %x_next is the x-coordinate of the proceeding point on the boundary
                y_first = y_SH(idx);   %same as x_first but for the y-coordinate
                y_next = y_SH(idx+1);  %same as x_next but for the y-coordinate
                
                   m_line = (y_next-y_first)/(x_next-x_first); %finds the slope of the line between the two points
                   mid_x_test = x_first+0.5*(x_next-x_first); %finds the x-coordinate of the midpoint between the two points
                   mid_y_test = y_first+0.5*(y_next-y_first); %finds the y-coordinate of the midpoint between the two points
                   
                   new_m = (-1)/(m_line);    %finds the slope of the perpendicular line to the line segment connecting the two points
                   intercept = (-(new_m)*(mid_x_test))+(mid_y_test); %finds the y-intercept of the perpendcular line (b = -mx+y)
                   
                   norm_x_1 = floor(mid_x_test-10); %x-coordinate of the cyan endpoint (NOTE: 10 is just an arbitrary value)
                   norm_x_2 = floor(mid_x_test+10); %x-coordinate of the dark blue endpoint
                   norm_y_1 = floor((new_m)*(norm_x_1)+intercept); %y-coordinate of the cyan endpoint using equation of the perpendicular line (y = mx+b)
                   norm_y_2 = floor((new_m)*(norm_x_2)+intercept); %y-coordinate of the dark blue endpoint using equation of perpendicular line
                   
                   Dist_test_1 = sqrt((norm_x_2-x_SH(idx+1)).^2+(norm_y_2-y_SH(idx+1)).^2); %method to determine the distance from the border point to the endpoint of the normal line
                   Dist_test_2 = sqrt((norm_x_2-x_SH(idx+1)).^2+(norm_y_2-y_SH(idx+1)).^2);
                   
%                    dist_limiter = sqrt((x_SH(idx+1) - norm_x_1).^2 + (y_SH(idx+1)-norm_y_1).^2);
%                    if m_line == 0
%                        plot(mid_x_test,mid_y_test, 'o') %plots a circle at the midpoint
%                        hold on
%                        plot([mid_x_test, mid_x_test], [mid_y_test, mid_y_test+10],'c','LineWidth',1) %plots vertical cyan line +10 from midpoint
%                        plot([mid_x_test,mid_x_test], [mid_y_test, mid_y_test-10],'b','LineWidth',1) %plots vertical dark blue line -10 from midpoint
%                        plot(mid_x_test,mid_y_test+10, 'o', 'Color', 'g');
%                        plot(mid_x_test,mid_y_test-10, 'o', 'Color', 'm');
%                        hold on

                   for idx_comp = 1:length(y_SH)
                       dist_limiter = sqrt((x_SH(idx_comp) - norm_x_1).^2 + (y_SH(idx_comp)-norm_y_1).^2);
                       dist_limiter_2 = sqrt((x_SH(idx_comp) - norm_x_2).^2 + (y_SH(idx_comp) - norm_y_2).^2);
                       
                       if y_SH(idx_comp) == norm_y_1&&dist_limiter<12.5
                            norm_y_1 = NaN;
                            norm_x_1 = NaN;
                       end
                       if y_SH(idx_comp) == norm_y_2&&dist_limiter_2<12.0
                           norm_y_2 = NaN;
                           norm_x_1 = NaN;
                       end
                       
                       in = inpolygon(norm_x_1,norm_y_1,x_SH,y_SH);
                       in_2 = inpolygon(norm_x_2,norm_y_2,x_SH,y_SH);
                       if ~in
                           norm_x_1 = NaN;
                           norm_y_1 = NaN;
                       end
                       if ~in_2
                           norm_x_2 = NaN;
                           norm_y_2 = NaN;
                       end
                   end
                   
                   if Dist_test_1>13.0&&Dist_test_2>12.5 %thresholds the line plotting to reduce noise
                       plot([x_first,x_next],[y_first,y_next],'r') %plots the red line segment connecting the two points
                       hold on
                       plot([mid_x_test,norm_x_1],[mid_y_test,norm_y_1],'c','LineWidth',1) %plots the cyan line
                       plot([mid_x_test,norm_x_2],[mid_y_test,norm_y_2],'b','LineWidth',1) %plots the dark blue line
                       plot(mid_x_test,mid_y_test, 'o') %plots a circle at the midpoint
                       plot(norm_x_1, norm_y_1, 'o', 'Color','g'); %green endpoint of the cyan line
                       plot(norm_x_2, norm_y_2, 'o', 'Color', 'm'); %magenta endpoint of the dark blue line
                       line(norm_x_1, norm_y_1)
                       hold on  
                   end
                   
%                   for idx_plot = 1:length(x_SH)-1
%                       normplot_x = norm_x_
%                       normplot_y = 
%                    
            end
        end
end
            
         
            plot(x_SH, y_SH, 'y-', 'LineWidth', 1);
            % For this blob, put a line between the points farthest away from each other.
            line([x_SH(index1), x_SH(index2)], [y_SH(index1), y_SH(index2)], 'Color', 'r', 'LineWidth', 1);
            plot(xMidPoint, yMidPoint, 'r*', 'MarkerSize', 15, 'LineWidth', 1);
            % Plot perpendicular line.  Make it green across the whole image but magenta inside the blob.
            line([1, columns], [y1, y2], 'Color', 'g', 'LineWidth', 1);	
            line([cx(firstIndex), cx(lastIndex)], [cy(firstIndex), cy(lastIndex)], 'Color', 'm', 'LineWidth', 2);
            impixelinfo
   % end
    hold off;
%     cd(Primary_Folder)
 %end
%end   