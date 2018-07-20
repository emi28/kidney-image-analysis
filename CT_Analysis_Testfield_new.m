clear, clc, format short g

fontsize = 20;

% Make some overall designations first

Number_Slices = input('Enter the number of slices you want to analyze per time point: ');
Number_Medulla_Cortex_Regions = input('Enter the number of medullary and cortical sections you want to analyze per slice: ');
radius = input('Enter the radius size [pixels] for grayscale analysis in circle form (~5): ');

% Define the primary folder through directory designation

A = dir; 
Primary_Folder = A.folder;

disp('  ')
 
% for idx = 1:length(Primary_Folder)
%     name = Primary_Folder(idx);
%     number = sscanf(name, '%02d.jpg', 1); 
%     movefile(name, [num2str(number) '.jpg']); 
% end

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
            
            I2 = adapthisteq(I);

            IC = uint8(I2);
            IC_size = size (IC);
            IC_rows = uint16(IC_size(1));
            IC_columns = uint16(IC_size(2));
            
            IC_2 = imcrop(IC, [30 115 220 280]);
            thresh = graythresh(IC)
            
            numberofcolorbands = 1
            [rows columns numberofcolorbands] = size(IC)
            
            BW = im2bw(IC)
            BW_2 = imclearborder(BW)
            BW_3 = bwareaopen(BW_2,150)
            figure
            imshow(BW_3)
            
            imshow(IC);
            hold on;
            boundaries = bwboundaries(BW_3);	
            numberOfBoundaries = size(boundaries);
            for i = 1 : numberOfBoundaries
                    thisBoundary = boundaries{i};
                    plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 1);
            end
            hold off;
            
            maskedImage = zeros(rows, columns, 'uint8');
            maskedImage(BW_3) = IC(BW_3); 
            imshow(maskedImage);
            
            BW_4 = im2bw(maskedImage);
            imshow(maskedImage);
            hold on;
            boundaries_2 = bwboundaries(maskedImage);	
            numberOfBoundaries_2 = size(boundaries_2);
            for j = 1 : numberOfBoundaries_2
                    thisBoundary_2 = boundaries_2{j};
                    plot(thisBoundary_2(:,2), thisBoundary_2(:,1), 'g', 'LineWidth', 1);
            end
            hold off;
            
            [r c] = find(bwperim(maskedImage))
            
            Hypotonic_region_pixels = zeros(I_rows,I_columns,'uint16');
            Loop_count_HR = uint16(0);

            for dd = 1:IC_rows
                for ee = 1:IC_columns
                    Hypotonic_region_pixels(dd) = r(dd);
                    Hypotonic_region_pixels(ee) = c(ee);
                    if Hypotonic_region_pixels(dd,ee) == 1
                        Loop_count_HR = Loop_count_HR+1;
                        I_hold_HR(Loop_count_HR) = I(dd,ee);
                    end
                end
            end
            Hypotonic_region_I(aa,bb,cc) = mean(I_hold_HR)  
        end
end