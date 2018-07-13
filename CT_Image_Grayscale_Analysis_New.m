clear, clc, format short g

% Make some overall designations first

Number_Slices = input('Enter the number of slices you want to analyze per time point: ');
Number_Medulla_Cortex_Regions = input('Enter the number of medullary and cortical sections you want to analyze per slice: ');
radius = input('Enter the radius size [pixels] for grayscale analysis in circle form (~5): ');

% Define the primary folder through directory designation

A = dir; 
Primary_Folder = A.folder;

disp('  ')

% Primary loop to obtain average grayscale values

for aa = 1:(length(dir)-3)
    
    % Selecting the folder for a given time
    
    Folder = A(aa+2).name;
    disp(strcat('You are analyzing the current folder:',{' '},Folder)) 
    cd(Folder)
    
        % Selecting the hypotonic image and isotonic image you want to analyze
        
        for bb = 1:uint16(Number_Slices)
            
            % Hypotonic image
            
            dirFolder = dir('*.jpg');                 %get the selected file data
            fileNames = {dirFolder.name};             %create a cell array of file names
            for iFile = 1:numel(fileNames)            %loop over the file names
                newName = sprintf('%02d.jpg',iFile);  %make the new name
                movefile(fileNames{iFile},newName);   %rename the file
            end 
            
            filename_designation = char(strcat('Enter hypotonic image',{' '},num2str(bb),{' '},'of',{' '},num2str(Number_Slices),':',{' '}));
            filename = input(filename_designation,'s');
            
            disp('  ')
            
            I = double(imread(filename));
            I_size = size(I);
            I_rows = uint16(I_size(1));
            I_columns = uint16(I_size(2));
            figure
            imshow(filename)
            
                % Select the points of interest for the hypotonic kidney
            
                for cc = 1:uint16(Number_Medulla_Cortex_Regions)
                    
                    % Choose the points and draw the circles

                    disp(strcat('Click the',{' '},num2str(cc),{' '},'of',{' '},num2str(Number_Medulla_Cortex_Regions),{' '},'hypotonic medulla point on the image and then hit enter'))
                    [Hypotonic_Medulla_Point_X,Hypotonic_Medulla_Point_Y] = getpts;
                    Hypotonic_Medulla_Point_Matrix = [Hypotonic_Medulla_Point_X,Hypotonic_Medulla_Point_Y];
                    viscircles(Hypotonic_Medulla_Point_Matrix,radius,'Color','b','LineStyle','--');

                    disp(strcat('Click the',{' '},num2str(cc),{' '},'of',{' '},num2str(Number_Medulla_Cortex_Regions),{' '},'hypotonic cortex point on the image and then hit enter'))
                    [Hypotonic_Cortex_Point_X,Hypotonic_Cortex_Point_Y] = getpts;
                    Hypotonic_Cortex_Point_Matrix = [Hypotonic_Cortex_Point_X,Hypotonic_Cortex_Point_Y];
                    viscircles(Hypotonic_Cortex_Point_Matrix,radius,'Color','b');
                
                    % Define the logical matrix for pixels to keep

                        % Preallocate

                            Hypotonic_Medulla_Pixels = zeros(I_rows,I_columns,'uint16');
                            Hypotonic_Cortex_Pixels = zeros(I_rows,I_columns,'uint16');

                        % Initialize the if statement count variable

                            Loop_Count_Hypotonic_Medulla = uint16(0);
                            Loop_Count_Hypotonic_Cortex = uint16(0);

                        % Loop to define logical values and keep grayscale values for averaging

                            for dd = 1:I_rows

                                for ee = 1:I_columns

                                    Hypotonic_Medulla_Pixels(dd,ee) = ((double(dd)-Hypotonic_Medulla_Point_Y)^2)+((double(ee)-Hypotonic_Medulla_Point_X)^2) <= (radius^2);
                                    Hypotonic_Cortex_Pixels(dd,ee) = ((double(dd)-Hypotonic_Cortex_Point_Y)^2)+((double(ee)-Hypotonic_Cortex_Point_X)^2) <= (radius^2);

                                        if Hypotonic_Medulla_Pixels(dd,ee) == 1

                                            Loop_Count_Hypotonic_Medulla = Loop_Count_Hypotonic_Medulla+1;

                                            I_hold_Hypotonic_Medulla(Loop_Count_Hypotonic_Medulla) = I(dd,ee);

                                        end

                                        if Hypotonic_Cortex_Pixels(dd,ee) == 1

                                            Loop_Count_Hypotonic_Cortex = Loop_Count_Hypotonic_Cortex+1;

                                            I_hold_Hypotonic_Cortex(Loop_Count_Hypotonic_Cortex) = I(dd,ee);

                                        end

                                end

                            end

                            Hypotonic_Medulla_I(aa,bb,cc) = mean(I_hold_Hypotonic_Medulla);
                            Hypotonic_Cortex_I(aa,bb,cc) = mean(I_hold_Hypotonic_Cortex);
                             
                end
                
                % Save figure as jpg and clear figure
                    
                saveas(gcf,char(strcat('Hypotonic',{' '},Folder,{' '},'Slice',{' '},num2str(bb))),'jpg')
                close
                            
            % Isotonic image
            
            dirFolder = dir('*.jpg');                 %get the selected file data
            fileNames = {dirFolder.name};             %create a cell array of file names
            for iFile = 1:numel(fileNames)            %loop over the file names
                newName = sprintf('%02d.jpg',iFile);  %make the new name
                movefile(fileNames{iFile},newName);   %rename the file
            end 
            
            filename_designation = char(strcat('Enter isotonic image',{' '},num2str(bb),{' '},'of',{' '},num2str(Number_Slices),':',{' '}));
            filename = input(filename_designation,'s');
            
            disp('  ')
            
            I = double(imread(filename));
            I_size = size(I);
            I_rows = uint16(I_size(1));
            I_columns = uint16(I_size(2));
            figure
            imshow(filename)
            
                % Select the points of interest for the isotonic kidney
            
                for cc = 1:uint16(Number_Medulla_Cortex_Regions)
                    
                    % Choose the points and draw the circles

                    disp(strcat('Click the',{' '},num2str(cc),{' '},'of',{' '},num2str(Number_Medulla_Cortex_Regions),{' '},'isotonic medulla point on the image and then hit enter'))
                    [Isotonic_Medulla_Point_X,Isotonic_Medulla_Point_Y] = getpts;
                    Isotonic_Medulla_Point_Matrix = [Isotonic_Medulla_Point_X,Isotonic_Medulla_Point_Y];
                    viscircles(Isotonic_Medulla_Point_Matrix,radius,'Color','g','LineStyle','--');

                    disp(strcat('Click the',{' '},num2str(cc),{' '},'of',{' '},num2str(Number_Medulla_Cortex_Regions),{' '},'isotonic cortex point on the image and then hit enter'))
                    [Isotonic_Cortex_Point_X,Isotonic_Cortex_Point_Y] = getpts;
                    Isotonic_Cortex_Point_Matrix = [Isotonic_Cortex_Point_X,Isotonic_Cortex_Point_Y];
                    viscircles(Isotonic_Cortex_Point_Matrix,radius,'Color','g');

                    % Define the logical matrix for pixels to keep

                        % Preallocate

                            Isotonic_Medulla_Pixels = zeros(I_rows,I_columns,'uint16');
                            Isotonic_Cortex_Pixels = zeros(I_rows,I_columns,'uint16');

                        % Initialize the if statement count variable

                            Loop_Count_Isotonic_Medulla = uint16(0);
                            Loop_Count_Isotonic_Cortex = uint16(0);

                        % Loop to define logical values and keep grayscale values for averaging

                            for dd = 1:I_rows

                                for ee = 1:I_columns

                                    Isotonic_Medulla_Pixels(dd,ee) = ((double(dd)-Isotonic_Medulla_Point_Y)^2)+((double(ee)-Isotonic_Medulla_Point_X)^2) <= (radius^2);
                                    Isotonic_Cortex_Pixels(dd,ee) = ((double(dd)-Isotonic_Cortex_Point_Y)^2)+((double(ee)-Isotonic_Cortex_Point_X)^2) <= (radius^2);

                                        if Isotonic_Medulla_Pixels(dd,ee) == 1

                                            Loop_Count_Isotonic_Medulla = Loop_Count_Isotonic_Medulla+1;

                                            I_hold_Isotonic_Medulla(Loop_Count_Isotonic_Medulla) = I(dd,ee);

                                        end

                                        if Isotonic_Cortex_Pixels(dd,ee) == 1

                                            Loop_Count_Isotonic_Cortex = Loop_Count_Isotonic_Cortex+1;

                                            I_hold_Isotonic_Cortex(Loop_Count_Isotonic_Cortex) = I(dd,ee);

                                        end

                                end

                            end

                            Isotonic_Medulla_I(aa,bb,cc) = mean(I_hold_Isotonic_Medulla);
                            Isotonic_Cortex_I(aa,bb,cc) = mean(I_hold_Isotonic_Cortex);
                    
                end
                
                % Save figure as jpg and close figure
                    
                saveas(gcf,char(strcat('Isotonic',{' '},Folder,{' '},'Slice',{' '},num2str(bb))),'jpg')
                close
                
        end
            
    cd(Primary_Folder)
    
end

% Results Presentation

Time_Plot = [0,30,31,32,33,34,35,40,45,50];

    % Average across slices and number of points for left medulla, right medulla, left cortex, and right cortex
    
    for ff = 1:aa
    
        Hypotonic_Medulla_I_Plot(ff) = mean(Hypotonic_Medulla_I(ff,:,:));
        Isotonic_Medulla_I_Plot(ff) = mean(Isotonic_Medulla_I(ff,:,:));
        Hypotonic_Cortex_I_Plot(ff) = mean(Hypotonic_Cortex_I(ff,:,:));
        Isotonic_Cortex_I_Plot(ff) = mean(Isotonic_Cortex_I(ff,:,:));

    end

figure

plot(Time_Plot,Hypotonic_Medulla_I_Plot,'bx--',Time_Plot,Isotonic_Medulla_I_Plot,'gx--',Time_Plot,Hypotonic_Cortex_I_Plot,'bo-',Time_Plot,Isotonic_Cortex_I_Plot,'go-')
xlabel('Time [min]')
ylabel('Grayscale Average')
Label = {'Hypotonic Medullary Region','Isotonic Medullary Region','Hypotonic Cortical Region','Isotonic Cortical Region'};
legend(Label)

% Display Results and write to an Excel spreadsheet for easy export

Results = [Time_Plot',Hypotonic_Medulla_I_Plot',Isotonic_Medulla_I_Plot',Hypotonic_Cortex_I_Plot',Isotonic_Cortex_I_Plot'];

fprintf('%30s %30s %30s %30s %30s\n','Time',Label{1,1},Label{1,2},Label{1,3},Label{1,4})
fprintf('%30.f %30.4f %30.4f %30.4f %30.4f\n',Results')

Results_Pres = cell((aa+1),5);

Results_Pres{1,gg} = 'Time';

    for gg = 1:5
        
        if gg ==1
        
            Results_Pres{1,gg} = 'Time';

        else
            
            Results_Pres{1,gg} = Label{1,(gg-1)};
        
        end

    end
    
    for hh = 1:5
        
        for ii = 2:(aa+1)
            
            Results_Pres{ii,hh} = Results((ii-1),hh);
            
        end
        
    end

xlswrite('CT_Image_Data.xlsx',Results_Pres)
