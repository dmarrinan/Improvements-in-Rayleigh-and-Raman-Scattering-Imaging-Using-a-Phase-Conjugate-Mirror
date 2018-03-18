clc
close all
clear all

% p = path;
save_info_cell = cell(18,5);
dark_folder = 'C:\Users\Dan\Desktop\Research\Ohio State\Prof Sutton 2013\11_7_13_dark';
image_folder = 'C:\Users\Dan\Desktop\Research\Ohio State\Prof Sutton 2013\11_7_13_air_singlepass';
resolution_folder = 'C:\Users\Dan\Desktop\Research\Ohio State\Prof Sutton 2013\11_7_13_resolution';

average_dark_image_title = [dark_folder '\dark_average_image.xlsx'];
average_resolution_image_title = [resolution_folder '\resolution_average_image.xlsx'];
resolution_value_title = [resolution_folder '\resolution_calculated_values.xlsx'];
image_title = '11 7 13 air single pass';
average_image_title = [image_folder '\Average Dark Image ' image_title '.jpg'];


calculate_dark_flag=1;
plot_dark_flag=0;
plot_image_flag=0;
analyze_resolution_images_flag=1;
calculate_resolution_flag=1;

num_pix_x = 1376; %688;
num_pix_y = 1040; %520;
images_to_skip = []; %[5,8,10,11,19,27,30,41,52,63,68,72,74,85,96]; %[15,17,24,35,37,46,47,55,57,65,68,73,78,79,85,90]; %[11,14,16,22,28,29,30,33,44,51,55,66,67,69,74,77,88,96,99]; %[6,10,17,28,31,39,41,50,58,61,72,76,83,91,92];
if ~isempty(images_to_skip)
    xlswrite([image_folder '\Images_to_skip_' image_title '.xlsx'],images_to_skip);
end

save_info_cell(1,:) = {'Description of Information', 'Value of Information','Min','Max','Step Size'};
save_info_cell(2,1:2) = {'dark average title',average_dark_image_title};
save_info_cell(3,1:2) = {'resolution title',average_resolution_image_title};
save_info_cell(4,1:2) = {'image title',average_image_title};
save_info_cell(5,1:2) = {'number of x pixels', num_pix_x};
save_info_cell(6,1:2) = {'number of y pixels', num_pix_y};

%find average background
number_of_dark_images = 100;
if calculate_dark_flag == 1
%     path(p,dark_folder);
    fprintf('Analyzing dark images to subtract background noise...\n');
    
    dark_sum = zeros(num_pix_y,num_pix_x);
    if plot_dark_flag==1
        figure
    end
    for i = 1:number_of_dark_images;
        if i<10
            filename = sprintf('B0000%i.tif',i);
        elseif i<100
            filename = sprintf('B000%i.tif',i);
        elseif i<1000
            filename = sprintf('B00%i.tif',i);
        end
        filename = [dark_folder '\' filename];
        dark_sum = dark_sum + double(imread(filename));
        if plot_dark_flag==1
            caxis([0 1500])
            image(imread(filename),'CDataMapping','scaled');
            title(filename);
            pause(0.2)
        end
    end
    
    dark_average = dark_sum/number_of_dark_images;
    
    xlswrite(average_dark_image_title,dark_average)
%     path(p);
else
    fprintf('Loading background noise to be subtracted...\n');
    dark_average = xlsread(average_dark_image_title);
end

figure
caxis([0 1500])
image(dark_average,'CDataMapping','scaled');
title('Average Dark Image')
colorbar
saveas(gcf,average_image_title)

% path(p,image_folder);
number_of_images = 500;
fprintf('Analyzing Images...\n');
image_sum=zeros(num_pix_y,num_pix_x);
all_images = zeros(num_pix_y,num_pix_x,number_of_images);
if plot_image_flag==1
    figure
end
for i=1:number_of_images
    if i<10
        filename=sprintf('B0000%i.tif',i);
    elseif i<100
        filename=sprintf('B000%i.tif',i);
    elseif i<1000
        filename=sprintf('B00%i.tif',i);
    end
    filename = [image_folder '\' filename];
    if(sum(i==images_to_skip)<1)
        all_images(:,:,i) = double(imread(filename))-dark_average;
        image_sum = image_sum+double(imread(filename));
        if plot_image_flag==1
            caxis([0 1500])
            image(imread(filename),'CDataMapping','scaled')
            title(filename)
            pause(0.2);
        end
    end
end
number_of_images_used = number_of_images-length(images_to_skip);
image_average = image_sum/number_of_images_used-dark_average;
figure
caxis([0 1500])
image(image_average,'CDataMapping','scaled')
title('Average Image')
colorbar
saveas(gcf,[image_folder '\Average Image ' image_title '.jpg'])
saveas(gcf,[image_folder '\Average Image ' image_title '.fig'])
% path(p);

%background subtraction
avg_top = mean(image_average(100:300,:),1);
avg_bottom = mean(image_average(800:1000,:),1);
background_slope = (avg_bottom - avg_top)/(900-200);
background_subtraction = zeros(size(image_average));
for i = 1:1040
    background_subtraction(i,:) = background_slope*(i-200)+avg_top;
end

background_subtracted_image_average = image_average-background_subtraction;
figure
caxis([0 1500])
image(background_subtracted_image_average,'CDataMapping','scaled')
title('Background Subtracted Average Image')
colorbar
saveas(gcf,[image_folder '\Background Subtracted Average Image ' image_title '.jpg'])
saveas(gcf,[image_folder '\Background Subtracted Average Image ' image_title '.fig'])

image_1D = sum(background_subtracted_image_average,1);
figure
plot(image_1D)
title('1D Image')
xlabel('X Pixel')
ylabel('Sum of Counts')
saveas(gcf,[image_folder '\Image 1D ' image_title '.fig'])
save([image_folder '\1D_image.mat'],'image_1D')

%for each column find number of rows that are at least 1/e^2 (~13.5%) of
    %max intensity
fprintf('Calculating Beam Width...\n')
avg_beam_width= zeros(1,num_pix_x);
%find value and index of maximum intensity
[avg_max_row, avg_ind_row] = max(background_subtracted_image_average);
[avg_max_val, avg_ind_col] = max(avg_max_row);
avg_ind_max = [avg_ind_row(avg_ind_col) avg_ind_col];
%find how many pixels in each column are larger than 1/e^2
for i = 1:num_pix_x
    avg_beam_width(i) = sum((background_subtracted_image_average(:,i)/avg_max_val)>1/exp(2));
end
%plot beam width for each column of pixel array
figure
plot(avg_beam_width)
title('Beam Width per Column for Average Image')
ylabel('Beam Width (pixels)')
xlabel('Column Number')
saveas(gcf,[image_folder '\Beam Width ' image_title '.jpg']);

%find beam profile of column with maximum intensity
beam_profile = background_subtracted_image_average(:,avg_ind_col);
figure
plot(beam_profile,'x')
title('Beam Profile of Max Intensity Column for Average Image')
ylabel('Intensity')
xlabel('Row Number')
saveas(gcf,[image_folder '\Beam Profile ' image_title '.jpg']);
saveas(gcf,[image_folder '\Beam Profile ' image_title '.fig']);

avg_width_val = beam_profile/avg_max_val>1/exp(2);
[avg_first_beam_val,avg_first_beam_ind] = max(avg_width_val);
avg_last_beam_ind = avg_first_beam_ind+length(beam_profile(avg_width_val))-1;
avg_last_beam_val = beam_profile(avg_last_beam_ind);
hold on
plot(avg_first_beam_ind,beam_profile(avg_first_beam_ind),'or')
plot(avg_last_beam_ind,beam_profile(avg_last_beam_ind),'ok')
hold off
saveas(gcf,[image_folder '\Beam Profile ' image_title '.jpg']);

%compare average image to individual images
beam_width = zeros(number_of_images_used,1);
peak_intensity = zeros(number_of_images_used,1);
peak_intensity_normalized = zeros(number_of_images_used,1);
displacement_peak_row = zeros(number_of_images_used,1);
displacement_peak_column = zeros(number_of_images_used,1);

max_row = zeros(number_of_images_used,num_pix_x);
ind_row = zeros(number_of_images_used,num_pix_x);
max_val = zeros(number_of_images_used,1);
ind_col = zeros(number_of_images_used,1);
ind_max = zeros(number_of_images_used,1);
rms_sum = zeros(num_pix_y,num_pix_x);

for i=1:number_of_images
        %find how many pixels in each column are larger than 1/e^2
        [peak_intensity(i),ind_max(i)] = max(all_images(:,avg_ind_col,i));
        peak_intensity_normalized(i) = peak_intensity(i)/avg_max_val;
        beam_width(i) = sum((all_images(:,avg_ind_col,i)/peak_intensity(i))>1/exp(2));
        displacement_peak_row(i) = avg_ind_max(1)-ind_max(i);
        rms_sum = rms_sum + (image_average - all_images(:,:,i)).^2;
end

rms = sqrt(rms_sum/(number_of_images-1));
figure
caxis([0 1500])
image(rms,'CDataMapping','scaled')
colorbar
title('RMS')

signal_to_noise = image_average./rms;
figure
image(signal_to_noise,'CDataMapping','scaled')
title('Image Average/RMS')

x_peak_intensity = 200:50:850;
n_peak_intensity = hist(peak_intensity,x_peak_intensity);
figure
bar(x_peak_intensity,n_peak_intensity)
% set(gca,'XTick',[x_peak_intensity-25, max(x_peak_intensity)+25])
peak_intensity_title = ['Histogram of Peak Intensity ' image_title];
title(peak_intensity_title)
xlabel('Peak Intensity')
saveas(gcf,[image_folder, '\' peak_intensity_title '.jpg'])

x_peak_intensity_normalized = 0.45:0.1:1.55;
n_peak_intensity_normalized = hist(peak_intensity_normalized,x_peak_intensity_normalized);
figure
bar(x_peak_intensity_normalized,n_peak_intensity_normalized);
peak_intensity_normalized_title = ['Histogram of Normalized Peak Intensity ' image_title];
title(peak_intensity_normalized_title)
xlabel('Normalized Peak Intensity')
saveas(gcf,[image_folder, '\' peak_intensity_normalized_title '.jpg'])

x_beam_width = 12:1:24;
n_beam_width = hist(beam_width,x_beam_width);
figure
bar(x_beam_width,n_beam_width)
beam_width_title = ['Histogram of Beam Width ' image_title];
title(beam_width_title)
xlabel('Beam Width (pixels)')
saveas(gcf,[image_folder '\' beam_width_title '.jpg'])
saveas(gcf,[image_folder '\' beam_width_title '.fig'])

x_displacement_peak_row = -3.5:2.5;
n_displacement_peak_row = hist(displacement_peak_row,x_displacement_peak_row);
figure
bar(x_displacement_peak_row,n_displacement_peak_row)
displacement_peak_row_title = ['Histogram of Displacement of Peak Row Location ' image_title];
title(displacement_peak_row_title)
xlabel('Displacement of Peak Location (pixels in y direction)')
saveas(gcf,[image_folder '\' displacement_peak_row_title '.jpg'])
saveas(gcf,[image_folder '\' displacement_peak_row_title '.fig'])
% x_displacement_peak_column = -450:100:350;
% n_displacement_peak_column = hist(displacement_peak_column,x_displacement_peak_column);
% figure
% bar(x_displacement_peak_column,n_displacement_peak_column)
% set(gca,'XTick',[x_displacement_peak_column-50,max(x_displacement_peak_column)+50])
% displacement_peak_column_title = ['Histogram of Displacement of Peak Column Location ' image_title];
% title(displacement_peak_column_title)
% xlabel('Displacement of Peak Location (pixels in x direction)')
% saveas(gcf,[image_folder '\' displacement_peak_column_title '.jpg'])

number_of_resolution_images = 10;
if analyze_resolution_images_flag == 1
%     path(p,resolution_folder);
    resolution_sum = zeros(num_pix_y,num_pix_x);
    for i=1:number_of_resolution_images
        if i<10
            filename = sprintf('B0000%i.tif',i);
        elseif i<100
            filename = sprintf('B000%i.tif',i);
        elseif i <1000
            filename = sprintf('B00%i.tif',i);
        end
        filename = [resolution_folder '\' filename];
        resolution_sum = resolution_sum + double(imread(filename));
    end
    resolution_avg = resolution_sum/number_of_resolution_images;
    xlswrite(average_resolution_image_title,resolution_avg);
else
    resolution_avg = xlsread(average_resolution_image_title);
end
if calculate_resolution_flag == 1
    figure
    image(resolution_avg)
    [resolution_x,resolution_y] = ginput(2);
    resolution = 1/4/abs(resolution_x(2)-resolution_x(1))*2.54/100; %mm/cross*(m/mm)/(pixels/cross) = m/pixel %4/1000/(resolution_x(2)-resolution_x(1)); %mm/cross*(m/mm)/(pixels/cross) = m/pixel    %1/4/abs(resolution_x(2)-resolution_x(1))*2.54/100; %(inch/marking)*(marking/pixels)*(cm/inch)*(m/cm)
    magnification = 6.675*10^-6/resolution;
    resolution_cell{1,1}='resolution';
    resolution_cell{1,2}=resolution;
    resolution_cell{2,1}='magnification';
    resolution_cell{2,2}=magnification;
    xlswrite(resolution_value_title,resolution_cell)
%     path(p);    
else
    num_resolution = xlsread(resolution_value_title);
    resolution = num_resolution(1);
    magnification = num_resolution(2);
end

beam_width = beam_width*resolution;
x_beam_width = x_beam_width*resolution;
n_beam_width = hist(beam_width,x_beam_width);
figure
bar(x_beam_width,n_beam_width)
beam_width_title = ['Histogram of Beam Width ' image_title];
title(beam_width_title)
xlabel('Beam Width (m)')
saveas(gcf,[image_folder '\' beam_width_title '_m.jpg'])

beam = background_subtracted_image_average/avg_max_val>1/exp(2);
sum_beam = sum(background_subtracted_image_average(background_subtracted_image_average/avg_max_val>1/exp(2)));

save_info_cell(7,1:2) = {'number of dark images', number_of_dark_images};
save_info_cell(8,1:2) = {'number of images', number_of_images};
save_info_cell(9,1:2) = {'number of images analyzed', number_of_images_used};
save_info_cell(10,1:2) = {'number of resolution images', number_of_resolution_images};
save_info_cell(11,:) = {'histogram range for beam width','', min(x_beam_width),max(x_beam_width),x_beam_width(2)-x_beam_width(1)};
save_info_cell(12,:) = {'histogram range for peak intensity','', min(x_peak_intensity),max(x_peak_intensity),x_peak_intensity(2)-x_peak_intensity(1)};
save_info_cell(13,:) = {'histogram range for normalized peak intensity','', min(x_peak_intensity_normalized),max(x_peak_intensity_normalized),x_peak_intensity_normalized(2)-x_peak_intensity_normalized(1)};
save_info_cell(14,:) = {'histogram range for relative row location of maximum value','', min(x_displacement_peak_row),max(x_displacement_peak_row),x_displacement_peak_row(2)-x_displacement_peak_row(1)};
%save_info_cell(15,:) = {'histogram range for relative column location of maximum value','', min(x_displacement_peak_column),max(x_displacement_peak_column),x_displacement_peak_column(2)-x_displacement_peak_column(1)};
save_info_cell(16,1:2) = {'resolution',resolution};
save_info_cell(17,1:2) = {'magnification',magnification};
save_info_cell(18,1:2) = {'binned signal',sum_beam};
xlswrite([image_folder '\Important Info' image_title '.xlsx'],save_info_cell);