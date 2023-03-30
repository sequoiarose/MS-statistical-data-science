%experiment set up for Phase III Final Project:
%Idea: smaller mult is consistently better, usually smaller deviation gain
%is also better. See if there is interaction between these results and
%other vals 

%Parameter tuning for ideal image set
path = fileparts(mfilename('fullpath'));
images_path = path + "\edge detection data";
exp_nums = [1,2,3,4,5]; %only looking at closly touching objects for each of the 5 experiments
nscale_vals = [3, 4, 5];
minWaveLength_vals = [2,3,4];
mult_vals = [1.5, 2.1, 2.5];
sigma_vals = [0.45, 0.55, 0.65];
k_vals = [1, 3, 10, 15];
cutoff_vals = [0.4, 0.5, 0.6];
deviationGain_vals = [1, 1.5, 2];
g_vals = [5, 10, 15];
noise_method_vals = [-2, -1, 1];
params = {'nscale', 'minWaveLength', 'sigmaOnf', 'k', 'cutOff', 'g', 'noiseMethod'};
num_exp = length(exp_nums);
data_dict = containers.Map();
for exp_num = exp_nums
    dev_count = 1;
    %read in images
    images = read_images(images_path, exp_num);
    %take only the third image to look at closly touching objects
    images = images(3);
    %run experiment, get output images, phase congrunecy, average phase,
    for deviationGain = deviationGain_vals
        mult_count = 1;
        for mult = mult_vals
            col = strings(length(params), 1);
            row_count = 1; 
            for param = params
                param_clean = param(1);
                if strcmp(param_clean,'nscale') == 1
                    param_list = nscale_vals;
                elseif strcmp(param_clean,'minWaveLength') == 1
                    param_list = minWaveLength_vals;
                elseif strcmp(param_clean,'sigmaOnf') == 1
                    param_list = sigma_vals;
                elseif strcmp(param_clean, 'k') == 1
                    param_list = k_vals;
                elseif strcmp(param_clean, 'g') == 1
                    param_list = g_vals;
                elseif strcmp(param_clean, 'cutOff') == 1
                    param_list = cutoff_vals;
                elseif strcmp(param_clean, 'noiseMethod') == 1
                    param_list = noise_method_vals;
                end
                %calculate mse
                mse_list = phase_congruency_parameter_tuning(images, mult, deviationGain, param_clean, param_list);
                %store mse_list
                tf = isKey( data_dict , "mult"+string(mult_count)+"dev"+string(dev_count)+param );
                if tf == true
                    data_dict("mult"+string(mult_count)+"dev"+string(dev_count)+param) = cat(2, data_dict("mult"+string(mult_count)+"dev"+string(dev_count)+param), mse_list);
                elseif tf == false
                    data_dict("mult"+string(mult_count)+"dev"+string(dev_count)+param) = mse_list;
                end
            end
            mult_count = mult_count+1;
        end
        dev_count = dev_count +1;
    end
end
mult1dev1 = strings(length(params), 1);
mult2dev1 = strings(length(params), 1);
mult3dev1 = strings(length(params), 1);
mult1dev2 = strings(length(params), 1);
mult2dev2 = strings(length(params), 1);
mult3dev2 = strings(length(params), 1);
mult1dev3 = strings(length(params), 1);
mult2dev3 = strings(length(params), 1);
mult3dev3 = strings(length(params), 1);
i = 1;
for param=params
    mult1dev1(i) = string(min(data_dict("mult1"+"dev1"+param)))+"-"+string(max(data_dict("mult1"+"dev1"+param)))+", "+string(max(data_dict("mult1"+"dev1"+param))-min(data_dict("mult1"+"dev1"+param)));
    mult2dev1(i) = string(min(data_dict("mult2"+"dev1"+param)))+"-"+string(max(data_dict("mult2"+"dev1"+param)))+", "+string(max(data_dict("mult2"+"dev1"+param))-min(data_dict("mult2"+"dev1"+param)));
    mult3dev1(i) = string(min(data_dict("mult3"+"dev1"+param)))+"-"+string(max(data_dict("mult3"+"dev1"+param)))+", "+string(max(data_dict("mult3"+"dev1"+param))-min(data_dict("mult3"+"dev1"+param)));
    mult1dev2(i) = string(min(data_dict("mult1"+"dev2"+param)))+"-"+string(max(data_dict("mult1"+"dev2"+param)))+", "+string(max(data_dict("mult1"+"dev2"+param))-min(data_dict("mult1"+"dev2"+param)));
    mult2dev2(i) = string(min(data_dict("mult2"+"dev2"+param)))+"-"+string(max(data_dict("mult2"+"dev2"+param)))+", "+string(max(data_dict("mult2"+"dev2"+param))-min(data_dict("mult2"+"dev2"+param)));
    mult3dev2(i) = string(min(data_dict("mult3"+"dev2"+param)))+"-"+string(max(data_dict("mult3"+"dev2"+param)))+", "+string(max(data_dict("mult3"+"dev2"+param))-min(data_dict("mult3"+"dev2"+param)));
    mult1dev3(i) = string(min(data_dict("mult1"+"dev3"+param)))+"-"+string(max(data_dict("mult1"+"dev3"+param)))+", "+string(max(data_dict("mult1"+"dev3"+param))-min(data_dict("mult1"+"dev3"+param)));
    mult2dev3(i) = string(min(data_dict("mult2"+"dev3"+param)))+"-"+string(max(data_dict("mult2"+"dev3"+param)))+", "+string(max(data_dict("mult2"+"dev3"+param))-min(data_dict("mult2"+"dev3"+param)));
    mult3dev3(i) = string(min(data_dict("mult3"+"dev3"+param)))+"-"+string(max(data_dict("mult3"+"dev3"+param)))+", "+string(max(data_dict("mult3"+"dev3"+param))-min(data_dict("mult3"+"dev3"+param)));
    i = i+1;
end

row_names = {'nscale', 'minWaveLength', 'sigmaOnf', 'k', 'cutOff', 'g', 'noiseMethod'};
t = table(table(mult1dev1, mult2dev1, mult3dev1, 'VariableNames', {'Mult=1.5', 'Mult=2.1', 'Mult=2.5'}), ...
table(mult1dev2, mult2dev2, mult3dev2, 'VariableNames', {'Mult=1.5', 'Mult=2.1', 'Mult=2.5'}), ...
table(mult1dev3, mult2dev3, mult3dev3, 'VariableNames', {'Mult=1.5', 'Mult=2.1', 'Mult=2.5'}), ...
'VariableNames', {'DeviationGain=1', 'DeviationGain=1.5', 'DeviationGain=2'}, ...
'RowNames', {'nscale', 'minWaveLength', 'sigmaOnf', 'k', 'cutOff', 'g', 'noiseMethod'});
%save table
myFolder = images_path;
file = fullfile(myFolder, 'table.csv');
t2 = splitvars(t);
writetable(t2, file, 'WriteRowNames', true);

%Parasite data:

images_path = path+"\juvenile"; %set to correct path, it should work as is
treatments = ["010312-B8-1-2", "010312-control-000-2", "032612-2CPT-10-3-B", "032612-Prava-1-3", ...
               "032911-PZQ-10-1", "040212-Simva-10-1", "061112-Ator-01-1-c", "082112-B8-10-4", ...
               "082112-Prom-1-4", "092611-A2-1-2", "092611-D7-10-4"];
i = 1;
for treatment = treatments
    [images, binary] = read_parasite_images(images_path, treatment);
    total_images{i} = images;
    total_binary{i} = binary;
    i = i + 1;
end

%flatten  images
k=1;
for i = 1:11
    im_set = total_images{i};
    len_im_set = length(total_images{i});
    for j = 1:len_im_set
        im = im_set{j};
        images{k} = im;
        k = k+1;
    end
end
k=1;
for i = 1:11
    im_set = total_binary{i};
    len_im_set = length(total_binary{i});
    for j = 1:len_im_set
        im = im_set{j};
        binarys{k} = im;
        k = k+1;
    end
end

%Parameter Tuning: WARNING VERY LONG RUN TIME (~0.5 day)
dev_count = 1;
data_dict = dictionary;
for deviationGain = deviationGain_vals
    mult_count = 1;
    for mult = mult_vals
        col = strings(length(params), 1);
        row_count = 1; 
        for param = params
            param_clean = param(1);
            if strcmp(param_clean,'nscale') == 1
                param_list = nscale_vals;
            elseif strcmp(param_clean,'minWaveLength') == 1
                param_list = minWaveLength_vals;
            elseif strcmp(param_clean,'sigmaOnf') == 1
                param_list = sigma_vals;
            elseif strcmp(param_clean, 'k') == 1
                param_list = k_vals;
            elseif strcmp(param_clean, 'g') == 1
                param_list = g_vals;
            elseif strcmp(param_clean, 'cutOff') == 1
                param_list = cutoff_vals;
            elseif strcmp(param_clean, 'noiseMethod') == 1
                param_list = noise_method_vals;
            end
            %calculate mse
            mse_list = phase_congruency_parameter_tuning_parasites(images, binarys, mult, deviationGain, param_clean, param_list);
            %calculate mse range
            min_mse = min(mse_list);
            max_mse = max(mse_list);
            mse_range = string(round(min_mse,10))+"-"+string(round(max_mse,10));
            %store values
            col(row_count) = mse_range;
            row_count = row_count + 1;
        end
        col = {col};
        data_dict("mult"+string(mult_count)+"dev"+string(dev_count)) = col;
        mult_count = mult_count+1;
    end
    dev_count = dev_count +1;
end
mult1dev1 = data_dict('mult1dev1');
mult1dev1 = mult1dev1{1};
mult2dev1 = data_dict('mult2dev1');
mult2dev1 = mult2dev1{1};
mult3dev1 = data_dict('mult3dev1');
mult3dev1 = mult3dev1{1};
mult1dev2 = data_dict('mult1dev2');
mult1dev2 = mult1dev2{1};
mult2dev2 = data_dict('mult2dev2');
mult2dev2 = mult2dev2{1};
mult3dev2 = data_dict('mult3dev2');
mult3dev2 = mult3dev2{1};
mult1dev3 = data_dict('mult1dev3');
mult1dev3 = mult1dev3{1};
mult2dev3 = data_dict('mult2dev3');
mult2dev3 = mult2dev3{1};
mult3dev3 = data_dict('mult3dev3');
mult3dev3 = mult3dev3{1};
row_names = {'nscale', 'minWaveLength', 'sigmaOnf', 'k', 'cutOff', 'g', 'noiseMethod'};
t = table(table(mult1dev1, mult2dev1, mult3dev1, 'VariableNames', {'Mult=1.5', 'Mult=2.1', 'Mult=2.5'}), ...
table(mult1dev2, mult2dev2, mult3dev2, 'VariableNames', {'Mult=1.5', 'Mult=2.1', 'Mult=2.5'}), ...
table(mult1dev3, mult2dev3, mult3dev3, 'VariableNames', {'Mult=1.5', 'Mult=2.1', 'Mult=2.5'}), ...
'VariableNames', {'DeviationGain=1', 'DeviationGain=1.5', 'DeviationGain=2'}, ...
'RowNames', {'nscale', 'minWaveLength', 'sigmaOnf', 'k', 'cutOff', 'g', 'noiseMethod'});
%save table
myFolder = images_path;
file = fullfile(myFolder, 'paramer_table.csv');
t2 = splitvars(t);
writetable(t2, file, 'WriteRowNames', true);


% Sharpening method experiment
mult =1.5;
deviationGain=1;
%test method 1
[fts1, pc_images1, Ts1, ors1, mse1, bw_images1] = phase_congruency_method1(images, binarys, mult, deviationGain);

%test method 2
[fts2, pc_images2, Ts2, ors2, mse2, bw_images2] = phase_congruency_method2(images, binarys, mult, deviationGain);

%control
[fts3, pc_images3, Ts3, ors3, mse3, bw_images3] = phase_congruency_control(images, binarys, mult, deviationGain);

%compare results:
T = table(mse1.', mse2.', mse3.', 'VariableNames', {'Method 1', 'Method 2', 'Control'});
format compact

summary(T)

"Control"
mean(mse3.')
"Method 1"
mean(mse1.')
"Method 2"
mean(mse2.')
"Control"
std(mse3.')
"Method 1"
std(mse1.')
"Method 2"
std(mse2.')

%example images:
input_path = path+"\juvenile\grey\010312-B8-1-2\org\frame0001.png";
binary_path = path+"\juvenile\binary\010312-B8-1-2\segmented\frame0001.png";
true_image = imread(binary_path);
image = imread(input_path);
images = {image};
figure()
[fts, pc_images, Ts, ors, mse, bw_images] = phase_congruency_method1(images, binarys, mult, deviationGain);
imshow(bw_images{1})
imwrite(bw_images{1}, path+"\output\method1.png")
figure()
[fts, pc_images, Ts, ors, mse, bw_images] = phase_congruency_method2(images, binarys, mult, deviationGain);
imshow(bw_images{1})
imwrite(bw_images{1}, path+"\output\method2.png")
figure()
[fts, pc_images, Ts, ors, mse, bw_images] = phase_congruency_control(images, binarys, 2.1, 1.5);
imshow(bw_images{1})
imwrite(bw_images{1}, path+"\output\control.png")
figure()
[fts, pc_images, Ts, ors, mse, bw_images] = phase_congruency_control(images, binarys, 2.5, 2);
imshow(bw_images{1})
imwrite(bw_images{1}, path+"\output\worst_param.png")
figure()
[fts, pc_images, Ts, ors, mse, bw_images] = phase_congruency_control(images, binarys, 1.5, 1);
imshow(bw_images{1})
imwrite(bw_images{1}, path+"\output\best_param.png")

%functions
function mse = phase_congruency_parameter_tuning_parasites(images, binarys, mult, deviationGain, param, param_list)
    num_images = length(images);
    ii=1;
    for i = 1:num_images
        image = rgb2gray(images{i});
        param = string(param(1));
        for param_val = param_list
            [PC, or, ft, T] = phasecongmono(image, 'mult', mult, 'deviationGain', deviationGain, char(param), double(param_val));
            nm = nonmaxsup(PC, or, 1.5);   % nonmaxima suppression
            bw = hysthresh(nm, 0.1, 0.3);
            mse(ii) = immse(double(binarys{i}-1), double(bw-1));
            ii = ii + 1; 
        end
    end
end

function mse = phase_congruency_parameter_tuning(images, mult, deviationGain, param, param_list)
    image = images{1};
    i = 1;
    param = string(param(1));
    for param_val = param_list
        [PC, or, ft, T] = phasecongmono(image, 'mult', mult, 'deviationGain', deviationGain, char(param), double(param_val));
        nm = nonmaxsup(PC, or, 1.5);   % nonmaxima suppression
        bw = hysthresh(nm, 0.1, 0.3);
        mse(i) = immse(double(image), double(bw-1));
        i = i + 1; 
    end
end

function [fts, pc_images, Ts, ors, mse, bw_images] = phase_congruency_method2(images, binarys, mult, deviationGain)
    num_images = length(images);
    for i = 1:num_images
        image = rgb2gray(images{i});
        image_sharp = imsharpen(image);
        [PC, or, ft, T] = phasecongmono(image_sharp, 'mult', mult, 'deviationGain', deviationGain);%, 4, 3, mult, 0.55, 3.0, -1, 0.5, 10, deviationGain);
        pc_images{i} = PC;
        fts{i} = ft;
        Ts{i} = T;
        ors{i} = or;
        nm = nonmaxsup(PC, or, 1.5);   % nonmaxima suppression
        bw = hysthresh(nm, 0.1, 0.3);
        bw_images{i} = bw;
        mse(i) = immse(double(binarys{i}-1), double(bw-1));
    end
end

function [fts, pc_images, Ts, ors, mse, bw_images] = phase_congruency_method1(images, binarys, mult, deviationGain)
    num_images = length(images);
    for i = 1:num_images
        image = rgb2gray(images{i});
        [PC, or, ft, T] = phasecongmono(image, 'mult', mult, 'deviationGain', deviationGain);%, 4, 3, mult, 0.55, 3.0, -1, 0.5, 10, deviationGain);
        pc_images{i} = PC;
        fts{i} = ft;
        Ts{i} = T;
        ors{i} = or;
        %sharpen after
        pc_sharp = imsharpen(PC);
        nm = nonmaxsup(pc_sharp, or, 1.5);   % nonmaxima suppression
        bw = hysthresh(nm, 0.1, 0.3);
        bw_images{i} = bw;
        mse(i) = immse(double(binarys{i}-1), double(bw-1));
    end
end

function [fts, pc_images, Ts, ors, mse, bw_images] = phase_congruency_control(images, binarys, mult, deviationGain)
    num_images = length(images);
    for i = 1:num_images
        image = rgb2gray(images{i});
        [PC, or, ft, T] = phasecongmono(image, 'mult', mult, 'deviationGain', deviationGain);%, 4, 3, mult, 0.55, 3.0, -1, 0.5, 10, deviationGain);
        pc_images{i} = PC;
        fts{i} = ft;
        Ts{i} = T;
        ors{i} = or;
        nm = nonmaxsup(PC, or, 1.5);   % nonmaxima suppression
        bw = hysthresh(nm, 0.1, 0.3);
        bw_images{i} = bw;
        mse(i) = immse(double(binarys{i}-1), double(bw-1));
    end
end

function images = read_images(images_path, exp)
    myFolder = images_path+"\exp_"+string(exp);
    filePattern = fullfile(myFolder, '*.png');
    imagefiles = dir(filePattern);  
    nfiles = length(imagefiles);    % Number of files found
    for ii=1:nfiles
       baseFileName = imagefiles(ii).name;
       fullFileName = fullfile(myFolder, baseFileName);
       fprintf(1, 'Now reading %s\n', fullFileName);
       currentimage = imread(fullFileName);
       images{ii} = currentimage;
    end
end

function [images, binary] = read_parasite_images(images_path, treatment)
    images_folder = images_path+"\grey\"+string(treatment)+"\org";
    binary_folder = images_path+"\binary\"+string(treatment)+"\segmented";
    folder = fullfile(images_folder, '*.png');
    image_files = dir(folder);
    bfolder = fullfile(binary_folder, '*.png');
    binary_files = dir(bfolder);
    nfiles = length(image_files);
    for ii=1:nfiles
       baseFileName = image_files(ii).name;
       fullFileName = fullfile(images_folder, baseFileName);
       fprintf(1, 'Now reading %s\n', fullFileName);
       currentimage = imread(fullFileName);
       images{ii} = currentimage;
       base_binary_file = binary_files(ii).name;
       full_binary_file = fullfile(binary_folder, base_binary_file);
       fprintf(1, 'Now reading %s\n', full_binary_file);
       currentimage = imread(full_binary_file);
       binary{ii} = currentimage;
    end

end

function animate(images, filename)
    nImages = length(images);
    for idx = 1:nImages
        I = images{idx};
        %Check if I is Grayscale
        if (size(I, 3) == 1)
            %Convert Grayscale to RGB by duplicating I three times.
            I = cat(3, I, I, I);
        end
        [A,map] = rgb2ind(I,256);
        if idx == 1
            imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",1);
        else
            imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",1);
        end
    end
end

