%experiment set up for Phase II Project:
path = fileparts(mfilename('fullpath'));
images_path = path + "\edge detection data"; %replace with correct path
exp_nums = [1, 2, 3, 4, 5]; %define the number of experiments
mult_vals = [1.5, 2.1, 2.5]; %define the mult values to tru
deviationGain_vals = [1, 1.5, 2]; %define the deviation gain values to try
num_exp = length(exp_nums);
disp(num_exp) %check the number of experiments
for exp_num = 1:num_exp %for each exp
    for mult = mult_vals %for each value
        for deviationGain = deviationGain_vals %for each value
            run_exp(images_path, exp_num, mult, deviationGain, true) %run experiment
        end
    end
end

%functions:

function run_exp(images_path, exp_num, mult, deviationGain, save_gif)
    %read in images
    images = read_images(images_path, exp_num);
    %run experiment, get output images, phase congrunecy, average phase
    [fts, pc_images, Ts, ors] = phase_congruency(images, mult, deviationGain);
    %animate output images
    myFolder = images_path+"\exp_"+string(exp_num);
    if save_gif == true
        output_filename = fullfile(myFolder, string(mult)+"_"+string(deviationGain)+"PC_output_test.gif");
    end
    animate(pc_images, output_filename);
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

function [fts, pc_images, Ts, ors] = phase_congruency(images, mult, deviationGain)
    num_images = length(images);
    for i = 1:num_images
        image = images{i};
        [PC, or, ft, T] = phasecongmono(image, 'mult', mult, 'deviationGain', deviationGain);
        pc_images{i} = PC;
        fts{i} = ft;
        Ts{i} = T;
        ors{i} = or;
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
