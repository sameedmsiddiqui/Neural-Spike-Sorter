%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sameed Siddiqui
% EMG decomposition Project
% BE 260, Fall 2015, UCLA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Instructions:
%    This code is devided into many sections of operation.
%    Each section has its own specific instructions and/or options. 
%
%    Thank you!

%% Input file
% This is where the user inputs his or her file. Please note that your file
% must include in column 1 the time series and in columns 2 through n+1 the
% signal from n channels
%
% Please note that every subsection of this code has its own options and
% parameters as well! I did not put all of these parameters and options in
% this input file because I felt that it would be easier for you to test
% each section if each section's options were with that section. Where you
% are given the option to chose between various techniques, the code
% currently labels the default settings. These are those techniques which
% seem to work the best and the most easily, but this may vary depending on
% the signal.
%
% Please select the channel of your input file you'd like to use:
channel_select= 1; % channel_select<= channel_number

M= csvread('EMG_example_2_fs_2k.csv'); %read in csv file
time= M(:,1); % first column is the time series
fs= (time(2)-time(1))^-1; % calculate the sample frequecy
channel_number= size(M,2)-1; % num of channels in the database
for i=1:channel_number,
    figure('Color',[1 1 1]);plot(time,M(:,i+1)); %plot each channel
    str= sprintf('Channel %d',i);
    xlabel('seconds');title(str);xlim([time(1) time(size(time,1))]); % label and title each plots
end
original_data= M(:,channel_select+1); % test_input will go through all the individual sections

%% Filter data
% 
%  By default, this system uses a high-pass filter with a cut-off of
%  200Hz.This was a "good" result from "emg-data-segment" file.
%
%  Options:
%    (1) For more sophisticated filtering, set smart_filtering = 1 and 
%        input the time of the "pure noise" and the time where your signal 
%        is present. This filtering technique compares the frequency 
%        content of each respective signal and determines what frequencies 
%        should be filtered. It is useful when your neural signal is pulsed 
%        and repetitive in nature. 
%        This technique uses more computational power. It is also
%        somewhat experimental - I came up with it myself to make the
%        project as thorough as possible! However, initial test results are
%        exceptionally positive. 
%        Finally, please make sure that the input ranges for noise and
%        pure signal have at least 300 indices. This is to make sure that
%        we are using enough data to actually come up with an accurate
%        estimate of frequency content of the signals.

sampling_rate = fs;
smart_filtering = 1; % 0 for high pass, 1 (default) for experimental smart filtering. 
% The below values are added for your convenience, and correspond to Channel 1 of the EMG_example_2_fs_2k file.
time_pure_noise_start = 0.4405;
time_pure_noise_end = 0.7435;
time_pure_signal_start = 0.802;
time_pure_signal_end = 1.285;


if smart_filtering == 0
    fmax = sampling_rate/2-1; 
    fmin= 200; % Even at higher frequencies, you probably shouldn't change this!

    % create an FIR filter
    bpfilt = designfilt('bandpassfir','FilterOrder',20, ...
             'CutoffFrequency1',fmin,'CutoffFrequency2',fmax, ...
             'SampleRate',sampling_rate);

    % use the filter
    filtered_data = filtfilt(bpfilt,original_data);
    
   
elseif smart_filtering == 1
    % Find the indices of your signal entries.
    [~, index_pure_noise_start] = min(abs(time - time_pure_noise_start));
    [~, index_pure_noise_end] = min(abs(time - time_pure_noise_end));
    [~, index_pure_signal_start] = min(abs(time - time_pure_signal_start));
    [~, index_pure_signal_end] = min(abs(time - time_pure_signal_end));
    
    % Check for any errors from the user's input: 
    if (isempty(index_pure_noise_start) || isempty(index_pure_noise_end) ||...
            isempty(index_pure_signal_start) || isempty(index_pure_signal_end) ||...
            index_pure_noise_start > index_pure_noise_end + 300 || ...
            index_pure_signal_start > index_pure_signal_end + 300)
        warning('There was some error in your smart_filtering parameters. No filtering used. A common mistake is that the time values for signals need to be exact')
    end
    noise_signal = original_data(index_pure_noise_start:index_pure_noise_end);
    pure_signal = original_data(index_pure_signal_start:index_pure_signal_end);
    
    if length(noise_signal) > length(pure_signal)
        % We need to interpolate the signals to make sure they have the same #
        % of indices. This way, we can compare their FFT's
        temp_signal_time_axis = linspace(time(index_pure_signal_start), time(index_pure_signal_end), length(noise_signal))';
        pure_signal = interp1(time(index_pure_signal_start:index_pure_signal_end), pure_signal, temp_signal_time_axis, 'spline');
        
    else % if length(pure_signal) >length(noise_signal), we have to do our interpolation in the other way.
        temp_noise_time_axis = linspace(time(index_pure_noise_start), time(index_pure_noise_end), length(pure_signal))';
        noise_signal = interp1(time(index_pure_noise_start:index_pure_noise_end), noise_signal, temp_noise_time_axis, 'spline'); 
    end
    
    fft_signal = fft(pure_signal);
    fft_noise = fft(noise_signal);

    dF = sampling_rate/length(noise_signal);
    frequency_axis = -sampling_rate/2:dF:sampling_rate/2-dF;
    
    % Let's FFTshift and get only the positive frequency components. 
    frequency_axis = fftshift(frequency_axis);
    frequency_axis = frequency_axis(1:floor(end/2));
    fft_signal = fftshift(fft_signal);
    fft_signal = fft_signal(1:floor(end/2));
    fft_noise = fftshift(fft_noise);
    fft_noise = fft_noise(1:floor(end/2));
    
    figure; hold on; plot(frequency_axis, 10*log10(abs(fft_signal))); plot(frequency_axis, 10*log10(abs(fft_noise)), 'color', 'r');
    title('FFT of noise versus FFT of signal')
    xlabel('frequency (Hz)')
    ylabel('10log(abs(signal))')
    legend('signal', 'noise')
    
    % Analyzing the signals you've provided, it seems that the "optimal"
    % cut-off frequency seems to vary somewhere betweeen 50Hz and 300Hz.
    % Thus by using the ffts of the signal, we can then determine where the
    % cut-off should be on a signal-by-signal basis.
    %
    % Future versions of this program would be even more complex in
    % filtering, wherein they would have selective filtering (e.g. 0-50Hz
    % could be filtered, 50-100Hz might not be filtered, 100-200Hz would be
    % filtered, etc., depending on where noise is high and where signal is
    % high).
    
    % Key: in 50Hz incremets, we'll see if the noise signal is higher
    % than the data signal. If it is, then we'll create a highpass filter
    % that get rids of that. 
    fmin = 50;
    cutoff_frequency = 50;
    for i = 100:50:300
        current_index =  find(frequency_axis > i, 1);
        previous_index = find(frequency_axis> i-50, 1);
        if mean(10*log10(abs(fft_noise(previous_index:current_index))))-...
                mean(10*log10(abs(fft_signal(previous_index:current_index)))) > 0 ...
                % If you change the 0 above, you can change the criteria for 
                % cutoff frequency (e.g. how the allowable maximum difference 
                % between noise and signal before the filter cuts off.
            cutoff_frequency = i;
        else
            break
        end
    end
    
    fmax = sampling_rate/2 -1; 
    bpfilt = designfilt('bandpassfir','FilterOrder',20, ...
         'CutoffFrequency1',cutoff_frequency,'CutoffFrequency2',fmax, ...
         'SampleRate',sampling_rate);
     
    % use the filter
    filtered_data = filtfilt(bpfilt,original_data);
    
else % if smart_filtering parameter was inputted incorrectly.
    warning('The smart_filtering paramater was neither 0 nor 1. No filtering used.')
end % end smart filtering!

figure; plot(time, original_data); hold on; plot(time, filtered_data, 'r');
title('Section 2: data filtereing');
legend('Original data', 'Filtered data');
ylabel('Signal amplitude')
xlabel('Time (s)'); xlim([time(1), time(end)*1.05]);


% Renaming some variables for ease-of-use throughout the program.
data = filtered_data;

%% Spike Detection
%  Instructions/Options
%    (1) Depending on sampling size, you should set the bin size of each
%        spike. For EMG_example_2_fs_2k channel 1, this bin-size has been
%        preset to 8. 
%    (2) Chose spike separating method: by approximate energy or by
%    findpeaks. For approximate energy, set energy_method = 1. Else for
%    findpeaks, set energy_method = 0; By default, we use the findpeaks
%    method. 

threshold = 0.3;
bin_size = 8;
energy_method = 0; % default = 0 -> uses findpeaks instead


% In case the user forgets to add a threshold value, this formula presented
% by Quiroga is an acceptable place to start. 
if threshold == 0
    threshold = 4*std(data);
    fprintf('Looks like you have not entered a threshold. One was autogenerated')    
end

% Now, there are two potential methods to finding individual spikes: (1)
% findpeaks and (2) Find the center point of a continuous region with
% spikes above the threshold. 

if energy_method == 1
    % We only look for positive threshold values because or signals tend to
    % be highly symmetric. If we include negative threshold, we would
    % likely double-count signals. 
    v_thresholded = find(data > threshold);
    
    % Visualize what points have made it through our thresholding.
    figure; hold on; plot(time,data); plot(time(v_thresholded), data(v_thresholded), 'color', 'r');
    xlabel('time (s)'); ylabel('Amplitude');
    title('Neural signals vs. background'); xlim([time(1), time(end)*1.05]);
    legend('Signal', 'Neural signal above threshold minimum')

    beginning_of_region = v_thresholded(1);
    end_of_region = [];
    for i = 1:length(v_thresholded)-1
        % If the threshold criteria isn't met for three contiguous indices
        % in row, we can conclude that we're at the end of a particular 
        % signal region.
        if v_thresholded(i) < v_thresholded(i+1) - 3
            end_of_region = [end_of_region, v_thresholded(i)];
            if i ~= length(v_thresholded)-1
                beginning_of_region = [beginning_of_region, v_thresholded(i+1)];
            end
        elseif i == length(v_thresholded)-1
            end_of_region = [end_of_region, v_thresholded(i+1)];
        end
    end
    
% FindPeaks method
elseif energy_method == 0
    
    % We want the peaks to have both minimum height of
    % threshold. Note that this method only looks for positive peaks. It
    % also seeks to minimize overlapping signals by restricting the space
    % between detected signals.
    [peak_height, peak_idx] = findpeaks(data, 'MinPeakHeight', threshold, 'MinPeakDistance', bin_size/2);
    figure; plot(time, data); hold on; scatter(time(peak_idx), peak_height);
    title('Peaks detected using FindPeaks'); xlim([time(1), time(end)*1.05]);
    xlabel('time (s)'); ylabel('Amplitude');
    
else
    error('Error! Your "energy_method" paratmeter is wrong! Must be 0 or 1')
end



%% Spike alignment
% This section takes the previous section's results and uses them to create 
% matrices of aligned neural spikes. No new input is required from the user
% - the user's previous input of spike detection method carries through to
% this section. 


if energy_method == 1
    % From spike detection, we have coordinate points for when each neural
    % signal is above the threshold value. Thus, the center point between
    % the first point above the threshold and the last point is
    % approximately the center of energy of a particular spike. Thus we can
    % align the spikes with respect to their approximate energy centers. 
    region_center_index = (end_of_region + beginning_of_region)/2;

    % This actually separates the individual signals and makes sure they
    % are aligned.
    regions = [];
    accepted_regions_index = [];
    for i = 1:length(region_center_index)
        if rem(region_center_index(i),1) ~= 0
            region_center_index(i) = region_center_index(i) +0.5;
        end
        if region_center_index(i) - bin_size/2 > 0 && region_center_index(i) + bin_size/2 + 1 < length(data)
            accepted_regions_index = [accepted_regions_index, region_center_index(i)];
            regions = cat(1, regions,( data(region_center_index(i)-bin_size/2:region_center_index(i)+bin_size/2-1))');
        end 
    end

elseif energy_method == 0 % findpeaks method
    regions = [];
    accepted_regions_index = [];
    % Separate the individual signals, making sure they are aligned:
    for i = 1:length(peak_idx)
        if peak_idx(i) - bin_size/2 > 0 && peak_idx(i) + bin_size/2 + 1 < length(data)
            accepted_regions_index = [accepted_regions_index, peak_idx(i)];
            regions = cat(1, regions, (data(peak_idx(i)-bin_size/2:peak_idx(i)+bin_size/2-1))');
        end 
    end
end


% Test to see if our spikes were picked out correctly. 
figure; hold on;
for i = 1: length(regions)
    plot(regions(i,:))
end
title('Detected and Aligned Spikes')
xlabel('index')
ylabel('Amplitude')


%% Feature Extraction
% This section uses PCA for feature extraction. By default, the program
% takes the first four PCA components. The user does not need to input any
% instructions to this section of the code.

[coeff, score, latent] = pca(regions);

% Show that PCA has "saved" us some variance.
figure; plot(100*cumsum(latent)./sum(latent));
title('Percent of total data variance maintained with n PCA components')
xlabel('Number of PCA components used')
ylabel('% of total variance contained within PCA components')


reduced_regions = coeff(1:4, :)*regions';

% Show our PCA reduction results, plotted in 3D with the color
% representing the 4th dimension
figure; hold on;
map = jet(100);
% To use the color map as our fourth dimension, we must normalize the 4th
% PCA dimension:
normalized_4d = reduced_regions(4,:) - min(reduced_regions(4,:));
norm_factor = 99/max(normalized_4d);
normalized_4d = normalized_4d*norm_factor;
scatter3(reduced_regions(1,:), reduced_regions(2,:), reduced_regions(3,:), 20, map(1+floor(normalized_4d)))
title('PCA-based compression of spike signals')
xlabel('PCA component 1')
ylabel('PCA component 2')
zlabel('PCA Component 3')
c = colorbar;
ylabel(c,'PCA component 4');

%% Clustering
% Instructions/Options:
%   Like most of our sections, the user has an option of how to go about
%   doing things
%   (1) Simply enter the number of clusters desired from k-means cluster
%       analysis (supervised sorting).
%   (2) Allow the system to repeatedly run k-means until it reaches a
%       number of clusters such that the smallest cluster contains 
%       minimum_percent of the total neural spikes. By default, this
%       minimum percent is equal to 7.75% - this seems to be reasonable. If
%       desired, the user can change this minimum percent value as well.
%       (semi-unsupervised). Set semi_unsupervised = 1 for this iterative
%       functionality. Set = 0 for fully supervised operation. 
%
%       Note that k-means cluster visualization method we use can only plot 
%       up to the first 26 spikes. However, so many spikes is unrealistic 
%       anyways, so we believe this is very acceptable. 

number_of_clusters = 4; % For supervised sorting
semi_unsupervised = 1; % Set 0 for supervised, 1 (default) for semi-unsupervised
minimum_percent = 7.75; % For semi-unsupervised sorting.
nspikes = size(reduced_regions, 2);

if semi_unsupervised == 0
    idx = kmeans(reduced_regions', number_of_clusters);
    
elseif semi_unsupervised == 1
    above_min_percent = true;
    
    number_of_clusters = 0;
    previous_idx = zeros(1, size(reduced_regions, 2)); % This will help us reduce computation time. 
    idx = 0;
    while(above_min_percent)
        number_of_clusters = number_of_clusters + 1;
        previous_idx = idx;
        idx = kmeans(reduced_regions', number_of_clusters);
        cluster_occurences = zeros(1, max(idx));
        
        for i = 1: max(idx)
            % Count how many spikes fell into each respective cluster.
            for j = 1: nspikes
                if idx(j) == i
                    cluster_occurences(i) =  cluster_occurences(i)+1;
                end
            end
            % Check if the idx we just counted meets the minimum percentage
            % criteria. 
            if cluster_occurences(i)/nspikes < minimum_percent/100
                above_min_percent = false;
                break;
            end
        end
    end
    
    % Because the previous loops ends when we reach too many clusters, that
    % means the actual number of clusters we should have is one less than
    % our previous calculation.
    number_of_clusters = number_of_clusters - 1;
    idx = previous_idx;
    
else % if semi_unsupervised wasn't 1 or 0
    warning('Your parameter entry for "semi_unsupervised" is incorrect. It should be either 1 or 0');
end % end kmeans indexing

% This code is a sanity check for me. Please ignore this professor Liu. 
if max(idx) ~= number_of_clusters
   error('Sameed, something went wrong!') 
end


% Plot k-means clusters. As before, we plot our 4-d PCA components with
% color determining the 4th dimension. The shape of the marker represents
% the cluster. Note that if there are more than 13 clusters, then there may
% be two clusters assigned to one marker due to the fact that Matlab
% supports only 13 markers.
marker_types = ['+', 'o', 's', 'x', '>', 'v', '<', '^', 'p', 'd', '*', '.', 'h'];
figure; hold on;
map = jet(100);

for i = 1:number_of_clusters
    if i > 26
        warning('Above 26 clusters, our visualization cuts off! However, this does not otherwise effect the program')
        break;
    end
    scatter3(reduced_regions(1,idx == i), reduced_regions(2,idx == i), reduced_regions(3,idx == i), 20, map(1+floor(normalized_4d(idx == i))),  marker_types(floor(i/13)+mod(i,13)))
end
    title(sprintf('K-means clustering results: each marker shape is one cluster\n %d clusters total', number_of_clusters))
xlabel('PCA component 1')
ylabel('PCA component 2')
zlabel('PCA Component 3')
c = colorbar;
ylabel(c,'PCA component 4');



%% Classification
% No instructions, just outputs!

% Count how many spikes fell into each respective cluster.
cluster_occurences = zeros(1, max(idx));
for i = 1: max(idx) % e.g. for the number of different clusters
    for j = 1: nspikes
        if idx(j) == i
            cluster_occurences(i) =  cluster_occurences(i)+1;
        end
    end
end

% Create cluster templates
templates = zeros(number_of_clusters, bin_size);
for i = 1:number_of_clusters
    template_match = find(idx == i); 
    templates(i,:) = median(regions(template_match',:));
end

% Visualize the templates
figure; hold on;
map = jet(number_of_clusters);
for i = 1:number_of_clusters
    plot(1/sampling_rate*(1:bin_size), templates(i,:), 'Color', map(i,:))
end
ylabel('Amplitude');
xlabel('time');
title('Median templates');
% Creating a good legend is tricky!: 
legend_string = cellstr(sprintf('Template 1, %d matches', cluster_occurences(1)));
for i = 2:number_of_clusters
   legend_string = cat(1, legend_string, sprintf('Template %d, %d matches', i, cluster_occurences(i))); 
end
legend(legend_string, 'location', 'best');



% Plot all the waveforms on top of each other.
figure; hold on;
h = zeros(1, number_of_clusters); % h array is needed so that the legend comes out okay.
for i = 1:number_of_clusters
    template_match = find(idx == i); 
    for j = 1:length(template_match)
        h(i) = plot(1/sampling_rate*(1:bin_size), regions(template_match(j),:), 'Color', map(i, :));
    end
end
title('All waveforms plotted, labeled with template identity')
ylabel('Amplitude')
xlabel('time')
legend(h, legend_string, 'location', 'best')

% Overlay all of the different spikes with each other. 
for i = 1:number_of_clusters
    template_match = find(idx == i);
    figure; hold on; title(sprintf('Cluster %d, with %d occurences', i, cluster_occurences(i)));
    ylabel('Amplitude')
    xlabel('time')
    for j = 1:length(template_match)
        plot(1/sampling_rate*(1:bin_size), regions(template_match(j),:))
    end
end

%% Analysis
% Instructions:
%   (1) it's nice for one to be able to know when in the signal there
%       were a lot of spikes for a given cluster. Thus, we split the total 
%       signal into some number of "bins" and determine how many spikes 
%       are in which bins. By default, this is set to 100 bins, but you can
%       control that.

number_of_bins = 100;

% Create a bar graph to show the relative amounts of the templates
figure; bar(cluster_occurences); title('Instances of each template');
xlabel('Template number')
ylabel('Number of instances')

% Label on the original data series each peak's corresponding template
figure; plot(time, data, 'color', 'black'); hold on;
h = zeros(1, number_of_clusters); % h array is needed so that the legend comes out okay.
for i = 1:number_of_clusters
    indices_in_this_cluster = accepted_regions_index(idx == i);
    h(i) = scatter(time(indices_in_this_cluster), data(indices_in_this_cluster), [],  map(i,:), 'filled');
end
legend(h, legend_string, 'location', 'best');
xlabel('time (s)'); ylabel('Amplitude'); xlim([0, time(end)*1.05]);
title('Original data superimposed with cluster identities')


% Same as above except with one cluster at a time. This is done because it
% might be simpler for someone to visualize a large number of spikes this
% way. 
for i = 1: number_of_clusters
    indices_in_this_cluster = accepted_regions_index(idx == i);
    figure; plot(time, data); hold on; scatter(time(indices_in_this_cluster), data(indices_in_this_cluster));
    xlabel('time (s)'); ylabel('Amplitude'); xlim([0, time(end)*1.05]);
    title(sprintf('Original data superimposed with cluster %d', i));
end

% Finally, it's nice for one to be able to know when in the signal there
% were a lot of spikes for a given cluster. Thus, we split the total signal
% into 100 "bins" and determine how many spikes were in which bins.
bin_ends = linspace(length(data)/number_of_bins, length(data), number_of_bins);
spikes_in_bin = zeros(number_of_clusters, number_of_bins);
for i = 1: number_of_clusters
    indices_in_this_cluster = accepted_regions_index(idx == i);
    spikes_in_bin(i, 1) = length(find(indices_in_this_cluster <= bin_ends(1)));
    for j = 2:number_of_bins
        spikes_in_bin(i, j) = length(find(indices_in_this_cluster <= bin_ends(j)...
            & indices_in_this_cluster > bin_ends(j-1)));
    end
end

figure; bar(spikes_in_bin');
title(sprintf('Instances of each template within intervals of \n 1/%dth of total signal time', number_of_bins));
ylabel('Number of instances')
xlabel(sprintf('Time interval: \n 0 is start of signal, %d is end', number_of_bins))
legend(legend_string, 'location', 'best')
xlim([-5,number_of_bins+5])

%%%%%%%% END PROGRAM %%%%%%%%
    
    
    
