% Get the ratio of oscillatory points for every participant x voxel x
% frequency averaging across time.

% sBOSC path
p.sBOSC = 'Z:\Toolbox\sourceBOSC';
addpath(genpath(p.sBOSC))

% Fieldtrip path
addpath('Z:\Toolbox\fieldtrip-20230118')  % Fieldtrip path
ft_defaults

% OMEGA data path
p.data = 'Z:\OMEGA\OMEGA_data';
p.raw = 'Z:\OMEGA\OMEGA_raw';

% Participant and session list
sub = [1 2 3 4 5 6 7 8 9 11 12 14 15 16 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 37 39 40 41 42 44 45 46 47 48 49 50 51 52 55 56 57 58 59 60 61 62 63 64 65 67 68 69 70 71 72 73 74 75 76 77 78 79 80 84 85 87 88 89 90 91 92 94 95 96 97 98 99 101 102 103 104 105 106 134 145 146 148 149 150 151 152 154 155 156 157 158 159 160 161 165 166 167 168 169 170 171 175 176 177 179 181 184 185 195 197 200 207 208 210 212]';
ses = [1 1 1 1 1 1 1 2 1  2  1  2  1  1  1  2  3  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  3  1  2  1  2  1  2  1  1  2  1  1  1  1  1  1  1  1  2  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 ]';

subs = {} ; sess = {};
for s = 1:length(sub)
    if sub(s) < 10
        subs{s} = ['000' num2str(sub(s))];
    elseif sub(s) >= 10 & sub(s) < 100
        subs{s} = ['00' num2str(sub(s))];
    else
        subs{s} = ['0' num2str(sub(s))];
    end
    sess{s} = ['000' num2str(ses(s))];
end

Nsub  = length(sub);

srate = 128;

load source_template_10mm_1925.mat

% Frex parameters 
frex = exp(0.6:0.1:3.7); % 1.8 Hz to 40 Hz
Nfrex = length(frex);

freqbands = [frex(1) frex(9) frex(16) frex(21) frex(25) frex(30)];
freqnames = {'delta', 'theta', 'alpha', 'lowbeta', 'highbeta'};

Nvoxin = length(find(source.inverse.inside));
Ntp = 33270;

total_osc_time = zeros(Nsub,Nvoxin); % in all freqs
anyvoxosc = zeros(Nsub,1); % If any voxel and freq is oscillating
anyfreqosc = zeros(Nsub,Nvoxin); % If any freq and for each voxel
fb_osc = zeros(Nsub,Nvoxin, 5); % for each voxel and freqband
freq_osc = zeros(Nsub,Nvoxin, length(frex)-2); % for each freq
overlap_subjects =  zeros(Nsub,Nvoxin,length(freqnames),length(freqnames));

%%
for s=1:Nsub
    s

%% Transform  epis from cells into a matrix

    episodes = false(Nvoxin,Nfrex,Ntp);
    for vx = 1:Nvoxin
            load([p.data '\sub-' subs{s} '\ses-' sess{s} '\episodes\epis_v' num2str(vx) '.mat'])

        for ep = 1:length(epis)
            fm = epis(ep).freq;
            fbin = dsearchn(frex',fm);
            tpts = epis(ep).timeps;
            episodes(vx,fbin(1),tpts) = 1;
        end
    end
    episodes = episodes(:,1:30,:);

%% Proportion of oscillatory time 

        % Any voxel and freq osclillating
        tmp = squeeze(sum(sum(episodes,2))>0);
        anyvoxosc(s) = sum(tmp) ./ size(episodes,3) * 100;

        % Any freq oscillating
        tmp = sum(squeeze(sum(episodes,2)>0),2); 
        anyfreqosc(s,:) = tmp ./ size(episodes,3) * 100;  

        % Oscillation for each bandfreq
        ct = 1;
             for fb = 1:length(freqnames)
                 if fb==length(freqnames) % last frequency
                     ct=0;
                 end
                 fband = [freqbands(fb) freqbands(fb+1)];
                 fbidx = dsearchn(frex',fband');
                 tmp = sum(squeeze(sum(episodes(:,fbidx(1):fbidx(2)-ct,:),2))>0,2);
                 fb_osc(s,:,fb) = tmp ./ size(episodes,3) * 100;
             end

         % Oscillation for each freq
        freq_osc(s,:,:) = sum(episodes,3) ./  size(episodes,3) * 100; % 

        % Overlaps
        overlap_voxels = zeros(Nvoxin, length(freqnames), length(freqnames));
        ct = 1;
        ct2 = 1;
        for fb = 1:length(freqnames)
            if fb==length(freqnames) % last frequency
                ct=0;
            end
            fband = [freqbands(fb) freqbands(fb+1)];
            fbidx = dsearchn(frex',fband');
            freq_node = squeeze(sum(episodes(:,fbidx(1):fbidx(2)-ct,:),2)>0);

            for fb2 = 1:length(freqnames)
                if fb2==length(freqnames) % last frequency
                    ct2=0;
                end
                fband2 = [freqbands(fb2) freqbands(fb2+1)];
                fbidx2 = dsearchn(frex',fband2');
                freq_overlap = squeeze(sum(episodes(:,fbidx2(1):fbidx2(2)-ct2,:),2)>0);

                % Calculate overlap
                overlap = freq_node .* freq_overlap;
                overlap_voxels(:,fb,fb2) = sum(overlap,2) ./ size(overlap,2) * 100;
            end
        end
        overlap_subjects(s,:,:,:) = overlap_voxels;


end

save([p.data '\oscillatory_results'], 'anyvoxosc', 'anyfreqosc', 'fb_osc', 'freq_osc', 'overlap_subjects')

