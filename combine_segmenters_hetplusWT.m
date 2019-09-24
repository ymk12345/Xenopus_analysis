%% Combine output of two segmenters: male/female/overlap/nosong and female/nosong
function combine_segmenters_hetplusWT(fname_song,fname_4c,fname_2c)
tic
% clear all
% close all
%% Training data set
% mo_fname='/tigress/christab/virilis_segmenter_cnn/segmented_data/4001_trainingplusnoise_NNsegmented.mat';
% fn_fname='/tigress/christab/virilis_segmenter_cnn/segmented_data/4001_femaleNoise_talmojan_trainingplusnoise_NNsegmented.mat';
% fname_tosave='/tigress/christab//virilis_segmenter_cnn/combined_data/4001_trainingplusnoise_NNsegmented_combined.mat';

% mo_fname='T:\virilis_segmenter_cnn\predicted_data\4001_trainingplusnoise_outp_envelopes.mat';
% fn_fname='T:\virilis_segmenter_cnn\predicted_data\4001_femaleNoise_talmojan_trainingplusnoise_outp_envelopes.mat';
% fname_tosave='T:\virilis_segmenter_cnn\combined_data\4001_femaleNoise_talmojan_trainingplusnoise_outp_envelopes.mat';

%% WT files
% mo_fname='T:/virilis_segmenter_cnn/predicted_data/4001_talmojan_160910_0913ch30_outp_envelopes.mat';
% fn_fname='T:/virilis_segmenter_cnn/predicted_data/2001_femaleNoise_talmojan_160910_0913ch30_outp_envelopes.mat';
% fname_tosave='T:/virilis_segmenter_cnn/combined_data/160910_0913ch30_2segs40012001_fmindur_segmented.mat';

% mo_fname='T:/virilis_segmenter_cnn/predicted_data/4001_talmojan_160913_1015ch30_outp_envelopes.mat';
% fn_fname='T:/virilis_segmenter_cnn/predicted_data/2001_femaleNoise_talmojan_160913_1015ch30_outp_envelopes.mat';
% fname_tosave='T:/virilis_segmenter_cnn/combined_data/160913_1015ch30_2segs40012001_fmindur_segmented.mat';
% 
% mo_fname='T:/virilis_segmenter_cnn/predicted_data/4001_talmojan_160920_0912ch26_outp_envelopes.mat';
% fn_fname='T:/virilis_segmenter_cnn/predicted_data/2001_femaleNoise_talmojan_160920_0912ch26_outp_envelopes.mat';
% fname_tosave='T:/virilis_segmenter_cnn/combined_data/160920_0912ch26_2segs40012001_fmindur_segmented.mat';

% mo_fname='T:/virilis_segmenter_cnn/predicted_data/4001_talmojan_160920_0912ch32_outp_envelopes.mat';
% fn_fname='T:/virilis_segmenter_cnn/predicted_data/2001_femaleNoise_talmojan_160920_0912ch32_outp_envelopes.mat';
% fname_tosave='T:/virilis_segmenter_cnn/combined_data/160920_0912ch32_2segs40012001_fmindur_segmented.mat';

% mo_fname='T:/virilis_segmenter_cnn/predicted_data/4001_talmojan_160920_1017ch28_outp_envelopes.mat';
% fn_fname='T:/virilis_segmenter_cnn/predicted_data/2001_femaleNoise_talmojan_160920_1017ch28_outp_envelopes.mat';
% fname_tosave='T:/virilis_segmenter_cnn/combined_data/160920_1017ch28_2segs40012001_fmindur_segmented.mat';
%% Load data (tigress)
disp(fname_song)
disp('Loading data.')
song_path='/tigress/christab/virilis_segmenter_cnn/data_to_segment/';
predicted_path='/tigress/christab/virilis_segmenter_cnn/predicted_data/';
combined_path='/tigress/christab/virilis_segmenter_cnn/combined_data/';
% % mo_fname=strcat(segmented_path,maleOverlap_fname);
% fn_fname=strcat(segmented_path,femaleNoise_fname);
% mo_fname=maleOverlap_fname;
% fn_fname=femaleNoise_fname;
path_tosave='/tigress/christab/virilis_segmenter_cnn/combined_data/';
fname_tosave=strcat(path_tosave,fname_song(1:end-4),'_combined_hetplusWT.mat');
%% Load data (local)
% disp(fname_song)
% disp('Loading data.')
% song_path='T:/virilis_segmenter_cnn/data_to_segment/';
% predicted_path='T:/virilis_segmenter_cnn/predicted_data/';
% combined_path='T:/virilis_segmenter_cnn/combined_data/';
% % % mo_fname=strcat(segmented_path,maleOverlap_fname);
% % fn_fname=strcat(segmented_path,femaleNoise_fname);
% % mo_fname=maleOverlap_fname;
% % fn_fname=femaleNoise_fname;
% path_tosave='T:/virilis_segmenter_cnn/combined_data/';
% fname_tosave=strcat(path_tosave,fname_song(1:end-4),'_combined_hetplusWT.mat');
%% Load data
disp('Loading data.')
temp_dir=cd('../utility/');
% 
out_p_4c = readNPY(strcat(predicted_path,fname_4c));
out_p_2c = readNPY(strcat(predicted_path,fname_2c));
d = load(strcat(song_path,fname_song),'data');
b = load(strcat(song_path,fname_song),'begNoiseEndSec');

data = d.data;
begNoiseEndSec=b.begNoiseEndSec;

L_4c=4001;
L_2c=2001;

out_p_4c=out_p_4c(((L_4c-1)/2+1):end-((L_4c-1)/2),:);
out_p_2c=out_p_2c(((L_2c-1)/2+1):end-((L_2c-1)/2),:);
% 
% dmo=load(strcat(mo_fname,'data');
% dfn=load(fn_fname,'data');
% data_mo=dmo.data;
% data_fn=dfn.data;
% if (sum(data_mo==data_fn)<length(data_mo))
%     error('Data do not match!')
% else
%     data=data_mo;
% end
% % out_p_filt_mp=load(strcat(segmented_path,mo_fname),'maleBoutInfo');
% mo=load(mo_fname,'out_p_filt');
% out_p_filt_mo=mo.out_p_filt;
% fn=load(fn_fname,'out_p_filt');
% out_p_filt_fn=fn.out_p_filt;
% 
% le=load(mo_fname,'lower_envelope');
% ue=load(mo_fname,'upper_envelope');
% lower_envelope=le.lower_envelope;
% upper_envelope=ue.upper_envelope;
% 
% clear mo fn le ue data_mp data_fn
%% Filter output probabilities
disp('Filtering and averaging output probabilities.')
for i=1:4
   out_p_4c_filt(:,i)=medfilt1(out_p_4c(:,i),100);  %median filter output probabilities over 10 ms window
end

for i=1:2
   out_p_2c_filt(:,i)=medfilt1(out_p_2c(:,i),100);  %median filter output probabilities over 10 ms window
end
%% Average no song and female output probabilities
allf=[out_p_4c_filt(:,2) out_p_2c_filt(:,2)];
meanf_out_p_filt=mean(allf,2);

allns=[out_p_4c_filt(:,1) out_p_2c_filt(:,1)];
meanns_out_p_filt=mean(allns,2);
%% Combine output probabilities
new_out_p=[meanns_out_p_filt meanf_out_p_filt out_p_4c_filt(:,3) out_p_4c_filt(:,4)];
%% Require probability of female song > 1.25* probability of male song
compareMF=(new_out_p(:,2)>(1.25*new_out_p(:,3)));
new_out_p(:,2)=compareMF.*new_out_p(:,2);
%% Require overlap prob >=0.85
compareO=(new_out_p(:,4)>=0.85);
new_out_p(:,4)=compareO.*new_out_p(:,4);
%% Make predictions
[~,predY]=max(new_out_p,[],2);
predY=predY-1;
%% Segment song 
disp('Segmenting song.')
 maleBoutInfo.w0=[];
 maleBoutInfo.w1=[];
 femaleBoutInfo.w0=[];
 femaleBoutInfo.w1=[];
 overlap.w0=[];
 overlap.w1=[];
 nosong.w0=[];
 nosong.w1=[];
 havemale=0;
 havefemale=0;
 haveoverlap=0;
 havenosong=0;
 
if (sum(predY==2)>1)
    [mInds ~]=find(predY==2);
%     [~, mInds]=find(predY==2);
    boutNum=1;
    maleBoutInfo.w0(boutNum)=mInds(1);
    trans=find(diff(mInds)~=1);
    maleBoutInfo.w0=[maleBoutInfo.w0 mInds(trans+1)'];
    maleBoutInfo.w1=mInds(trans);
    maleBoutInfo.w1=[maleBoutInfo.w1' mInds(end)]; 
    havemale=1;
else
    maleBoutInfo.w0=NaN;
    maleBoutInfo.w1=NaN;
end
    clear trans
% find female segments
if (sum(predY==1)>1)
    [fInds ~]=find(predY==1);
%     [~, fInds]=find(predY==1);
    boutNum=1;
    femaleBoutInfo.w0(boutNum)=fInds(1);
    trans=find(diff(fInds)~=1);
    femaleBoutInfo.w0=[femaleBoutInfo.w0 fInds(trans+1)'];
    femaleBoutInfo.w1=fInds(trans);
    femaleBoutInfo.w1=[femaleBoutInfo.w1' fInds(end)]; 
    havefemale=1;
else
    femaleBoutInfo.w0=NaN;
    femaleBoutInfo.w1=NaN;
end
clear trans
% find no song segments
if (sum(predY==0))
%     [nsInds ~]=find(predY==0);
    [~,nsInds]=find(predY==0);
    boutNum=1;
    nosong.w0(boutNum)=nsInds(1);
    trans=find(diff(nsInds)~=1);
    nosong.w0=[nosong.w0 nsInds(trans+1)'];
    nosong.w1=nsInds(trans);
    nosong.w1=[nosong.w1' nsInds(end)]; 
    havenosong=1;
else
    nosong.w0=NaN;
    nosong.w1=NaN;
end
clear trans
% find overlap segments
clear oInds
if (sum(predY==3)>1)
    [oInds ~]=find(predY==3);
    boutNum=1;
    overlap.w0(boutNum)=oInds(1);
    trans=find(diff(oInds)~=1);
    overlap.w0=[overlap.w0 oInds(trans+1)'];
    overlap.w1=oInds(trans);
    overlap.w1=[overlap.w1' oInds(end)]; 
    haveoverlap=1;
else 
    overlap.w0=NaN;
    overlap.w1=NaN;
end
    clear trans
%% Eliminate very short female pulses
% assign to second-place out_p if female duration is too short 
disp('Throwing out short female pulses.')
fdur=femaleBoutInfo.w1-femaleBoutInfo.w0;
fdur_min=4; %ms 
f_toreassign=fdur<(fdur_min/1000*10000);

% figure
% a(1)=subplot(2,1,1);
% plot(data)
% hold on
f_reassigned=0;
for i=1:length(f_toreassign)
    if (f_toreassign(i))
%         i
%         plot(femaleBoutInfo.w0(i):femaleBoutInfo.w1(i),data(femaleBoutInfo.w0(i):femaleBoutInfo.w1(i)),'g')
%         xlim([femaleBoutInfo.w0(i)-500 femaleBoutInfo.w1(i)+500])
%         pause
%           new_out_p(femaleBoutInfo.w0(i):femaleBoutInfo.w1(i),2)=0;
           [~, temp]=max(new_out_p(femaleBoutInfo.w0(i):femaleBoutInfo.w1(i),[1 3 4]),[],2);
           temp(temp==1)=0;
           predY(femaleBoutInfo.w0(i):femaleBoutInfo.w1(i))=temp;
           clear temp
           f_reassigned=1;
    end
end
%% Re-segment song 
if f_reassigned
    disp('Re-segmenting after throwing out short female pulses.')
    maleBoutInfo.w0=[];
     maleBoutInfo.w1=[];
     femaleBoutInfo.w0=[];
     femaleBoutInfo.w1=[];
     overlap.w0=[];
     overlap.w1=[];
     nosong.w0=[];
     nosong.w1=[];

if (sum(predY==2)>1)
    [mInds ~]=find(predY==2);
    boutNum=1;
    maleBoutInfo.w0(boutNum)=mInds(1);
    trans=find(diff(mInds)~=1);
    maleBoutInfo.w0=[maleBoutInfo.w0 mInds(trans+1)'];
    maleBoutInfo.w1=mInds(trans);
    maleBoutInfo.w1=[maleBoutInfo.w1' mInds(end)]; 
    havemale=1;
else
    maleBoutInfo.w0=NaN;
    maleBoutInfo.w1=NaN;
    havemale=0;
end
    clear trans
% find female segments
if (sum(predY==1)>1)
    [fInds ~]=find(predY==1);
    boutNum=1;
    femaleBoutInfo.w0(boutNum)=fInds(1);
    trans=find(diff(fInds)~=1);
    femaleBoutInfo.w0=[femaleBoutInfo.w0 fInds(trans+1)'];
    femaleBoutInfo.w1=fInds(trans);
    femaleBoutInfo.w1=[femaleBoutInfo.w1' fInds(end)]; 
    havefemale=1;
else
    femaleBoutInfo.w0=NaN;
    femaleBoutInfo.w1=NaN;
    havefemale=0;
end
clear trans
% find no song segments
if (sum(predY==0))
    [nsInds ~]=find(predY==0);
    boutNum=1;
    nosong.w0(boutNum)=nsInds(1);
    trans=find(diff(nsInds)~=1);
    nosong.w0=[nosong.w0 nsInds(trans+1)'];
    nosong.w1=nsInds(trans);
    nosong.w1=[nosong.w1' nsInds(end)]; 
    havenosong=1;
else
    nosong.w0=NaN;
    nosong.w1=NaN;
    havenosong=0;
end
clear trans
% find overlap segments
clear oInds
if (sum(predY==3)>1)
    [oInds ~]=find(predY==3);
    boutNum=1;
    overlap.w0(boutNum)=oInds(1);
    trans=find(diff(oInds)~=1);
    overlap.w0=[overlap.w0 oInds(trans+1)'];
    overlap.w1=oInds(trans);
    overlap.w1=[overlap.w1' oInds(end)]; 
    haveoverlap=1;
else 
    overlap.w0=NaN;
    overlap.w1=NaN;
    haveoverlap=0;
end
    clear trans
end

%% At edges of male song, require female song to have classification probability of 1 
% find female song occuring within 40 ms before male song
f_reassigned=0;
for i=1:length(maleBoutInfo.w0)
    if (~isnan(maleBoutInfo.w0(i)))
        tempdiff(i,:)=maleBoutInfo.w0(i)-femaleBoutInfo.w1;
        tempdiff_abs(i,:)=abs(tempdiff(i,:));
    %     tempdiff(i,:)=abs(tempdiff(i,:));
    %     f_before_m(i)=find(tempdiff(i,:)==min(tempdiff(i,:)));
        f_before_m_temp=find(tempdiff_abs(i,:)==min(tempdiff_abs(i,:)));
        if (length(f_before_m_temp)==1)
            f_before_m(i)=f_before_m_temp;
        else
            f_before_m(i)=min(f_before_m_temp);
        end
        clear f_before_m_temp
        f_before_m(i)=find(tempdiff(i,:)==min(tempdiff(i,:)));
        if ((maleBoutInfo.w0(i)-femaleBoutInfo.w1(f_before_m(i)))<0)
            f_before_m(i)=NaN;
        elseif (maleBoutInfo.w0(i)-femaleBoutInfo.w1(f_before_m(i))>400)
            f_before_m(i)=NaN;
        end
    else
        f_before_m(i)=NaN;
    end
end
f_before_m(isnan(f_before_m))=[];

% find female song occuring within 10 ms after male song
clear tempdiff
for i=1:length(maleBoutInfo.w1)
    if (~isnan(maleBoutInfo.w1(i)))
        tempdiff(i,:)=maleBoutInfo.w1(i)-femaleBoutInfo.w0;
    %     tempdiff(i,:)=abs(tempdiff(i,:));
    %     f_after_m(i)=find(tempdiff(i,:)==min(tempdiff(i,:)));
        tempdiff_abs(i,:)=abs(tempdiff(i,:));
        f_after_m_temp=find(tempdiff_abs(i,:)==min(tempdiff_abs(i,:)));
        if (length(f_after_m_temp)==1)
            f_after_m(i)=f_after_m_temp;
        else
            f_after_m(i)=max(f_after_m_temp);
        end
        if ((maleBoutInfo.w1(i)-femaleBoutInfo.w0(f_after_m(i)))>0)
            f_after_m(i)=NaN;
        elseif (maleBoutInfo.w1(i)-femaleBoutInfo.w0(f_after_m(i))<-50)
            f_after_m(i)=NaN;
        end
    else
        f_after_m(i)=NaN;
    end
end
f_after_m(isnan(f_after_m))=[];
% if probability of female song within these two segments is less than 1,
% call portion male
for i=1:length(f_before_m)
    checkprob(i)=max(new_out_p(femaleBoutInfo.w0(f_before_m(i)):femaleBoutInfo.w1(f_before_m(i)),2));
    if (~(checkprob(i)>0.99))
        predY(femaleBoutInfo.w0(f_before_m(i)):femaleBoutInfo.w1(f_before_m(i)))=2;
        f_reassigned=1;
    end
end
clear checkprob

for i=1:length(f_after_m)
    checkprob(i)=max(new_out_p(femaleBoutInfo.w0(f_after_m(i)):femaleBoutInfo.w1(f_after_m(i)),2));
    if (~(checkprob(i)>0.99))
        predY(femaleBoutInfo.w0(f_after_m(i)):femaleBoutInfo.w1(f_after_m(i)))=2;
        f_reassigned=1;
    end
end
clear checkprob
%% Re-segment
if f_reassigned
    disp('Re-segmenting after throwing out short female pulses.')
    maleBoutInfo.w0=[];
     maleBoutInfo.w1=[];
     femaleBoutInfo.w0=[];
     femaleBoutInfo.w1=[];
     overlap.w0=[];
     overlap.w1=[];
     nosong.w0=[];
     nosong.w1=[];

    if (sum(predY==2)>1)
        [mInds ~]=find(predY==2);
        boutNum=1;
        maleBoutInfo.w0(boutNum)=mInds(1);
        trans=find(diff(mInds)~=1);
        maleBoutInfo.w0=[maleBoutInfo.w0 mInds(trans+1)'];
        maleBoutInfo.w1=mInds(trans);
        maleBoutInfo.w1=[maleBoutInfo.w1' mInds(end)]; 
        havemale=1;
    else
        maleBoutInfo.w0=NaN;
        maleBoutInfo.w1=NaN;
        havemale=0;
    end
        clear trans
    % find female segments
    if (sum(predY==1)>1)
        [fInds ~]=find(predY==1);
        boutNum=1;
        femaleBoutInfo.w0(boutNum)=fInds(1);
        trans=find(diff(fInds)~=1);
        femaleBoutInfo.w0=[femaleBoutInfo.w0 fInds(trans+1)'];
        femaleBoutInfo.w1=fInds(trans);
        femaleBoutInfo.w1=[femaleBoutInfo.w1' fInds(end)]; 
        havefemale=1;
    else
        femaleBoutInfo.w0=NaN;
        femaleBoutInfo.w1=NaN;
        havefemale=0;
    end
    clear trans
    % find no song segments
    if (sum(predY==0))
        [nsInds ~]=find(predY==0);
        boutNum=1;
        nosong.w0(boutNum)=nsInds(1);
        trans=find(diff(nsInds)~=1);
        nosong.w0=[nosong.w0 nsInds(trans+1)'];
        nosong.w1=nsInds(trans);
        nosong.w1=[nosong.w1' nsInds(end)]; 
        havenosong=1;
    else
        nosong.w0=NaN;
        nosong.w1=NaN;
        havenosong=0;
    end
    clear trans
    % find overlap segments
    clear oInds
    if (sum(predY==3)>1)
        [oInds ~]=find(predY==3);
        boutNum=1;
        overlap.w0(boutNum)=oInds(1);
        trans=find(diff(oInds)~=1);
        overlap.w0=[overlap.w0 oInds(trans+1)'];
        overlap.w1=oInds(trans);
        overlap.w1=[overlap.w1' oInds(end)]; 
        haveoverlap=1;
    else 
        overlap.w0=NaN;
        overlap.w1=NaN;
        haveoverlap=0;
    end
        clear trans
end
%% Stitch together beginning male pulses
m_reassigned=0;
for i=1:length(maleBoutInfo.w0)-1
    following_diff(i)=maleBoutInfo.w0(i+1)-maleBoutInfo.w1(i);
    following_diff(i)=following_diff(i)/10;
    if (following_diff(i)<40)
        maleBoutInfo.w0(i+1)=maleBoutInfo.w0(i);
        maleBoutInfo.w0(i)=NaN;
        maleBoutInfo.w1(i)=NaN;
        m_reassigned=1;
    end
end
maleBoutInfo.w0(isnan(maleBoutInfo.w0))=[];
maleBoutInfo.w1(isnan(maleBoutInfo.w1))=[];
%% Collect pulses
disp('Collecting all pulses.')
if havemale
    for i=1:length(maleBoutInfo.w0)
        maleBoutInfo.x{i}=data(maleBoutInfo.w0(i):maleBoutInfo.w1(i));
    end
end

if havefemale
    for i=1:length(femaleBoutInfo.w0)
        femaleBoutInfo.x{i}=data(femaleBoutInfo.w0(i):femaleBoutInfo.w1(i));
    end
end

if haveoverlap
    for i=1:length(overlap.w0)
        overlap.x{i}=data(overlap.w0(i):overlap.w1(i));
    end
end

% save(fname_tosave,'data','maleBoutInfo','femaleBoutInfo','overlap','nosong','lower_envelope','upper_envelope','-v7.3');
% save(fname_tosave,'data','maleBoutInfo','femaleBoutInfo','overlap','nosong','-v7.3');
% disp('Saved segmented data.')
  %% 3.Calculate envelopes
disp('Calculating envelopes.')
if ~exist('upper_envelope','var')
   [upper_envelope lower_envelope]=envelope(data,50,'peak');
end
diffenv=upper_envelope-lower_envelope;
disp('Envelopes calculated.')
save(fname_tosave,'data','maleBoutInfo','femaleBoutInfo','overlap','nosong','diffenv','-v7.3');
disp('Saved segmented data in progress.')
%% Morlet transform
% morlet_freqs=100:10:900;
% omega0=2*pi;
% dt=1/10000;
% [morlet_amps,~] = fastWavelet_morlet_convolution_parallel(data,morlet_freqs,omega0,dt);

%% IPIs
% fname='T:\virilis_segmenter_cnn\combined_data\4001_femaleNoise_talmojan_trainingplusnoise_combined_withminfdur.mat';
% fname='T:\virilis_segmenter_cnn\combined_data\160910_0913ch30_2segs40012001_fmindur_segmented.mat'; % 1
% fname='T:\virilis_segmenter_cnn\combined_data\160913_1015ch30_2segs40012001_fmindur_segmented.mat'; % 2
% fname='T:\virilis_segmenter_cnn\combined_data\160920_0912ch26_2segs40012001_fmindur_segmented.mat'; % 3
% fname='T:\virilis_segmenter_cnn\combined_data\160920_0912ch32_2segs40012001_fmindur_segmented.mat'; % 4
% fname='T:\virilis_segmenter_cnn\combined_data\160920_1017ch28_2segs40012001_fmindur_segmented.mat'; % 4
% load(fname)
cd('../combine_segmenters/');
disp('Calculating IPIs.')
deltaF=0.2;
deltaM=0.03;
havemale=1;
havefemale=1;
haveoverlap=0;
havenosong=1;

if havemale
    disp('calculating male ipis')
    maleBoutInfo.ipis={};
    for i=1:length(maleBoutInfo.w0)
        env_temp=diffenv(maleBoutInfo.w0(i):maleBoutInfo.w1(i));
        [maxes mins]=peakdet(env_temp,deltaM);
        if isempty(maxes)
            [~, maxes]=max(env_temp);
        end
        max_ints=maxes(:,1)+maleBoutInfo.w0(i);
        maleBoutInfo.envpks{i}=max_ints;
        if (length(max_ints)>1)
            min_ints=mins(:,1)+maleBoutInfo.w1(i);
            maleBoutInfo.ipis{i}=diff(max_ints(:,1))./10; %ipis in ms
        else
            maleBoutInfo.ipis{i}=NaN;
        end
        clear maxes mins max_ints min_ints env_temp
    end
end

if haveoverlap
    disp('calculating overlap ipis')
    overlap.ipis={};
    for i=1:length(overlap.w0)
        env_temp=diffenv(overlap.w0(i):overlap.w1(i));
        [maxes mins]=peakdet(env_temp,deltaM);
        if isempty(maxes)
            [~, maxes]=max(env_temp);
        end
        max_ints=maxes(:,1)+overlap.w0(i);
        overlap.envpks{i}=max_ints;
        max_ints=maxes(:,1)+overlap.w0(i);
        if (length(max_ints)>1)
            min_ints=mins(:,1)+overlap.w1(i);
            overlap.ipis{i}=diff(max_ints(:,1))./10; %ipis in ms
        else
            overlap.ipis{i}=NaN;
        end
        clear maxes mins max_ints min_ints env_temp
    end
end

if havefemale
    disp('calculating female ipis')
    femaleBoutInfo.envpks={};
    femaleBoutInfo.ipis={};
    for i=1:length(femaleBoutInfo.w0)
        env_temp=diffenv(femaleBoutInfo.w0(i):femaleBoutInfo.w1(i));
        [maxes mins]=peakdet(env_temp,deltaF);
        if isempty(maxes)
            [~, maxes]=max(env_temp);
        end
        max_ints=maxes(:,1)+femaleBoutInfo.w0(i);
        femaleBoutInfo.envpks{i}=max_ints;
        if (length(max_ints)>1)
            min_ints=mins(:,1)+femaleBoutInfo.w1(i);
            femaleBoutInfo.ipis{i}=diff(max_ints(:,1))./10; %ipis in ms
        else
            femaleBoutInfo.ipis{i}=NaN;
        end
        clear maxes mins max_ints min_ints env_temp
    end
end

% figure
% a(1)=subplot(2,1,1);
% hold on
% plot(data)
% if ~isempty(nosong.w0)
%     if ~isnan(nosong.w0(1))
%         for i=1:length(nosong.w0)
%             if ~isnan(nosong.w0(i))
%                 plot(nosong.w0(i):nosong.w1(i),data(nosong.w0(i):nosong.w1(i)),'k')
%             end
%         end
%     end
% end
% if ~isempty(femaleBoutInfo.w0)
%     if ~isnan(femaleBoutInfo.w0(1))
%         for i=1:length(femaleBoutInfo.w0)
%             plot(femaleBoutInfo.w0(i):femaleBoutInfo.w1(i),data(femaleBoutInfo.w0(i):femaleBoutInfo.w1(i)),'c')
%         end
%     end
% end
% if ~isempty(maleBoutInfo.w0) 
%     if ~isnan(maleBoutInfo.w0(1))
%         for i=1:length(maleBoutInfo.w0)
%             plot(maleBoutInfo.w0(i):maleBoutInfo.w1(i),data(maleBoutInfo.w0(i):maleBoutInfo.w1(i)),'r')
%         end
%     end
% end
% if ~isempty(overlap.w0)
%     if ~isnan(overlap.w0(1))
%         for i=1:length(overlap.w0)
%             plot(overlap.w0(i):overlap.w1(i),data(overlap.w0(i):overlap.w1(i)),'y')
%         end
%     end
% end
% 
% a(2)=subplot(2,1,2);
% plot(diffenv,'k')
% hold on
% for i=1:length(maleBoutInfo.w0)
%     scatter(maleBoutInfo.envpks{i},diffenv(maleBoutInfo.envpks{i}),'ro')
% end
% for i=1:length(femaleBoutInfo.w0)
%     scatter(femaleBoutInfo.envpks{i},diffenv(femaleBoutInfo.envpks{i}),'co')
% end
% if haveoverlap
%     for i=1:length(overlap.w0)
%         scatter(overlap.envpks{i},diffenv(overlap.envpks{i}),'yo')
%     end
% end
% linkaxes(a,'x')
%% Separate female pulses
new.w0=[];
new.w1=[];
new.envpks=[];
for i=1:length(femaleBoutInfo.w0)
    pulsenum=length(femaleBoutInfo.envpks{i});
    if (pulsenum>1)
        temp.w0=[];
        temp.w0(1)=femaleBoutInfo.w0(i);
        temp.w1(pulsenum)=femaleBoutInfo.w1(i);
        temp.envpks=femaleBoutInfo.envpks{i};
        for j=2:pulsenum
             [a, b]=min(diffenv(femaleBoutInfo.envpks{i}(j-1):femaleBoutInfo.envpks{i}(j)));
             temp.w1(j-1)=b-1+femaleBoutInfo.envpks{i}(j-1);
             temp.w0(j)=b+femaleBoutInfo.envpks{i}(j-1);
             clear a b
        end
        new.w0=[new.w0 temp.w0];
        new.w1=[new.w1 temp.w1];
        new.envpks=[new.envpks femaleBoutInfo.envpks{i}'];
        clear temp
    else
        new.w0=[new.w0 femaleBoutInfo.w0(i)];
        new.w1=[new.w1 femaleBoutInfo.w1(i)];
        new.envpks=[new.envpks femaleBoutInfo.envpks{i}];
    end
end
for i=1:length(new.w0)
    new.x{i}=data(new.w0(i):new.w1(i));
end
%% Switch femaleBoutInfo (so that now femaleBoutInfo has pulses separated)
femaleBoutInfo_notseparated=femaleBoutInfo;
femaleBoutInfo=[];
femaleBoutInfo=new;
clear new
%% Measure SNR
signal=[];
noise=[]; 
if havemale
    for i=1:length(maleBoutInfo.w0)
        signal=[signal; data(maleBoutInfo.w0(i):maleBoutInfo.w1(i))];
    end
end
if havefemale
    for i=1:length(femaleBoutInfo.w0)
        signal=[signal; data(femaleBoutInfo.w0(i):femaleBoutInfo.w1(i))];
    end
end
if haveoverlap
    for i=1:length(overlap.w0)
        signal=[signal; data(overlap.w0(i):overlap.w1(i))];
    end
end
if havenosong
    for i=1:length(nosong.w0)
        noise=[noise; data(nosong.w0(i):nosong.w1(i))];
    end
end
snr.meas = mean( signal .^ 2 ) / mean( noise .^ 2 );
snr.meas_db = 10 * log10( snr.meas );
toc
%% Save
% save(fname_tosave,'data','maleBoutInfo','femaleBoutInfo','overlap','nosong','lower_envelope','upper_envelope','-v7.3');
% disp('Saved segmented data and ipis.')
%% FFT
% disp('Computing power spectra.')
% fname_tosave=strcat(fname_tosave(1:end-4),'_combined_EnvIpiPs.mat');
% fname_tosave=strcat(fname_tosave(1:end-4),'_EnvIpi.mat');
% Fs=10000;
% PrwSpec.info.length=10000;
% PwrSpec.info.Fs=Fs;
% % Male pulses
% disp('Calculating male power spectra')
% if havemale
%     mlength=max(maleBoutInfo.w1-maleBoutInfo.w0);
% %     PwrSpec.info.mlength=mlength;
%     for i=1:length(maleBoutInfo.w0)
%         if (maleBoutInfo.w1(i)-maleBoutInfo.w0(i)>1)
%             [Pm(i,:), Fm(i,:)]=pwelch(data(maleBoutInfo.w0(i):maleBoutInfo.w1(i)),maleBoutInfo.w1(i)-maleBoutInfo.w0(i),[],PrwSpec.info.length,Fs);
%             Pm_norm(i,:)=Pm(i,:)/max(Pm(i,:));
%             maleBoutInfo.pwrSpecFreqs=Fm(i,:);
%         else
%            Pm(i,:)=NaN*ones(1,length(Fm(i-1)));
%            Fm(i,:)=NaN*ones(1,length(Fm(i-1)));
%            Pm_norm(i,:)=NaN*ones(1,length(Fm(i-1)));
%         end
%     end
%     mean_Pm=nanmean(Pm,1);
%     mean_Pm_norm=nanmean(Pm_norm,1);
% 
%     maleBoutInfo.pwrSpec=Pm;
%     maleBoutInfo.NormPwrSpec=Pm_norm;
% else
%     PwrSpec.info.mlength=NaN;
%     maleBoutInfo.pwrSpecFreqs=NaN;
% end
% % Female pulses
% disp('Calculating female power spectra')
% if havefemale
%     flength_notseparated=max(femaleBoutInfo_notseparated.w1-femaleBoutInfo_notseparated.w0);
% %     PwrSpec.info.flength=flength;
%     for i=1:length(femaleBoutInfo.w0)
%         if ((femaleBoutInfo.w1(i)-femaleBoutInfo.w0(i))>1)
%             [Pf(i,:), Ff(i,:)]=pwelch(data(femaleBoutInfo.w0(i):femaleBoutInfo.w1(i)),femaleBoutInfo.w1(i)-femaleBoutInfo.w0(i),[],PrwSpec.info.length,Fs);
%             Pf_norm(i,:)=Pf(i,:)/max(Pf(i,:));
%             femaleBoutInfo.pwrSpecFreqs=Ff(i,:);
%         else 
%            Pf(i,:)=NaN*ones(1,length(Ff(i-1)));
%            Ff(i,:)=NaN*ones(1,length(Ff(i-1)));
%            Pf_norm(i,:)=NaN*ones(1,length(Ff(i-1)));
%         end
%     end
%     mean_Pf=nanmean(Pf,1);
%     mean_Pf_norm=nanmean(Pf_norm,1);
%     
%     femaleBoutInfo.pwrSpec=Pf;
%     femaleBoutInfo.NormPwrSpec=Pf_norm;
% else
%     PwrSpec.info.flength=NaN;
%     femaleBoutInfo.pwrSpecFreqs=NaN;
% end
% 
% disp('Calculating female power spectra')
% if havefemale
%     flength=max(femaleBoutInfo_notseparated.w1-femaleBoutInfo_notseparated.w0);
% %     PwrSpec.info.flength=flength;
%     for i=1:length(femaleBoutInfo_notseparated.w0)
%         if ((femaleBoutInfo_notseparated.w1(i)-femaleBoutInfo_notseparated.w0(i))>1)
%             [Pf(i,:), Ff(i,:)]=pwelch(data(femaleBoutInfo_notseparated.w0(i):femaleBoutInfo_notseparated.w1(i)),femaleBoutInfo_notseparated.w1(i)-femaleBoutInfo_notseparated.w0(i),[],PrwSpec.info.length,Fs);
%             Pf_norm(i,:)=Pf(i,:)/max(Pf(i,:));
%             femaleBoutInfo_notseparated.pwrSpecFreqs=Ff(i,:);
%         else 
%            Pf(i,:)=NaN*ones(1,length(Ff(i-1)));
%            Ff(i,:)=NaN*ones(1,length(Ff(i-1)));
%            Pf_norm(i,:)=NaN*ones(1,length(Ff(i-1)));
%         end
%     end
%     mean_Pf=nanmean(Pf,1);
%     mean_Pf_norm=nanmean(Pf_norm,1);
%     
%     femaleBoutInfo_notseparated.pwrSpec=Pf;
%     femaleBoutInfo_notseparated.NormPwrSpec=Pf_norm;
% else
%     PwrSpec.info.flength=NaN;
%     femaleBoutInfo_notseparated.pwrSpecFreqs=NaN;
% end
% % Overlap
% disp('Calculating overlap power spectra')
% if haveoverlap
%     olength=max(overlap.w1-overlap.w0);
% %     PwrSpec.info.olength=olength;
%     for i=1:length(overlap.w0)
%         if ((overlap.w1(i)-overlap.w0(i))>1)
%             [Po(i,:), Fo(i,:)]=pwelch(data(overlap.w0(i):overlap.w1(i)),overlap.w1(i)-overlap.w0(i),[],PrwSpec.info.length,Fs);
%             Po_norm(i,:)=Po(i,:)/max(Po(i,:));
%             overlap.pwrSpecFreqs=Fo(i,:);
%         else 
%            Po(i,:)=NaN*ones(1,length(Fo(i-1)));
%            Fo(i,:)=NaN*ones(1,length(Fo(i-1)));
%            Po_norm(i,:)=NaN*ones(1,length(Fo(i-1)));
%         end
%     end
%     mean_Po=nanmean(Po,1);
%     mean_Po_norm=nanmean(Po_norm,1);
%     
%     overlap.pwrSpec=Po;
%     overlap.NormPwrSpec=Po_norm;
% else
%     PwrSpec.info.olength=NaN;
%     overlap.pwrSpecFreqs=NaN;
% end

% Save file
% save(fname_tosave,'femaleBoutInfo','femaleBoutInfo_notseparated','maleBoutInfo','nosong','overlap','PwrSpec','data','snr','diffenv','-v7.3');
save(fname_tosave,'femaleBoutInfo','femaleBoutInfo_notseparated','maleBoutInfo','nosong','overlap','data','snr','diffenv','begNoiseEndSec','-v7.3');
disp('file with envelope and ipis saved.')
toc
% end