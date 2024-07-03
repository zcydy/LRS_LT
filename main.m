%ZCY LRS all in one 
% the file is created for GRL paper submission
% structre of this file:
% LRS data read 
% Step1~step4(pls refer to the submitted paper)

% note that the colormap , colorbar or makers of plots may vary from the
% ones in submitted paper
%% LRS read
clc
clear all
% load dB and linear scale bscan, which are read from DARTS
load("echo_db.mat")
load("echo_ln.mat")
% Display the image
% imagesc(echoPower_dB);
range=linspace(-25*210,25*790,1000);
figure(1) % subplot of figure 1 in paper
imagesc(1:45196,range,echoPower_dB-max(max(echoPower_dB)));
colormap('Parula');
colorbar;
xlabel('Latitude[degree]');
ylabel('Range[m]');
ylim([-5000,5000])
title('SAR Processed Echo Power [dB]');
xticks([-1e-5 0 1e-5])
xticklabels({'-0.1',' 0','0.1'})
%%

time=range*2/3e8; %figure 2 in paper

figure(2)
subplot(2,5,1)%anomaly
% subplot(2,1,1)%anomaly
N= 19444;
Fs=12e6;
[sst,f]=wsst(echoPower_ln(:,N),Fs);
hp = pcolor(time,f+3.9e6,abs(sst));
hp.EdgeColor = "none";
ylim([4e6,6e6])
xlim([-1e-5,1.5e-5])
xticks([-1e-5 0 1e-5])
xticklabels({'-0.1',' 0','0.1'})
ylabel("Frequency (Hz)")
xlabel("time[\mus]")
title('(a)')


subplot(2,5,6)%anomaly
% subplot(2,1,2)%anomaly
sst_edge=edge(abs(sst(160:190,200:225)),"sobel",3);
imshow(flipud(sst_edge))
% [A_hat1,E_hat1,iter1] = inexact_alm_rpca(double(sst_edge));
% A_hat1(A_hat1<0.8*max(max(A_hat1)))=0;
% imagesc(A_hat1)

% ylabel("Frequency (Hz)")
% xlabel("time[\mus]")
% set(gca, 'Visible', 'off')
title('(b)')


subplot(2,5,2) %flat
N= 17761;
[sst,f]=wsst(echoPower_ln(:,N),Fs);
hp = pcolor(time,f+3.9e6,abs(sst));
hp.EdgeColor = "none";
ylim([4e6,6e6])
xlim([-1e-5,1.5e-5])
ylabel("Frequency (Hz)")
xlabel("time[\mus]")
title('(c)')


subplot(2,5,7)%anomaly
sst_edge=edge(abs(sst(160:190,200:225)),"sobel",3);
imshow(flipud(sst_edge))
title('(d)')


subplot(2,5,3)% complicated 
N= 13332;
[sst,f]=wsst(echoPower_ln(:,N),Fs);
hp = pcolor(time,f+3.9e6,abs(sst));
hp.EdgeColor = "none";
ylim([4e6,6e6])
xlim([-1e-5,1.5e-5])
ylabel("Frequency (Hz)")
xlabel("time[\mus]")
title('(e)')


subplot(2,5,8)%anomaly
sst_edge=edge(abs(sst(160:190,200:225)),"sobel",3);
imshow(flipud(sst_edge))
title('(f)')

load("sim_rough_lt.mat") %simulation data of lava tube
subplot(2,5,4)% simulation_lt
s=bscan_img(:,15);
down_s=downsample(s,16);
down_s(205:end)=0;
t=linspace(0e-5,4e-5,length(down_s));
% t=linspace(time(1),time(end),length(down_s));
Fs=1/(4e-5/length(down_s));
[sst,fs]=wsst(down_s,Fs);
hp=pcolor(t,fs,abs(sst));
hp.EdgeColor = "none";
ylim([4e6,10e6])
xlim([0.6e-5,0.85e-5])
xticks([0.6e-5 0.7e-5 0.8e-5])
xticklabels({'-0.1',' 0','0.1'})
ylabel("Frequency (Hz)")
xlabel("time[\mus]")
title('(g)')

load("sim_rough.mat") %simulation data of flat surface
subplot(2,5,9)% simulation_lt
sst_edge=edge(abs(sst(170:200,170:195)),"sobel",0.05);
% sst_edge=edge(abs(sst),"sobel",0.1);
imshow(flipud(sst_edge))
[A_hat1,E_hat1,iter1] = inexact_alm_rpca(double(sst_edge));
% A_hat1(A_hat1<0.8*max(max(A_hat1)))=0;
s1=rank(A_hat1);
title('(h)')

subplot(2,5,5)% simulation_flat
s=bscan_img(:,3);
down_s=downsample(s,16);
down_s(205:end)=0;
t=linspace(0e-5,4e-5,length(down_s));
% t=linspace(time(1),time(end),length(down_s));
Fs=1/(4e-5/length(down_s));
[sst,fs]=wsst(down_s,Fs);
hp=pcolor(t,fs,abs(sst));
hp.EdgeColor = "none";
ylim([4e6,10e6])
xlim([0.6e-5,0.85e-5])
xticks([0.6e-5 0.7e-5 0.8e-5])
xticklabels({'-0.1',' 0','0.1'})
ylabel("Frequency (Hz)")
xlabel("time[\mus]")
title('(i)')


subplot(2,5,10)%anomaly
sst_edge=edge(abs(sst(170:200,170:195)),"sobel",0.05);
% sst_edge=edge(abs(sst),"sobel",0.1);
imshow(flipud(sst_edge))
[A_hat1,E_hat1,iter1] = inexact_alm_rpca(double(sst_edge));
% A_hat1(A_hat1<0.8*max(max(A_hat1)))=0;
s2=rank(A_hat1);
title('(j)')




%%
%step 1~4 note that they should be inside one loop, but for reaseach
%purposes , it is disassembled to 4, which causes high computational load.
%for time save, the results are also unloadded as sparsity.mat ranky.mat
%and ano.mat
start=13861;
stop=36770;
sparsity1=zeros(1,stop-start+1);
ranky1=zeros(1,stop-start+1);
ano1=zeros(1,stop-start+1);

for C=start:stop %step 1  spasity
    sst=wsst(echoPower_ln(:,C));%WSSTs
    sst_edge=edge(abs(sst(160:190,200:225)),"sobel",0.5); %edge detecion
    numNonZero =nnz(sst_edge);
    numZero=numel(sst_edge)-numNonZero;% Sparsity

    if numZero == 0
        ratio = Inf;% avoid dividing zero
    else
        ratio =numNonZero /numZero;
    end
        sparsity1(C-(start-1))=ratio;

    if sparsity1(C-(start-1))>=0.25
        [A_hat1,E_hatl,iter1]= inexact_alm_rpca(double(sst_edge));
         %A_hatl(A_hat1<0.8*max(max(A_hat1)))=0;
         A_hat1(A_hat1<0.8*max(max(A_hat1)))=0;
        ranky1(C-(start-1))=rank(A_hat1);
    else
        ranky1(C-(start-1))=0;
    end


    if ranky1(C-(start-1))<=2
        ano1(C-(start-1))=ranky1(C-(start-1));
    else
        ano1(C-(start-1))=0;
    end

end

n=5;% number of CAPs
sequenceLength =length(ano1);
indices =[];

for i=1:(sequenceLength-n+1)
    if all(ano1(i:i+n-1)>= 2)
        indices =[indices,i];%记录起始索引
    end
end

CAPs=zeros(1,sequenceLength );
CAPs(indices)=1.5;


%%
%lat. and lon
load lat.mat
load lon.mat

figure(3)
set(gcf, 'Color', 'white')

subplot(4,1,1)
plot(sparsity1,'color','[0.00,0.00,0.00]')
title('(a)')
setupPlot(start, stop)


subplot(4,1,2)
plot(ranky1,'color','[0.00,0.00,0.00]')
title('(b)')
setupPlot(start, stop)

subplot(4,1,3)
plot(ano1,'color','[0.00,0.00,0.00]')
title('(c)')
ylim([0,2.5])
setupPlot(start, stop)

subplot(4,1,4)
plot(CAPs,'color','[0.00,0.00,0.00]')
title('(d)')
ylim([0,3])

last_annotated_idx = 0;
for i=1:length(CAPs)
   if CAPs(i)==1.5
        if abs(i - last_annotated_idx) >= 300
            lat_val = latitudes(i+start);
            lon_val = longitudes(i+start);
            annotation = sprintf('%.2f°, %.2f°', lat_val, lon_val);
            text(i,1.8,annotation, 'Rotation', 45,'FontSize',8)
            last_annotated_idx = i;  % Update the last annotated index
        end
   end
end

setupPlot(start, stop)
    % xlim([1, 23000]);
    % xticks([1, stop - start]);
    % xticklabels({'-4.2739°', '53.7312°'});
    % ylabel('Sparsity');
    % xlabel('latitude[°]');


function setupPlot(start, stop)
    xlim([1, 23000]);
    xticks([1, stop - start]);
    xticklabels({'-4.2739°', '53.7312°'});
    ylabel('Sparsity');
    xlabel('latitude[°]');
    box off;  % Remove the box

    % Add a secondary axes
    % ax1 = axes('Position', get(gca, 'Position'), 'XAxisLocation', 'top', ...
    %     'YAxisLocation', 'right', 'Color', 'none', 'XColor', 'k', 'YColor', 'k');
    % set(ax1, 'XTick', [], 'YTick', []);  % Remove ticks and labels on the secondary axes
end