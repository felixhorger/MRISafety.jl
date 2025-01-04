%% Calculate forbidden frequencies amplitude
clearvars; close all; clc

%% Input filenames

test_dsvfile = 'DataForKCL/SimulationHighRes/SimulationProtocol_GRX.dsv';

reference_dsvfile = 'CoC_ForbiddenFrequencies/7T_EPI_REFERENCE_GRX.dsv';

%% Resonance frequency values 

gamma_H = 42.576; % [MHz/T]
gamma_P = 17.235; % [MHz/T]

ratio = gamma_H/gamma_P;

%%% Terra values
% Resonance frequency #1
freqres{1}(1) = 550; %595;
freqres{1}(2) = 100; %130;
% Resonance frequency #2
freqres{2}(1) = 1100; %1030;
freqres{2}(2) = 300; %250;

%% Compute MSA for test an reference

[TEST_MSA,TEST_freq] = calc_gradient_MSA(test_dsvfile);
[REF_MSA,REF_freq]   = calc_gradient_MSA(reference_dsvfile);

REF_df  = REF_freq(2)-REF_freq(1);
TEST_df = TEST_freq(2)-TEST_freq(1);
for n=1:2
    idx_REFfreqres{n} = REF_freq>=(freqres{n}(1)-freqres{n}(2)/2) & REF_freq<=(freqres{n}(1)+freqres{n}(2)/2);
    idx_TESTfreqres{n} = TEST_freq>=(freqres{n}(1)-freqres{n}(2)/2) & TEST_freq<=(freqres{n}(1)+freqres{n}(2)/2);
    
    Total_REF_MSA{n} = sum(REF_MSA(idx_REFfreqres{n}))*REF_df; %[mT/m/s]
    Total_TEST_MSA{n} = sum(TEST_MSA(idx_TESTfreqres{n}))*TEST_df; %[mT/m/s]
end

%% Output comparison

cm = lines(4);
figure('position',[100 100 900 500],'color','w')
plot(REF_freq,abs(REF_MSA),'linewidth',2); hold on
plot(TEST_freq,abs(TEST_MSA)*ratio,'linewidth',2); 
for n=1:2
    area(freqres{n}(1)+[-freqres{n}(2)/2, freqres{n}(2)/2],ones(2,1)*max([REF_MSA(:);TEST_MSA(:)]),'FaceColor',cm(n+2,:),'FaceAlpha',0.25,'EdgeAlpha',0)
end
ylim([0 max([REF_MSA(:);TEST_MSA(:)])]); xlim([0 2e3])
xlabel('Frequency (Hz)'); ylabel('Gradient Fourier Transform (mT/m)')
legend('Reference','Test','Res. Freq. #1','Res. Freq. #2')

for n=1:2
    MSA_ratio(n) = Total_TEST_MSA{n}/Total_REF_MSA{n};
end

if all(MSA_ratio<=1)
    fprintf(1,'Mean Squared Amplitude is below reference protocol - test protocol should be safe to run.\n')
else
    idx_fail = find(MSA_ratio>1);
    if numel(idx_fail)==numel(MSA_ratio)
        fprintf(1,'Mean Squared Amplitude is above reference protocol for ALL resonant frequencies - DO NOT run provided test protocol.\n')    
    else
        fprintf(1,'Mean Squared Amplitude is above reference protocol for AT LEAST ONE resonant frequency - NOT RECOMMENDED to run provided test protocol.\n')
    end
    
    for ii=idx_fail
        fprintf(1,'\tMSA is %.1f%% above the reference for resonant frequency %dHz\n',100*(MSA_ratio(ii)-1),freqres{ii}(1))
    end 
end
    