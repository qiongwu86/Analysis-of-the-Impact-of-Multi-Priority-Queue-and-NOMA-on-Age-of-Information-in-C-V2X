function [phyParams,varargin] = initiatePhyParameters(simParams,appParams,fileCfg,varargin)
% function [phyParams,varargin]= initiatePhyParameters(fileCfg,varargin)
%
% Settings of the PHY layer
% It takes in input the name of the (possible) file config and the inputs
% of the main function
% It returns the structure "phyParams"

fprintf('Physical layer settings\n');

%% Common PHY parameters

K = 1.38e-23;           % Avogadro's constant (J/K)
T0 = 290;               % Reference temperature (�K)

% [BwMHz]
% Bandwidth (MHz)
[phyParams,varargin]= addNewParam([],'BwMHz',10,'Bandwidth (MHz)','double',fileCfg,varargin{1});
if phyParams.BwMHz~=1.4 && phyParams.BwMHz~=5 && phyParams.BwMHz~=10 && phyParams.BwMHz~=20
    error('Invalid Bandwidth. Possible values: 1.4, 5, 10, 20 MHz');
end
if simParams.technology ~= 1 && phyParams.BwMHz~=10 % not only LTE and not 10MHz
    error('Invalid Bandwidth. Only with LTE it can be different from 10 MHz');
end

% [Raw]
% Sets the awareness range in meters
% The awareness range can be a number or a string of values
[phyParams,varargin]= addNewParam(phyParams,'Raw',150,'Awareness range (m)','integerOrArrayString',fileCfg,varargin{1});
for iRaw=1:length(phyParams.Raw)
    if phyParams.Raw(iRaw)<=0
        error('Error: "appParams.Raw" cannot be <= 0');
    end
    if iRaw>1 && phyParams.Raw(iRaw) <= phyParams.Raw(iRaw-1)
        error('Error: "appParams.Raw" must include numbers in ascending order');
    end
end

% [Ptx_dBm]
% Transmitted power (dBm)
[phyParams,varargin]= addNewParam(phyParams,'Ptx_dBm',23,'Transmitted power (dBm)','double',fileCfg,varargin{1});
 phyParams.Ptx = 10^((phyParams.Ptx_dBm-30)/10);

% [FixedPdensity]
[phyParams,varargin]= addNewParam(phyParams,'FixedPdensity',false,'Fixed power density (instead of fixed power)','bool',fileCfg,varargin{1});

% [Gt_dB]
% Transmitter antenna gain (dB)
[phyParams,varargin]= addNewParam(phyParams,'Gt_dB',3,'Transmitter antenna gain (dB)','double',fileCfg,varargin{1});
phyParams.Gt = 10^(phyParams.Gt_dB/10);

% [PtxERP_dBm]
% Effective radiated power (dBm)
%phyParams.PtxERP_dBm = phyParams.Ptx_dBm + phyParams.Gt_dB;
%phyParams.PtxERP = 10^((phyParams.PtxERP_dBm-30)/10);

% [Gr_dB]
% Receiver antenna gain (dB)
[phyParams,varargin]= addNewParam(phyParams,'Gr_dB',3,'Receiver antenna gain (dB)','double',fileCfg,varargin{1});
phyParams.Gr = 10^(phyParams.Gr_dB/10);

% [F_dB]
% Noise figure of the receiver (dB)
[phyParams,varargin]= addNewParam(phyParams,'F_dB',9,'Noise figure of the receiver (dB)','double',fileCfg,varargin{1});

% Noise power per MHz
phyParams.Pnoise_MHz_dBm = 10*log10(K*T0*1e6) + phyParams.F_dB + 30;
phyParams.Pnoise_MHz = 10^((phyParams.Pnoise_MHz_dBm-30)/10);

% [folderPERcurves]
% Select if using a BLER curve
[phyParams,varargin]= addNewParam(phyParams,'folderPERcurves','null','Name of folder with PER vs. SINR curves - Null if thresholds are used','string',fileCfg,varargin{1});
%phyParams.folderPERcurves
if ~strcmpi(phyParams.folderPERcurves,'null')
    if ~exist(phyParams.folderPERcurves, 'dir')
        error('Folder with PER vs. SINR curves not existing');
    end
    phyParams.PERcurves = true;
else
    phyParams.PERcurves = false;
end

if (simParams.typeOfScenario == 4) %Urban-ETSI has NLOS PER curve
    [phyParams,varargin]= addNewParam(phyParams,'folderPERcurvesNLOS','null','Name of folder with PER vs. SINR curves - Null if thresholds are used','string',fileCfg,varargin{1});
    %phyParams.folderPERcurves
    if ~strcmpi(phyParams.folderPERcurvesNLOS,'null')
        if ~exist(phyParams.folderPERcurvesNLOS, 'dir')
            error('Folder with PER vs. SINR curves not existing');
        end
    end
else
    [phyParams,varargin]= addNewParam(phyParams,'folderPERcurvesNLOS','null','Name of folder with PER vs. SINR curves - Null if thresholds are used','string',fileCfg,varargin{1});
    %phyParams.folderPERcurves
    if ~strcmpi(phyParams.folderPERcurvesNLOS,'null')
        if ~exist(phyParams.folderPERcurvesNLOS, 'dir')
            phyParams.folderPERcurvesNLOS=phyParams.folderPERcurvesLOS;
        end
    end
    
    
end

if simParams.technology ~= 1 % not only LTE
    
    %% PHY parameters related to 802.11p
    
    % power in 11p - all the bandwidth is used
    phyParams.P_ERP_MHz_11p_dBm = (phyParams.Ptx_dBm + phyParams.Gt_dB) - 10*log10(phyParams.BwMHz);
    phyParams.P_ERP_MHz_11p = 10^((phyParams.P_ERP_MHz_11p_dBm-30)/10);
    
    % [pWithLTEPHY]
    % Boolean to simulate LTE PHY in 11p
    [phyParams,varargin]= addNewParam(phyParams,'pWithLTEPHY',false,'Boolean to simulate LTE PHY in 11p','bool',fileCfg,varargin{1});
    
    if ~phyParams.pWithLTEPHY
        % [Mode]
        % 802.11p TX Mode, from 0 to 7 (was 1-8 before version 5.0.11)
        % Mode 2 is normally considered as the best trade-off
        [phyParams,varargin]= addNewParam(phyParams,'MCS_11p',2,'TX Mode','integer',fileCfg,varargin{1});
        if phyParams.MCS_11p<0 || phyParams.MCS_11p>7
            error('Error: "phyParams.MCS11p" must be within [0,7]');
        end
    else
        % [MCS]
        % 802.11p with LTE-PHY MCS
        % Sets the modulation and coding scheme between 0 and 28
        [phyParams,varargin]= addNewParam(phyParams,'MCS_pWithLTEphy',3,'Modulation and coding scheme','integer',fileCfg,varargin{1});
        if phyParams.MCS_pWithLTEphy<0 || phyParams.MCS_pWithLTEphy>28
            error('Error: "phyParams.MCS_pWithLTEphy" must be within [0,28]');
        end
    end
    
    % [rilModel11p]
    % Boolean to use the relative interference level (RIL) model in 11p
    [phyParams,varargin]= addNewParam(phyParams,'rilModel11p',false,'Boolean to use the relative interference level (RIL) model in 11p','bool',fileCfg,varargin{1});
    
    % [CW]
    % Contention Window
    [phyParams,varargin]= addNewParam(phyParams,'CW',15,'Contention Window','integer',fileCfg,varargin{1});
    if phyParams.CW~=3 && phyParams.CW~=7 && phyParams.CW~=15 && phyParams.CW~=31 && phyParams.CW~=63
        error('Error: "phyParams.CW" must be equal to either 3, 7, 15, 31, or 63');
    end
    
    % Slot time (s)
    phyParams.tSlot = 13e-6;
    
    % [AIFS]
    % Arbitration Inter-frame Space (s)
    % By specifications, should be AifsN=2 for high priority DENM, AifsN=3
    % for DENM, AifsN=6 for CAM
    [phyParams,varargin]= addNewParam(phyParams,'AifsN',6,'Arbitration inter-frame space','integer',fileCfg,varargin{1});
    if phyParams.AifsN<0
        error('Error: "phyParams.AifsN" cannot be negative');
    end
    tSIFS = 32e-6;
    phyParams.tAifs = tSIFS + phyParams.tSlot * phyParams.AifsN;
    
    % [CCAthr11p_notsync]
    % Sensitivity when the preamble is not decoded (clear channel assessment, CCA, threshold)
    % In IEEE 802.11-2007, pag. 618, the value indicated is -65 dBm
    [phyParams,varargin]= addNewParam(phyParams,'CCAthr11p_notsync',-65,'CCA threshold to set busy if not decodable [dBm]','double',fileCfg,varargin{1});
    phyParams.PrxSensNotSynch = 10^((phyParams.CCAthr11p_notsync-30)/10);
    
    % [CCAthr11p_sync]
    % Sensitivity when the preamble is decoded (clear channel assessment, CCA, threshold)
    % In IEEE 802.11-2007, pag. 618, the value indicated is -85 dBm
    % Used also for the CBR assessment
    [phyParams,varargin]= addNewParam(phyParams,'CCAthr11p_sync',-85,'CCA threshold to set busy if decodable [dBm]','double',fileCfg,varargin{1});
    phyParams.PrxSensWhenSynch = 10^((phyParams.CCAthr11p_sync-30)/10);
    
    if phyParams.PERcurves
        % If PER vs. SINR curves are used, the vector of SINR thresholds must be set
        % Set the file name
        fileName = sprintf('%s/PER_%s_MCS%d_%dB.txt',phyParams.folderPERcurves,'11p',phyParams.MCS_11p,appParams.beaconSizeBytes);
        [phyParams.sinrVector11p_LOS_dB,~,~,~,sinrInterp,perInterp] = readPERtable(fileName,10000);
        
        if (simParams.typeOfScenario == 4)
            fileName = sprintf('%s/PER_%s_MCS%d_%dB.txt',phyParams.folderPERcurvesNLOS,'11p',phyParams.MCS_11p,appParams.beaconSizeBytes);
            [phyParams.sinrVector11p_NLOS_dB,~,~,~,sinrInterpNLOS,perInterpNLOS] = readPERtable(fileName,10000);
        else
            phyParams.sinrVector11p_NLOS_dB=phyParams.sinrVector11p_LOS_dB;
            sinrInterpNLOS=sinrInterp;
            perInterpNLOS=perInterp;
        end
        %[phyParams.sinrVector11p_LOS_dB,perVector11p,~,~,sinrInterp,perInterp] = readPERtable(fileName,10000);
        %% DEBUG PER vs. SINR curves
        %plot(phyParams.sinrVector11p_LOS_dB,perVector11p,'pc');
        %hold on
        %plot(10*log10(sinrInterp),perInterp,'.k');
        %
        % Calculate the sinrThreshold, corresponding to 90% PER
        phyParams.sinrThreshold11p_NLOS_dB = 10*log10( sinrInterpNLOS(find(perInterpNLOS>=0.9,1)) );
        phyParams.sinrThreshold11p_LOS_dB = 10*log10( sinrInterp(find(perInterp>=0.9,1)) );
        % Find the duration of a packet, derived from the packet payload and Mode
        phyParams.tPck11p = packetDuration11p(appParams.beaconSizeBytes,phyParams.MCS_11p,-1,-1,false);
    else
        % If PER vs. SINR curves are not used, the SINR threshold must be set
        
        % [SINRthreshold11p]
        % SINR threshold for IEEE 802.11p
        [phyParams,varargin]= addNewParam(phyParams,'sinrThreshold11p_LOS',-1000,'SINR threshold for error assessment [dB]','double',fileCfg,varargin{1});
        phyParams.sinrThreshold11p_LOS_dB = phyParams.sinrThreshold11p_LOS;
        
        if ~phyParams.pWithLTEPHY
            % Find minimum SINR
            % -1000 means automatic setting
            if phyParams.sinrThreshold11p_LOS == -1000
                phyParams.sinrThreshold11p_LOS_dB = SINRmin11p(phyParams.MCS_11p);
            end
            %phyParams.gammaMin11p_dB = 10*log10(phyParams.gammaMin11p);
            
            % Find the duration of a packet, derived from the packet payload and Mode
            phyParams.tPck11p = packetDuration11p(appParams.beaconSizeBytes,phyParams.MCS_11p,-1,-1,phyParams.pWithLTEPHY);
        else
            % Find minimum SINR when using 11p with LTE PHY
            if phyParams.sinrThreshold11p_LOS == -1000
                [~,phyParams.sinrThreshold11p_LOS_dB,phyParams.NbitsHz] = findRBsBeaconSINRmin(phyParams.MCS_pWithLTEphy,appParams.beaconSizeBytes);
            else
                [~,~,phyParams.NbitsHz] = findRBsBeaconSINRmin(phyParams.MCS_pWithLTEphy,appParams.beaconSizeBytes);
            end
            %phyParams.gammaMin11p = 10^(phyParams.gammaMin11p_dB/10);
            
            % Find the duration of a packet when using 11p with LTE PHY
            phyParams.tPck11p = packetDuration11p(appParams.beaconSizeBytes,-1,phyParams.NbitsHz,phyParams.BwMHz,phyParams.pWithLTEPHY);
        end
        
        phyParams.sinrVector11p_LOS_dB = phyParams.sinrThreshold11p_LOS_dB;


        phyParams.sinrVector11p_NLOS_dB = phyParams.sinrVector11p_LOS_dB;
        phyParams.sinrThreshold11p_NLOS_dB = phyParams.sinrThreshold11p_LOS_dB;
    end
    phyParams.sinrThreshold11p_LOS = 10.^(phyParams.sinrThreshold11p_LOS_dB/10);
    phyParams.sinrVector11p_LOS = 10.^(phyParams.sinrVector11p_LOS_dB/10);


    phyParams.sinrThreshold11p_NLOS = 10.^(phyParams.sinrThreshold11p_NLOS_dB/10);
    phyParams.sinrVector11p_NLOS = 10.^(phyParams.sinrVector11p_NLOS_dB/10);
    
    if simParams.technology==2 && appParams.variableBeaconSize
        error('This part needs revision in version 5');
        %         if ~phyParams.pWithLTEPHY
        %             % If variable beacon size is selected, derive small packet
        %             % duration (11p standard PHY)
        %             phyParams.tPck11pSmall = packetDuration11p(appParams.beaconSizeSmallBytes,phyParams.MCS_11p,-1,-1,phyParams.pWithLTEPHY);
        %         else
        %            % If variable beacon size is selected, derive small packet
        %             % duration (11p with LTE PHY)
        %             phyParams.tPck11pSmall = packetDuration11p(appParams.beaconSizeSmallBytes,-1,phyParams.NbitsHz,phyParams.BwMHz,phyParams.pWithLTEPHY);
        %         end
    end
    
end

if simParams.technology ~= 2 % not only 11p
    
    %% PHY parameters related to LTE-V2V
    
    phyParams.Tfr = 0.01;                        % LTE frame period (s)
    phyParams.Tsf = phyParams.Tfr/10;            % LTE subframe period (s)
    phyParams.Tslot = phyParams.Tsf/2;           % LTE time slot period (s)
    phyParams.TsfGap = phyParams.Tsf * (1/14);   % LTE gap at the end of the subframe
    
    % [PnRB_dBm]
    % Noise power over a RB (dBm)
    phyParams.RBbandwidth = 180e3;    % Bandwidth of a RB (Hz)
    %phyParams.PnRB_dBm = 10*log10(K*T0*phyParams.RBbandwidth) + phyParams.F_dB + 30;
    %phyParams.PnRB = 10^((phyParams.PnRB_dBm-30)/10);
    
    % [MCS]
    % Sets the modulation and coding scheme between 0 and 28
    [phyParams,varargin]= addNewParam(phyParams,'MCS_LTE',3,'Modulation and coding scheme','integer',fileCfg,varargin{1});
    if phyParams.MCS_LTE<0 || phyParams.MCS_LTE>28
        error('Error: "phyParams.MCS_LTE" must be within [0,28]');
    end
    
    % [duplex]
    % Sets the duplexing: HD or FD
    [phyParams,varargin]= addNewParam(phyParams,'duplexLTE','HD','Duplexing type','string',fileCfg,varargin{1});
    if ~strcmp(phyParams.duplexLTE,'HD') && ~strcmp(phyParams.duplexLTE,'FD')
        error('Error: "phyParams.duplexLTE" must be equal to HD (half duplex) or FD (full duplex)');
    end
    
    % [Ksi_dB]
    % Self-interference cancellation coefficient (dB)
    if strcmp(phyParams.duplexLTE,'HD')
        phyParams.Ksi_dB = Inf;
        phyParams.Ksi = Inf;
    else
        [phyParams,varargin]= addNewParam(phyParams,'Ksi_dB',-110,'Self-interference cancellation coefficient (dB)','double',fileCfg,varargin{1});
        phyParams.Ksi = 10^(phyParams.Ksi_dB/10);
    end
    
    % [NumBeaconsFrequency]
    % Specify the number of BRs in the frequency domain
    % -1 -> default value: exploitation of all the available BRs
    [phyParams,varargin]= addNewParam(phyParams,'NumBeaconsFrequency',-1,'Specify the number of BRs in the frequency domain','integer',fileCfg,varargin{1});
    if phyParams.NumBeaconsFrequency<=0 && phyParams.NumBeaconsFrequency~=-1
        error('Error: "phyParams.NumBeaconsFrequency" must be -1 (all) or larger than 0');
    end
    
    % [ifAdjacent]
    % Select the subchannelization scheme: adjacent or non-adjacent PSCCH and PSSCH
    [phyParams,varargin]= addNewParam(phyParams,'ifAdjacent',true,'If using adjacent PSCCH and PSSCH','bool',fileCfg,varargin{1});
    if phyParams.ifAdjacent<0
        error('Error: "phyParams.ifAdjacent" must be equal to false or true');
    end
    
    % [BRoverlapAllowed]
    % Allow partial overlap of the allocations in the frequency domain
    [phyParams,varargin]= addNewParam(phyParams,'BRoverlapAllowed',false,'If a pratial overlap in frequency is allowed','bool',fileCfg,varargin{1});
    
    % [sizeSubchannel]
    % Select the subchannel size according to 3GPP TS 36.331
    % -1 -> default value: automatically select the best value
    [phyParams,varargin]= addNewParam(phyParams,'sizeSubchannel',-1,'Subchannel size','integer',fileCfg,varargin{1});
    if phyParams.sizeSubchannel<=0 && phyParams.sizeSubchannel~=-1
        error('Error: "phyParams.sizeSubchannel" must be -1 (best choice) or larger than 0');
    end
end

%% Channel Model parameters

% % [winnerModel]
% % Boolean to activate the use of WINNER+ B1 channel model (3GPP specifications)
% [phyParams,varargin]= addNewParam(phyParams,'winnerModel',true,'If using Winner+ channel model','bool',fileCfg,varargin{1});
% if phyParams.winnerModel < 0
%     error('Error: "phyParams.winnerModel" must be equal to false or true');
% end

% [channelModel]
% Integer to selecct the channel model 
% 0: WINNER+ B1 (3GPP specifications)
% 1: single slope model
% 2: two slopes model
% 3: three slopes model
[phyParams,varargin]= addNewParam(phyParams,'channelModel',0,'Channel model (0=WINNER+ B1, N=N slopes with N=1,2,3)','integer',fileCfg,varargin{1});
if phyParams.channelModel < 0 || phyParams.channelModel>3
    error('Error: "phyParams.channelModel" must be between 0 and 3');
end

% Non WINNER
if phyParams.channelModel>0 %~phyParams.winnerModel
    
    %% LOS non-Winner
    

    % [L0]
    % Path loss at 1m (dB)
    [phyParams,varargin]= addNewParam(phyParams,'L0_dB',47.86,'Path loss at 1m (dB)','double',fileCfg,varargin{1});
    if phyParams.L0_dB<=0
        error('L0_dB = %f not accceptable',phyParams.L0_dB);
    end
    
    % [beta]
    % Path loss exponent
    [phyParams,varargin]= addNewParam(phyParams,'beta',2.20,'Path loss exponent','double',fileCfg,varargin{1});
    if phyParams.beta<=0
        error('beta = %f not acceptable',phyParams.beta);
    end

    if phyParams.channelModel==1 % Single Slope
        
        phyParams.L0_near = 1000;
        phyParams.b_near = 0;
        phyParams.d_threshold1 = 0;
        phyParams.L0_mid = 1000;
        phyParams.b_mid = 0;
        phyParams.d_threshold2 = 0;
        phyParams.L0_far = 10^(phyParams.L0_dB/10);
        phyParams.b_far = phyParams.beta;        
    
    else
        
        % [dthreshold1]
        % threshold distance
        [phyParams,varargin]= addNewParam(phyParams,'d_threshold1',10,'First distance threshold','double',fileCfg,varargin{1});
        if phyParams.d_threshold1 < 1
            error('d_threshold1 = %f not larger than 1',phyParams.d_threshold1);
        end

        % [beta2]
        % Path loss exponent
        [phyParams,varargin]= addNewParam(phyParams,'beta2',2.20,'Path loss exponent of the second slope','double',fileCfg,varargin{1});
        if phyParams.beta2<=0
            error('beta2 = %f not acceptable',phyParams.beta2);
        end

        if phyParams.channelModel==2 % Two Slopes
        
            phyParams.d_threshold2 = phyParams.d_threshold1;
            phyParams.d_threshold1 = 0;
            phyParams.L0_near = 1000;
            phyParams.b_near = 0;
            phyParams.L0_mid = 10^(phyParams.L0_dB/10);
            phyParams.b_mid = phyParams.beta;
            phyParams.L0_far = phyParams.L0_mid * phyParams.d_threshold1^(phyParams.beta-phyParams.beta2);            
            phyParams.b_far = phyParams.beta2;
        
        elseif phyParams.channelModel==3 % three Slopes

            % [dthreshold2]
            % threshold distance
            [phyParams,varargin]= addNewParam(phyParams,'d_threshold2',10,'Second distance threshold','double',fileCfg,varargin{1});
            if phyParams.d_threshold2 < phyParams.d_threshold1
                error('d_threshold2 = %f cannot be smaller than d_threshold1 = %f',phyParams.d_threshold2,phyParams.d_threshold1);
            end

            % [beta3]
            % Path loss exponent
            [phyParams,varargin]= addNewParam(phyParams,'beta3',2.20,'Path loss exponent of the second slope','double',fileCfg,varargin{1});
            if phyParams.beta3<=0
                error('beta3 = %f not acceptable',phyParams.beta3);
            end

            phyParams.L0_near = 10^(phyParams.L0_dB/10);
            phyParams.b_near = phyParams.beta;
            phyParams.L0_mid = phyParams.L0_near * phyParams.d_threshold1^(phyParams.beta-phyParams.beta2);            
            phyParams.b_mid = phyParams.beta2;
            phyParams.L0_far = phyParams.L0_mid * phyParams.d_threshold2^(phyParams.beta2-phyParams.beta3);            
            phyParams.b_far = phyParams.beta3;
            
        else
            error('Error: channel model %d not implemented',phyParams.channelModel);
        end
    end
        
    %% NLOS non-Winner
    
    % [Abuild]
    % Attenuation every meter inside buildings (dB)
    [phyParams,varargin]= addNewParam(phyParams,'Abuild_dB',0.4,'Attenuation every meter inside buildings (dB)','double',fileCfg,varargin{1});
    phyParams.Abuild = 10^(phyParams.Abuild_dB/10);
    
    % [Awall]
    % Attenuation for each wall crossed (dB)
    [phyParams,varargin]= addNewParam(phyParams,'Awall_dB',6,'Attenuation for each wall crossed (dB)','double',fileCfg,varargin{1});
    phyParams.Awall = 10^(phyParams.Awall_dB/10);

    phyParams.L0_NLOS = 1000;
    phyParams.b_NLOS = 0;

else % Initialization of parameters for path loss calculation
%if phyParams.channelModel==0 %phyParams.winnerModel

    %% WINNER
    h = 1.5;                            % Antenna height [m]
    h_eff = h - 1;                      % Effective antenna height [m]
    f = 5.9;                            % Central frequency [GHz]
    c = 3e8;                            % Speed of light [m/s]
    Dbp = 4*h_eff*h_eff*f*10^9/c;       % Breakpoint distance [m]
    
    phyParams.d_threshold1 = 0;
    phyParams.d_threshold2 = Dbp;
    phyParams.L0_near = 1000;
    phyParams.b_near = 0;
    phyParams.L0_mid = 10^((27.0 + 20.0*log10(f))/10);
    phyParams.b_mid = 2.27;
    phyParams.L0_far = 10^((7.56-34.6*log10(h_eff)+2.7*log10(f))/10);
    phyParams.b_far = 4;    
    
    phyParams.L0_NLOS = 10^((5.83*log10(h)+18.38+23*log10(f))/10);
    phyParams.b_NLOS = (44.9-6.55*log10(h));

    % NLOS
    if simParams.typeOfScenario == 4 %ETSI Manhattan Layout        
        phyParams.L0_NLOS_0=10^((17.3+3*log10(f))/10);
        %n_j
        phyParams.b_nj0=-12.5;
        phyParams.n_threshold=1.84;
        phyParams.n_j0=2.8;
        phyParams.a_bj0=0.0024; %n_j=max(n_j0-a_bj0*d_k,n_threshold) %N_0=10
        phyParams.L0_NLOS=10^((17.3+3*log10(f)- phyParams.b_nj0.*phyParams.n_threshold)/10);
        phyParams.b_NLOS=1.84;
    else
        %Hexagonal Layout
        phyParams.L0_NLOS = 10^((5.83*log10(h)+18.38+23*log10(f))/10);
        phyParams.b_NLOS = (44.9-6.55*log10(h));        
    end
    
end

% [stdDevShadowLOS_dB]
% Standard deviation of shadowing in LOS (dB)
[phyParams,varargin]= addNewParam(phyParams,'stdDevShadowLOS_dB',3,'Standard deviation of shadowing in LOS (dB)','integer',fileCfg,varargin{1});

% [stdDevShadowNLOS_dB]
% Standard deviation of shadowing in NLOS (dB)
[phyParams,varargin]= addNewParam(phyParams,'stdDevShadowNLOS_dB',4,'Standard deviation of shadowing in NLOS (dB)','integer',fileCfg,varargin{1});

% Temporary parameters
% %% MCO additional settings
% if simParams.mco_nVehInterf>0
%     % [mco_interfERP]
%     % Power of MCO interferers
%     [phyParams,varargin]= addNewParam(phyParams,'mco_interfERP',0,'Power of MCO interferers (dBm)','double',fileCfg,varargin{1});
%     % From dBm to linear
%     phyParams.mco_interfERP = 10^((phyParams.mco_interfERP-30)/10);
%     % [mco_interfERP]
%     % Residual power from adjacent channel
%     [phyParams,varargin]= addNewParam(phyParams,'mco_resPowerFromAdjacent',-33,'Residual power from main to adjacent channel  (dB)','double',fileCfg,varargin{1});
%     if phyParams.mco_resPowerFromAdjacent > 0
%         error('phyParams.mco_resPowerFromAdjacent must be <= 0');
%     end
%     % From dBm to linear
%     phyParams.mco_resPowerFromAdjacent = 10^(phyParams.mco_resPowerFromAdjacent/10);
%     
%     [phyParams,varargin]= addNewParam(phyParams,'mco_interfNeglectedToMainChannel',false,'If interference from adjacent to main should be neglected','bool',fileCfg,varargin{1});
% end

fprintf('\n');

end
