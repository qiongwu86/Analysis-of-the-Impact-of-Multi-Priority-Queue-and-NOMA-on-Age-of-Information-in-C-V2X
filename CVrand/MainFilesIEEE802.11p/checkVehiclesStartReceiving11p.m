function [timeManagement,stationManagement,sinrManagement,outputValues] = checkVehiclesStartReceiving11p(idEvent,indexEventInActiveIDs,timeManagement,stationManagement,sinrManagement,simParams,phyParams,outParams,outputValues,simValues)

% Variable used for easier reading
activeIDs = stationManagement.activeIDs;

if idEvent==-1
    sameChannel = 1;
else
    sameChannel = (stationManagement.vehicleChannel==stationManagement.vehicleChannel(idEvent));
end

%% The nodes that start receiving are identified
% They are those that:
% A.  are in idle or in backoff (do not transmit and are not
% already receiving)
% B. do not end the backoff in the next time slot (a 1e-10
% margin is added due to problems with the representation of
% floating point numbers)
% C+D. receive this signal with sufficient quality (= (C) are able to
% decode the preamble, since SINR>SINR_min) OR (D) do not receive the
% signal with sufficient quality, but perceive the channel as
% busy
rxPowerTotNow_MHz = sinrManagement.P_RX_MHz * (stationManagement.vehicleState(activeIDs)==3);

idleOBU = ( (stationManagement.vehicleState(activeIDs)==1) + (stationManagement.vehicleState(activeIDs)==2) );
notEndingBackoff = (timeManagement.timeNextTxRx11p(activeIDs) >= (timeManagement.timeNow+phyParams.tSlot-1e-10));

alreadyReceingAnotherPreamble = ( stationManagement.vehicleState(activeIDs)==9 .* ...
    ( (timeManagement.timeNow-sinrManagement.instantThisSINRstarted11p(activeIDs)) < 4e-6 ...
    | sinrManagement.idFromWhichRx11p(activeIDs) == activeIDs ) );

if indexEventInActiveIDs>0
    % Normal case - one node is transmitting
    %
    % Approach 1: I check the SINR and if above a threshold I start decoding
    % The SINR corresponding to PER = 0.9 is used
    sinrThr=  phyParams.LOS(:,indexEventInActiveIDs)*phyParams.sinrThreshold11p_LOS+... %LOS
            (1-phyParams.LOS(:,indexEventInActiveIDs))*phyParams.sinrThreshold11p_LOS;

    decodingThisPreample = sinrManagement.P_RX_MHz(:,indexEventInActiveIDs) ./ (phyParams.Pnoise_MHz + (rxPowerTotNow_MHz-sinrManagement.P_RX_MHz(:,indexEventInActiveIDs)) + sinrManagement.coex_InterfFromLTEto11p(activeIDs)/phyParams.BwMHz) >  sinrThr;
    % % Approach 2: I check that the received power is above the sensing
    % % threshold for decodable signals (-85 dBm in specs)
    % C = (sinrManagement.P_RX_MHz(:,indexEventInActiveIDs)*phyParams.BwMHz) >= phyParams.PrxSensWhenSynch;
    % - Approach 1 is used - gives significantly better performance
    %
    % From version 5.3.1, the channel check is added
    % If the channel is not the same, then C is set to 0
    decodingThisPreample = decodingThisPreample .* sameChannel(stationManagement.activeIDs);
else
    % Only interference, with no real useful signal
    % In this case C is always 'false'
    decodingThisPreample = zeros(length(sinrManagement.P_RX_MHz(:,1)),1);
end
% in D, an additional factor 'coex_virtualInterference' is used to block all 11p transmissions (used
% in coexistence method A) - it could be added also in C, but adding it in D
% is sufficient
receivingHighPower = (rxPowerTotNow_MHz*phyParams.BwMHz +...
    sinrManagement.coex_InterfFromLTEto11p(activeIDs) + ...
    sinrManagement.coex_virtualInterference(activeIDs)) >= phyParams.PrxSensNotSynch;

% Revised From version 5.3.2
ifStartReceiving = logical(alreadyReceingAnotherPreamble .* decodingThisPreample + ...
        idleOBU .* notEndingBackoff .* (decodingThisPreample+receivingHighPower) );
%%

%%
% Focusing on those that start receiving
% The backoff is freezed if the node was in vState==2 (backoff)
% State is set to 9 (receiving)
% SINR is reset and initial instant is set to now
% 'timeNextTxRx' is set to infinity
% The node from which receiving is set
vehiclesFreezingList=activeIDs(logical(ifStartReceiving.*(stationManagement.vehicleState(activeIDs)==2)));
for idVehicle = vehiclesFreezingList'
    stationManagement.nSlotBackoff11p(idVehicle) = freezeBackoff11p(timeManagement.timeNow,timeManagement.timeNextTxRx11p(idVehicle),phyParams.tSlot,stationManagement.nSlotBackoff11p(idVehicle));
    % DEBUG BACKOFF
    printDebugBackoff11p(timeManagement.timeNow,'11p backoff freeze',idVehicle,stationManagement,outParams)
end
stationManagement.vehicleState(activeIDs(ifStartReceiving)) = 9;
sinrManagement.sinrAverage11p(activeIDs(ifStartReceiving)) = 0;
sinrManagement.interfAverage11p(activeIDs(ifStartReceiving)) = 0;
sinrManagement.instantThisSINRavStarted11p(activeIDs(ifStartReceiving)) = timeManagement.timeNow;
sinrManagement.instantThisSINRstarted11p(activeIDs(ifStartReceiving)) = timeManagement.timeNow;
timeManagement.timeNextTxRx11p(activeIDs(ifStartReceiving)) = Inf;
%%

%% Update of idFromWhichRx11p
if idEvent>0
    sinrManagement.idFromWhichRx11p(activeIDs(ifStartReceiving)) = idEvent * sameChannel(activeIDs(ifStartReceiving)) + activeIDs(ifStartReceiving) .* (~sameChannel(activeIDs(ifStartReceiving)));
    %    
    % Coexistence   
    % Save 11p as detected by LTE if the SINR is above the threshold
    if simParams.technology == 4 && simParams.coexMethod==3 
        sinrManagement.coex_lteDetecting11pTx(activeIDs) = logical( stationManagement.vehicleState(activeIDs)==100 .* decodingThisPreample );
    end

    % Coexistence, dynamic slots, calculation of CBR_11p: LTE nodes check if
    % starting detecting an 11p message
    if simParams.technology == 4 && simParams.coexMethod~=0 && simParams.coex_slotManagement==2 && simParams.coex_cbrTotVariant==2
        % In this case, nodes must be LTE, not already
        % detecting, receiving with sufficient quality
        % 
        % The inteference is saved in sensedPowerByLteNo11p per beacon
        % resource and needs to be converted into "per subchannel"
        if ~isempty(sinrManagement.sensedPowerByLteNo11p)
            interferenceFromLTEnodesPerSubframe = (sum(sinrManagement.sensedPowerByLteNo11p)/length(sinrManagement.sensedPowerByLteNo11p(:,1)))';
        else
            interferenceFromLTEnodesPerSubframe = 0; %zeros(length(),1);
        end
        sinrThr=phyParams.LOS(stationManagement.indexInActiveIDs_ofLTEnodes,indexEventInActiveIDs).*phyParams.sinrThreshold11p_LOS+...
                (1-phyParams.LOS(stationManagement.indexInActiveIDs_ofLTEnodes,indexEventInActiveIDs)).*phyParams.sinrThreshold11p_NLOS;
        decodingThePreamble = (sinrManagement.P_RX_MHz(stationManagement.indexInActiveIDs_ofLTEnodes,indexEventInActiveIDs) ./ (phyParams.Pnoise_MHz + (rxPowerTotNow_MHz(stationManagement.activeIDsLTE)-sinrManagement.P_RX_MHz(stationManagement.indexInActiveIDs_ofLTEnodes,indexEventInActiveIDs)) + interferenceFromLTEnodesPerSubframe) > sinrThr);
        ifStartDetecting11p = logical( (stationManagement.vehicleState(stationManagement.activeIDsLTE)==100)...
            .* ~sinrManagement.coex_detecting11p(stationManagement.activeIDsLTE)...
            .* decodingThePreamble );
        sinrManagement.coex_detecting11p(stationManagement.activeIDsLTE(ifStartDetecting11p)) = true; 
        timeManagement.cbr11p_timeStartBusy(stationManagement.activeIDsLTE(ifStartDetecting11p)) = timeManagement.timeNow;
        % DEBUG
        %fp = fopen('_Debug_LTE_CBR11p.xls','a');
        %fprintf(fp,'%f\t%d\t\n',timeManagement.timeNow,sum(decodingThePreamble));
        %fclose(fp);
    end    
else
    % if State 9 is due to interference, the idFromWhichRx must be set to
    % 'self'
    sinrManagement.idFromWhichRx11p(activeIDs(ifStartReceiving)) = activeIDs(ifStartReceiving);
end
%%

%% The channel busy ratio is updated 
% A different threshold (-85 dBm in ETSI EN 302 571) needs to be considered in this case
% 
if simParams.cbrActive && ~isempty(stationManagement.channelSensedBusyMatrix11p)
    % The cbr11p_timeStartBusy must be set to those that were not sensing
    % the channel as busy (cbr11p_timeStartBusy=-1) and are now sensing it
    % busy

    % First of all, E is calculated, which are those nodes transmitting
    % or perceiving interference above -85
    nodesThatMightStart = (stationManagement.vehicleState(activeIDs)==3) | ...
        ((rxPowerTotNow_MHz*phyParams.BwMHz + sinrManagement.coex_InterfFromLTEto11p(activeIDs) +...
         sinrManagement.coex_virtualInterference(activeIDs)) >= phyParams.PrxSensWhenSynch);
    
    % Then the nodes starting the channel busy are calculated
    ifStartCBR = (timeManagement.cbr11p_timeStartBusy(activeIDs)==-1 & nodesThatMightStart);
    
    % Time is updated
    timeManagement.cbr11p_timeStartBusy(activeIDs(ifStartCBR)) = timeManagement.timeNow;

end
%%

% Temporary
% if simParams.mco_nVehInterf>0 && outParams.mco_printInterfStatistic
%     outputValues = mco_updateInterfFromAdjacent(timeManagement,stationManagement,sinrManagement,phyParams,simParams,outputValues);
% end

