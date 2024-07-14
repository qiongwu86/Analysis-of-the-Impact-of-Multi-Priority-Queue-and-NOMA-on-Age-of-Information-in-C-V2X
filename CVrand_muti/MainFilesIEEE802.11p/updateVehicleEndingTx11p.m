function [timeManagement,stationManagement,sinrManagement] = updateVehicleEndingTx11p(idEvent,indexEvent,timeManagement,stationManagement,sinrManagement,phyParams)
% A transmission is concluded in IEEE 802.11p

% If the vehicle is exiting the scenario, indexEvent is set to -1
% thus shouldn't call this function
if indexEvent<=0
    error('Call to updateVehicleEndingTx11p with indexEvent<=0');
end

% The number of packets in the queue is reduced
stationManagement.pckBuffer(idEvent) = stationManagement.pckBuffer(idEvent)-1;
% The medium is sensed to check if it is free
% (note: 'vState(idEvent)' is set to 9 in order not to contribute
% to the sensed power)
stationManagement.vehicleState(idEvent)=9; % rx
if ((stationManagement.vehicleState(stationManagement.activeIDs)==3) * sinrManagement.P_RX_MHz(indexEvent,:)) > (phyParams.PrxSensNotSynch/phyParams.BwMHz)
    % If it is busy, State 9 with error is entered
    stationManagement.vehicleState(idEvent)=9; % rx
    sinrManagement.idFromWhichRx11p(idEvent) = idEvent;
    sinrManagement.interfAverage11p(idEvent) = 0;
    sinrManagement.instantThisSINRavStarted11p(idEvent) = timeManagement.timeNow;
    timeManagement.timeNextTxRx11p(idEvent) = Inf;
else
    % If it is free, then: the idle state is entered if the queue is empty
    % otherwise a new backoff is started
    if stationManagement.pckBuffer(idEvent)==0
        % If no other packets are in the queue, the node goes
        % in idle state
        stationManagement.vehicleState(idEvent)=1; % idle
        timeManagement.timeNextTxRx11p(idEvent) = Inf;
    elseif stationManagement.pckBuffer(idEvent)>=1
        % If there are other packets in the queue, a new
        % backoff is initialized and started
        stationManagement.vehicleState(idEvent)=2; % backoff
        if simParams.technology~=4 || simParams.coexMethod~=3 || ~simParams.coexCmodifiedCW
            [stationManagement.nSlotBackoff11p(idEvent), timeManagement.timeNextTxRx11p(idEvent)] = startNewBackoff11p(timeManagement.timeNow,stationManagement.CW_11p(idEvent),stationManagement.tAifs_11p(idEvent),phyParams.tSlot);
        else
            relativeTime = timeManagement.timeNow-simParams.coex_superFlength*floor(timeManagement.timeNow/simParams.coex_superFlength);
            subframeIndex = floor(relativeTime/phyParams.Tsf);
            [stationManagement.nSlotBackoff11p(idEvent), timeManagement.timeNextTxRx11p(idEvent)] = coexistenceStartNewBackoff11pModified(timeManagement.timeNow,stationManagement.CW_11p(idEvent),stationManagement.tAifs_11p(idEvent),phyParams.tSlot,subframeIndex,simParams.coex_superframeSF);
        end        
    else
        error('Error: pckBuffer<0');
    end
    %% Not anymore used since version 5.2.9
    %% In both cases, the channel busy ratio is updated 
    %if ~isempty(stationManagement.channelSensedBusyMatrix11p)
    %    stationManagement.channelSensedBusyMatrix11p(1,idEvent) = stationManagement.channelSensedBusyMatrix11p(1,idEvent) + (timeManagement.timeNow-timeManagement.cbr11p_timeStartBusy(idEvent));
    %    timeManagement.cbr11p_timeStartBusy(idEvent) = -1;
    %end
end

