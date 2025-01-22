% ==========================================================================
% -- Jammer-Resilient Synchronization in the MIMO Uplink
% --------------------------------------------------------------------------
% -- (c) 2025 Flurin Arquint and Gian Marti
% -- e-mail: farquin@student.ethz.ch and marti@iis.ee.ethz.ch
% --------------------------------------------------------------------------
% -- If you use this simulator or parts of it, then you must cite our paper:
%
% -- Flurin Arquint, Gian Marti, and Christoph Studer
% -- "Jammer-Resilient Synchronization in the MIMO Uplink"
% -- under review (preprint available)
%
% =========================================================================

function res = multi_ant_sync_ROC_simulator(id, par)

rng(id); % Set random seed for reproducibility

par.Es = 1;                                 % signal power
par.N0 = par.Es*10.^(-par.SNRdB/10);        % noise power
par.J0 = par.Es*10.^(par.rhoEdB/10);        % jammer power 

% Initialize result vectors
res.false_positive = zeros(1, length(par.thresholds));
res.false_negative = zeros(1, length(par.thresholds));

if par.isQuadriga
  par.U = 1; % argument needed
  par.J = 1; % argument needed
  if par.use_3GPP 
    loaded_channels = load("channels/3GPP_38.901_UMa_U1_J1_I4_B16_W2048-1_angSepMin1_trials1000_OG_0.mat");
  elseif par.use_mmWave
    loaded_channels = load("channels/mmMAGIC_UMi_LOS_U1_J1_I4_B128_W2048-1_angSepMin1_trials1000_OG_0.mat");
  else
    error("channels not properly specified")
  end
  UE_channels = loaded_channels.UE_channels;
  jammer_channels = loaded_channels.jammer_channels;
  [UE_channels, jammer_channels] = power_control(par, UE_channels, jammer_channels);
end

% Simulation loop
for nn=1:par.trials 
    par.s = (randn(1, par.seq_length)+1i*randn(1, par.seq_length))/sqrt(2); % draw Gaussian synchronization sequence

    % Draw the time index L at which the user equipment (UE) sends the
    % synchronization sequence, from a geometric probability distribution 
    % with expectation E[L]= par.seq_length^2
    L = geornd(1/(par.seq_length^2));
 
    signal_length = L+length(par.s); % Total signal length only as long as necessary
 
    % Generate channel coefficients using an iid Rayleigh fading model
    if ~par.isQuadriga
      h = sqrt(0.5)*(randn(par.B, 1)+1i*randn(par.B, 1));
      J = sqrt(0.5)*(randn(par.B, par.I)+1i*randn(par.B, par.I));
    else
      h = squeeze(UE_channels(mod(nn-1,1000)+1,:,1)).';
      J = reshape(jammer_channels(mod(nn-1,1000)+1,:,:,:), [par.B, par.I]);
    end
    
    n = sqrt(0.5)*(randn(par.B, signal_length)+1i*randn(par.B, signal_length)); % Generate additive white Gaussian noise
    s = zeros(1, signal_length);
    s(1, L+1:signal_length) = par.s; % Generate transmit signal
 
    % Generate jammer signal
    switch (par.jammer)
        case 'None' % No jammer
            w = zeros(par.I, signal_length); 
        case 'Barrage' % Jammer that jams perpetually on all antennas
            w = sqrt(0.5)*(randn(par.I, signal_length)+1i*randn(par.I, signal_length));              
        case 'Off-Sync' % Jammer that jams perpetually on all antennas, except when the synchronization sequence is transmitted
            w = sqrt(0.5)*(randn(par.I, signal_length)+1i*randn(par.I, signal_length)); 
            w(:,s~=0) = 0;                                                       
        case 'On-Sync' % Jammer that jams on all antennas only when the synchronization sequence is transmitted
            w = sqrt(0.5)*(randn(par.I, signal_length)+1i*randn(par.I, signal_length)); 
            w(:,s==0) = 0;                                                        
        case 'Spoofer' % Jammer that spoofs the sync sequence on all antennas as soon as it is transmitted
            w = zeros(par.I, signal_length);
            w(:,s~=0) = ones(par.I,1)*par.s;
        case 'Delayed-Spoofer' % Jammer that spoofs the sync sequence on all antennas with one sample delay
            w = zeros(par.I, signal_length);
            w(:,signal_length-par.seq_length+2:end) = ones(par.I,1)*par.s(1:end-1);
        case 'Erratic' % Jammer that jams/is silent in bursts of random length between 1 and par.seq_length 
            burst_lengths = randi([1,par.seq_length],1,signal_length);
            change_indices = [0,cumsum(burst_lengths)];
            jam = 0;
            w = zeros(par.I, signal_length);
            for ii=1:length(change_indices)-1
              burst_idxs = change_indices(ii)+1:min(change_indices(ii+1),signal_length);
              if jam
                w(:,burst_idxs) = sqrt(0.5)*(randn(par.I, length(burst_idxs))+1i*randn(par.I, length(burst_idxs)));
              end
              jam = ~jam;
              if burst_idxs(end)==signal_length
                break;
              end
            end
      case 'Dynamic-Beamformer' % Jammer that uses a random subset of its antennas for bursts of random length between 1 and par.seq_length 
            assert(par.I>1)
            burst_lengths = randi([1,par.seq_length],1,signal_length);
            change_indices = [0,cumsum(burst_lengths)];
            w = zeros(par.I, signal_length);
            for ii=1:length(change_indices)-1
              burst_idxs = change_indices(ii)+1:min(change_indices(ii+1),signal_length);
              dimension = randi([1, par.I-1]);
              A = [eye(dimension), zeros(dimension, par.I-dimension); zeros(par.I-dimension, par.I)];
              A = A(randperm(par.I),:);
              w(:,burst_idxs) = A*sqrt(0.5)*(randn(par.I, length(burst_idxs))+1i*randn(par.I, length(burst_idxs)));
            end
        otherwise
            error('par.jammer not defined')
    end

    % Generate receive signal
    y = h*s + sqrt(par.J0)*J*w + sqrt(par.N0)*n;

    % Call synchronization algorithm
    switch(par.mitigation)
      case 'None'
        scores = correlation_sync(par,y,L); % No mitigation, just (normalized) correlation with synchronization sequence
      case 'BAJASS' 
        scores = BAJASS(par,y,L); % Null Strongest Dimension using running window before performing (normalized) correlation
      case 'JASS'
        scores = JASS(par,y,L);  % Proposed algorithm
      otherwise
        error('mitigation method not defined')
    end
    for mm=1:length(par.thresholds) % loop over thresholds
      l_hat= min(find(scores>=par.thresholds(mm)*norm(par.s)^2))-1;
      if isempty(l_hat) % Synchronization sequence not detected ("false negative")
        res.false_negative(mm) = res.false_negative(mm)+1.;
      elseif l_hat~=L % Synchronization sequence detected before it was there ("false positive")
        res.false_positive(mm) = res.false_positive(mm)+1.;
      end
    end
end
end


function [UE_channels, jammer_channels] = power_control(par, UE_channels, jammer_channels)
  range_scale = 2; %gives +/- 3dB power control, i.e., the strongest user has 6dB more receive power than the weakest one

  % -- power control of user channels
  powers = vecnorm(UE_channels, 2, 2).^2;
  minpower = min(powers, [], 3);
  user_channels = zeros(size(UE_channels));
  for u=1:par.U
    user_channels(:, :,u) = min(squeeze(powers(:,:,u)), range_scale*minpower)./sqrt(squeeze(powers(:,:,u))).*UE_channels(:,:,u);
  end
  UE_channels = user_channels;
  avg_powers = mean(vecnorm(UE_channels, 2, 2).^2, 3);
  UE_channels = sqrt(par.B)*(UE_channels./sqrt(avg_powers));
  
  % -- power control of jammer channels
  if par.I==1 && par.J==1
    % one-single-antenna jammer
    jammer_channels = sqrt(par.B)*normalize(jammer_channels, 2, 'norm');
  elseif par.I==1 && par.J>1
    % multiple single-antenna jammers
    jammer_channels = sqrt(par.B)*normalize(jammer_channels, 2, 'norm');
  elseif par.I>1 && par.J==1
    %one multi-antenna jammer
    powers = vecnorm(jammer_channels,2,2).^2;
    avg_power_per_jammer = mean(powers,4);
    jammer_channels = sqrt(par.B)*jammer_channels./sqrt(avg_power_per_jammer);    
  else
    % multiple multi-antenna jammers
    powers = vecnorm(jammer_channels,2,2).^2;
    avg_power_per_jammer = mean(powers,4);
    jammer_channels = sqrt(par.B)*jammer_channels./sqrt(avg_power_per_jammer);
  end
end