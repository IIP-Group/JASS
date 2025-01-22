addpath('./funcs')
addpath('./channels')

id = 0;

par.B = 16;                 % number of antennas at BS
par.I = 4;                  % number of antennas at the jammer
par.I_est = 4;              % assumed number of antennas at the jammer
par.seq_length = 16;        %length of synchronization sequence
par.trials = 1e2;           % number of Monte Carlo trials
par.num_power_its = 4;      % number of power iterations when calculating eigenvectors
par.thresholds = 0:0.005:1; % list of thresholds
par.isQuadriga = 0;         % use Rayleigh fading channel model
synchronizer_dict = {"None","BAJASS","JASS"};

par.jammer = "Barrage";
par.rhoEdB = 30; % receive energy ratio

% start with SNR=-10dB
par.SNRdB = -10; % SNR (dB)

fprs_at_minus10dB = zeros(length(synchronizer_dict),length(par.thresholds));
fnrs_at_minus10dB = zeros(length(synchronizer_dict),length(par.thresholds));
for ii=1:length(synchronizer_dict)
  par.mitigation = synchronizer_dict{ii};
  res = multi_ant_sync_ROC_simulator(id,par);
  fprs_at_minus10dB(ii,:) = res.false_positive/par.trials;
  fnrs_at_minus10dB(ii,:) = res.false_negative/par.trials;
end

% now with SNR=0dB
par.SNRdB = 0; % SNR (dB)
fprs_at_0dB = zeros(length(synchronizer_dict),length(par.thresholds));
fnrs_at_0dB = zeros(length(synchronizer_dict),length(par.thresholds));
for ii=1:length(synchronizer_dict)
  par.mitigation = synchronizer_dict{ii};
  res = multi_ant_sync_ROC_simulator(id,par);
  fprs_at_0dB(ii,:) = res.false_positive/par.trials;
  fnrs_at_0dB(ii,:) = res.false_negative/par.trials;
end

% plot results
title_str = "B="+num2str(par.B)+", I="+num2str(par.I)+", Iest="+num2str(par.I_est)+", seqlength="+num2str(par.seq_length)+", numpowerits="+num2str(par.num_power_its)+", jammer="+num2str(par.jammer)+", rhoEdB="+num2str(par.rhoEdB)+", trials="+num2str(par.trials);
figure(1)
xlim([0 1])
ylim([0 1])
colors = get(gca,'colororder');
for ii=1:length(synchronizer_dict)
  plot(fprs_at_minus10dB(ii,:),fnrs_at_minus10dB(ii,:),'LineWidth',2,'LineStyle','-','Color',colors(ii,:),'DisplayName',synchronizer_dict{ii}+" @SNR=-10dB");
  if ii==1
    hold on
  end
  legend
end
for ii=1:length(synchronizer_dict)
  plot(fprs_at_0dB(ii,:),fnrs_at_0dB(ii,:),'LineWidth',2,'LineStyle','none','Marker','*','MarkerEdgeColor',colors(ii,:),'DisplayName',synchronizer_dict{ii}+" @SNR=0dB");
end
xlabel("false positive ratio (FPR)")
ylabel("false negative ratio (FNR)")
title(title_str)
grid on 
hold off
set(gcf,'position',[10,10,400,400])


figure(2)
xlim([0 1])
ylim([1e-4 1])
colors = get(gca,'colororder');
for ii=1:length(synchronizer_dict)
  semilogy(par.thresholds,fprs_at_minus10dB(ii,:)+fnrs_at_minus10dB(ii,:),'LineWidth',2,'LineStyle','-','Color',colors(ii,:),'DisplayName',synchronizer_dict{ii}+" @SNR=-10dB");
  if ii==1
    hold on
  end
  legend
end
for ii=1:length(synchronizer_dict)
  semilogy(par.thresholds,fprs_at_0dB(ii,:)+fnrs_at_0dB(ii,:),'LineWidth',2,'LineStyle','none','Marker','*','MarkerEdgeColor',colors(ii,:),'DisplayName',synchronizer_dict{ii}+" @SNR=0dB");
end
xlabel("Treshold $\tau$",Interpreter="latex")
ylabel("total error rate (TER)")
title(title_str)
grid on 
hold off
set(gcf,'position',[10,10,400,400])
