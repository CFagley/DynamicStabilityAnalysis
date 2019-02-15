function [data o] = run_multiple_datasets(in)


%% Input parameters
% in.base_dir   %% Base directory location
% in.fprefix    %% file name prefix
% in.files      %% Array of files
% in.aoff       %% Aoff file indec
% in.dcfile     %% DC File index



%% Default parameters
d.filt.scale = .2;
d.cond.Temp = 20;    
d.cond.Patm = 11.3;
d.devices.torque.cal = 413.73;  % Nm / Volt
% d.devices.torque.cal = 613.73;
% d.devices.torque.cal = in.torque;
d.devices.torque.DC  = 0;
d.sysID.Aoff_coeffs = [0 0 0];
d.phase.nbins = 101;
d.params.C = in.diam*.0254;

%% Find DC Values

if isfield(in,'dcfile')
    d.file = [in.base_dir filesep in.dcfile '.tdms'];
    Baseline = Capsule_Exp(d);
    d.devices.torque.DC = mean(Baseline.data.Torque_Cell);
    display([' Offset Torque ' num2str(mean(Baseline.data.Torque_Cell)) ] )
end


%% Find Aoff coeffs

if isfield(in,'aoff')
    d.file = [in.base_dir filesep in.aoff '.tdms'];
    Aoff = Capsule_Exp(d);
    d.sysID.Aoff_coeffs = Aoff.sysID.Total_coeffs;
    display([' Air Off coefficients ' num2str(d.sysID.Aoff_coeffs) ] )
end



%%  Loop over AOn files

for i = 1:length(in.files)
    d.file = [in.base_dir filesep in.fprefix num2str(in.files(i),'%03d') '.tdms'];
    display(['Processing run ' num2str((i)) ' of ' num2str(length(in.files))])
    data(i) = Capsule_Exp(d);
    sysid(i) = data(i).sysID;
end

o.sysid=sysid;
o.sum.a0=[sysid.a0];
o.sum.Cma=[sysid.Cma];
o.sum.Cmq=[sysid.Cmq];
o.sum.Coeffs = [sysid.Aon_coeffs];

o.sum.error=[sysid.error];
o.sum.error_freq=[sysid.error_freq];

% o.sum.a0_resamp = -75:1:75;
% [~,ind]=unique(o.sum.a0);
% o.sum.Cmq_resamp = interp1(o.sum.a0(ind),o.sum.Cmq(ind),o.sum.a0_resamp,'pchip','extrap');
% o.sum.Cma_resamp = interp1(o.sum.a0(ind),o.sum.Cma(ind),o.sum.a0_resamp,'pchip','extrap');
% o.sum.error_resamp = interp1(o.sum.a0(ind),o.sum.error(ind),o.sum.a0_resamp,'pchip','extrap');
% o.sum.error_freq_resamp = interp1(o.sum.a0(ind),o.sum.error_freq(ind),o.sum.a0_resamp,'pchip','extrap');
% 




 
