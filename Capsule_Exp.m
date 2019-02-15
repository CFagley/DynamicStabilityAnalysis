function d = Capsule_Exp(d)

d=d;

[~,~,ext]=fileparts(d.file);

if(strcmp(ext,'.tdms'))
    d.raw = readTDMSFile(d.file);
else
    d.raw = load(d.file);
end
vars = fields(d.raw);

% Double everything
for i = 1:length(vars)
    d.data.(vars{i}) = double(d.raw.(vars{i}).Data);
end

if isfield(d.raw,'Position')
    d,d.time = extract_params(d,d.raw.Position);
    pos = d.raw.Position.Data;
else
     d = extract_params(d,d.raw.Actual_Position);
    pos = d.raw.Actual_Position.Data;
end

% Set Geometric parameters
 
% Filter everything based on forcing frequency
 d = filter_data(d);
 
 
 % System Identification /// ADD SYS ID FOR FORCED OSCILLATION TESTS
 if isfield(d.raw,'Disturbance_Torque')
    d = system_ID(d);
 end

  
% Compute tunnel velocity
[d.tunnel, d.data.Q]= Compute_Velocity(d);
% d.params.C = 5*.0254;
d.params.S = d.params.C^2/4*pi;
d.params.k = 1/(d.tunnel.U_Infinity/2/pi/d.params.f/d.params.C);

 
 % Compute inertial loads
% d = compute_inert_loads(d);

% Compute aerodynamic torque
d = Compute_Aero_Coeffs(d);
d = Frequency_Info(d);
d = Bode_info(d);

%Phase Average
var = {'Cm','aero_torque','Q','Velocity'};
for i  = 1:length(var)
    temp = [var{i} '_p'];
    temp2 = [var{i} '_n'];
%     [d.phase.theta_p, d.phase.(temp), d.phase.theta_n, d.phase.(temp2)]=Bin_Phase_Average(d.data.Actual_Position,d.data.(var{i}),d.phase.nbins);
    [d.phase.theta_p, d.phase.(temp), ~,~,d.phase.([temp '_std']),~]=Bin_Phase_Average(pos,d.data.(var{i}),d.phase.nbins);
    [~,~,d.phase.theta_n, d.phase.(temp2),~,d.phase.([temp2 '_std'])]=Bin_Phase_Average(pos,d.data.(var{i}),d.phase.nbins);
    d.phase.theta = [d.phase.theta_n d.phase.theta_p];
    d.phase.(var{i}) = [d.phase.(temp2) d.phase.(temp)];
    d.phase.([var{i} '_std']) = [d.phase.([temp2 '_std']) d.phase.([temp '_std'])];
end
% [d.phase.theta_p, d.phase.Cm_p, d.phase.theta_n, d.phase.Cm_n]=Bin_Phase_Average(d.data.Actual_Position,d.data.Cm,d.phase.nbins);

 function [r, Q] = Compute_Velocity(in)

volt_to_torr = 20;
torr_to_pa = 133.32;
c_to_kel = 273.2;
psi_to_pa = 6894.75729;
Q = volt_to_torr*torr_to_pa*in.data.Tunnel_Pressure;
R = 287.058; %J/kg/k
r.Rho = in.cond.Patm*psi_to_pa/R/(in.cond.Temp+c_to_kel);
if mean(Q)<0
    Q=-Q;
end
r.U_Infinity = sqrt(2/r.Rho*mean(Q));
r.Dyn_Pressure = mean(Q);
r.c = sqrt(1.4*8.414/.0289645*(in.cond.Temp+c_to_kel));
r.Mach = r.U_Infinity/r.c;


function d = filter_data(d)

scale = d.filt.scale;
d.filt.Wn = scale;
% d.filt.Wn = scale * d.params.f/d.time.fs/2;

[d.filt.B,d.filt.A] = butter(7,d.filt.Wn,'low');
vars = fields(d.data);
for i = 1:length(vars)
    d.data.([vars{i} '_filtered']) = filtfilt(d.filt.B,d.filt.A,d.data.(vars{i}));
end

% function d = compute_inert_loads(d)
% if isfield(d.data,'Actual_Position_filtered')
%     d.struct.dpdt = gradient(d.data.Actual_Position_filtered,d.time.dt);
% else
%     d.struct.dpdt = gradient(d.data.Position_filtered,d.time.dt);
% end
% d.struct.dp2dt2 = gradient(d.struct.dpdt,d.time.dt);
% d.struct.interial = d.struct.G*d.struct.dp2dt2*pi/180;
% d.struct.friction = d.struct.D*d.struct.dpdt*pi/180;
% d.struct.total = d.struct.interial+d.struct.friction;


function d = Compute_Aero_Coeffs(d)

% Calibrate torque value
% d.data.struct_torque = d.struct.total;
d.data.total_torque = (d.data.Torque_Cell_filtered-d.devices.torque.DC)*d.devices.torque.cal;

alpha = deg2rad(d.data.Actual_Position_filtered);

adot = gradient(alpha,d.time.dt);
addot = gradient(adot,d.time.dt);
X = [alpha;adot;addot];
d.data.Velocity = adot;
d.data.Accelleration = addot;

if d.sysID.Aoff_coeffs(2)==0
    X = [0*alpha;adot;addot];
end

if d.sysID.Aoff_coeffs(3)~=0
     X = [alpha;adot;addot];
end

ind = randi(length(d.data.total_torque),1,length(d.data.total_torque));
% ind = 1:length(d.data.total_torque);

d.sysID.Total_coeffs = d.data.total_torque(ind)/X(:,ind);

% d.sysID.Total_coeffs = Total_coeffs;
d.sysID.Aon_coeffs = d.sysID.Total_coeffs-d.sysID.Aoff_coeffs;

d.data.str_torque = d.sysID.Aoff_coeffs*X;
d.data.total_torque_est = d.sysID.Total_coeffs*X;
d.data.aero_torque = d.sysID.Aon_coeffs*X;
% d.sysID.error = norm(d.data.total_torque-d.data.total_torque_est)/norm(d.data.total_torque);
d.sysID.error = rms(d.data.total_torque-d.data.total_torque_est)/rms(d.data.total_torque);


if d.sysID.Aoff_coeffs(2)~=0
    d.data.Cm = d.data.aero_torque/d.tunnel.Dyn_Pressure/d.params.S/d.params.C;     
else
    d.data.Cm = d.data.total_torque/d.tunnel.Dyn_Pressure/d.params.S/d.params.C;     
end

d.sysID.a0 =  rad2deg(mean(alpha));
d.sysID.Cma = d.sysID.Aon_coeffs(1)*deg2rad(d.sysID.a0)/(d.tunnel.Dyn_Pressure*d.params.C^3*pi/4);
d.sysID.Cmq = 2*d.sysID.Aon_coeffs(2)*d.tunnel.U_Infinity/(d.tunnel.Dyn_Pressure*d.params.C^4*pi/4);
% d.sysID.Cmq = -2*d.tunnel.U_Infinity/max(d.time.time)/d.tunnel.Dyn_Pressure/(d.params.C^4*pi/4)*trapz(d.time.time,d.data.aero_torque.*adot);
% d.sysID.Cmq = -2*d.tunnel.U_Infinity/d.tunnel.Dyn_Pressure/(d.params.C^4*pi/4)*trapz(d.time.time,d.data.aero_torque.*adot);
% d.sysID.Cmq = -d.data.Cm/(adot*d.params.C/2/d.tunnel.U_Infinity);

% d.sysID.sys = tf(1,d.sysID.Total_coeffs([3,2,1]));
% Find Aero Torque
% Figure alpha
% plot(adot*d.params.C/2/d.tunnel.U_Infinity,d.data.Cm)
function d = Frequency_Info(d)

fields = fieldnames(d.data);
for i=1:length(fields)
    n=length(d.data.(fields{i}));
    d.FFT.([fields{i} '_cmplx'])=fft(d.data.(fields{i}))/n;
    temp = 2*abs(d.FFT.([fields{i} '_cmplx']));
    d.FFT.(fields{i})=temp(1:n/2+1);
    d.FFT.([fields{i} '_ind']) = d.FFT.(fields{i})(d.params.ind);
   
%     blah = (d.data.(fields{i})-mean(d.data.(fields{i})));
%     d.Hilb.([fields{i} 'unwrap']) = unwrap(angle(hilbert(blah)));
%     d.Hilb.([fields{i} ]) = (angle(hilbert(blah)));

end
d.FFT.Omega = d.time.fs/2*linspace(0,1,n/2+1);
d.FFT.K = 1./(d.tunnel.U_Infinity./2/pi./d.FFT.Omega./d.params.C);

function d = Bode_info(d)
    
    d.bode.Gain = 20*log10(abs(deg2rad(d.FFT.Actual_Position_filtered)./d.FFT.total_torque));
    d.bode.Gain_est = 20*log10(abs(deg2rad(d.FFT.Actual_Position_filtered)./d.FFT.total_torque_est));
    d.bode.Omega = d.FFT.Omega;
    d.bode.error = norm(d.FFT.total_torque-d.FFT.total_torque_est)/norm(d.FFT.total_torque);
    d.sysID.error_freq = d.bode.error;


function    d = extract_params(d,Position)

    d=d;
pos = double(Position.Data);
time = double(Position.Axis);
dt = time(2) - time(1);
fs = 1/dt;

d.params.a0 = mean(pos);
% a1 = (max(pos)-min(pos))/2;
temp = findpeaks(pos,'MinPeakHeight',d.params.a0);
d.params.a_max = mean(temp);
[temp] = -findpeaks(-pos,'MinPeakHeight',-d.params.a0);
d.params.a_min = mean(temp);
d.params.a1 = (d.params.a_max-d.params.a_min)/2;

n=length(pos);
pos_fft = 2*abs(fft(pos)/n);
pos_fft = pos_fft(1:end/2+1);
freq = fs/2*linspace(0,1,n/2+1);

[~,ind]=findpeaks(pos_fft(1:end),'NPeaks',1,'MinPeakHeight',.1);
if isempty(ind)
    ind = 5;
end
d.params.ind=ind+1;
d.params.f = freq(ind);
d.time.time = time;
d.time.dt = dt;
d.time.fs = fs;

function [xp,yp,xn,yn,ypstd,ynstd] = Bin_Phase_Average(position,input,N_bins)


position = double(position);
Min = min(position);
Max = max(position);
delta = (Max-Min)/(N_bins-2);
x = Min-delta/2:delta:Max+delta/2;
win = 100;
posfilt = filtfilt(1/win*ones(1,win),1,position);
dpdt = gradient(posfilt);
n = length(x)-1;

xn = zeros(1,n);
y = xn;
ystd = xn;
% Figure blah
% plot(position,input)
for i=1:n
    locs = position>x(i)&position<x(i+1)&dpdt>0;
%     Figure blah
%     plot(position(locs),input(locs),'ro')
    xn(i) = (x(i)+x(i+1))/2;
    y(i) = mean(input(locs));
    ystd(i) = std(input(locs));
end

xp = xn;
yp = y;
ypstd = ystd;

xn = zeros(1,n);

for i=1:n
    x2 = fliplr(x);
    locs = position<x2(i)&position>x2(i+1)&dpdt<=0;
    xn(i) = (x2(i)+x2(i+1))/2;
    y(i) = mean(input(locs));
    ystd(i) = std(input(locs));
end

yn = y;
ynstd = ystd;

function d = system_ID(d)

% Motor_Kt      =27;  %ozin/amp %%%% ENSURE THIS NUMBER IS CORRECT %%%
% Conv_Oz_Nm    =0.007061552;  % Nm/ozin
% d.sys.Input_Data=d.data.Disturbance_Torque'*Motor_Kt*Conv_Oz_Nm;


Input_Data=d.data.Disturbance_Torque'*0.38792/2;
Output_Data=double(-d.data.Position_filtered')*pi/180;

ind = find(diff(d.data.Disturbance_Torque)~=0);

% data = iddata(d.sys.Output_Data,d.sys.Input_Data,d.time.dt);
% count = 1;
% count2 = 1;
for i = 1:length(ind)-1
    if i==1
        shift = ind(1);
     else
         shift = ceil((ind(2)-ind(1))/3);
    end
    input{i} = Input_Data(ind(i)-shift+1:ind(i+1));
    output{i} = Output_Data(ind(i)-shift+1:ind(i+1));    
end

% insert case selection
switch d.sys.type
    case 'Beginning'
        display('Computing Sys ID for Beginning transient')
        sub_input = input(1:2:length(input));
        sub_output = output(1:2:length(output));
    case 'Ending'
        display('Computing Sys ID for Ending transient')
        sub_input = input(2:2:length(input));
        sub_output = output(2:2:length(output));
    case 'All'
        display('Computing Sys ID for All transient')
        sub_input = input(1:length(input));
        sub_output = output(1:length(output));
end

switch d.sys.sets
    case 'Total'
        
        data = iddata(sub_output,sub_input,d.time.dt);
        init_sys=idproc('p2U');
%         disp('Computing System Dynamics by Step Input')
        sys = pem(data,init_sys);
        
        for i=1:length(sub_input)
            d.sys.time{i} = d.data.Time(1:length(sub_input{i}));
            d.sys.est_output{i} = lsim(sys,sub_input{i},d.sys.time{i},[mean(sub_output{i}(1:10))/2/pi,0]);
        end
        
        d.sys.input = sub_input;
        d.sys.output = sub_output;
        
        % d.sys.bt_est = lsim(d.sys.bt_sys,bt_input{1},d.data.Time(1:length(bt_input{1})));
        % d.sys.et_est = lsim(d.sys.et_sys,et_input{1},d.data.Time(1:length(et_input{1})));
        
        d.sys.G = sys.Tw^2/sys.Kp;  %kg-m^2/rad
        d.sys.D = 2*sys.Tw*sys.Zeta/sys.Kp;  %Nm*s/rad
        d.sys.K = 1/sys.Kp;
        d.sys.Zeta = d.sys.D/(2*sqrt(d.sys.K*d.sys.G));
        d.sys.system=sys;
        
    case 'Each'
        display('Computing Sys ID for EACH transient')
        for i =1:length(sub_input)
            display(['Computing transient # ' num2str(i) ' of ' num2str(length(sub_input))])
            data = iddata(sub_output{i},sub_input{i},d.time.dt);
            init_sys=idproc('p2U');
            sys = pem(data,init_sys);
            d.sys.time{i} = d.data.Time(1:length(sub_input{i}));
            d.sys.est_output{i} = lsim(sys,sub_input{i},d.sys.time{i},[mean(sub_output{i}(1:100))/2/pi,0]);
            
          
            d.sys.G{i} = sys.Tw^2/sys.Kp;  %kg-m^2/rad
            d.sys.D{i} = 2*sys.Tw*sys.Zeta/sys.Kp;  %Nm*s/rad
            d.sys.K{i} = 1/sys.Kp;
            d.sys.Zeta{i} = d.sys.D{i}/(2*sqrt(d.sys.K{i}*d.sys.G{i}));
            d.sys.system{i}=sys;
        end
        d.sys.input = sub_input;
        d.sys.output=sub_output;
end







