function [refl,std,wl,c,include,outdirname,leafbio,angles,canopy,soilpar,wlrange,Esundata,Eskydata,std_provided,iparams,pcf,atmfile,sensor] = input_data

[d, txt] = xlsread('input_data.xlsx');
d = d(~isnan(d));

PCflu = xlsread('..\data\input\PC_flu.xlsx');
pcf   = PCflu(2:end,2:5);

%% Tell me, where do I find your measured soil reflectance spectrum and measured incident light
reflfilename = [char(txt(2,2)),'\',char(txt(3,2))];
stdfilename  = char(txt(4,2));
wlfilename   = [char(txt(2,2)),'\',char(txt(5,2))];
outdirname   = [char(txt(6,2)),'\',char(txt(7,2))];

%% atmospheric file
atmfile      = char(txt(10,2));
Esunfile     = char(txt(8,2));
Eskyfile     = char(txt(9,2));

wl           = load(wlfilename);
refl         = load(reflfilename);%load('..\data\measured\Reflectance.txt');

if ~isempty(Esunfile)
    Esundata     = load(Esunfile);
else
    Esundata = [];
end

if ~isempty(Eskyfile)
    Eskydata     = load(Eskyfile);
else
    Eskydata = [];
end
%rsdata       = load(soilfile);

if ~isempty(stdfilename),   
    measured.stdmeas = load([char(txt(2,2)),'\' stdfilename]); 
    std_provided = 1; 
    std = measured.stdmeas;
else std = ones(length(wl),size(refl,2)); std_provided = 0; 
end


%% which parameters should I tune?
c               = d(1);

%%
include.B       = d(2);
include.lat     = d(3);
include.lon     = d(4);
include.SMp     = d(5);
include.Cab     = d(6);
include.Cdm     = d(7);
include.Cw      = d(8);
include.Cs      = d(9);
include.Cca     = d(10);
include.Cant    = d(11);
include.N       = d(12);
include.LAI     = d(13);
include.LIDF    = d(14);
include.ChlF    = d(15);

ip              = d(2:15);%find(d(2:14)>0);
ip1             = [ip(1:13); ip(13)>0*14; (ip(14)>0)*15; (ip(14)>0)*16; (ip(14)>0)*17; (ip(14))>0*18];
iparams         = find(ip1>0);

%% which spectral region should I calibrate?
wlmin           = d(16);          % starting wavelength (nm)
wlmax           = d(17);         % ending wavelength (nm)

%% initialize parameters for retrieval 
% these will be calibrated to your reflectance data if you said so above
soilpar.B       = d(18);
soilpar.lat     = d(19);
soilpar.lon     = d(20);
soilpar.SMp     = d(21);
leafbio.Cab     = d(22);           % chlorophyll content               [ug cm-2]
leafbio.Cdm     = d(23);        % dry matter content                [g cm-2]
leafbio.Cw      = d(24);        % leaf water thickness equivalent   [cm]
leafbio.Cs      = d(25);          % senescent material                [fraction]
leafbio.Cca     = d(26);            % carotenoids                       [?]
leafbio.Cant    = d(27);            % anthocyanins                       [?]
leafbio.N       = d(28);          % leaf structure parameter (affects the ratio of refl: transmittance) []
canopy.LAI      = d(29);
canopy.LIDFa    = d(30);
canopy.LIDFb    = d(31);

%% The following parameters will not be tuned to the data
% but I need to know what their values were
angles.tts      = d(32);      % solar zenith angle (deg)
angles.tto      = d(33);           % observation zenith angle (deg)
angles.psi      = d(34);            % absolute value (deg) of the difference between solar and zenith azimuth angle (arbitrary if angles.tto = 0)
canopy.hot      = d(35);         % ratio between leaf width: canopy height. This is the hot spot parameter. 
                                % you were probably avoiding the hotspot,
                                % so in that case the model is insensitive
                                % to canopy.hot
wlrange.wlmin     = wlmin;
wlrange.wlmax     = wlmax;

%% Characteristics of the instrument
sensor.FWHM  = d(36);
