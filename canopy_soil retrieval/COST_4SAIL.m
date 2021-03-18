function [er,rad,refl,rmse,soil] = COST_4SAIL(p,measurement,input)

soilpar = input{1};
leafbio = input{2};
canopy  = input{3};
angles  = input{5};
optipar = input{6};
t        = input{7};
spectral = input{8};

include = input{10};
%wlrange = input{11};
iparams = input{12};
soilemp = input{13};
soilspec = input{14};
%params0  = input{15};
pcf      = input{16};
prior    = input{17};
I        = input{18};
Esundata = input{19};
Eskydata = input{20};

params = zeros(18,1);
if ~isempty(iparams)
    params(iparams) = p;
end

leafbio.Cab = params(5) * include.Cab + (1-include.Cab)*leafbio.Cab;
leafbio.Cdm = params(6) * include.Cdm + (1-include.Cdm)*leafbio.Cdm;
leafbio.Cw = params(7) * include.Cw + (1-include.Cw)*leafbio.Cw;
leafbio.Cs = params(8) * include.Cs + (1-include.Cs)*leafbio.Cs;
leafbio.Cca = params(9) * include.Cca + (1-include.Cca)*leafbio.Cca;
leafbio.Cant = params(10) * include.Cca + (1-include.Cant)*leafbio.Cant;
leafbio.N = params(11) * include.N + (1-include.N)*leafbio.N;
SIF = zeros(length(spectral.wlS),1);
SIF(640-399:850-399)   = pcf*(params(15:18) * include.ChlF);

if ~include.LAI, params(12) = 1-exp(-0.2*canopy.LAI); end
if ~include.LIDF
    params(13)      = canopy.LIDFa+canopy.LIDFb;
    params(14)      = canopy.LIDFa-canopy.LIDFb;
end

soilpar.B = params(1) * include.B + (1-include.B)*soilpar.B;
soilpar.lat = params(2) * include.lat + (1-include.B)*soilpar.lat;
soilpar.lon = params(3) * include.lon + (1-include.lat)*soilpar.lon;
soilpar.SMp = params(4) * include.SMp + (1-include.lon)*soilpar.SMp;
canopy.LAI      = -5*log(1-params(12));
canopy.LIDFa    = (params(14)+params(13))/2;
canopy.LIDFb    = (params(13)-params(14))/2;
canopy.lidf     = leafangles(canopy.LIDFa,canopy.LIDFb); 

%[leafopt]       = fluspect_bcar(spectral,leafbio,optipar);
leafbio.V2Z = 0;
[leafopt]       = fluspect_B_CX_PSI_PSII_combined(spectral,leafbio,optipar);

IwlP            = spectral.IwlP;
IwlT            = spectral.IwlT;
nwlP            = length(IwlP);
nwlT            = length(IwlT);

[rho,tau,rs]    = deal(zeros(nwlP + nwlT,1));
rho(IwlP)       = leafopt.refl;
tau(IwlP)       = leafopt.tran;
rs(IwlP)        = BSM(soilpar,soilspec,soilemp);
rs(IwlT)        = rs(nwlP) * ones(nwlT,1);
tau(IwlT)       = leafbio.tau_thermal;
rho(IwlT)       = leafbio.rho_thermal;

leafopt.refl    = rho;     % extended wavelength ranges are stored in structures
leafopt.tran    = tau; 
soil.refl       = rs;

%% run the model, calculate the difference between measured and modeled data
rad   = RTMo_lite(soil,leafopt,canopy,angles);

rso     = interp1(spectral.wlS, rad.rso ,spectral.wlM,'splines',1E-4);
rdo     = interp1(spectral.wlS, rad.rdo ,spectral.wlM,'splines',1E-4);
rdd     = interp1(spectral.wlS, rad.rdd ,spectral.wlM,'splines',1E-4);
rsd     = interp1(spectral.wlS, rad.rsd ,spectral.wlM,'splines',1E-4);
SIFs    = interp1(spectral.wlS,SIF ,spectral.wlM,'splines');

t1      = t(:,1);
t3      = t(:,2);
t4      = t(:,3);
t5      = t(:,4);
t12     = t(:,5);

if isempty(Esundata)
    Esun_   = pi*t1.*t4;
else
    Esun_   = interp1(Esundata(:,1), Esundata(:,2),spectral.wlM,'splines');
end
if isempty(Eskydata)
    Esky_   = pi./(1-t3.*rdd).*(t1.*(t5+t12.*rsd));
else
    Esky_   = interp1(Eskydata(:,1), Eskydata(:,2),spectral.wlM,'splines');
end
 

Esun_   = max(nanmean(Esun_)*1E-3,Esun_);      % replace zeros by small values
Esky_   = max(nanmean(Esky_)*1E-3,Esky_);      % replace zeros by small values
piL_      = rso .* Esun_ + rdo .* Esky_ + SIFs;
refl    = piL_./(Esun_ + Esky_);

%reflP   = interp1(spectral.wlM,refl, spectral.wlP,'splines');

refl(spectral.wlM<spectral.wlS(1)) = NaN;
refl(spectral.wlM>spectral.wlS(end)) = NaN;
soil.refl  = interp1(spectral.wlS,soil.refl,spectral.wlM,'splines');
rad.SIF = SIF(640-399:850-399);

%i = find(spectral.wlM>=wlrange.wlmin & spectral.wlM<=wlrange.wlmax & ~isnan(measurement.refl));

er1 = (refl(I) - measurement.refl(I));
er2 = (params - prior.Apm) ./ prior.Aps; 
er = [er1 ; 3E-2* er2];

rmse = sqrt((refl(I) - measurement.refl(I))'*(refl(I) - measurement.refl(I))/length(I));

