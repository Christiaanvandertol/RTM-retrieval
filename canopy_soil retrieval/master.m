%function [er,params0,params1,spectral,leafbio,canopy,rad,reflSAIL_all,J,fluorescence]= master
%% start fresh
close all
progressbar(0)

%% get constants
global constants
[constants] = define_constants();

%% get the data we need
[refl,std,wl,c,include,outdirname,leafbio,angles,canopy,soilpar,wlrange,Esundata,Eskydata,std_provided,iparams,pcf,atmfile,sensor] = input_data;

if c==-999; c = (1:size(refl,2));
else
    if min(c)<0, fprintf('%s \r', 'negative incides are not allowed (except -999 for all spectra), program ends'), return, end
    if max(c)>size(refl,2), fprintf('%s \r', ['you only have ' num2str(size(refl,2)) ' spectra, not ' num2str(max(c)) ' in your measurement file, program ends']), return, end
end

measured.refl = refl;
measured.std = std;

%% other inputs, not important for reflectance, but necessary to run the model
leafbio.rho_thermal   = 0.01;               % reflectance in the thermal range
leafbio.tau_thermal  = 0.01;                % transmittance in the thermal range
leafbio.fqe(2)  = 0.02;                     % quantum yield
leafbio.fqe(1)  = 0.02/5;
meteo.Rin       = -999;
meteo.Ta        = 20;
options.calc_vert_profiles  = 0;    % calculate vertical profiles of fluxes

%% the input for FLUSPECT and BSM
load '../data/input/fluspect_data/Optipar2017_ProspectD'
% O = xlsread('../data/input/fluspect_data/Optipar_2015.xlsx');
% optipar.nr     = O(:,2);
% optipar.Kab    = O(:,3);
% optipar.Kca    = O(:,4);
% optipar.Ks     = O(:,5);
% optipar.Kw     = O(:,6);
% optipar.Kdm    = O(:,7);
% optipar.nw     = O(:,8);
% optipar.phiI   = O(:,9);
% optipar.phiII  = O(:,10);
% optipar.GSV    = O(:,11:13);

soilemp.SMC   = 25;        % soil moisture content
soilemp.film  = 0.015;     % water film optical thickness

soilspec.wl  = wl;
soilspec.GSV = optipar.GSV;
soilspec.kw  = optipar.Kw;
soilspec.nw  = optipar.nw;

%% Define canopy structure
canopy.nlayers  = 60;
nl              = canopy.nlayers;
canopy.x        = (-1/nl : -1/nl : -1)';         % a column vector
canopy.xl       = [0; canopy.x];                 % add top level
canopy.nlincl   = 13;
canopy.nlazi    = 36;
canopy.litab    = [ 5:10:75 81:2:89 ]';   % a column, never change the angles unless 'ladgen' is also adapted
canopy.lazitab  = ( 5:10:355 );           % a row

%% Define spectral regions
spectral = define_bands;
spectral.wlM = wl;%refl(:,1);
nwlP = length(spectral.wlP);
nwlT = length(spectral.wlT);
spectral.IwlP = 1 : nwlP;
spectral.IwlT = nwlP+1 : nwlP+nwlT;
%spectral.wlRs = rsdata(:,1);
%spectral.wlEsun = Esundata(:,1);
%spectral.wlEsky = Eskydata(:,1);
measured.stdP = interp1(spectral.wlM,measured.std,spectral.wlP);
measured.stdP(isnan(measured.stdP)) = 0;
I1 = find(spectral.wlM>=wlrange.wlmin & spectral.wlM<=wlrange.wlmax );

%% define output structures
if length(c)>1
    [rmse_all,leafbio_all.Cab,leafbio_all.Cw,leafbio_all.Cdm, leafbio_all.Cs,leafbio_all.Cca,leafbio_all.N,leafopt_all.er2,canopy_all.LAI,canopy_all.LIDFa,canopy_all.LIDFb,soilpar_all.B,soilpar_all.lat,soilpar_all.lon,soilpar_all.SMp] = deal(zeros(length(c),1));
    leafbio_all.std = zeros(length(c),17);
    [reflSAIL_all,rsoil_all] = deal(zeros(length(spectral.wlM),length(c)));
    [SIF_all,SIFnorm_all] = deal(zeros(length(spectral.wlF),length(c)));
    fluorescence_all.wpcf = zeros(4,length(c));
    J_all = zeros(length(I1),17,length(c));
end

%% atmospheric data
if ~isempty(atmfile)
    t       = specbin(atmfile,spectral.wlM,sensor.FWHM);
    tS      = specbin(atmfile,spectral.wlS,1);
else 
    [t,tS ] = deal(ones(1,5));  
    %soilfile     = char(txt(11,2)); 
end

%% initial parameter values, and boundary boxes
params0         = ones(18,1);
params0(5)      = leafbio.Cab;
params0(6)      = leafbio.Cdm;
params0(7)      = leafbio.Cw;
params0(8)      = leafbio.Cs;
params0(9)      = leafbio.Cca;
params0(10)      = leafbio.Cant;
params0(11)      = leafbio.N;
params0(12)      = 1-exp(-0.2*canopy.LAI);
params0(13)      = canopy.LIDFa+canopy.LIDFb;
params0(14)      = canopy.LIDFa-canopy.LIDFb;

params0(1)      = soilpar.B;
params0(2)      = soilpar.lat;
params0(3)      = soilpar.lon;
params0(4)      = soilpar.SMp;

wpcf            = zeros(4,1);
params0(15:18)     = wpcf;

%params0(10)     = fluor.PSI;
%params0(11)     = fluor.PSII;

LB              = [0   20 40 5    0    0.00   0     0     0 0 1     0    -1  -1   0 -20 -10 -2]'; % lower boundaries
UB              = [0.9 40 60 55  100   0.02   0.2  0.4  25 40 3.5  .7     1   1  40  20  10  2]'; % upper boundaries

%Apm             = [0.5 25 45 30  40  0.005   0.02   0.1  8 0  1.5  .4     0  0  0   0   0     0    ]; % A priori estimates
Apm             = params0;
Aps             = [0.3 12  9 12  30  0.006   0.02   0.2  4  12 0.75 .2    .6 .6  1E9  1E9  1E9   1E9  ]'; % Uncertainties in stdev

prior.Apm = Apm;
prior.Aps = Aps;

%% do the job: fit the model to the data
progressbar(20/1000)

for j = c
    measurement.refl            = measured.refl(:,j);
    measurement.std             = measured.std(:,j);
    measurement.stdP            = measured.stdP(:,j);
    I2                          = find(spectral.wlM>=wlrange.wlmin & spectral.wlM<=wlrange.wlmax & ~isnan(measurement.refl));
    input           = {soilpar,leafbio,canopy,meteo,angles,optipar,t,spectral,options,include,wlrange,iparams,soilemp,soilspec,params0,pcf,prior,I2,Esundata,Eskydata};
    stoptol = 1E-3;
    opt = optimset('MaxIter',30,'TolFun',stoptol);
    params1 = params0;
    if sum(structfun(@(x)x, include)) >0
        f               = @(params)COST_4SAIL(params,measurement,input);
        tic
        
        [paramsout,~,~,~,~,~,J]= lsqnonlin(f,params0(iparams),LB(iparams),UB(iparams),opt);
        params1(iparams) = paramsout;
        toc

    else
        paramsout = params0;
    end
    
    %% done, this is what comes out
    
    [er,rad,reflSAIL,rmse,soil] = COST_4SAIL(paramsout,measurement,input);
    
    include_all.Cab = 1;
    include_all.Cdm = 1;
    include_all.Cw = 1;
    include_all.Cs = 1;
    include_all.Cca = 1;
    include_all.Cant = 1;
    include_all.N = 1;
    include_all.LAI = 1;
    include_all.LIDF = 1;
    include_all.B = 1;
    include_all.lat = 1;
    include_all.lon = 1;
    include_all.SMp = 1;
    include_all.ChlF = 1;
    
    input{10} = include_all;
    input{12} = (1:18);
    if ~std_provided
        measurement.stdP = 0.01*ones(length(spectral.wlP),1);
    end
    J2=numjacobian(params1,measurement,input);
    J2 = J2(I1,:);
    leafbio.std = abs((inv(J2.'*J2)) * J2.' * measurement.std(spectral.wlM>=wlrange.wlmin & spectral.wlM<=wlrange.wlmax & ~isnan(measurement.refl)));
        
    params1(12)     = -5*log(1-params1(12));
    soilpar.B       = params1(1);
    soilpar.lat     = params1(2);
    soilpar.lon     = params1(3);
    soilpar.SMp     = params1(4);
    leafbio.Cab     = params1(5);
    leafbio.Cdm     = params1(6);
    leafbio.Cw      = params1(7);
    leafbio.Cs      = params1(8);
    leafbio.Cca     = params1(9);
    leafbio.Cant    = params1(10);
    leafbio.N       = params1(11);
    canopy.LAI      = params1(12);
    canopy.LIDFa    = (params1(13)+params1(14))/2;
    canopy.LIDFb    = (params1(13)-params1(14))/2;
    %params1(12)      = canopy.LIDFa;
    %params1(13)      = canopy.LIDFb;
    fluorescence.wpcf= params1(15:18);
    
    t1      = tS(:,1);
    t3      = tS(:,2);
    t4      = tS(:,3);
    t5      = tS(:,4);
    t12     = tS(:,5);
    
    if isempty(Esundata)
        Esun_   = pi*t1.*t4;
    else
        Esun_   = interp1(Esundata(:,1), Esundata(:,2),spectral.wlP,'splines');
    end
    if isempty(Eskydata)
        Esky_   = pi./(1-t3.*rad.rdd).*(t1.*(t5+t12.*rad.rsd));
    else
        Esky_   = interp1(Eskydata(:,1), Eskydata(:,2),spectral.wlP,'splines');
    end
    
    %     Esun_   = pi*t1.*t4;
    %     Esky_   = pi./(1-t3.*rad.rdd).*(t1.*(t5+t12.*rad.rsd));
    E       = 1E-3*Sint(Esun_(spectral.IwlP)+Esky_(spectral.IwlP),spectral.wlP);
    SIF     = rad.SIF;
    SIFnorm = SIF/E;
    %rmse            = sqrt(er'*er./length(er));
    
    if length(c)>1
        soilpar_all.B(j)    = soilpar.B;
        soilpar_all.lat(j)    = soilpar.lat;
        soilpar_all.lon(j)    = soilpar.lon;
        soilpar_all.SMp(j)    = soilpar.SMp;
        
        leafbio_all.Cab(j) = leafbio.Cab;
        leafbio_all.Cca(j) = leafbio.Cca;
        leafbio_all.Cant(j) = leafbio.Cant;
        leafbio_all.Cdm(j) = leafbio.Cdm;
        leafbio_all.Cw(j)  = leafbio.Cw;
        leafbio_all.Cs(j)  = leafbio.Cs;
        leafbio_all.N(j)   = leafbio.N;
        leafbio_all.std(j,:) = leafbio.std';
        
        rmse_all(j) = rmse;
        
        canopy_all.LAI(j)       = canopy.LAI;
        canopy_all.LIDFa(j)     = canopy.LIDFa;
        canopy_all.LIDFb(j)     = canopy.LIDFb;
        reflSAIL_all(:,j)       = reflSAIL;
        SIF_all(:,j)            = SIF;
        SIFnorm_all(:,j)        = SIFnorm;
        rsoil_all(:,j)          = soil.refl;
        fluorescence_all.wpcf(:,j) = fluorescence.wpcf;
        J_all(:,:,j)            = J2;
    else
        reflSAIL_all = reflSAIL;
        SIF_all = SIF;
        SIFnorm_all = SIFnorm;
        leafbio_all = leafbio;
        soilpar_all = soilpar;
        leafbio_all.std = leafbio_all.std';
        canopy_all = canopy;
        rsoil_all = soil.refl;
        rmse_all  = rmse;
        fluorescence_all.wpcf = fluorescence.wpcf;
        J_all = J2;
    end
    progressbar(1/100  + j/length(c)*.8)
    
end
if ~isempty(iparams)
    
    [U,S,V] = svd(J2(:,iparams).*((ones(length(I1),1)*(UB(iparams)-LB(iparams))')),0);
    
%     for x = 1:length(iparams)
%         fprintf('%9.2f',V(x,:))
%         fprintf('\r')
%     end
    
    %%
    
    F11 = figure(11); clf
    %set(F11,'Position',[215 552 978 370])
    set(F11,'Position',[215 64 388 574])
    [~,I] = max(abs(V),[],1);
    diagS = diag(S);
     
    parnames = {'GSVB','GSVlat','GSVlon','SMp','C_{ab}', 'C_{dm}', 'C_w', 'C_{s}', 'C_{ca}','Cant','N', 'L','LIDFsum','LIDFdif','C_1','C_2', 'C_3', 'C_4'};
    for k = 1:length(iparams)
        s(k) = subplot(4,ceil(length(iparams)/4),k);
        plot(spectral.wlM(I1),U(:,k),'k')
        set(gca, 'xlim' ,[400 900])
        if k > 9, xlabel('wl (nm)'), end
        
        str1 = ['U_{' num2str(k) '}'];
        
       str2 = [' (S_{' num2str(k) '}=' num2str(round(diagS(k)*1E3)/1E3) ', \uparrow' parnames{(iparams(I(k)))} ')'];
       str = {str1,str2};
       title(str,'FontSize',9)
        
        %if k==1 || k == 6, ylabel('U'), end
    end
    %resizefigure(s,3,4,.07,.05,.05,.1,.97,.9)
end

%%
% figure(12), clf
%
% J_norm = J2(:,iparams).*(ones(length(I1),1)*(UB(iparams)-LB(iparams))');
% for k = 1:12
%     subplot(2,6,k)
%
%     plot(spectral.wlM(I1),J_norm(:,k))
%     title(parnames(k))
% end


%% save the output
string          = clock;
Output_dir      =fullfile(outdirname, sprintf('%4.0f-%02.0f-%02.0f-%02.0f%02.0f/',[string(1) string(2) string(3) string(4) string(5)]));
%Output_dir      = Output_dir{1};
mkdir(Output_dir)
copyfile('input_data.xlsx',[Output_dir, 'output_data.xlsx'],'f');
outfile = [Output_dir, 'output_data.xlsx'];
xlswrite(outfile,measured.refl(:,c),'Rmeas','B2'  );

progressbar(3/4)
xlswrite(outfile,spectral.wlM,'Rmod','A2'  );
xlswrite(outfile,spectral.wlM,'Rmeas','A2'  );
xlswrite(outfile,spectral.wlM,'Rsoilmod','A2'  );
xlswrite(outfile,reflSAIL_all,'Rmod','B2'  );
xlswrite(outfile,rsoil_all,'Rsoilmod','B2'  );
xlswrite(outfile,soilpar_all.B','output','B2'  );
xlswrite(outfile,soilpar_all.lat','output','B3'  );
xlswrite(outfile,soilpar_all.lon','output','B4'  );
xlswrite(outfile,soilpar_all.SMp','output','B5'  );
xlswrite(outfile,leafbio_all.Cab','output','B6'  );
xlswrite(outfile,leafbio_all.Cw','output','B7'  );
xlswrite(outfile,leafbio_all.Cdm','output','B8'  );
xlswrite(outfile,leafbio_all.Cs','output','B9'  );
xlswrite(outfile,leafbio_all.Cca','output','B10'  );
xlswrite(outfile,leafbio_all.Cant','output','B11'  );
xlswrite(outfile,leafbio_all.N','output','B12'  );
xlswrite(outfile,canopy_all.LAI','output','B13'  );
xlswrite(outfile,canopy_all.LIDFa','output','B14'  );
xlswrite(outfile,canopy_all.LIDFb','output','B15'  );
xlswrite(outfile,fluorescence_all.wpcf,'output','B16'  );
xlswrite(outfile,SIF_all,'Fluorescence','B2'  );
xlswrite(outfile,SIFnorm_all,'Fluorescence_norm','B2'  );
xlswrite(outfile,rmse_all','output','B20'  );
xlswrite(outfile,leafbio_all.std','output','B21');

save([Output_dir 'J'],'J_all')
progressbar(1)

%% lets plot the result
figure(1), clf
plot(spectral.wlM,[reflSAIL,measurement.refl])
legend('RTMo','Measured')
xlabel('wavelength (nm)')
ylabel('reflectance')
string          = clock;
Output_dir      =fullfile(outdirname, sprintf('%4.0f-%02.0f-%02.0f-%02.0f%02.0f/',[string(1) string(2) string(3) string(4) string(5)]));
%Output_dir      = Output_dir{1};
mkdir(Output_dir)
copyfile('input_data.xlsx',[Output_dir, 'output_data.xlsx'],'f');
outfile = [Output_dir, 'output_data.xlsx'];
xlswrite(outfile,measured.refl(:,c),'Rmeas','B2'  );

progressbar(3/4)
xlswrite(outfile,spectral.wlM,'Rmod','A2'  );
xlswrite(outfile,spectral.wlM,'Rmeas','A2'  );
xlswrite(outfile,spectral.wlM,'Rsoilmod','A2'  );
xlswrite(outfile,reflSAIL_all,'Rmod','B2'  );
xlswrite(outfile,rsoil_all,'Rsoilmod','B2'  );
xlswrite(outfile,soilpar_all.B','output','B2'  );
xlswrite(outfile,soilpar_all.lat','output','B3'  );
xlswrite(outfile,soilpar_all.lon','output','B4'  );
xlswrite(outfile,soilpar_all.SMp','output','B5'  );
xlswrite(outfile,leafbio_all.Cab','output','B6'  );
xlswrite(outfile,leafbio_all.Cdm','output','B7'  );
xlswrite(outfile,leafbio_all.Cw','output','B8'  );
xlswrite(outfile,leafbio_all.Cs','output','B9'  );
xlswrite(outfile,leafbio_all.Cca','output','B10'  );
xlswrite(outfile,leafbio_all.Cant','output','B11'  );
xlswrite(outfile,leafbio_all.N','output','B12'  );
xlswrite(outfile,canopy_all.LAI','output','B13'  );
xlswrite(outfile,canopy_all.LIDFa','output','B14'  );
xlswrite(outfile,canopy_all.LIDFb','output','B15'  );
xlswrite(outfile,fluorescence_all.wpcf,'output','B16'  );
xlswrite(outfile,SIF_all,'Fluorescence','B2'  );
xlswrite(outfile,SIFnorm_all,'Fluorescence_norm','B2'  );
xlswrite(outfile,rmse_all','output','B20'  );
xlswrite(outfile,leafbio_all.std','output','B21');

save([Output_dir 'J'],'J_all')
progressbar(1)
