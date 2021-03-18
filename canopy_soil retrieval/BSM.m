function rwet = BSM(soilpar,spec,emp)
    
    % Spectral parameters
    
    %wl  = spec.wl;          % wavelengths
    GSV = spec.GSV;         % Global Soil Vectors spectra (nwl * 3)
    kw  = spec.kw;          % water absorption spectrum
    nw  = spec.nw;          % water refraction index spectrum
    
    % Soil parameters
    
    B   = soilpar.B;        % soil brightness
    lat = soilpar.lat;      % spectral shape latitude (range = 20 - 40 deg)
    lon = soilpar.lon;      % spectral shape longitude (range = 45 - 65 deg)
    SMp = soilpar.SMp;      % soil moisture volume percentage (5 - 55)
    
    % Empirical parameters
    
    SMC  = emp.SMC;         % soil moisture capacity parameter
    film = emp.film;        % single water film optical thickness
    
    f1 = B * sind(lat);
    f2 = B * cosd(lat) * sind(lon);
    f3 = B * cosd(lat) * cosd(lon);
    
    rdry = f1 * GSV(:,1) + f2 * GSV(:,2) + f3 * GSV(:,3);
    
    % Soil moisture effect
    
    rwet = soilwat(rdry,nw,kw,SMp,SMC,film);
    
end

