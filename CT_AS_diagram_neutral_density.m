function [SA, CT]=CT_AS_diagram_neutral_density(P, SP, T, Lon, Lat, p_ref, isopycs_neutral, title_string)

% CT_AS_diagram_neutral_density 
% calculates the TEOS-10 derived Conservative Temperature (CT) and Absolute Salinity (AS) 
% together with the corresponding neutral density and plots the CT-AS diagram including the 
% freezing temperature line, standard potential density contours and selected neutral density contours.
%
%=====================================================================================================================
%
% USAGE: 
% CT_AS_diagram_neutral_density(P, SP, T, Lon, Lat, p_ref, isopycs_neutral, title_string)
%
% DESCRIPTION:
% This function is based on the TEOS-10 Gibbs SeaWater (GSW) Oceanographic Toolbox 
% and uses the eos80_legacy_gamma_n function for calculations of neutral density.
% It plots the Conservative Temperature-Absolute Salinity diagram including
% the freezing temperature line, potential density contours, and user-selected neutral density contours.
%
% INPUTS:
%  P   =  sea pressure (absolute pressure - 10.1325 dbar) [dbar]
%  SP  =  Practical Salinity (PSS-78)                     [unitless]
%  T   =  in-situ temperature (ITS-90)                    [deg C]
%  Lon =  longitude in decimal degrees                    [0 ... 360] or [-180 ... 180]
%  Lat =  latitude in decimal degrees north               [-90 ... 90]
%  
% Optional: 
%  p_ref = reference sea pressure for the isopycnals      [dbar] (default: 0 dbar)
%  isopycs_neutral = array of desired neutral density isopycnals (default: 28, 28.27)
%  title_string = title text for the plot
%
% OUTPUTS:
%  SA  =  Absolute Salinity                               [g/kg]
%  CT  =  Conservative Temperature                        [deg C]
%
% 
% REFERENCES:
% Jackett, D.R., and T.J. McDougall, 1997: A Neutral Density Variable for the World's Oceans. 
% J. Phys. Oceanogr., 27, 237-263. doi: 10.1175/1520-0485(1997)027<0237:ANDVFT>2.0.CO;2
%
% Unesco, 1983: Algorithms for computation of fundamental properties of seawater. 
% Unesco technical papers in marine science 44, 53pp.
%
% 
% % % EXAMPLE:
% % % Generate 100 random samples of temperature, salinity, and pressure
% % % within typical oceanographic limits of the Ross Sea and use them to 
% % % compute and plot the Conservative Temperature-Absolute Salinity diagram.
% % 
% % % Define the number of samples 
% % num_samples = 100;
% % 
% % % Generate random temperature from -2 to 2 degrees Celsius
% % T = -2 + 4 * rand(num_samples, 1);  % Scaling from 0-1 to -2 to 2
% % 
% % % Generate random salinity from 34 to 35 (PSS-78)
% % SP = 34 + (35 - 34) * rand(num_samples, 1);  % Scaling from 0-1 to 34 to 35
% % 
% % % Generate random pressure from 0 to 1000 dbar
% % P = 0 + 1000 * rand(num_samples, 1);  % Scaling from 0-1 to 0 to 1000
% % 
% % % Longitude and Latitude for a specific region (Antarctic)
% % Lon = 170 + 10 * rand(num_samples, 1);  % Longitude from 170 to 180 degrees
% % Lat = -75 + 5 * rand(num_samples, 1);   % Latitude from -75 to -70 degrees
% % 
% % % Reference pressure for the isopycnals
% % p_ref = 0;  % Surface reference
% % 
% % % Array of desired neutral density isopycnals 
% % isopycs_neutral = [27.5, 27.8, 28.1];
% % 
% % % Title for the plot
% % title_string = 'Example CT-AS Diagram';
% % 
% % % Call the function
% % [SA, CT] = CT_AS_diagram_neutral_density(P, SP, T, Lon, Lat, p_ref, isopycs_neutral, title_string);
% % %
% 
%
% AUTHOR: 
% Naomi Krauzig, Politechnic University of Marche, Ancona, 2024.
% https://github.com/naomikrauzig/CT_AS_diagram_neutral_density
%=====================================================================================================================


%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------


% Checking number of input arguments
if nargin < 5
    error('CT_AS_diagram_neutral_density: Requires at least five inputs');
elseif nargin > 9
    error('CT_AS_diagram_neutral_density: Too many inputs');
end

% Dimension checks and resizing arrays to match dimensions if necessary
[ms, ns] = size(SP);
[mt, nt] = size(T);
[mp, np] = size(P);

% Ensure SP and T have the same dimensions
if (mt ~= ms || nt ~= ns)
    error('T and SP must have the same dimensions');
end

% Resizing P to match the dimensions of SP and T if it is a scalar or a vector
if (mp == 1 && np == 1)              % P scalar - fill to size of SP
    P = P * ones(size(SP));
elseif (ns == np && mp == 1)         % P is row vector, copy down each column
    P = repmat(P, ms, 1);
elseif (ms == mp && np == 1)         % P is column vector, copy across each row
    P = repmat(P, 1, ns);
elseif (ns == mp && np == 1)         % P is a transposed row vector, transpose then copy down each column
    P = repmat(P.', ms, 1);
elseif (ms ~= mp || ns ~= np)
    error('Dimensions of P do not match dimensions of SP and T');
end

% Handling out-of-range values for T based on pressure
T(P < 100 & (T > 80 | T < -12)) = NaN;
T(P >= 100 & (T > 40 | T < -12)) = NaN;

% Transpose vectors if necessary (if all vectors are row vectors)
if ms == 1
    SP = SP.';
    T = T.';
    P = P.';
    transposed = true;
else
    transposed = false;
end

% Ensure dimensions of Lon and Lat match SP
[mla, nla] = size(Lat);
[mlo, nlo] = size(Lon);

% Latitude and longitude resizing
if any([mla nla mlo nlo] == 1)  % Check if any are scalars or vectors and resize
    if (mla == 1 && nla == 1)
        Lat = repmat(Lat, size(SP));
    elseif (ms == mla && nla == 1)
        Lat = repmat(Lat, 1, ns);
    elseif (ns == nla && mla == 1)
        Lat = repmat(Lat, ms, 1);
    elseif (ms ~= mla || ns ~= nla)
        error('Dimensions of Lat do not match dimensions of SP and T');
    end

    if (mlo == 1 && nlo == 1)
        Lon = repmat(Lon, size(SP));
    elseif (ms == mlo && nlo == 1)
        Lon = repmat(Lon, 1, ns);
    elseif (ns == nlo && mlo == 1)
        Lon = repmat(Lon, ms, 1);
    elseif (ms ~= mlo || ns ~= nlo)
        error('Dimensions of Lon do not match dimensions of SP and T');
    end
end

Lon(Lon < 0) = Lon(Lon < 0) + 360;  % Adjust longitude values

% Remove out-of-range values for SP
SP(P < 100 & SP > 120) = NaN;
SP(P >= 100 & SP > 42) = NaN;

% Change standard blank fill values to NaNs
SP(abs(SP) == 99999 | abs(SP) == 999999) = NaN;
P(abs(P) == 99999 | abs(P) == 999999) = NaN;
Lon(abs(Lon) == 9999 | abs(Lon) == 99999) = NaN;
Lat(abs(Lat) == 9999 | abs(Lat) == 99999) = NaN;

% Additional error checks for data ranges
if any(P < -1.5 | P > 12000)
    error('P is out of range');
end
if any(Lon < 0 | Lon > 360)
    error('Longitude is out of range');
end
if any(abs(Lat) > 90)
    error('Latitude is out of range');
end


%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% Handle optional inputs
if nargin < 8 || isempty(p_ref)
    p_ref = 0; % Default reference pressure
end

if nargin < 9 || isempty(isopycs_neutral)
    isopycs_neutral = [28, 28.27]; % Default neutral densities
end


SA=gsw_SA_from_SP(SP,P,Lon,Lat);
CT=gsw_CT_from_t(SA,T,P);

min_SA_data = min(min(SA(:)));
max_SA_data = max(max(SA(:)));
min_CT_data = min(min(CT(:)));
max_CT_data = max(max(CT(:)));

SA_min = min_SA_data - 0.1*(max_SA_data - min_SA_data);
SA_min(SA_min < 0) = 0;
SA_max = max_SA_data + 0.1*(max_SA_data - min_SA_data);
SA_axis = [SA_min:(SA_max-SA_min)/200:SA_max];

CT_freezing = gsw_CT_freezing(SA_axis,p_ref); 
CT_min = min_CT_data - 0.1*(max_CT_data - min_CT_data);
CT_max = max_CT_data + 0.1*(max_CT_data - min_CT_data);
if CT_min > min(CT_freezing) 
    CT_min = min_CT_data - 0.1*(max_CT_data - min(CT_freezing));
end
CT_axis = [CT_min:(CT_max-CT_min)/200:CT_max];

SA_gridded = meshgrid(SA_axis,1:length(CT_axis));
CT_gridded = meshgrid(CT_axis,1:length(SA_axis))';

isopycs_gridded = gsw_rho(SA_gridded,CT_gridded,p_ref)-1000;

Lat_mean=mean(Lat,'all','omitnan');
Lon_mean=mean(Lon,'all','omitnan');


%grid gamma_n
[gamma_n_gridded] = eos80_legacy_gamma_n(SA_gridded,CT_gridded,p_ref,Lon_mean,Lat_mean);

if transposed
    CT = CT.';
end

if transposed
    SA = SA.';
end

%--------------------------------------------------------------------------
% End of the calculation
%--------------------------------------------------------------------------


figure
axis('square');
xlim([SA_min-0.01 SA_max+0.01])
ylim([CT_min-0.05 CT_max+0.05])
[cd,ch] = contour(SA_gridded,CT_gridded,isopycs_gridded,':','Color',[.5 .5 .5]);
clabel(cd,ch,'LabelSpacing',100,'FontSize',6,'Color',[.5 .5 .5]);
hold on 
[cd2,ch2] = contour(SA_gridded,CT_gridded,gamma_n_gridded,isopycs_neutral,'k');
clabel(cd2,ch2,'LabelSpacing',125,'FontSize',10);
hold on
scatter(SA(:),CT(:),10,'filled')
ylabel('Conservative Temperature (^oC)','FontWeight','bold')
xlabel('Absolute Salinity (g/kg)','FontWeight','bold')
box on
hold on 
plot(SA_axis,CT_freezing,'b','Linewidth',2)
hold on
if exist('title_string','var')
    title([title_string],'fontsize',12);
else
    title('CT-AS diagram','fontsize',12);
end
set(findobj(gcf,'type','axes'),'FontSize',10,'FontWeight','bold');


end

