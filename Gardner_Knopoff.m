function[declustered_catalogue] = Gardner_Knopoff(catalogue,mainshock,type)
%==========================================================================
%
%   Input:
%    catalogue := earthquake catalogue (filtered for min. and max.
%    magnitudes), where column
%
%   (1) := month
%   (2) := day
%   (3) := year
%   (4) := hour
%   (5) := minute
%   (6) := second
%   (7) := magnitude
%   (8) := latitude
%   (9) := longitude
%   (10 := depth (km)
%
%    mainshock := earthquake catalogue of "mainshocks" formatted as above 
%
%    type := 1 (filter for foreshocks), 2 (filter for aftershocks), or 3
%    (filter for foreshocks and aftershocks)
%
%   Output:
%    filtered_catalogue := result of Gardner and Knopoff's declustering
%    algorithm
%   
%==========================================================================

reference = referenceEllipsoid('WGS84');

declustered_catalogue = catalogue;

for ii = 1:size(mainshock,1)
    
    t_mainshock = datetime(mainshock(ii,3),mainshock(ii,1),...
        mainshock(ii,2),mainshock(ii,4),mainshock(ii,5),mainshock(ii,6));
    
    if mainshock(ii,7) >= 6.5
        t_radius = 10^(0.0320*mainshock(ii,7) + 2.7389);
        %^days
    elseif mainshock(ii,7) < 6.5
        t_radius = 10^(0.5409*mainshock(ii,7) - 0.5470);
        %^days
    end
    
    if type == 1
        %^declustering foreshocks
        t_lower = t_mainshock - t_radius;
        t_upper = t_mainshock;
    elseif type == 2
        %^declustering aftershocks
        t_lower = t_mainshock;
        t_upper = t_mainshock + t_radius;
    elseif type == 3
        %^declustering foreshocks and aftershocks
        t_lower = t_mainshock - t_radius;
        t_upper = t_mainshock + t_radius;
    end
    
    lat1 = mainshock(ii,8);
    lon1 = mainshock(ii,9);
    
    d_radius = 10^(0.1238*mainshock(ii,7) + 0.9873);
    %^km 
    
    index = zeros(size(declustered_catalogue,1),1);
    
    for jj = 1:size(declustered_catalogue,1)
        
        check = 0;
        
        t = datetime(declustered_catalogue(jj,3),...
        declustered_catalogue(jj,1),declustered_catalogue(jj,2),...
        declustered_catalogue(jj,4),declustered_catalogue(jj,5),...
        declustered_catalogue(jj,6));
        %^date of earthquake jj
        
        if isbetween(t,t_lower,t_upper)
            check = check + 0.5;
        end
        %^time window check
        
        lat2 = declustered_catalogue(jj,8);
        lon2 = declustered_catalogue(jj,9);
        
        d = distance(lat1,lon1,lat2,lon2,reference)/1000;
        %^distance (km) of earthquake jj from mainshock ii
        
        if d <= d_radius
            check = check + 0.5;
        end
        %^space window check
        
        if check == 1
            index(jj) = 1;
        end
    end 
    
    declustered_catalogue(logical(index),:) = [];
    %^declustering with respect to mainshock ii

end

end