function BAND = adjustBand(TT,band,crit)
    if (TT >= 2.9) && (TT < 4) BAND = [band(1)*crit, band(2), band(3)*crit]; % Adjust BAND
    else
        BAND = band; 
    end
end