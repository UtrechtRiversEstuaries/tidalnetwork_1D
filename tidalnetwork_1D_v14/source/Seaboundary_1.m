% this script is to add more than one tidal constituent at sea boundary(ies). 
% You just need to addanother sinusoidal function in water level
% calculation. This script works for Boyo V4

    for channel=1:Nc

        if any(channel == seabound)
            Z{channel}(Nx{channel},1)= (Tide_amp(find(seabound==channel),1)*sin((2*pi*t(j)/T_tide(1,1))-Tide_phase(find(seabound==channel),1)))...
                                        +(Tide_amp(find(seabound==channel),2)*sin((2*pi*t(j)/T_tide(2,1))-Tide_phase(find(seabound==channel),2)));
            Z{channel}(Nx{channel},2)= (Tide_amp(find(seabound==channel),1)*sin((2*pi*t(j+1)/T_tide(1,1))-Tide_phase(find(seabound==channel),1)))...
                                        +(Tide_amp(find(seabound==channel),2)*sin((2*pi*t(j+1)/T_tide(2,1))-Tide_phase(find(seabound==channel),2)));
        
        end
    end