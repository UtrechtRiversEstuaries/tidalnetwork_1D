%sediment transport and bed update w/ tide-averaged morfac

for channel = 1:Nc
 %store sediment transport for one tidal cycle
        %for first time step
        if t(j)==0
           Q_storetide{channel}(:,1) = Qtot{channel}(:,1);  
        end
        %reboot the matrix to store the next tidal cycle
        if mod(j,jTide)==0
                 Q_storetide{channel} = zeros(Nx{channel},jTide);
        end
        
        %for first tidal cycle
        if j+1<=jTide
            Q_storetide{channel}(:,j+1) = Qtot{channel}(:,2);
            if mod(j,jTide)==0
                 Q_storetide{channel} = zeros(Nx{channel},jTide);
            end
        else
            checkJold = floor((j+1)/jTide)-1;
            checkJ = floor((j+1)/jTide);
            
            % store transpport
            %at the end of cycle
            if mod(j+1,jTide)==0
                Q_storetide{channel}(:,(j+1)-(checkJold*jTide)) = Qtot{channel}(:,2);
            %in tidal cycle
            else 
                Q_storetide{channel}(:,(j+1)-(checkJ*jTide)) = Qtot{channel}(:,2);
            end
        end
        
        %average transport for one tidal cycle
        if t(j+1)>=T_tide(1) && mod(j+1,jTide)==0
            for i=1:Nx{channel}
                    Qtot_tideAve{channel}(i) = morfac*mean(Q_storetide{channel}(i,:)); 
            end
         
        %bed update
            if t(j)>=morstt        
                for i = 2:Nx{channel}-1
                    
                            if Qtot{channel}(i,2)>=0
                                eta{channel}(2,i) = (-(((1-theta)*(Qtot_tideAve{channel}(i+1)-Qtot_tideAve{channel}(i)))+((1+theta)*(Qtot_tideAve{channel}(i)-Qtot_tideAve{channel}(i-1))))*(T_tide(1)/(2*dx))*(1/B{channel}(i))*(1/(1-lamp)))+eta{channel}(1,i);

                            else
                                eta{channel}(2,i) = (-(((1+theta)*(Qtot_tideAve{channel}(i+1)-Qtot_tideAve{channel}(i)))+((1-theta)*(Qtot_tideAve{channel}(i)-Qtot_tideAve{channel}(i-1))))*(T_tide(1)/(2*dx))*(1/B{channel}(i))*(1/(1-lamp)))+eta{channel}(1,i);
                            end

                end
                
                    % at in-channel and open boundaries
                    eta{channel}(2,1) = (-((-(Qtot_tideAve{channel}(3)))+(4*(Qtot_tideAve{channel}(2)))-(3*(Qtot_tideAve{channel}(1))))*(T_tide(1)/(2*dx))*(1/B{channel}(1))*(1/(1-lamp)))+eta{channel}(1,1);
                    eta{channel}(2,end) = (-(((Qtot_tideAve{channel}(end-2)))-(4*(Qtot_tideAve{channel}(end-1)))+(3*(Qtot_tideAve{channel}(end))))*(T_tide(1)/(2*dx))*(1/B{channel}(end))*(1/(1-lamp)))+eta{channel}(1,end);
            end
        else
            
            eta{channel}(2,:)= eta{channel}(1,:);
            
        end
end
