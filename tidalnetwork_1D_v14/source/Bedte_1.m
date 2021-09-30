if t(j)>=morstt

    for channel = 1:Nc
            % at in-channel and open boundaries

            for i = 2:Nx{channel}-1
                    if Qtot{channel}(i,2)>=0
                        eta{channel}(2,i) = (-(((1-theta)*(Qtot{channel}(i+1,2)-Qtot{channel}(i,2)))+((1+theta)*(Qtot{channel}(i,2)-Qtot{channel}(i-1,2))))*(dt/(2*dx))*(1/B{channel}(i))*(1/(1-lamp)))+eta{channel}(1,i);
    %                     eta{channel}(2,i) = (-(((Qtot{channel}(i+1,2)+Qtot{channel}(i,2)))-((Qtot{channel}(i,2)+Qtot{channel}(i-1,2))))*(dt/(2*dx))*(1/B{channel}(i))*(1/(1-lamp)))+eta{channel}(1,i);
                    else
                        eta{channel}(2,i) = (-(((1+theta)*(Qtot{channel}(i+1,2)-Qtot{channel}(i,2)))+((1-theta)*(Qtot{channel}(i,2)-Qtot{channel}(i-1,2))))*(dt/(2*dx))*(1/B{channel}(i))*(1/(1-lamp)))+eta{channel}(1,i);
                    end
    %             eta{channel}(2,i) = (-(((theta)*((Qtot{channel}(i+1,2)+Qtot{channel}(i,2))-(Qtot{channel}(i-1,2)+Qtot{channel}(i,2))))+((1-theta)*((Qtot{channel}(i+1,1)+Qtot{channel}(i,1))-(Qtot{channel}(i-1,1)+Qtot{channel}(i,1)))))*(dt/(2*dx))*(1/B{channel}(i))*(1/(1-lamp)))+eta{channel}(1,i);
            end
            % at in-channel and open boundaries
    %         eta{channel}(2,1) = eta{channel}(1,1);
    %         eta{channel}(2,end) = eta{channel}(1,end);
                eta{channel}(2,1) = (-((-(Qtot{channel}(3,2)))+(4*(Qtot{channel}(2,2)))-(3*(Qtot{channel}(1,2))))*(dt/(2*dx))*(1/B{channel}(1))*(1/(1-lamp)))+eta{channel}(1,1);
                eta{channel}(2,end) = (-(((Qtot{channel}(end-2,2)))-(4*(Qtot{channel}(end-1,2)))+(3*(Qtot{channel}(end,2))))*(dt/(2*dx))*(1/B{channel}(end))*(1/(1-lamp)))+eta{channel}(1,end);    

    %             eta{channel}(2,1) = (((-((theta*(Qtot{channel}(2,2)-Qtot{channel}(1,2)))+((1-theta)*(Qtot{channel}(2,1)-Qtot{channel}(1,1))))*(1/dx)*(1/(0.5*(B{channel}(1)+B{channel}(2))))*(1/(1-lamp)))...
    %                                 +((eta{channel}(1,1)+eta{channel}(1,2))/(2*dt)))*2*dt)-eta{channel}(2,2);
    %                  
    %             eta{channel}(2,Nx{channel}) = (((-((theta*(Qtot{channel}(Nx{channel},2)-Qtot{channel}(Nx{channel}-1,2)))+((1-theta)*(Qtot{channel}(Nx{channel},1)-Qtot{channel}(Nx{channel}-1,1))))*(1/dx)*(1/(0.5*(B{channel}(Nx{channel}-1)+B{channel}(Nx{channel}))))*(1/(1-lamp)))...
    %                             +((eta{channel}(1,Nx{channel}-1)+eta{channel}(1,Nx{channel}))/(2*dt)))*2*dt)-eta{channel}(2,Nx{channel}-1);


    end

 end