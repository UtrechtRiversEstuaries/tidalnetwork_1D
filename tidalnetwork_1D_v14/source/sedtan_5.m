%sediment transport calculation to test
for channel = 1:Nc
     C_salt = sqrt(g/Cd);
     perc_salt = 1/((10^(C_salt/18))/12);
     for i = 1:Nx{channel} 
         if switchCp ==1 
            za{channel}(i) = perc_salt*(Z{channel}(i,2)-eta{channel}(1,i));
            za_susp{channel}(i) = 0.01*(Z{channel}(i,2)-eta{channel}(1,i));
         elseif switchCp == 0
            za{channel}(i) = perc_salt*(Z{channel}(i,2)-eta{channel}(1,i));
            za_susp{channel}(i) = 2.5*Sizes(channel,4);
         elseif switchCp == 2
            za{channel}(i) = perc_salt*(Z{channel}(i,2)-eta{channel}(1,i));
            za_susp{channel}(i) = perc_salt*(Z{channel}(i,2)-eta{channel}(1,i));
         end
         qf{channel}(i) = Q{channel}(i,2)/B{channel}(i);


         if switchtransp == 1
            Cf{channel}(i) = ( alc.*log10(12.2.*P{channel}(i,1)./ks) ).^(-2); %friction
                tau{channel}(i) = Cf{channel}(i)*(qf{channel}(i)^2)...
                    /(Rr*g*Sizes(channel,4)*((-eta{channel}(1,i)+Z{channel}(i,2))^2));
         elseif switchtransp == 2
             Cf{channel}(i) = 0.24/((log10(12*(-eta{channel}(1,i)+Z{channel}(i,2))/za{channel}(i)))^2);
         
             tau{channel}(i) = Cd*(qf{channel}(i)^2)...
                    /(Rr*g*Sizes(channel,4)*((-eta{channel}(1,i)+Z{channel}(i,2))^2));
                
             miu{channel}(i) = ((18*log10(12*(-eta{channel}(1,i)+Z{channel}(i,2))/za{channel}(i)))...
                    /(18*log10(12*(-eta{channel}(1,i)+Z{channel}(i,2))/(2.5*Sizes(channel,4)))))^2; %D90 = 2.5D50
                
           

        end
        %EH
        if switchtransp == 1
            tau_prime{channel}(i)=(tau{channel}(i));
            if Q{channel}(i,1)>=0
                qb_star{channel}(i,1) = (al_EH/Cf{channel}(i))*(tau{channel}(i)^nt_EH);
            elseif Q{channel}(i,1)<0
                qb_star{channel}(i,1) = -(al_EH/Cf{channel}(i))*(tau{channel}(i)^nt_EH);
            end
            %total transport m3/s
            if calc_TiMor == 0
                Qtot{channel}(i,2) = morfac*qb_star{channel}(i,1)*(sqrt(Rr*g*Sizes(channel,4))*Sizes(channel,4))*B{channel}(i);
            elseif calc_TiMor == 1
                Qtot{channel}(i,2) = qb_star{channel}(i,1)*(sqrt(Rr*g*Sizes(channel,4))*Sizes(channel,4))*B{channel}(i);
            end
        %vRijn 1984       
        elseif switchtransp==2
%     if channel == Upbound
            tau_prime{channel}(i)=(miu{channel}(i)*tau{channel}(i));       
            if (miu{channel}(i)*tau{channel}(i))<tau_cr
                T{channel}(i) = 0;
            else
                T{channel}(i) = ((miu{channel}(i)*tau{channel}(i))-tau_cr)/tau_cr;
                
            end
%             if T{channel}(i)>=3
                
                alb = alb1;
                nt = nt1;
%             elseif T{channel}(i)<3
%                 alb = alb2;
%                 nt = nt2;
%                 
%             end
%             

            %bedload calculation 
            if Q{channel}(i,2)>=0
                qb_star{channel}(i) = alb*(T{channel}(i)^nt)*(Dstar(channel)^-0.3);
            elseif Q{channel}(i,2)<0
                qb_star{channel}(i) = -alb*(T{channel}(i)^nt)*(Dstar(channel)^-0.3);
            end

            %bedload transport flux  m3/s
            qb{channel}(i) = qb_star{channel}(i)*(sqrt(Rr*g*Sizes(channel,4))*Sizes(channel,4));
            Qb{channel}(i) = qb{channel}(i)*B{channel}(i);

            % %suspended load calculation
            Ca{channel}(i) = 0.015*alI*(Sizes(channel,4)/za_susp{channel}(i))*(T{channel}(i)^1.5)*(Dstar(channel)^-0.3);
            ustar{channel}(i) = abs((qf{channel}(i)/(-eta{channel}(1,i)+Z{channel}(i,2)))*sqrt(Cf{channel}(i)/8));
            if (miu{channel}(i)*tau{channel}(i))<tau_cr
              beta{channel}(i) = 0;
              phi{channel}(i) = 0;
              Zc{channel}(i) = 0;
              F{channel}(i) = 0;
            else
                beta{channel}(i) = 1;
%                 if beta{channel}(i)>1.5
%                     beta{channel}(i)=1.5;
%                 end
%                 phi{channel}(i) = 2.5*((ws(channel)/ustar{channel}(i))^0.8)*((Ca{channel}(i)/0.65)^0.4);
                Zc{channel}(i) = (ws(channel)/(beta{channel}(i)*vk*ustar{channel}(i)));%+phi{channel}(i);
                if Zc{channel}(i)>20
                    Zc{channel}(i)=20;
                end
                if Zc{channel}(i) == 1.2
                    F{channel}(i) = (((za_susp{channel}(i)/(-eta{channel}(1,i)+Z{channel}(i,2)))...
                        /(1-(za_susp{channel}(i)/(-eta{channel}(1,i)+Z{channel}(i,2)))))^1.2)*log(za_susp{channel}(i)...
                        /(-eta{channel}(1,i)+Z{channel}(i,2)));
                else
                    F{channel}(i) = ((((za_susp{channel}(i)/(-eta{channel}(1,i)+Z{channel}(i,2))))^Zc{channel}(i))...
                        -(((za_susp{channel}(i)/(-eta{channel}(1,i)+Z{channel}(i,2))))^1.2))/...
                        ((((1-(za_susp{channel}(i)/(-eta{channel}(1,i)+Z{channel}(i,2)))))^1.2)*(1.2-Zc{channel}(i)));
                end
            end
            qs_star{channel}(i) = F{channel}(i)*(Q{channel}(i,2)/A{channel}(i,1))*((-eta{channel}(1,i)+Z{channel}(i,2)))*Ca{channel}(i)...
                /(sqrt(Rr*g*Sizes(channel,4))*Sizes(channel,4));

            %suspended transport flux  m3/s
            qs{channel}(i) = qs_star{channel}(i)*(sqrt(Rr*g*Sizes(channel,4))*Sizes(channel,4));
            Qs{channel}(i) = (qs{channel}(i)*B{channel}(i));

            if calc_TiMor == 0
                Qtot{channel}(i,2) = morfac*(Qs{channel}(i)+Qb{channel}(i));
            elseif calc_TiMor == 1
                Qtot{channel}(i,2) = (Qs{channel}(i)+Qb{channel}(i));
            end
        end

     end
     

 end