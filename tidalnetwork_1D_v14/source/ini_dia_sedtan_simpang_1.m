%initial transport just in case the simulation is hot start
for channel = 1:Nc
    for i = 1:Nx{channel}
            za{channel}(i,1) = 0.01*(Z{channel}(i,1)-eta{channel}(1,i)); % reference depth to distinguish suspended load and bedload according to van rijn
            qf{channel}(i,1) = Q{channel}(i,1)/B{channel}(i);
            
            % in case za~= ks condition must be defined
            if switchtransp == 1
                Cf{channel}(i) = ( alc.*log10(12.2.*P{channel}(i,1)./ks) ).^(-2); %friction
                tau{channel}(i) = Cf{channel}(i)*(qf{channel}(i)^2)...
                    /(Rr*g*Sizes(channel,4)*((-eta{channel}(1,i)+Z{channel}(i,1))^2));
            elseif switchtransp == 2
                Cf{channel}(i) = 0.24/((log10(12*(-eta{channel}(1,i)+Z{channel}(i,1))/za{channel}(i)))^2);
                tau{channel}(i) = 0.125*Cf{channel}(i)*(qf{channel}(i)^2)...
                            /(Rr*g*Sizes(channel,4)*((-eta{channel}(1,i)+Z{channel}(i,1))^2));
                miu{channel}(i) = ((18*log10(12*(-eta{channel}(1,i)+Z{channel}(i,1))/za{channel}(i)))...
                            /(18*log10(12*(-eta{channel}(1,i)+Z{channel}(i,1))/(3*2*Sizes(channel,4)))))^2;
                end
       %EH
       if switchtransp == 1;
            if Q{channel}(i,1)>=0
                qb_star{channel}(i,1) = (al_EH/Cf{channel}(i))*tau{channel}(i)^nt_EH;
            elseif Q{channel}(i,1)<0
                qb_star{channel}(i,1) = -(al_EH/Cf{channel}(i))*tau{channel}(i)^nt_EH;
            end
            %total transport
            if calc_TiMor ==0 
                Qtot{channel}(i,1) = morfac*qb_star{channel}(i,1)*(sqrt(Rr*g*Sizes(channel,4))*Sizes(channel,4))*B{channel}(i);
            elseif calc_TiMor ==1
                Qtot{channel}(i,1) = qb_star{channel}(i,1)*(sqrt(Rr*g*Sizes(channel,4))*Sizes(channel,4))*B{channel}(i);
            end
            
        %vanRijn 1984
        elseif switchtransp == 2;
            %bedload calculation
            tau_prime{channel}(i)=(miu{channel}(i)*tau{channel}(i));
            if (miu{channel}(i)*tau{channel}(i,1))<=tau_cr
                T{channel}(i,1) = 0;
            else
                T{channel}(i,1) = (miu{channel}(i)*tau{channel}(i,1))-tau_cr/tau_cr;
            end
            
%             if T{channel}(i)>=3
                
                alb = alb1;
                nt = nt1;
%             elseif T{channel}(i)<3
%                 alb = alb2;
%                 nt = nt2;
%                 
%             end

            if Q{channel}(i,1)>=0
                qb_star{channel}(i,1) = alb*(T{channel}(i,1)^nt)*(Dstar(channel)^-0.3);
            elseif Q{channel}(i,1)<0
                qb_star{channel}(i,1) = -alb*(T{channel}(i,1)^nt)*(Dstar(channel)^-0.3);
            end

            %sediment transport flux
            qb{channel}(i) = qb_star{channel}(i)*(sqrt(Rr*g*Sizes(channel,4))*Sizes(channel,4));
            Qb{channel}(i) = qb{channel}(i)*B{channel}(i);

            %suspended load calculation
            Ca{channel}(i,1) = 0.015*alI*(Sizes(channel,4)/za{channel}(i,1))*(T{channel}(i,1)^1.5)*(Dstar(channel)^-0.3);
            ustar{channel}(i,1) = abs((qf{channel}(i)/(-eta{channel}(1,i)+Z{channel}(i,1)))*sqrt(Cf{channel}(i)/8));
            if (miu{channel}(i)*tau{channel}(i,1))<=tau_cr
              beta{channel}(i,1) = 0;
              phi{channel}(i,1) = 0;
              Zc{channel}(i,1) = 0;
              F{channel}(i,1) = 0;
            else
                beta{channel}(i,1) = 1+(2*((ws(channel)/ustar{channel}(i,1))^2));
                if beta{channel}(i,1)>1.5
                    beta{channel}(i,1)=1.5;
                end
                phi{channel}(i,1) = 2.5*((ws(channel)/ustar{channel}(i,1))^0.8)*((Ca{channel}(i,1)/0.65)^0.4);
                Zc{channel}(i,1) = (ws(channel)/(beta{channel}(i,1)*vk*ustar{channel}(i,1)))+phi{channel}(i,1);
                if Zc{channel}(i,1)>20
                    Zc{channel}(i,1)=20;
                end
                if Zc{channel}(i,1) == 1.2
                    F{channel}(i,1) = (((za{channel}(i,1)/(-eta{channel}(i)+Z{channel}(i,1)))...
                        /(1-(za{channel}(i,1)/(-eta{channel}(i)+Z{channel}(i,1)))))^1.2)*log(za{channel}(i,1)...
                        /(-eta{channel}(i)+Z{channel}(i,1)));
                else
                    F{channel}(i,1) = ((((za{channel}(i,1)/(-eta{channel}(i)+Z{channel}(i,1))))^Zc{channel}(i,1))...
                        -(((za{channel}(i,1)/(-eta{channel}(i)+Z{channel}(i,1))))^1.2))/...
                        ((((1-(za{channel}(i,1)/(-eta{channel}(i)+Z{channel}(i,1)))))^1.2)*(1.2-Zc{channel}(i,1)));
                end
            end

            qs_star{channel}(i,1) = F{channel}(i,1)*(Q{channel}(i,1)...
                /A{channel}(i,1))*((-eta{channel}(i)+Z{channel}(i,1)))*Ca{channel}(i,1)...
                /(sqrt(Rr*g*Sizes(channel,4))*Sizes(channel,4));
            qs{channel}(i,1) = qs_star{channel}(i,1)*(sqrt(Rr*g*Sizes(channel,4))*Sizes(channel,4));


                Qs{channel}(i,1) = (qs{channel}(i,1)*B{channel}(i));

            %total transport
            if calc_TiMor ==0 
                Qtot{channel}(i,1) = morfac*(Qs{channel}(i,1)+Qb{channel}(i,1)); 
            elseif calc_TiMor ==1
                Qtot{channel}(i,1) = (Qs{channel}(i,1)+Qb{channel}(i,1));
            end
        end
    end
    %for nodal point
    %EH
% for nodal point
%      %EH
%      if switchtransp == 1
%          for bifur = 1:Nbif
%              qvRsusnode(bifur,1) = 1;
%          end
%      %vRijn 1984    
%      elseif switchtransp == 2
%         for bifur = 1:Nbif
%            if Q{topob{bifur}(1)}(end,1)>=0 && Q{topob{bifur}(2)}(1,1)>=0 && Q{topob{bifur}(3)}(1,1)>=0
%                 if qb_star{topob{bifur}(1)}(end) == 0 &&  qs_star{topob{bifur}(1)}(end) == 0
%                     qvRsusnode(bifur,1) = 0;
%                 else
%                     qvRsusnode(bifur,1) = 1;%qb_star{topob{bifur}(1)}(end) / (qb_star{topob{bifur}(1)}(end) + qs_star{topob{bifur}(1)}(end));
%                 end
%                 
%            elseif Q{topob{bifur}(1)}(end,1)<0 && Q{topob{bifur}(2)}(1,1)<0 && Q{topob{bifur}(3)}(1,1)>=0
%                 if qb_star{topob{bifur}(2)}(1) == 0 &&  qs_star{topob{bifur}(2)}(1) == 0
%                     qvRsusnode(bifur,1) = 0;
%                 else
%                     qvRsusnode(bifur,1) = 1;%qb_star{topob{bifur}(2)}(1) / (qb_star{topob{bifur}(2)}(1) + qs_star{topob{bifur}(2)}(1));
%                 end
%                
%            elseif Q{topob{bifur}(1)}(end,1)<0 && Q{topob{bifur}(2)}(1,1)>=0 && Q{topob{bifur}(3)}(1,1)<0
%                if qb_star{topob{bifur}(3)}(1) == 0 &&  qs_star{topob{bifur}(3)}(1) == 0
%                     qvRsusnode(bifur,1) = 0;
%                else
%                     qvRsusnode(bifur,1) = 1;%qb_star{topob{bifur}(3)}(1) / (qb_star{topob{bifur}(3)}(1) + qs_star{topob{bifur}(3)}(1));
%                end
%                            
%            else
%                qvRsusnode(bifur,1) = 0;
%            end
%         end
%      end

 %% Shields at junctions
 if Nbif >0
  for bifur = 1:Nbif
  tau_bolla1(bifur,1) = Cd*(qf{topob{bifur}(1)}(end)^2)...
                    /(Rr*g*Sizes(channel,4)*((-eta{topob{bifur}(1)}(1,end) +Z{topob{bifur}(1)}(end,1))^2));
  tau_bolla2(bifur,1) = Cd*(qf{topob{bifur}(2)}(1)^2)...
                    /(Rr*g*Sizes(channel,4)*((-eta{topob{bifur}(2)}(1,1) +Z{topob{bifur}(2)}(1,1))^2));
  tau_bolla3(bifur,1) = Cd*(qf{topob{bifur}(3)}(1)^2)...
                    /(Rr*g*Sizes(channel,4)*((-eta{topob{bifur}(3)}(1,1) +Z{topob{bifur}(3)}(1,1))^2));
  end
 end
 
   if Nconf>0
  for conflu = 1:Nconf
  tau_bolla1_conf(conflu,1) = Cd*(qf{topoc{conflu}(1)}(end)^2)...
                    /(Rr*g*Sizes(channel,4)*((-eta{topoc{conflu}(1)}(1,end) +Z{topoc{conflu}(1)}(end,1))^2));
  tau_bolla2_conf(conflu,1) = Cd*(qf{topoc{conflu}(2)}(end)^2)...
                    /(Rr*g*Sizes(channel,4)*((-eta{topoc{conflu}(2)}(1,1) +Z{topoc{conflu}(2)}(1,1))^2));
  tau_bolla3_conf(conflu,1) = Cd*(qf{topoc{conflu}(3)}(1)^2)...
                    /(Rr*g*Sizes(channel,4)*((-eta{topoc{conflu}(3)}(1,1) +Z{topoc{conflu}(3)}(1,1))^2));
  end
 end
 %% nodal point calculation
 if Nbif >0
for bifur = 1:Nbif
     
     
     %% first condition - bifurcation
        if Q{topob{bifur}(1)}(end,1)>=0 && Q{topob{bifur}(2)}(1,1)>=0 && Q{topob{bifur}(3)}(1,1)>=0
          if switchndpt == 0 || switchndpt == 1
            Qy(bifur,1) = ( Q{topob{bifur}(2)}(1,1)-Q{topob{bifur}(3)}(1,1)-(Q{topob{bifur}(1)}(end,1)*...
                ((B{topob{bifur}(2)}(1)-B{topob{bifur}(3)}(1))/(B{topob{bifur}(2)}(1)+B{topob{bifur}(3)}(1)))) )/2;
            H123(bifur,1) = ( (( (-eta{topob{bifur}(2)}(1,1)+Z{topob{bifur}(2)}(1,1))...
                +(-eta{topob{bifur}(3)}(1,1)+Z{topob{bifur}(3)}(1,1)) )/2)...
                + (-eta{topob{bifur}(1)}(1,end)+Z{topob{bifur}(1)}(end,1)) )/2;
  
            dzdy(bifur,1) = (eta{topob{bifur}(2)}(1,1) - eta{topob{bifur}(3)}(1,1))/(B{topob{bifur}(1)}(end)/2);

           if switchndpt == 0
               
               qsy(bifur,1) = ...
                   (Qtot{topob{bifur}(1)}(end,1)/B{topob{bifur}(1)}(end)) * ( ((Qy(bifur,1)*(-eta{topob{bifur}(1)}(1,end)+Z{topob{bifur}(1)}(end,1)))...
                   /(Q{topob{bifur}(1)}(end,1)*alw*H123(bifur,1))) - ...
                  (dzdy(bifur,1)*rbolla/sqrt(tau_bolla1(bifur,1))) );
              
           elseif  switchndpt == 1
               qsy(bifur,1) =  ...
                   ( qb{topob{bifur}(1)}(end) * ( ((Qy(bifur,1)*(-eta{topob{bifur}(1)}(1,end) +Z{topob{bifur}(1)}(end,1)))...
                   /(Q{topob{bifur}(1)}(end,1)*alw*H123(bifur,1))) - ...
                  (dzdy(bifur,1)*rbolla/sqrt(tau_bolla1(bifur,1))) ) ) + ( qs{topob{bifur}(1)}(end)...
                  * ((Qy(bifur,1)*(-eta{topob{bifur}(1)}(1,end) +Z{topob{bifur}(1)}(end,1)))...
                   /(Q{topob{bifur}(1)}(end,1)*alw*H123(bifur,1))) );
           end
           
           if calc_TiMor == 0
                    Qsy(bifur,1) = morfac*qsy(bifur,1)*alw*B{topob{bifur}(1)}(end);
           elseif calc_TiMor == 1
                    Qsy(bifur,1) = qsy(bifur,1)*alw*B{topob{bifur}(1)}(end);
           end
           
           
           Qtot{topob{bifur}(2)}(1,1) = Qsy(bifur,1) + ...
                 ((B{topob{bifur}(2)}(1)/(B{topob{bifur}(2)}(1)+B{topob{bifur}(3)}(1))) * Qtot{topob{bifur}(1)}(end,1));
          
          elseif switchndpt == 2
              rat_Qtot(bifur,1) = ((B{topob{bifur}(2)}(1)/B{topob{bifur}(3)}(1))^(1-k_wang)) * (Q{topob{bifur}(2)}(1,1)/Q{topob{bifur}(3)}(1,1))^k_wang;
              Qtot{topob{bifur}(2)}(1,1) = rat_Qtot(bifur,1) * (Qtot{topob{bifur}(1)}(end,1)) / (1+rat_Qtot(bifur,1));
                
          end
          
           Qtot{topob{bifur}(3)}(1,1) = Qtot{topob{bifur}(1)}(end,1)-Qtot{topob{bifur}(2)}(1,1);

         %% second condition - bifurcation    
         elseif Q{topob{bifur}(1)}(end,1)<0 && Q{topob{bifur}(2)}(1,1)<0 && Q{topob{bifur}(3)}(1,1)>=0
           if switchndpt == 0 || switchndpt == 1   
             Qy(bifur,1) = ( -Q{topob{bifur}(1)}(end,1)-Q{topob{bifur}(3)}(1,1)+(Q{topob{bifur}(2)}(1,1)*...
                ((B{topob{bifur}(1)}(end)-B{topob{bifur}(3)}(1))/(B{topob{bifur}(1)}(end)+B{topob{bifur}(3)}(1)))) )/2;
             
             H123(bifur,1) = ( (( (-eta{topob{bifur}(1)}(1,end)+Z{topob{bifur}(1)}(end,1))...
                +(-eta{topob{bifur}(3)}(1,1)+Z{topob{bifur}(3)}(1,1)) )/2)...
                + (-eta{topob{bifur}(2)}(1,1)+Z{topob{bifur}(2)}(1,1)) )/2;
             
             dzdy(bifur,1) = (eta{topob{bifur}(1)}(1,end) - eta{topob{bifur}(3)}(1,1))/(B{topob{bifur}(2)}(1)/2);
             
             if switchndpt == 0
                 qsy(bifur,1) = ...
                   (-Qtot{topob{bifur}(2)}(1,1)/B{topob{bifur}(2)}(1)) * ( ((Qy(bifur,1)*(-eta{topob{bifur}(2)}(1,1)+Z{topob{bifur}(2)}(1,1)))...
                   /(-Q{topob{bifur}(2)}(1,1)*alw*H123(bifur,1))) - ...
                  (dzdy(bifur,1)*rbolla/sqrt(tau_bolla2(bifur,1))) );
             elseif switchndpt == 1
                 qsy(bifur,1) = ...
                       ( -qb{topob{bifur}(2)}(1) * ( ((Qy(bifur,1)*(-eta{topob{bifur}(2)}(1,1) +Z{topob{bifur}(2)}(1,1)))...
                       /(-Q{topob{bifur}(2)}(1,1)*alw*H123(bifur,1))) - ...
                      (dzdy(bifur,1)*rbolla/sqrt(tau_bolla2(bifur,1))) ) )+ ( -qs{topob{bifur}(2)}(1) * ((Qy(bifur,1)*(-eta{topob{bifur}(2)}(1,1) +Z{topob{bifur}(2)}(1,1)))...
                       /(-Q{topob{bifur}(2)}(1,1)*alw*H123(bifur,1))) );
             end
             
            if calc_TiMor == 0
                            Qsy(bifur,1) = morfac*qsy(bifur,1)*alw*B{topob{bifur}(2)}(1);
            elseif calc_TiMor == 1
                            Qsy(bifur,1) = qsy(bifur,1)*alw*B{topob{bifur}(2)}(1);
            end
            
            Qtot{topob{bifur}(1)}(end,1) = (-Qsy(bifur,1)) + ...
                ((B{topob{bifur}(1)}(end)/(B{topob{bifur}(1)}(end)+B{topob{bifur}(3)}(1))) * Qtot{topob{bifur}(2)}(1,1));

          elseif switchndpt == 2
              rat_Qtot(bifur,1) = ((B{topob{bifur}(1)}(end)/B{topob{bifur}(3)}(1))^(1-k_wang)) * (-Q{topob{bifur}(1)}(end,1)/Q{topob{bifur}(3)}(1,1))^k_wang;
              Qtot{topob{bifur}(1)}(end,1) = rat_Qtot(bifur,1) * (Qtot{topob{bifur}(2)}(1,1)) / (1+rat_Qtot(bifur,1));
                
          end
            Qtot{topob{bifur}(3)}(1,1) = -Qtot{topob{bifur}(2)}(1,1)+Qtot{topob{bifur}(1)}(end,1);
         
         %% third condition - bifurcation 
         elseif Q{topob{bifur}(1)}(end,1)<0 && Q{topob{bifur}(2)}(1,1)>=0 && Q{topob{bifur}(3)}(1,1)<0
          if switchndpt == 0 || switchndpt == 1             
           Qy(bifur,1) = ( -Q{topob{bifur}(1)}(end,1)-Q{topob{bifur}(2)}(1,1)+(Q{topob{bifur}(3)}(1,1)*...
                ((B{topob{bifur}(1)}(end)-B{topob{bifur}(2)}(1))/(B{topob{bifur}(1)}(end)+B{topob{bifur}(2)}(1)))) )/2;
            
            H123(bifur,1) = ( (( (-eta{topob{bifur}(1)}(1,end)+Z{topob{bifur}(1)}(end,1))...
                +(-eta{topob{bifur}(2)}(1,1)+Z{topob{bifur}(2)}(1,1)) )/2)...
                + (-eta{topob{bifur}(3)}(1,1)+Z{topob{bifur}(3)}(1,1)) )/2;
    
            dzdy(bifur,1) = (eta{topob{bifur}(1)}(1,end) - eta{topob{bifur}(2)}(1,1))/(B{topob{bifur}(3)}(1)/2);
   
           if switchndpt == 0     
               qsy(bifur,1) = ...
                   (-Qtot{topob{bifur}(3)}(1,1)/B{topob{bifur}(3)}(1)) * ( ((Qy(bifur,1)*(-eta{topob{bifur}(3)}(1,1)+Z{topob{bifur}(3)}(1,1)))...
                   /(-Q{topob{bifur}(3)}(1,1)*alw*H123(bifur,1))) - ...
                  (dzdy(bifur,1)*rbolla/sqrt(tau_bolla3(bifur,1))) );
           elseif switchndpt == 1
               qsy(bifur,1) = ...
                   (-qb{topob{bifur}(3)}(1) * ( ((Qy(bifur,1)*(-eta{topob{bifur}(3)}(1,1) +Z{topob{bifur}(3)}(1,1)))...
                   /(-Q{topob{bifur}(3)}(1,1)*alw*H123(bifur,1))) - ...
                  (dzdy(bifur,1)*rbolla/sqrt(tau_bolla3(bifur,1))) ) ) + ( -qs{topob{bifur}(3)}(1) * ((Qy(bifur,1)*(-eta{topob{bifur}(3)}(1,1) +Z{topob{bifur}(3)}(1,1)))...
                   /(-Q{topob{bifur}(3)}(1,1)*alw*H123(bifur,1))) );
           end

           if calc_TiMor == 0
                        Qsy(bifur,1) = morfac*qsy(bifur,1)*alw*B{topob{bifur}(3)}(1);
                elseif calc_TiMor == 1
                        Qsy(bifur,1) = qsy(bifur,1)*alw*B{topob{bifur}(3)}(1);
                end
          
           Qtot{topob{bifur}(1)}(end,1) = (-Qsy(bifur,1)) + ...
                ((B{topob{bifur}(1)}(end)/(B{topob{bifur}(1)}(end)+B{topob{bifur}(2)}(1))) * Qtot{topob{bifur}(3)}(1,1));

          elseif switchndpt == 2
              rat_Qtot(bifur,1) = ((B{topob{bifur}(1)}(end)/B{topob{bifur}(2)}(1))^(1-k_wang)) * (-Q{topob{bifur}(1)}(end,1)/Q{topob{bifur}(2)}(1,1))^k_wang;
              Qtot{topob{bifur}(1)}(end,1) = rat_Qtot(bifur,1) * (Qtot{topob{bifur}(3)}(1,1)) / (1+rat_Qtot(bifur,1));
                
          end
            
           Qtot{topob{bifur}(2)}(1,1) = -Qtot{topob{bifur}(3)}(1,1)+Qtot{topob{bifur}(1)}(end,1);
         
         %% forth - sixth conditions -  confluence
         elseif Q{topob{bifur}(1)}(end,1)<0 && Q{topob{bifur}(2)}(1,1)<0 && Q{topob{bifur}(3)}(1,1)<0
             
           Qtot{topob{bifur}(1)}(end,1) = Qtot{topob{bifur}(3)}(1,1)+Qtot{topob{bifur}(2)}(1,1);
         
         elseif Q{topob{bifur}(1)}(end,1)>=0 && Q{topob{bifur}(2)}(1,1)<0 && Q{topob{bifur}(3)}(1,1)>=0
             
           Qtot{topob{bifur}(3)}(1,1) = Qtot{topob{bifur}(1)}(end,1)+(-Qtot{topob{bifur}(2)}(1,1));
             
         elseif Q{topob{bifur}(1)}(end,1)>=0 && Q{topob{bifur}(2)}(1,1)>=0 && Q{topob{bifur}(3)}(1,1)<0
             
           Qtot{topob{bifur}(2)}(1,1) = (-Qtot{topob{bifur}(3)}(1,1))+Qtot{topob{bifur}(1)}(end,1);
             
         end

        
        
 end
 end
 %% moving nodal point - for confluence
 if Nconf>0
 for conf = 1:Nconf
     
     
     %% first condition - bifurcation - for confluencing downstream - main channel is channel 1
        if Q{topoc{conf}(1)}(end,1)>=0 && Q{topoc{conf}(2)}(end,1)<0 && Q{topoc{conf}(3)}(1,1)>=0
          if switchndpt == 0 || switchndpt == 1   
            Qy_conf(conf,1) = ( -Q{topoc{conf}(2)}(end,1)-Q{topoc{conf}(3)}(1,1)-(Q{topoc{conf}(1)}(end,1)*...
                ((B{topoc{conf}(2)}(end)-B{topoc{conf}(3)}(1))/(B{topoc{conf}(2)}(end)+B{topoc{conf}(3)}(1)))) )/2;
            H123_conf(conf,1) = ( (( (-eta{topoc{conf}(2)}(1,end) +Z{topoc{conf}(2)}(end,1))...
                +(-eta{topoc{conf}(3)}(1,1) +Z{topoc{conf}(3)}(1,1)) )/2)...
                + (-eta{topoc{conf}(1)}(1,end) +Z{topoc{conf}(1)}(end,1)) )/2;
  
            dzdy_conf(conf,1) = (eta{topoc{conf}(2)}(1,end) - eta{topoc{conf}(3)}(1,1))/(B{topoc{conf}(1)}(end)/2);
               if switchndpt == 0
               
                   qsy_conf(conf,1) = ...
                       (Qtot{topoc{conf}(1)}(end,1)/B{topoc{conf}(1)}(end)) * ( ((Qy_conf(conf,1)*(-eta{topoc{conf}(1)}(1,end)+Z{topoc{conf}(1)}(end,1)))...
                       /(Q{topoc{conf}(1)}(end,1)*alw*H123_conf(conf,1))) - ...
                      (dzdy_conf(conf,1)*rbolla/sqrt(tau_bolla1_conf(conf,1))) );
              
                elseif  switchndpt == 1
                   qsy_conf(conf,1) =  ...
                       ( qb{topoc{conf}(1)}(end) * ( ((Qy_conf(conf,1)*(-eta{topoc{conf}(1)}(1,end) +Z{topoc{conf}(1)}(end,1)))...
                       /(Q{topoc{conf}(1)}(end,1)*alw*H123_conf(conf,1))) - ...
                      (dzdy_conf(conf,1)*rbolla/sqrt(tau_bolla1_conf(conf,1))) ) ) + ( qs{topoc{conf}(1)}(end)...
                      * ((Qy_conf(conf,1)*(-eta{topoc{conf}(1)}(1,end) +Z{topoc{conf}(1)}(end,1)))...
                       /(Q{topoc{conf}(1)}(end,1)*alw*H123_conf(conf,1))) );
               end
                
                if calc_TiMor == 0
                    Qsy_conf(conf,1) = morfac*qsy_conf(conf,1)*alw*B{topoc{conf}(1)}(end);
                elseif calc_TiMor == 1
                    Qsy_conf(conf,1) = qsy_conf(conf,1)*alw*B{topoc{conf}(1)}(end);
                end

                Qtot{topoc{conf}(2)}(end,1) = (- Qsy_conf(conf,1)) - ...
                     ((B{topoc{conf}(2)}(end)/(B{topoc{conf}(2)}(end)+B{topoc{conf}(3)}(1))) * Qtot{topoc{conf}(1)}(end,1));

          elseif switchndpt == 2
              rat_Qtot_conf(conf,1) = ((B{topoc{conf}(2)}(end)/B{topoc{conf}(3)}(1))^(1-k_wang)) * (-Q{topoc{conf}(2)}(end,1)/Q{topoc{conf}(3)}(1,1))^k_wang;
              Qtot{topoc{conf}(2)}(end,1) = - rat_Qtot_conf(bifur,1) * (Qtot{topoc{conf}(1)}(end,1)) / (1+rat_Qtot_conf(bifur,1));
                
          end
                 
          Qtot{topoc{conf}(3)}(1,1) = Qtot{topoc{conf}(1)}(end,1)+Qtot{topoc{conf}(2)}(end,1);

         %% second condition - bifurcation - for confluencing downstream - main channel is channel 2  
         elseif Q{topoc{conf}(1)}(end,1)<0 && Q{topoc{conf}(2)}(end,1)>=0 && Q{topoc{conf}(3)}(1,1)>=0
           if switchndpt == 0 || switchndpt == 1  
             Qy_conf(conf,1) = ( -Q{topoc{conf}(1)}(end,1)-Q{topoc{conf}(3)}(1,1)-(Q{topoc{conf}(2)}(end,1)*...
                ((B{topoc{conf}(1)}(end)-B{topoc{conf}(3)}(1))/(B{topoc{conf}(1)}(end)+B{topoc{conf}(3)}(1)))) )/2;
             
             H123_conf(conf,1) = ( (( (-eta{topoc{conf}(1)}(1,end) +Z{topoc{conf}(1)}(end,1))...
                +(-eta{topoc{conf}(3)}(1,1) +Z{topoc{conf}(3)}(1,1)) )/2)...
                + (-eta{topoc{conf}(2)}(1,end) +Z{topoc{conf}(2)}(end,1)) )/2;
             
             dzdy_conf(conf,1) = (eta{topoc{conf}(1)}(1,end) - eta{topoc{conf}(3)}(1,1))/(B{topoc{conf}(2)}(end)/2);
             
             if switchndpt == 0
                     qsy_conf(conf,1) = ...
                       (Qtot{topoc{conf}(2)}(end,1)/B{topoc{conf}(2)}(end)) * ( ((Qy_conf(conf,1)*(-eta{topoc{conf}(2)}(1,end)+Z{topoc{conf}(2)}(end,1)))...
                       /(Q{topoc{conf}(2)}(end,1)*alw*H123_conf(conf,1))) - ...
                      (dzdy_conf(conf,1)*rbolla/sqrt(tau_bolla2_conf(conf,1))) );
             elseif switchndpt == 1
                     qsy_conf(conf,1) = ...
                       ( qb{topoc{conf}(2)}(end) * ( ((Qy_conf(conf,1)*(-eta{topoc{conf}(2)}(1,end) +Z{topoc{conf}(2)}(end,1)))...
                       /(Q{topoc{conf}(2)}(end,1)*alw*H123_conf(conf,1))) - ...
                      (dzdy_conf(conf,1)*rbolla/sqrt(tau_bolla2_conf(conf,1))) ) )+ ( qs{topoc{conf}(2)}(end) * ((Qy(conf,1)*(-eta{topoc{conf}(2)}(1,end) +Z{topoc{conf}(2)}(end,1)))...
                       /(Q{topoc{conf}(2)}(end,1)*alw*H123_conf(conf,1))) );
             end
                    if calc_TiMor == 0
                            Qsy_conf(conf,1) = morfac*qsy_conf(conf,1)*alw*B{topoc{conf}(2)}(end);
                    elseif calc_TiMor == 1
                            Qsy_conf(conf,1) = qsy_conf(conf,1)*alw*B{topoc{conf}(2)}(end);
                    end

                    Qtot{topoc{conf}(1)}(end,1) = (-Qsy_conf(conf,1)) - ...
                        ((B{topoc{conf}(1)}(end)/(B{topoc{conf}(1)}(end)+B{topoc{conf}(3)}(1))) * Qtot{topoc{conf}(2)}(end,1));

          elseif switchndpt == 2
              rat_Qtot_conf(conf,1) = ((B{topoc{conf}(1)}(end)/B{topoc{conf}(3)}(1))^(1-k_wang)) * (-Q{topoc{conf}(1)}(end,1)/Q{topoc{conf}(3)}(1,1))^k_wang;
              Qtot{topoc{conf}(1)}(end,1) = - rat_Qtot_conf(bifur,1) * (Qtot{topoc{conf}(2)}(end,1)) / (1+rat_Qtot_conf(bifur,1));
                
          end
                    
          Qtot{topoc{conf}(3)}(1,1) = Qtot{topoc{conf}(2)}(end,1)+Qtot{topoc{conf}(1)}(end,1);
 
         %% third condition - bifurcation - for confluencing downstream - main channel is channel 3 
         elseif Q{topoc{conf}(1)}(end,1)<0 && Q{topoc{conf}(2)}(end,1)<0 && Q{topoc{conf}(3)}(1,1)<0
          if switchndpt == 0 || switchndpt == 1   
           Qy_conf(conf,1) = ( -Q{topoc{conf}(1)}(end,1)+Q{topoc{conf}(2)}(end,1)+(Q{topoc{conf}(3)}(1,1)*...
                ((B{topoc{conf}(1)}(end)-B{topoc{conf}(2)}(end))/(B{topoc{conf}(1)}(end)+B{topoc{conf}(2)}(end)))) )/2;
            
            H123_conf(conf,1) = ( (( (-eta{topoc{conf}(1)}(1,end) +Z{topoc{conf}(1)}(end,1))...
                +(-eta{topoc{conf}(2)}(1,end) +Z{topoc{conf}(2)}(end,1)) )/2)...
                + (-eta{topoc{conf}(3)}(1,1) +Z{topoc{conf}(3)}(1,1)) )/2;
    
            dzdy_conf(conf,1) = (eta{topoc{conf}(1)}(1,end) - eta{topoc{conf}(2)}(1,end))/(B{topoc{conf}(3)}(1)/2);
            
            
           if switchndpt == 0     
               qsy_conf(conf,1) = ...
                   (-Qtot{topoc{conf}(3)}(1,1)/B{topoc{conf}(3)}(1)) * ( ((Qy_conf(conf,1)*(-eta{topoc{conf}(3)}(1,1)+Z{topoc{conf}(3)}(1,1)))...
                   /(-Q{topoc{conf}(3)}(1,1)*alw*H123_conf(conf,1))) - ...
                  (dzdy_conf(conf,1)*rbolla/sqrt(tau_bolla3_conf(conf,1))) );
           elseif switchndpt == 1     
               qsy_conf(conf,1) = ...
                   (-qb{topoc{conf}(3)}(1) * ( ((Qy_conf(conf,1)*(-eta{topoc{conf}(3)}(1,1) +Z{topoc{conf}(3)}(1,1)))...
                   /(-Q{topoc{conf}(3)}(1,1)*alw*H123_conf(conf,1))) - ...
                  (dzdy_conf(conf,1)*rbolla/sqrt(tau_bolla3_conf(conf,1))) ) ) + ( -qs{topoc{conf}(3)}(1) * ((Qy_conf(conf,1)*(-eta{topoc{conf}(3)}(1,1) +Z{topoc{conf}(3)}(1,1)))...
                   /(-Q{topoc{conf}(3)}(1,1)*alw*H123_conf(conf,1))) );
           end
           
                if calc_TiMor == 0
                        Qsy_conf(conf,1) = morfac*qsy_conf(conf,1)*alw*B{topoc{conf}(3)}(1);
                elseif calc_TiMor == 1
                        Qsy_conf(conf,1) = qsy_conf(conf,1)*alw*B{topoc{conf}(3)}(1);
                end
               Qtot{topoc{conf}(1)}(end,1) = (-Qsy_conf(conf,1)) + ...
                    ((B{topoc{conf}(1)}(end)/(B{topoc{conf}(1)}(end)+B{topoc{conf}(2)}(end))) * Qtot{topoc{conf}(3)}(1,1));

         elseif switchndpt == 2
              rat_Qtot_conf(conf,1) = ((B{topoc{conf}(1)}(end)/B{topoc{conf}(2)}(end))^(1-k_wang)) * (Q{topoc{conf}(1)}(end,1)/Q{topoc{conf}(2)}(end,1))^k_wang;
              Qtot{topoc{conf}(1)}(end,1) = rat_Qtot_conf(bifur,1) * (Qtot{topoc{conf}(3)}(1,1)) / (1+rat_Qtot_conf(bifur,1));
                
         end
                
         Qtot{topoc{conf}(2)}(end,1) = Qtot{topoc{conf}(3)}(1,1)-Qtot{topoc{conf}(1)}(end,1);

         %% 4 - 6 condition -  confluence
         elseif Q{topoc{conf}(1)}(end,1)>=0 && Q{topoc{conf}(2)}(end,1)>=0 && Q{topoc{conf}(3)}(1,1)>=0
 
            Qtot{topoc{conf}(3)}(1,1) = Qtot{topoc{conf}(1)}(end,1)+Qtot{topoc{conf}(2)}(end,1);

         elseif Q{topoc{conf}(1)}(end,1)<0 && Q{topoc{conf}(2)}(end,1)>=0 && Q{topoc{conf}(3)}(1,1)<0
 
            Qtot{topoc{conf}(1)}(end,1) = Qtot{topoc{conf}(3)}(1,1)+(-Qtot{topoc{conf}(2)}(end,1));

         elseif Q{topoc{conf}(1)}(end,1)>=0 && Q{topoc{conf}(2)}(end,1)<0 && Q{topoc{conf}(3)}(1,1)<0

            Qtot{topoc{conf}(2)}(end,1) = (Qtot{topoc{conf}(3)}(1,1))+(-Qtot{topoc{conf}(1)}(end,1));

         end

        
        
 end
 end
 

end