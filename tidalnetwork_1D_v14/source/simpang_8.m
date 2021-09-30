%% nodal point calculation
  %% for ratio sediment taking into account at nodal point
     
 %% Shields at the junctions
 if Nbif>0
  for bifur = 1:Nbif
      if switchtransp==1
          tau_bolla1(bifur,1) = Cd*(qf{topob{bifur}(1)}(end)^2)...
                            /(Rr*g*Sizes(channel,4)*((-eta{topob{bifur}(1)}(1,end) +Z{topob{bifur}(1)}(end,2))^2));
          tau_bolla2(bifur,1) = Cd*(qf{topob{bifur}(2)}(1)^2)...
                            /(Rr*g*Sizes(channel,4)*((-eta{topob{bifur}(2)}(1,1) +Z{topob{bifur}(2)}(1,2))^2));
          tau_bolla3(bifur,1) = Cd*(qf{topob{bifur}(3)}(1)^2)...
                          /(Rr*g*Sizes(channel,4)*((-eta{topob{bifur}(3)}(1,1) +Z{topob{bifur}(3)}(1,2))^2));
      elseif switchtransp==2
          if switchnpShi == 0
              tau_bolla1(bifur,1) = tau_prime{1}(end);
              tau_bolla2(bifur,1) = tau_prime{2}(1);
              tau_bolla3(bifur,1) = tau_prime{3}(1);
          elseif switchnpShi == 1
                tau_bolla1(bifur,1) = tau{1}(end);
                tau_bolla2(bifur,1) = tau{2}(1);
                tau_bolla3(bifur,1) = tau{3}(1);
          elseif switchnpShi == 2
                  tau_bolla1(bifur,1) = 1;
                  tau_bolla2(bifur,1) = 1;
                  tau_bolla3(bifur,1) = 1;                              
          end
          
      end
  end
 end
 

  
 if Nconf>0
      for conflu = 1:Nconf
          if switchtransp==1
              tau_bolla1_conf(conflu,1) = Cd*(qf{topoc{conflu}(1)}(end)^2)...
                                /(Rr*g*Sizes(channel,4)*((-eta{topoc{conflu}(1)}(1,end) +Z{topoc{conflu}(1)}(end,2))^2));
              tau_bolla2_conf(conflu,1) = Cd*(qf{topoc{conflu}(2)}(end)^2)...
                                /(Rr*g*Sizes(channel,4)*((-eta{topoc{conflu}(2)}(1,end) +Z{topoc{conflu}(2)}(end,2))^2));
              tau_bolla3_conf(conflu,1) = Cd*(qf{topoc{conflu}(3)}(1)^2)...
                                /(Rr*g*Sizes(channel,4)*((-eta{topoc{conflu}(3)}(1,1) +Z{topoc{conflu}(3)}(1,2))^2));
          elseif switchtransp==2
              if switchnpShi == 0
                  tau_bolla1_conf(conflu,1) = tau_prime{1}(end);
                  tau_bolla2_conf(conflu,1) = tau_prime{2}(end);
                  tau_bolla3_conf(conflu,1) = tau_prime{3}(1);
              elseif switchnpShi == 1
                     tau_bolla1_conf(conflu,1) = tau{1}(end);
                    tau_bolla2_conf(conflu,1) = tau{2}(end);
                    tau_bolla3_conf(conflu,1) = tau{3}(1);
              elseif switchnpShi == 2
                      tau_bolla1_conf(conflu,1) = 1;
                      tau_bolla2_conf(conflu,1) = 1;
                      tau_bolla3_conf(conflu,1) = 1;                              
              end

          end
      end
 end
  
  
 %% moving nodal point - for bifurcation
 if Nbif>0
 for bifur = 1:Nbif
     
     
     %% first condition - bifurcation - for bifurcating downstream - main channel is channel 1
        if Q{topob{bifur}(1)}(end,2)>=0 && Q{topob{bifur}(2)}(1,2)>=0 && Q{topob{bifur}(3)}(1,2)>=0
          if switchndpt == 0 || switchndpt == 1
            Qy(bifur,1) = ( Q{topob{bifur}(2)}(1,2)-Q{topob{bifur}(3)}(1,2)-(Q{topob{bifur}(1)}(end,2)*...
                ((B{topob{bifur}(2)}(1)-B{topob{bifur}(3)}(1))/(B{topob{bifur}(2)}(1)+B{topob{bifur}(3)}(1)))) )/2;
            
            H123(bifur,1) = ( (B{topob{bifur}(2)}(1)*(-eta{topob{bifur}(2)}(1,1) +Z{topob{bifur}(2)}(1,2)))...
                + (B{topob{bifur}(3)}(1)*(-eta{topob{bifur}(3)}(1,1) +Z{topob{bifur}(3)}(1,2)))...
                + (B{topob{bifur}(1)}(end)*(-eta{topob{bifur}(1)}(1,end) +Z{topob{bifur}(1)}(end,2))) )...
                /(B{topob{bifur}(2)}(1)+B{topob{bifur}(3)}(1)+B{topob{bifur}(1)}(end));
            
%             H123(bifur,1) = ( (( (-eta{topob{bifur}(2)}(1,1) +Z{topob{bifur}(2)}(1,2))...
%                 +(-eta{topob{bifur}(3)}(1,1) +Z{topob{bifur}(3)}(1,2)) )/2)...
%                 + (-eta{topob{bifur}(1)}(1,end) +Z{topob{bifur}(1)}(end,2)) )/2;
  
            dzdy(bifur,1) = (eta{topob{bifur}(2)}(1,1) - eta{topob{bifur}(3)}(1,1))/(B{topob{bifur}(1)}(end)/2);
               if switchndpt == 0
               
                   qsy(bifur,1) = ...
                       (Qtot{topob{bifur}(1)}(end,2)/B{topob{bifur}(1)}(end))* ( ((Qy(bifur,1)*(-eta{topob{bifur}(1)}(1,end)+Z{topob{bifur}(1)}(end,2)))...
                       /(Q{topob{bifur}(1)}(end,2)*alw*H123(bifur,1))) - ...
                      (dzdy(bifur,1)*rbolla/sqrt(tau_bolla1(bifur,1))) );
              
                elseif  switchndpt == 1
%                    Qb_cry(bifur,1)  = qb{topob{bifur}(1)}(end)*alw*B{topob{bifur}(1)}(end) * ((Qy(bifur,1)*(-eta{topob{bifur}(1)}(1,end) +Z{topob{bifur}(1)}(end,2)))...
%                                         /(Q{topob{bifur}(1)}(end,2)*alw*H123(bifur,1)));                                %bedload cross flow effect
%                    Qb_slef(bifur,1) = qb{topob{bifur}(1)}(end)*alw*B{topob{bifur}(1)}(end) * (dzdy(bifur,1)*rbolla/sqrt(tau_bolla1(bifur,1)));      %bedslope effect
%                    Qs_cry(bifur,1)  = ( qs{topob{bifur}(1)}(end)*alw*B{topob{bifur}(1)}(end)...
%                                       * ((Qy(bifur,1)*(-eta{topob{bifur}(1)}(1,end) +Z{topob{bifur}(1)}(end,2)))...
%                                        /(Q{topob{bifur}(1)}(end,2)*alw*H123(bifur,1))) );                               %suspension cross flow
                                   
                   qsyB(bifur,1) =  ...
                       ( qb{topob{bifur}(1)}(end) * ( ((Qy(bifur,1)*(-eta{topob{bifur}(1)}(1,end) +Z{topob{bifur}(1)}(end,2)))...
                       /(Q{topob{bifur}(1)}(end,2)*alw*H123(bifur,1))) - ...
                      (dzdy(bifur,1)*rbolla/sqrt(tau_bolla1(bifur,1))) ) );
                  qsyS(bifur,1) =( qs{topob{bifur}(1)}(end)...
                      * ((Qy(bifur,1)*(-eta{topob{bifur}(1)}(1,end) +Z{topob{bifur}(1)}(end,2)))...
                       /(Q{topob{bifur}(1)}(end,2)*alw*H123(bifur,1))) );
                   qsy(bifur,1) =  qsyB(bifur,1) + qsyS(bifur,1);
               end
                
                if calc_TiMor == 0
                    if switchndpt == 0
                        Qsy(bifur,1) = qsy(bifur,1)*alw*B{topob{bifur}(1)}(end);
                    elseif  switchndpt == 1
                        QsyB(bifur,1) = morfac*qsyB(bifur,1)*alw*B{topob{bifur}(1)}(end);
                        QsyS(bifur,1) = morfac*qsyS(bifur,1)*alw*B{topob{bifur}(1)}(end);
                        Qsy(bifur,1) = QsyB(bifur,1)+QsyS(bifur,1);
                    end
                elseif calc_TiMor == 1
                    Qsy(bifur,1) = qsy(bifur,1)*alw*B{topob{bifur}(1)}(end);
                end

                Qtot{topob{bifur}(2)}(1,2) = Qsy(bifur,1) + ...
                     ((B{topob{bifur}(2)}(1)/(B{topob{bifur}(2)}(1)+B{topob{bifur}(3)}(1))) * Qtot{topob{bifur}(1)}(end,2));
          
          elseif switchndpt == 2
              rat_Qtot(bifur,1) = ((B{topob{bifur}(2)}(1)/B{topob{bifur}(3)}(1))^(1-k_wang)) * (Q{topob{bifur}(2)}(1,2)/Q{topob{bifur}(3)}(1,2))^k_wang;
              Qtot{topob{bifur}(2)}(1,2) = rat_Qtot(bifur,1) * (Qtot{topob{bifur}(1)}(end,2)) / (1+rat_Qtot(bifur,1));
                
          end

          Qtot{topob{bifur}(3)}(1,2) = Qtot{topob{bifur}(1)}(end,2)-Qtot{topob{bifur}(2)}(1,2);

         %% second condition - bifurcation - for bifurcating downstream - main channel is channel 2   
         elseif Q{topob{bifur}(1)}(end,2)<0 && Q{topob{bifur}(2)}(1,2)<0 && Q{topob{bifur}(3)}(1,2)>=0
          if switchndpt == 0 || switchndpt == 1             
             Qy(bifur,1) = ( -Q{topob{bifur}(1)}(end,2)-Q{topob{bifur}(3)}(1,2)+(Q{topob{bifur}(2)}(1,2)*...
                ((B{topob{bifur}(1)}(end)-B{topob{bifur}(3)}(1))/(B{topob{bifur}(1)}(end)+B{topob{bifur}(3)}(1)))) )/2;
             
            H123(bifur,1) = ( (B{topob{bifur}(2)}(1)*(-eta{topob{bifur}(2)}(1,1) +Z{topob{bifur}(2)}(1,2)))...
                + (B{topob{bifur}(3)}(1)*(-eta{topob{bifur}(3)}(1,1) +Z{topob{bifur}(3)}(1,2)))...
                + (B{topob{bifur}(1)}(end)*(-eta{topob{bifur}(1)}(1,end) +Z{topob{bifur}(1)}(end,2))) )...
                /(B{topob{bifur}(2)}(1)+B{topob{bifur}(3)}(1)+B{topob{bifur}(1)}(end));
            
            
%              H123(bifur,1) = ( (( (-eta{topob{bifur}(1)}(1,end) +Z{topob{bifur}(1)}(end,2))...
%                 +(-eta{topob{bifur}(3)}(1,1) +Z{topob{bifur}(3)}(1,2)) )/2)...
%                 + (-eta{topob{bifur}(2)}(1,1) +Z{topob{bifur}(2)}(1,2)) )/2;
             
             dzdy(bifur,1) = (eta{topob{bifur}(1)}(1,end) - eta{topob{bifur}(3)}(1,1))/(B{topob{bifur}(2)}(1)/2);
             
             if switchndpt == 0
                     qsy(bifur,1) = ...
                       (-Qtot{topob{bifur}(2)}(1,2)/B{topob{bifur}(2)}(1)) * ( ((Qy(bifur,1)*(-eta{topob{bifur}(2)}(1,1)+Z{topob{bifur}(2)}(1,2)))...
                       /(-Q{topob{bifur}(2)}(1,2)*alw*H123(bifur,1))) - ...
                      (dzdy(bifur,1)*rbolla/sqrt(tau_bolla2(bifur,1))) );
             elseif switchndpt == 1
                     qsy(bifur,1) = ...
                       ( -qb{topob{bifur}(2)}(1) * ( ((Qy(bifur,1)*(-eta{topob{bifur}(2)}(1,1) +Z{topob{bifur}(2)}(1,2)))...
                       /(-Q{topob{bifur}(2)}(1,2)*alw*H123(bifur,1))) - ...
                      (dzdy(bifur,1)*rbolla/sqrt(tau_bolla2(bifur,1))) ) )+ ( -qs{topob{bifur}(2)}(1) * ((Qy(bifur,1)*(-eta{topob{bifur}(2)}(1,1) +Z{topob{bifur}(2)}(1,2)))...
                       /(-Q{topob{bifur}(2)}(1,2)*alw*H123(bifur,1))) );
             end
                    if calc_TiMor == 0
                        if switchndpt == 0
                            Qsy(bifur,1) = qsy(bifur,1)*alw*B{topob{bifur}(2)}(1);
                        elseif  switchndpt == 1    
                            Qsy(bifur,1) = morfac*qsy(bifur,1)*alw*B{topob{bifur}(2)}(1);
                        end
                    elseif calc_TiMor == 1
                            Qsy(bifur,1) = qsy(bifur,1)*alw*B{topob{bifur}(2)}(1);
                    end

                    Qtot{topob{bifur}(1)}(end,2) = (-Qsy(bifur,1)) + ...
                        ((B{topob{bifur}(1)}(end)/(B{topob{bifur}(1)}(end)+B{topob{bifur}(3)}(1))) * Qtot{topob{bifur}(2)}(1,2));
                    
          elseif switchndpt == 2
              rat_Qtot(bifur,1) = ((B{topob{bifur}(1)}(end)/B{topob{bifur}(3)}(1))^(1-k_wang)) * (-Q{topob{bifur}(1)}(end,2)/Q{topob{bifur}(3)}(1,2))^k_wang;
              Qtot{topob{bifur}(1)}(end,2) = rat_Qtot(bifur,1) * (Qtot{topob{bifur}(2)}(1,2)) / (1+rat_Qtot(bifur,1));
                
          end                    

          Qtot{topob{bifur}(3)}(1,2) = -Qtot{topob{bifur}(2)}(1,2)+Qtot{topob{bifur}(1)}(end,2);
 
         %% third condition - bifurcation - for bifurcating downstream - main channel is channel 3
         elseif Q{topob{bifur}(1)}(end,2)<0 && Q{topob{bifur}(2)}(1,2)>=0 && Q{topob{bifur}(3)}(1,2)<0
          if switchndpt == 0 || switchndpt == 1             
           Qy(bifur,1) = ( -Q{topob{bifur}(1)}(end,2)-Q{topob{bifur}(2)}(1,2)+(Q{topob{bifur}(3)}(1,2)*...
                ((B{topob{bifur}(1)}(end)-B{topob{bifur}(2)}(1))/(B{topob{bifur}(1)}(end)+B{topob{bifur}(2)}(1)))) )/2;
            
            H123(bifur,1) = ( (B{topob{bifur}(2)}(1)*(-eta{topob{bifur}(2)}(1,1) +Z{topob{bifur}(2)}(1,2)))...
                + (B{topob{bifur}(3)}(1)*(-eta{topob{bifur}(3)}(1,1) +Z{topob{bifur}(3)}(1,2)))...
                + (B{topob{bifur}(1)}(end)*(-eta{topob{bifur}(1)}(1,end) +Z{topob{bifur}(1)}(end,2))) )...
                /(B{topob{bifur}(2)}(1)+B{topob{bifur}(3)}(1)+B{topob{bifur}(1)}(end));
            
            
%             H123(bifur,1) = ( (( (-eta{topob{bifur}(1)}(1,end) +Z{topob{bifur}(1)}(end,2))...
%                 +(-eta{topob{bifur}(2)}(1,1) +Z{topob{bifur}(2)}(1,2)) )/2)...
%                 + (-eta{topob{bifur}(3)}(1,1) +Z{topob{bifur}(3)}(1,2)) )/2;
    
            dzdy(bifur,1) = (eta{topob{bifur}(1)}(1,end) - eta{topob{bifur}(2)}(1,1))/(B{topob{bifur}(3)}(1)/2);
            
            
           if switchndpt == 0     
               qsy(bifur,1) = ...
                   (-Qtot{topob{bifur}(3)}(1,2)/B{topob{bifur}(3)}(1)) * ( ((Qy(bifur,1)*(-eta{topob{bifur}(3)}(1,1)+Z{topob{bifur}(3)}(1,2)))...
                   /(-Q{topob{bifur}(3)}(1,2)*alw*H123(bifur,1))) - ...
                  (dzdy(bifur,1)*rbolla/sqrt(tau_bolla3(bifur,1))) );
           elseif switchndpt == 1     
               qsy(bifur,1) = ...
                   (-qb{topob{bifur}(3)}(1) * ( ((Qy(bifur,1)*(-eta{topob{bifur}(3)}(1,1) +Z{topob{bifur}(3)}(1,2)))...
                   /(-Q{topob{bifur}(3)}(1,2)*alw*H123(bifur,1))) - ...
                  (dzdy(bifur,1)*rbolla/sqrt(tau_bolla3(bifur,1))) ) ) + ( -qs{topob{bifur}(3)}(1) * ((Qy(bifur,1)*(-eta{topob{bifur}(3)}(1,1) +Z{topob{bifur}(3)}(1,2)))...
                   /(-Q{topob{bifur}(3)}(1,2)*alw*H123(bifur,1))) );
           end
           
                if calc_TiMor == 0
                    if switchndpt == 0
                        Qsy(bifur,1) = qsy(bifur,1)*alw*B{topob{bifur}(3)}(1);
                    elseif  switchndpt == 1  
                        Qsy(bifur,1) = morfac*qsy(bifur,1)*alw*B{topob{bifur}(3)}(1);
                    end
                elseif calc_TiMor == 1
                        Qsy(bifur,1) = qsy(bifur,1)*alw*B{topob{bifur}(3)}(1);
                end
               Qtot{topob{bifur}(1)}(end,2) = (-Qsy(bifur,1)) + ...
                    ((B{topob{bifur}(1)}(end)/(B{topob{bifur}(1)}(end)+B{topob{bifur}(2)}(1))) * Qtot{topob{bifur}(3)}(1,2));
         
          elseif switchndpt == 2
              rat_Qtot(bifur,1) = ((B{topob{bifur}(1)}(end)/B{topob{bifur}(2)}(1))^(1-k_wang)) * (-Q{topob{bifur}(1)}(end,2)/Q{topob{bifur}(2)}(1,2))^k_wang;
              Qtot{topob{bifur}(1)}(end,2) = rat_Qtot(bifur,1) * (Qtot{topob{bifur}(3)}(1,2)) / (1+rat_Qtot(bifur,1));
                
          end
          
          Qtot{topob{bifur}(2)}(1,2) = -Qtot{topob{bifur}(3)}(1,2)+Qtot{topob{bifur}(1)}(end,2);

         %% 4 - 6 condition -  confluence
         elseif Q{topob{bifur}(1)}(end,2)<0 && Q{topob{bifur}(2)}(1,2)<0 && Q{topob{bifur}(3)}(1,2)<0
 
            Qtot{topob{bifur}(1)}(end,2) = Qtot{topob{bifur}(3)}(1,2)+Qtot{topob{bifur}(2)}(1,2);

         elseif Q{topob{bifur}(1)}(end,2)>=0 && Q{topob{bifur}(2)}(1,2)<0 && Q{topob{bifur}(3)}(1,2)>=0
 
            Qtot{topob{bifur}(3)}(1,2) = Qtot{topob{bifur}(1)}(end,2)+(-Qtot{topob{bifur}(2)}(1,2));

         elseif Q{topob{bifur}(1)}(end,2)>=0 && Q{topob{bifur}(2)}(1,2)>=0 && Q{topob{bifur}(3)}(1,2)<0

            Qtot{topob{bifur}(2)}(1,2) = (-Qtot{topob{bifur}(3)}(1,2))+Qtot{topob{bifur}(1)}(end,2);

         end

        
        
 end
 end
 
 %% moving nodal point - for confluence
 if Nconf>0
 for conf = 1:Nconf
     
     
     %% first condition - bifurcation - for confluencing downstream - main channel is channel 1
        if Q{topoc{conf}(1)}(end,2)>=0 && Q{topoc{conf}(2)}(end,2)<0 && Q{topoc{conf}(3)}(1,2)>=0
          if switchndpt == 0 || switchndpt == 1             
            Qy_conf(conf,1) = ( -Q{topoc{conf}(2)}(end,2)-Q{topoc{conf}(3)}(1,2)-(Q{topoc{conf}(1)}(end,2)*...
                ((B{topoc{conf}(2)}(end)-B{topoc{conf}(3)}(1))/(B{topoc{conf}(2)}(end)+B{topoc{conf}(3)}(1)))) )/2;
            
            H123_conf(conf,1) = ( (B{topoc{conf}(2)}(end)*(-eta{topoc{conf}(2)}(1,end) +Z{topoc{conf}(2)}(end,2)))...
                + (B{topoc{conf}(3)}(1)*(-eta{topoc{conf}(3)}(1,1) +Z{topoc{conf}(3)}(1,2)))...
                + B{topoc{conf}(1)}(end)*(-eta{topoc{conf}(1)}(1,end) +Z{topoc{conf}(1)}(end,2)) )...
                /(B{topoc{conf}(2)}(end)+B{topoc{conf}(3)}(1)+B{topoc{conf}(1)}(end));
            
%             H123_conf(conf,1) = ( (( (-eta{topoc{conf}(2)}(1,end) +Z{topoc{conf}(2)}(end,2))...
%                 +(-eta{topoc{conf}(3)}(1,1) +Z{topoc{conf}(3)}(1,2)) )/2)...
%                 + (-eta{topoc{conf}(1)}(1,end) +Z{topoc{conf}(1)}(end,2)) )/2;
  
            dzdy_conf(conf,1) = (eta{topoc{conf}(2)}(1,end) - eta{topoc{conf}(3)}(1,1))/(B{topoc{conf}(1)}(end)/2);
               if switchndpt == 0
               
                   qsy_conf(conf,1) = ...
                       (Qtot{topoc{conf}(1)}(end,2)/B{topoc{conf}(1)}(end)) * ( ((Qy_conf(conf,1)*(-eta{topoc{conf}(1)}(1,end)+Z{topoc{conf}(1)}(end,2)))...
                       /(Q{topoc{conf}(1)}(end,2)*alw*H123_conf(conf,1))) - ...
                      (dzdy_conf(conf,1)*rbolla/sqrt(tau_bolla1_conf(conf,1))) );
              
                elseif  switchndpt == 1
                   qsy_conf(conf,1) =  ...
                       ( qb{topoc{conf}(1)}(end) * ( ((Qy_conf(conf,1)*(-eta{topoc{conf}(1)}(1,end) +Z{topoc{conf}(1)}(end,2)))...
                       /(Q{topoc{conf}(1)}(end,2)*alw*H123_conf(conf,1))) - ...
                      (dzdy_conf(conf,1)*rbolla/sqrt(tau_bolla1_conf(conf,1))) ) ) + ( qs{topoc{conf}(1)}(end)...
                      * ((Qy_conf(conf,1)*(-eta{topoc{conf}(1)}(1,end) +Z{topoc{conf}(1)}(end,2)))...
                       /(Q{topoc{conf}(1)}(end,2)*alw*H123_conf(conf,1))) );
               end
                
                if calc_TiMor == 0
                    if switchndpt == 0
                        Qsy_conf(conf,1) = qsy_conf(conf,1)*alw*B{topoc{conf}(1)}(end);
                    elseif  switchndpt == 1
                        Qsy_conf(conf,1) = morfac*qsy_conf(conf,1)*alw*B{topoc{conf}(1)}(end);
                    end
                elseif calc_TiMor == 1
                    Qsy_conf(conf,1) = qsy_conf(conf,1)*alw*B{topoc{conf}(1)}(end);
                end

                Qtot{topoc{conf}(2)}(end,2) = (- Qsy_conf(conf,1)) - ...
                     ((B{topoc{conf}(2)}(end)/(B{topoc{conf}(2)}(end)+B{topoc{conf}(3)}(1))) * Qtot{topoc{conf}(1)}(end,2));

          elseif switchndpt == 2
              rat_Qtot_conf(conf,1) = ((B{topoc{conf}(2)}(end)/B{topoc{conf}(3)}(1))^(1-k_wang)) * (-Q{topoc{conf}(2)}(end,2)/Q{topoc{conf}(3)}(1,2))^k_wang;
              Qtot{topoc{conf}(2)}(end,2) = - rat_Qtot_conf(bifur,1) * (Qtot{topoc{conf}(1)}(end,2)) / (1+rat_Qtot_conf(bifur,1));
                
          end
          
          Qtot{topoc{conf}(3)}(1,2) = Qtot{topoc{conf}(1)}(end,2)+Qtot{topoc{conf}(2)}(end,2);

         %% second condition - bifurcation - for confluencing downstream - main channel is channel 2  
         elseif Q{topoc{conf}(1)}(end,2)<0 && Q{topoc{conf}(2)}(end,2)>=0 && Q{topoc{conf}(3)}(1,2)>=0
          if switchndpt == 0 || switchndpt == 1             
             Qy_conf(conf,1) = ( -Q{topoc{conf}(1)}(end,2)-Q{topoc{conf}(3)}(1,2)-(Q{topoc{conf}(2)}(end,2)*...
                ((B{topoc{conf}(1)}(end)-B{topoc{conf}(3)}(1))/(B{topoc{conf}(1)}(end)+B{topoc{conf}(3)}(1)))) )/2;
             
            H123_conf(conf,1) = ( (B{topoc{conf}(2)}(end)*(-eta{topoc{conf}(2)}(1,end) +Z{topoc{conf}(2)}(end,2)))...
                + (B{topoc{conf}(3)}(1)*(-eta{topoc{conf}(3)}(1,1) +Z{topoc{conf}(3)}(1,2)))...
                + B{topoc{conf}(1)}(end)*(-eta{topoc{conf}(1)}(1,end) +Z{topoc{conf}(1)}(end,2)) )...
                /(B{topoc{conf}(2)}(end)+B{topoc{conf}(3)}(1)+B{topoc{conf}(1)}(end));
            
%              H123_conf(conf,1) = ( (( (-eta{topoc{conf}(1)}(1,end) +Z{topoc{conf}(1)}(end,2))...
%                 +(-eta{topoc{conf}(3)}(1,1) +Z{topoc{conf}(3)}(1,2)) )/2)...
%                 + (-eta{topoc{conf}(2)}(1,end) +Z{topoc{conf}(2)}(end,2)) )/2;
             
             dzdy_conf(conf,1) = (eta{topoc{conf}(1)}(1,end) - eta{topoc{conf}(3)}(1,1))/(B{topoc{conf}(2)}(end)/2);
             
             if switchndpt == 0
                     qsy_conf(conf,1) = ...
                       (Qtot{topoc{conf}(2)}(end,2)/B{topoc{conf}(2)}(end)) * ( ((Qy_conf(conf,1)*(-eta{topoc{conf}(2)}(1,end)+Z{topoc{conf}(2)}(end,2)))...
                       /(Q{topoc{conf}(2)}(end,2)*alw*H123_conf(conf,1))) - ...
                      (dzdy_conf(conf,1)*rbolla/sqrt(tau_bolla2_conf(conf,1))) );
             elseif switchndpt == 1
                     qsy_conf(conf,1) = ...
                       ( qb{topoc{conf}(2)}(end) * ( ((Qy_conf(conf,1)*(-eta{topoc{conf}(2)}(1,end) +Z{topoc{conf}(2)}(end,2)))...
                       /(Q{topoc{conf}(2)}(end,2)*alw*H123_conf(conf,1))) - ...
                      (dzdy_conf(conf,1)*rbolla/sqrt(tau_bolla2_conf(conf,1))) ) )+ ( qs{topoc{conf}(2)}(end) * ((Qy_conf(conf,1)*(-eta{topoc{conf}(2)}(1,end) +Z{topoc{conf}(2)}(end,2)))...
                       /(Q{topoc{conf}(2)}(end,2)*alw*H123_conf(conf,1))) );
             end
                    if calc_TiMor == 0
                        if switchndpt == 0
                            Qsy_conf(conf,1) = qsy_conf(conf,1)*alw*B{topoc{conf}(2)}(end);
                        elseif  switchndpt == 1    
                            Qsy_conf(conf,1) = morfac*qsy_conf(conf,1)*alw*B{topoc{conf}(2)}(end);
                        end
                    elseif calc_TiMor == 1
                            Qsy_conf(conf,1) = qsy_conf(conf,1)*alw*B{topoc{conf}(2)}(end);
                    end

                    Qtot{topoc{conf}(1)}(end,2) = (-Qsy_conf(conf,1)) - ...
                        ((B{topoc{conf}(1)}(end)/(B{topoc{conf}(1)}(end)+B{topoc{conf}(3)}(1))) * Qtot{topoc{conf}(2)}(end,2));

          elseif switchndpt == 2
              rat_Qtot_conf(conf,1) = ((B{topoc{conf}(1)}(end)/B{topoc{conf}(3)}(1))^(1-k_wang)) * (-Q{topoc{conf}(1)}(end,2)/Q{topoc{conf}(3)}(1,2))^k_wang;
              Qtot{topoc{conf}(1)}(end,2) = - rat_Qtot_conf(bifur,1) * (Qtot{topoc{conf}(2)}(end,2)) / (1+rat_Qtot_conf(bifur,1));
                
          end
          
          Qtot{topoc{conf}(3)}(1,2) = Qtot{topoc{conf}(2)}(end,2)+Qtot{topoc{conf}(1)}(end,2);
 
         %% third condition - bifurcation - for confluencing downstream - main channel is channel 3 
         elseif Q{topoc{conf}(1)}(end,2)<0 && Q{topoc{conf}(2)}(end,2)<0 && Q{topoc{conf}(3)}(1,2)<0
          if switchndpt == 0 || switchndpt == 1   
           Qy_conf(conf,1) = ( -Q{topoc{conf}(1)}(end,2)+Q{topoc{conf}(2)}(end,2)+(Q{topoc{conf}(3)}(1,2)*...
                ((B{topoc{conf}(1)}(end)-B{topoc{conf}(2)}(end))/(B{topoc{conf}(1)}(end)+B{topoc{conf}(2)}(end)))) )/2;
            
           H123_conf(conf,1) = ( (B{topoc{conf}(2)}(end)*(-eta{topoc{conf}(2)}(1,end) +Z{topoc{conf}(2)}(end,2)))...
                + (B{topoc{conf}(3)}(1)*(-eta{topoc{conf}(3)}(1,1) +Z{topoc{conf}(3)}(1,2)))...
                + B{topoc{conf}(1)}(end)*(-eta{topoc{conf}(1)}(1,end) +Z{topoc{conf}(1)}(end,2)) )...
                /(B{topoc{conf}(2)}(end)+B{topoc{conf}(3)}(1)+B{topoc{conf}(1)}(end));
            
%             H123_conf(conf,1) = ( (( (-eta{topoc{conf}(1)}(1,end) +Z{topoc{conf}(1)}(end,2))...
%                 +(-eta{topoc{conf}(2)}(1,end) +Z{topoc{conf}(2)}(end,2)) )/2)...
%                 + (-eta{topoc{conf}(3)}(1,1) +Z{topoc{conf}(3)}(1,2)) )/2;
    
            dzdy_conf(conf,1) = (eta{topoc{conf}(1)}(1,end) - eta{topoc{conf}(2)}(1,end))/(B{topoc{conf}(3)}(1)/2);
            
            
           if switchndpt == 0     
               qsy_conf(conf,1) = ...
                   (-Qtot{topoc{conf}(3)}(1,2)/B{topoc{conf}(3)}(1)) * ( ((Qy_conf(conf,1)*(-eta{topoc{conf}(3)}(1,1)+Z{topoc{conf}(3)}(1,2)))...
                   /(-Q{topoc{conf}(3)}(1,2)*alw*H123_conf(conf,1))) - ...
                  (dzdy_conf(conf,1)*rbolla/sqrt(tau_bolla3_conf(conf,1))) );
           elseif switchndpt == 1     
               qsy_conf(conf,1) = ...
                   (-qb{topoc{conf}(3)}(1) * ( ((Qy_conf(conf,1)*(-eta{topoc{conf}(3)}(1,1) +Z{topoc{conf}(3)}(1,2)))...
                   /(-Q{topoc{conf}(3)}(1,2)*alw*H123_conf(conf,1))) - ...
                  (dzdy_conf(conf,1)*rbolla/sqrt(tau_bolla3_conf(conf,1))) ) ) + ( -qs{topoc{conf}(3)}(1) * ((Qy_conf(conf,1)*(-eta{topoc{conf}(3)}(1,1) +Z{topoc{conf}(3)}(1,2)))...
                   /(-Q{topoc{conf}(3)}(1,2)*alw*H123_conf(conf,1))) );
           end
           
                if calc_TiMor == 0
                    if switchndpt == 0
                        Qsy_conf(conf,1) = qsy_conf(conf,1)*alw*B{topoc{conf}(3)}(1);
                    elseif switchndpt == 1     
                        Qsy_conf(conf,1) = morfac*qsy_conf(conf,1)*alw*B{topoc{conf}(3)}(1);
                    end
                elseif calc_TiMor == 1
                        Qsy_conf(conf,1) = qsy_conf(conf,1)*alw*B{topoc{conf}(3)}(1);
                end
               Qtot{topoc{conf}(1)}(end,2) = (-Qsy_conf(conf,1)) + ...
                    ((B{topoc{conf}(1)}(end)/(B{topoc{conf}(1)}(end)+B{topoc{conf}(2)}(end))) * Qtot{topoc{conf}(3)}(1,2));

         elseif switchndpt == 2
              rat_Qtot_conf(conf,1) = ((B{topoc{conf}(1)}(end)/B{topoc{conf}(2)}(end))^(1-k_wang)) * (Q{topoc{conf}(1)}(end,2)/Q{topoc{conf}(2)}(end,2))^k_wang;
              Qtot{topoc{conf}(1)}(end,2) = rat_Qtot_conf(bifur,1) * (Qtot{topoc{conf}(3)}(1,2)) / (1+rat_Qtot_conf(bifur,1));
                
          end
                
               Qtot{topoc{conf}(2)}(end,2) = Qtot{topoc{conf}(3)}(1,2)-Qtot{topoc{conf}(1)}(end,2);

         %% 4 - 6 condition -  confluence
         elseif Q{topoc{conf}(1)}(end,2)>=0 && Q{topoc{conf}(2)}(end,2)>=0 && Q{topoc{conf}(3)}(1,2)>=0
 
            Qtot{topoc{conf}(3)}(1,2) = Qtot{topoc{conf}(1)}(end,2)+Qtot{topoc{conf}(2)}(end,2);

         elseif Q{topoc{conf}(1)}(end,2)<0 && Q{topoc{conf}(2)}(end,2)>=0 && Q{topoc{conf}(3)}(1,2)<0
 
            Qtot{topoc{conf}(1)}(end,2) = Qtot{topoc{conf}(3)}(1,2)+(-Qtot{topoc{conf}(2)}(end,2));

         elseif Q{topoc{conf}(1)}(end,2)>=0 && Q{topoc{conf}(2)}(end,2)<0 && Q{topoc{conf}(3)}(1,2)<0

            Qtot{topoc{conf}(2)}(1,2) = (Qtot{topoc{conf}(3)}(1,2))+(-Qtot{topoc{conf}(1)}(end,2));

         end

        
        
 end
 end