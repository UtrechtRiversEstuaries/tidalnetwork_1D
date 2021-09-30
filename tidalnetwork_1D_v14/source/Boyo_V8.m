

%% Main script
%   by AP Iwantoro

%% necessary file
% - ini_dia_1
% - ini_dia_matrix_1
% - ini_dia_sedtan_simpang_1
% - Aquae_2
% - sedtan_1
% - simpang_1
% - bedte_1
% - harmfit
%% initialisation %%

ini_dia_matrix_3


%% calculation
waittxt = [' please wait ;3 '];
h = waitbar(0,waittxt,'Name',name);
for j=1:Nt-1
    waitbar(j/(Nt-1),h);
    %prescribed downstream boundary
    Seaboundary_1

    %prescribed upstream boundary 
    for channel= 1:Nc
        if any(channel == Upbound) 
        Q{channel}(1,:) = discharge(find(Upbound==channel));                      % river discharge.
        end
    end
    
      % move the unknown column to the known column
      for channel=1:Nc
        if j~=1
            if any(channel == Upbound)
            Z{channel}(1:Nx{channel},1)=Z{channel}(1:Nx{channel},2);
            Q{channel}(2:Nx{channel},1)=Q{channel}(2:Nx{channel},2);
            elseif any(channel == seabound)
            Z{channel}(1:Nx{channel}-1,1)=Z{channel}(1:Nx{channel}-1,2);
            Q{channel}(1:Nx{channel},1)=Q{channel}(1:Nx{channel},2);
            else
            Z{channel}(1:Nx{channel},1)=Z{channel}(1:Nx{channel},2);
            Q{channel}(1:Nx{channel},1)=Q{channel}(1:Nx{channel},2);   
            end
            if calc_sedtan ==1
                Qtot{channel}(:,1) = Qtot{channel}(:,2);      %total transport
                if calc_bedte ==1
                 eta{channel}(1,:) = eta{channel}(2,:);      %bed change
                end
            end
        end
 
        %store data
        if any(t(j)==t_store)
              Z_store{channel}(:,find(t_store==t(j))) = Z{channel}(:,1);
              Q_store{channel}(:,find(t_store==t(j))) = Q{channel}(:,1);
              if calc_sedtan == 1
%                 tau_store{channel}(:,find(t_store==t(j))) = tau{channel};
                tau_prime_store{channel}(:,find(t_store==t(j))) = tau_prime{channel};
                if Nbif>0
                    tau_store1(:,find(t_store==t(j))) = tau_bolla1;
                    tau_store2(:,find(t_store==t(j))) = tau_bolla2;
                    tau_store3(:,find(t_store==t(j))) = tau_bolla3;
                end
                
                if Nconf>0
                    tau_store1_conf(:,find(t_store==t(j))) = tau_bolla1_conf;
                    tau_store2_conf(:,find(t_store==t(j))) = tau_bolla2_conf;
                    tau_store3_conf(:,find(t_store==t(j))) = tau_bolla3_conf;
                end
                
                ustar_store{channel}(:,find(t_store==t(j))) = ustar{channel};
                Qs_store{channel}(:,find(t_store==t(j))) = Qs{channel};
                Qb_store{channel}(:,find(t_store==t(j))) = Qb{channel};
                Qtot_store{channel}(:,find(t_store==t(j)))= Qtot{channel}(:,1);              
                if calc_bedte == 1
                    eta_store{channel}(:,find(t_store==t(j))) = eta{channel}(1,:)';
                end
              end
        end
      end
      
        %harmonic analysis at the end of morphodynamic simulation
        
     if harm_ini_final == 1        
        for channel = 1:Nc
            if any(t(j)==t_harm)
                     Z_harm{channel}(:,find(t_harm==t(j))) = Z{channel}(:,1);
                     Q_harm{channel}(:,find(t_harm==t(j))) = Q{channel}(:,1);
                     tau_harm{channel}(:,find(t_harm==t(j))) = tau_prime{channel}(:,1);
                     if calc_sedtan ==1
                        Qtot_harm{channel}(:,find(t_harm==t(j))) = Qtot{channel}(:,1);
                        Qs_harm{channel}(:,find(t_harm==t(j))) = Qs{channel};
                        Qb_harm{channel}(:,find(t_harm==t(j))) = Qb{channel};                        
                        za_harm{channel}(:,find(t_harm==t(j))) = za{channel};
                        T_harm{channel}(:,find(t_harm==t(j))) = T{channel};
                        Ca_harm{channel}(:,find(t_harm==t(j))) = Ca{channel};
                        F_harm{channel}(:,find(t_harm==t(j))) = F{channel};
                        if calc_bedte ==1
                            eta_harm{channel}(:,find(t_harm==t(j))) = eta{channel}(1,:)';
                        end
                    end
            end
        end
      
     end
    
    %% flow calculation
    Aquae_2
    
%    check if channel has been avulsed
    for channel =1:Nc
    
        if isreal(Z{channel}(:,2))==0
            day = t(j)/(3600*24)
            Fdl = findall(0,'type','figure','tag','TMWWaitbar');
            delete(Fdl)
            error('ENG: congratulation one of your channel has been abandoned.IDN:selokan anda ada yang mampet')
            return
        end
        
        
    end
    % storing stuff
     %store cross section
        for channel = 1:Nc
            if any(t(j)==t_store)
                  A_store{channel}(:,find(t_store==t(j))) = A{channel};
            end
            
            if harm_ini_final == 1
                if any(t(j)==t_harm)
                      A_harm{channel}(:,find(t_harm==t(j))) = A{channel};
                end
            end
            
        end
    %update channel properties for coupling with sedtan
    for channel = 1:Nc
        A{channel}(:,1)=B{channel}'.*(-eta{channel}(1,:)'+Z{channel}(:,2));
        
        if switchwetper == 1
            P{channel}(:,1)=B{channel}';
        elseif switchwetper == 0
            P{channel}(:,1)=B{channel}'+((2.*-eta{channel}(1,:)')+(2*Z{channel}(:,2)));
        end    
            
        if any(A{channel}(:,1)<0)
            day = t(j)/(3600*24)
            Fdl = findall(0,'type','figure','tag','TMWWaitbar');
            delete(Fdl)
            error('something goes wrong! check bed elevation and water level matrix')
            return;
        end
    end
    

 
    if calc_sedtan ==1    
    %% sediment transport

        sedtan_5

     %% nodal point

        simpang_8
           

        %% bed update
        if calc_bedte == 1
            %for normal morfac
            if calc_TiMor ==0  
                if t(j)>=morstt
                    for channel =1:Nc
                        courant{channel} = real(abs(Qtot{channel}./A{channel})).*dt./(dx);   %courant criteria
                        if any(courant{channel}(:,2)>=1)
                            courant_info = max(courant{channel}(:,2))

                            ['warning! your courant number is higher than its threshold. Suggestion: please reduce your dt, increase your dx, or reduce your morfac (if >1)']
                            Fdl = findall(0,'type','figure','tag','TMWWaitbar');
                            delete(Fdl);
                            return;
                        end


                    end
                end
                
                Bedte_1

            %tidal aaverage morfac
            elseif calc_TiMor == 1
                TiMor_1
            end
        end
    end




    
    %store the last time step if necessary
    for channel = 1:Nc
        if(t(j+1)==t(end)) && any(t(j+1)==t_store)
                 Z_store{channel}(:,find(t_store==t(j+1))) = Z{channel}(:,2);
                 Q_store{channel}(:,find(t_store==t(j+1))) = Q{channel}(:,2);
                 A_store{channel}(:,find(t_store==t(j+1))) = B{channel}'.*(-eta{channel}(2,:)'+Z_store{channel}(:,find(t_store==t(j+1))));
                 if calc_sedtan ==1
%                     tau_store{channel}(:,find(t_store==t(j+1))) = tau{channel};
                    tau_prime_store{channel}(:,find(t_store==t(j+1))) = tau_prime{channel};
                    if Nbif>0
                        tau_store1(:,find(t_store==t(j+1))) = tau_bolla1;
                        tau_store2(:,find(t_store==t(j+1))) = tau_bolla2;
                        tau_store3(:,find(t_store==t(j+1))) = tau_bolla3;
                    end
                    
                    if Nconf>0
                        tau_store1_conf(:,find(t_store==t(j+1))) = tau_bolla1_conf;
                        tau_store2_conf(:,find(t_store==t(j+1))) = tau_bolla2_conf;
                        tau_store3_conf(:,find(t_store==t(j+1))) = tau_bolla3_conf;
                    end
                    
                    ustar_store{channel}(:,find(t_store==t(j+1))) = ustar{channel};
                    Qs_store{channel}(:,find(t_store==t(j+1))) = Qs{channel};
                    Qb_store{channel}(:,find(t_store==t(j+1))) = Qb{channel};
                    Qtot_store{channel}(:,find(t_store==t(j+1))) = Qtot{channel}(:,2);
                    if calc_bedte ==1
                        eta_store{channel}(:,find(t_store==t(j+1))) = eta{channel}(2,:)';
                    end
                end
        end
    end
    
    
    %harmonic analysis at the end of morphodynamic simulation
    if harm_ini_final == 1
        for channel = 1:Nc
            if(t(j+1)==t(end)) && any(t(j+1)==t_harm)
                     Z_harm{channel}(:,find(t_harm==t(j+1))) = Z{channel}(:,2);
                     Q_harm{channel}(:,find(t_harm==t(j+1))) = Q{channel}(:,2);
                     A_harm{channel}(:,find(t_harm==t(j+1))) = B{channel}'.*(-eta{channel}(2,:)'+Z_harm{channel}(:,find(t_harm==t(j+1))));
                     tau_harm{channel}(:,find(t_harm==t(j+1))) = tau_prime{channel};
                     
                     if calc_sedtan ==1
                        Qtot_harm{channel}(:,find(t_harm==t(j+1))) = Qtot{channel}(:,2);
                        Qs_harm{channel}(:,find(t_harm==t(j+1))) = Qs{channel};
                        Qb_harm{channel}(:,find(t_harm==t(j+1))) = Qb{channel};
                        za_harm{channel}(:,find(t_harm==t(j+1))) = za{channel};
                        T_harm{channel}(:,find(t_harm==t(j+1))) = T{channel};
                        Ca_harm{channel}(:,find(t_harm==t(j+1))) = Ca{channel};
                        F_harm{channel}(:,find(t_harm==t(j+1))) = F{channel};
                        if calc_bedte ==1
                            eta_harm{channel}(:,find(t_harm==t(j+1))) = eta{channel}(2,:)';
                        end
                    end
            end
        end
    end
end

close(h);

%% calculate velocity
for channel =1:Nc
U{channel}=Q_store{channel}./A_store{channel};         % Flow velocity in m/s

if harm_ini_final == 1
    U_harm{channel}=Q_harm{channel}./A_harm{channel};         % Flow velocity in m/s
end
end
%% HD analysis
if harm_anal == 0
    ['no harmonic analysis']
elseif harm_anal == 1
    global wn

    wn(1)=2*pi/T_tide(1);
    wn(2)=2*wn(1);
    wn(3)=3*wn(1);
    Nsteps=floor(T_tide/store);
    
    if harm_ini_final == 0

        for channel = 1:Nc
            for px=1:Nx{channel}
                coefin{channel}=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
                coefout{channel}=nlinfit(t_store(end-Nsteps:end),Z_store{channel}(px,end-Nsteps:end),@harmfit,coefin{channel});
                Z0{channel}(px)=coefout{channel}(1);
                Z_tide{channel}(px)=sqrt(coefout{channel}(2)^2 + coefout{channel}(5)^2);
                Z_overtide2{channel}(px)=sqrt(coefout{channel}(3)^2 + coefout{channel}(6)^2);
                Z_overtide3{channel}(px)=sqrt(coefout{channel}(4)^2 + coefout{channel}(7)^2);
                phase_tide{channel}(px)=atan2(coefout{channel}(5),coefout{channel}(2));
                phase_overtide2{channel}(px)=atan2(coefout{channel}(6),coefout{channel}(3));
                phase_overtide3{channel}(px)=atan2(coefout{channel}(7),coefout{channel}(4));

            end

            for px=1:Nx{channel}
                coefinU{channel}=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
                coefoutU{channel}=nlinfit(t_store(end-Nsteps:end),U{channel}(px,end-Nsteps:end),@harmfit,coefinU{channel});
                U0{channel}(px)=coefoutU{channel}(1);
                U_tide{channel}(px)=sqrt(coefoutU{channel}(2)^2 + coefoutU{channel}(5)^2);
                U_overtide2{channel}(px)=sqrt(coefoutU{channel}(3)^2 + coefoutU{channel}(6)^2);
                U_overtide3{channel}(px)=sqrt(coefoutU{channel}(4)^2 + coefoutU{channel}(7)^2);
                phaseU_tide{channel}(px)=atan2(coefoutU{channel}(5),coefoutU{channel}(2));
                phaseU_overtide2{channel}(px)=atan2(coefoutU{channel}(6),coefoutU{channel}(3));
                phaseU_overtide3{channel}(px)=atan2(coefoutU{channel}(7),coefoutU{channel}(4));

            end
        end
    
    elseif harm_ini_final == 1
        
        for channel = 1:Nc
            for px=1:Nx{channel}
                for i = 1:length(harm_end)
                    coefin{channel}=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
                    coefout{channel}=nlinfit(t_harm,Z_harm{channel}(px,:),@harmfit,coefin{channel});
                    Z0{channel}(px,i)=coefout{channel}(1);
                    Z_tide{channel}(px,i)=sqrt(coefout{channel}(2)^2 + coefout{channel}(5)^2);
                    Z_overtide2{channel}(px,i)=sqrt(coefout{channel}(3)^2 + coefout{channel}(6)^2);
                    Z_overtide3{channel}(px,i)=sqrt(coefout{channel}(4)^2 + coefout{channel}(7)^2);
                    phase_tide{channel}(px,i)=atan2(coefout{channel}(5),coefout{channel}(2));
                    phase_overtide2{channel}(px,i)=atan2(coefout{channel}(6),coefout{channel}(3));
                    phase_overtide3{channel}(px,i)=atan2(coefout{channel}(7),coefout{channel}(4));
                end

            end

            for px=1:Nx{channel}
                for i = 1:length(harm_end)
                    coefinU{channel}=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
                    coefoutU{channel}=nlinfit(t_harm,U_harm{channel}(px,:),@harmfit,coefinU{channel});
                    U0{channel}(px,i)=coefoutU{channel}(1);
                    U_tide{channel}(px,i)=sqrt(coefoutU{channel}(2)^2 + coefoutU{channel}(5)^2);
                    U_overtide2{channel}(px,i)=sqrt(coefoutU{channel}(3)^2 + coefoutU{channel}(6)^2);
                    U_overtide3{channel}(px,i)=sqrt(coefoutU{channel}(4)^2 + coefoutU{channel}(7)^2);
                    phaseU_tide{channel}(px,i)=atan2(coefoutU{channel}(5),coefoutU{channel}(2));
                    phaseU_overtide2{channel}(px,i)=atan2(coefoutU{channel}(6),coefoutU{channel}(3));
                    phaseU_overtide3{channel}(px,i)=atan2(coefoutU{channel}(7),coefoutU{channel}(4));
                end
            end
        end        
    end
end

% 

