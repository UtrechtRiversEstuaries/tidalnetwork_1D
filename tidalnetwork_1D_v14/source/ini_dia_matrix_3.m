






%% variables matrices
for channel = 1:Nc
    x{channel}=0:dx:Sizes(channel,2); %spatial nodes
    Nx{channel} = length(x{channel});   %number of node including ghost nodes
    
 if switchstart == 0

    eta{channel}(1,1:Nx{channel}) = -H0(channel)+ (fliplr((x{channel}.*slope(channel))+eta_end(channel,1)));            %bed elevation
    eta{channel}(2,1:Nx{channel}) = -H0(channel)+ (fliplr((x{channel}.*slope(channel))+eta_end(channel,1)));            %bed elevation
    

                                                         
    if calc_Lb ==1
%         B{channel}(1:Nx{channel})=fliplr(Sizes(channel,1)*exp(-x{channel}/Sizes(channel,3))); %width
        if channel == 1
            B{channel}(1:Nx{channel})=fliplr(Sizes(channel,1)*exp(-x{channel}/Sizes(channel,3))); %width
            B{channel}(1:find(x{channel}==max(x{channel})-Lc_up-(dx)))=B{channel}(find(x{channel}==max(x{channel})-Lc_up-dx)); %width
        else
            B{channel}(1:Nx{channel})=(Sizes(channel,1)*exp(x{channel}/Sizes(channel,3))); %width
        end
    
     elseif calc_Lb ==0
          B{channel}(1:Nx{channel})=fliplr(Sizes(channel,1));

     end

    Z{channel}(:,1) =(fliplr((x{channel}.*slope(channel))+eta_end(channel,1)))';                    %WL
    Z{channel}(:,2) =(fliplr((x{channel}.*slope(channel))+eta_end(channel,1)))';                    %WL
    Q{channel}=zeros(Nx{channel},2);                                            %discharge        
    A{channel}=B{channel}'.*(-eta{channel}(1,:)'+Z{channel}(:,1));             % cross section area 
   
    % Wetted perimeter options
        if switchwetper == 1
            P{channel}=B{channel}';     
        elseif switchwetper == 0
            P{channel}=B{channel}'+(2.*(-eta{channel}(1,:)'+(2*Z{channel}(:,1)))); 
        end
    
  elseif switchstart == 1
        if channel == 2
            eta{channel}(1,1:Nx{channel}) = perturbdepth+coldmorf.eta_store{channel}(:,end)';            %bed elevation
            eta{channel}(2,1:Nx{channel}) = perturbdepth+coldmorf.eta_store{channel}(:,end)';            %bed elevation
        elseif channel ==3
            eta{channel}(1,1:Nx{channel}) = -perturbdepth+coldmorf.eta_store{channel}(:,end)';            %bed elevation
            eta{channel}(2,1:Nx{channel}) = -perturbdepth+coldmorf.eta_store{channel}(:,end)';            %bed elevation
        else
        eta{channel}(1,1:Nx{channel}) = coldmorf.eta_store{channel}(:,end)';                              %bed elevation
        eta{channel}(2,1:Nx{channel}) = coldmorf.eta_store{channel}(:,end)';                              %bed elevation
        end
    

                                                         
    if calc_Lb ==1
%         B{channel}(1:Nx{channel})=fliplr(Sizes(channel,1)*exp(-x{channel}/Sizes(channel,3))); %width
        if channel == 1
            B{channel}(1:Nx{channel})=fliplr(Sizes(channel,1)*exp(-x{channel}/Sizes(channel,3))); %width
            B{channel}(1:find(x{channel}==max(x{channel})-40000-(dx)))=B{channel}(find(x{channel}==max(x{channel})-40000-dx)); %width
        else
            B{channel}(1:Nx{channel})=(Sizes(channel,1)*exp(x{channel}/Sizes(channel,3))); %width
        end
    
     elseif calc_Lb ==0
          B{channel}(1:Nx{channel})=fliplr(Sizes(channel,1));

     end

    Z{channel}(:,1) =coldHD.Z0{channel}';                    %WL
    Z{channel}(:,2) =coldHD.Z0{channel}';                    %WL
    Q{channel}=zeros(Nx{channel},2);                    %discharge        
    A{channel}=B{channel}'.*(-eta{channel}(1,:)'+Z{channel}(:,1));             % cross section area                    
        
    % Wetted perimeter options
        if switchwetper == 1
            P{channel}=B{channel}';    
        elseif switchwetper == 0
            P{channel}=B{channel}'+(2.*(-eta{channel}(1,:)'+(2*Z{channel}(:,1))));     
        end
        
  end
    
    A12{channel} = zeros(Nx{channel},1);                %A average for 2 points
    %preissmann variable for flow
    C1{channel} = repmat(NaN,Nx{channel},1);
    D1{channel} = repmat(NaN,Nx{channel},1);
    E1{channel} = repmat(NaN,Nx{channel},1);
    CC{channel} = repmat(NaN,Nx{channel},1);
    DD{channel} = repmat(NaN,Nx{channel},1);
    EE{channel} = repmat(NaN,Nx{channel},1);
    FF{channel} = repmat(NaN,Nx{channel},1);
    MM{channel} = repmat(NaN,Nx{channel},1);
    NN{channel} = repmat(NaN,Nx{channel},1);
    MN{channel} = repmat(NaN,Nx{channel},1);
    ML{channel} = repmat(NaN,Nx{channel},1);
    LL{channel} = repmat(NaN,Nx{channel},1);
    NL{channel} = repmat(NaN,Nx{channel},1);
    K{channel} = repmat(NaN,Nx{channel},1);
    M{channel} = repmat(NaN,Nx{channel},1);
    N{channel} = repmat(NaN,Nx{channel},1);
    KK{channel} = repmat(NaN,Nx{channel},1);
    MK{channel} = repmat(NaN,Nx{channel},1);
    NK{channel} = repmat(NaN,Nx{channel},1);
    
end

%matrices to store data

t_store = zeros(1,floor(t(end)/store)+1);

for j=1:length(t_store)
    
    j_t = (j*store)-store;
        t_store(j) = j_t;

end

if harm_ini_final == 1
    Nsteps_harm=floor((2*T_tide)/harm_step);
    
    t_harm = harm_end(1)-(Nsteps_harm*harm_step):harm_step:harm_end(1);
    
  for channel =1:Nc  
    Z_harm{channel} = repmat(NaN,Nx{channel},length(t_harm));
    Q_harm{channel} = repmat(NaN,Nx{channel},length(t_harm));
    A_harm{channel} = repmat(NaN,Nx{channel},length(t_harm));
    tau_harm{channel} = repmat(NaN,Nx{channel},length(t_harm));
  end
  
end
    


for channel =1:Nc
    Z_store{channel} = repmat(NaN,Nx{channel},length(t_store));
    Q_store{channel} = repmat(NaN,Nx{channel},length(t_store));
    A_store{channel} = repmat(NaN,Nx{channel},length(t_store));
end


%% Boundary conditions
%downstream boundary for each channel
seabound = find(isnan(Topo(:,3))==1);    %channel conncected to the sea
bifurDowbound = find(~isnan(Topo(:,4))); %channel connected to other channel
confDowbound = find(~isnan(Topo(:,3))&isnan(Topo(:,4))); %channel connected to other channel
%upatream boundary for each channel
Upbound = find(isnan(Topo(:,1))==1);    %channel connected to the upstream boundary
bifurUpbound = find(~isnan(Topo(:,1))&isnan(Topo(:,2))==1); %channel connected to the other channel
confUpbound = find(~isnan(Topo(:,2))); %channel connected to the other channel



%continuity requirement
for tel=1:Nc
    topou{tel}=[ Topo(tel,find(~isnan(Topo(tel,1:2)))) ];
    topod{tel}=[ Topo(tel,2+find(~isnan(Topo(tel,3:4)))) ];
end

%channels connected at junction
Nbif = length(bifurDowbound);	%number of bifurcations
Nconf = length(confUpbound);	%number of bifurcations

% for tel=1:Nbif
%     topob{tel}=[ bifurDowbound(tel) Topo(bifurDowbound(tel),3:4) ];
% end

% member of of bifurcations
if Nbif>0
    for tel=1:Nbif
        topob{tel}=[ bifurDowbound(tel) Topo(bifurDowbound(tel),3:4) ];
        tobif(tel,:)=[ bifurDowbound(tel) Topo(bifurDowbound(tel),3:4) ];
    end
else
    topob=[];
end

% member of of confluences
if Nconf>0
    for tel=1:Nconf
        topoc{tel}=[ Topo(confUpbound(tel),1:2) confUpbound(tel) ];
        tocon(tel,:)=[ Topo(confUpbound(tel),1:2) confUpbound(tel) ];
    end
else
    topoc=[];
end

% initial condition at junction

if Nbif>0
    for bif = 1:length(bifurDowbound)
%                     bif = find(any(topob==channel));
                    Q{topob{bif}(1)}(Nx{topob{bif}(1)},1) = Q{topob{bif}(2)}(1,1)+Q{topob{bif}(3)}(1,1); 

    end
end

if Nconf>0
    for conf = 1:length(confUpbound)
%                     conf = find(any(topoc{:}==channel));
                    Q{topoc{conf}(3)}(1,1) = Q{topoc{conf}(1)}(Nx{topoc{conf}(1)},1)+Q{topoc{conf}(2)}(Nx{topoc{conf}(2)},1); 

    end
end
%matrices to calculate boundary conditions
M_bound = zeros(Nc*Nb,Nc*Nb);
F_bound =zeros(Nc*Nb,1);
ZQ_bound = zeros(Nc*Nb,1);

% M_bound(1,1) = 1;
chan_bound = zeros(Nb,Nb*Nc);

% for channel=1:Nc
%     if any (channel == seabound)
%         M_bound(channel*Nb,channel*Nb)=1;
%     end
% end
%% 
for i = 1:length(seabound)
Z_end(i)=Z{seabound(i)}(Nx{seabound(i)},1);
end
%% equilibrium discharge calculation based on Bolla 2015
if calc_Q == 1
    for i = 1:length(Upbound)
    discharge = B{Upbound(i)}(1)*(1/sqrt(Cd))*H0(Upbound(i))*sqrt(g*H0(Upbound(i))*slope(Upbound(i)));
    end
end

%% sediment transport 
% matrices
%sediment transport matrices
if calc_sedtan == 1
    for channel = 1:Nc
        qf{channel} = repmat(NaN,Nx{channel},1);            %flow velocity*width
        miu{channel} = repmat(NaN,Nx{channel},1);           %for effective shear stress in vRijn
        Cf{channel} = repmat(NaN,Nx{channel},1);            %friction term in sedtan
        tau{channel} = repmat(NaN,Nx{channel},1);           %shieds number for sediment transport
        Dstar(channel) = Sizes(channel,4)*(Rr*g/nu^2).^(1/3); %Bonnefille dimensionless grain size
%         if Sizes(channel,4)<=0.0001
%             ws(channel) = Rr*g*((Sizes(channel,4))^2)/(nu*18);
%         elseif Sizes(channel,4)>0.0001 && Sizes(channel,4)<=0.001
%             ws(channel) = (10*nu/Sizes(channel,4))*((sqrt(1+(0.01*Rr*g*((Sizes(channel,4))^3)/(nu^2))))-1);
%         elseif Sizes(channel,4)> 0.001
%             ws(channel) = 1.1*sqrt(Rr*g*Sizes(channel,4));
%         end
        ws(channel) = (nu/Sizes(channel,4))*(-10.36+sqrt((10.36^2)+(1.01*(Dstar(channel)^3)))); %soulsby(1997) calculation for settling velocity
        T{channel} = repmat(NaN,Nx{channel},1);                 %effective Shields
        tau_prime{channel} = repmat(NaN,Nx{channel},1);         %grain related shields
        Ca{channel} = repmat(NaN,Nx{channel},1);                %reference concentration for vRijn
        ustar{channel} = repmat(NaN,Nx{channel},1);             %nondimesionalised flow
        beta{channel} = repmat(NaN,Nx{channel},1);              %coef in vRijn
        phi{channel} = repmat(NaN,Nx{channel},1);               %coef in vRijn
        Zc{channel} = repmat(NaN,Nx{channel},1);                %coef in vRijn
        F{channel} = repmat(NaN,Nx{channel},1);                 %suspended sediment profile in vRijn
        qb_star{channel} = repmat(NaN,Nx{channel},1);           %nondimensional bedload
        qs_star{channel} = repmat(NaN,Nx{channel},1);           %nondimensional suspension
        Qb{channel} = repmat(NaN,Nx{channel},1);                %flux bedload
        Qs{channel} = repmat(NaN,Nx{channel},1);                %flux suspended load
        qs{channel} = repmat(NaN,Nx{channel},1);                %width-averaged suspended load
        qb{channel} = repmat(NaN,Nx{channel},1);                %width-averaged bedload
        Qtot{channel} = zeros(Nx{channel},2);                   %total sediment load
        if calc_TiMor == 1
           Q_storetide{channel} = zeros(Nx{channel},jTide); 
           Qtot_tideAve{channel}= zeros(Nx{channel},1);
        end
        za{channel} = zeros(Nx{channel},1);                     %reference level for vRijn

        %store data purpose 
        Qs_store{channel} = zeros(Nx{channel},length(t_store));
        Qb_store{channel} = zeros(Nx{channel},length(t_store));
        Qtot_store{channel} = zeros(Nx{channel},length(t_store));
        eta_store{channel} = zeros(Nx{channel},length(t_store));
        tau_store{channel} = zeros(Nx{channel},length(t_store));
        tau_prime_store{channel} = zeros(Nx{channel},length(t_store));
        ustar_store{channel} = zeros(Nx{channel},length(t_store));
    end
    %ndlpt relation for bifurcation
    if Nbif >0
        for bifur = 1:Nbif
          if switchndpt == 0 || switchndpt == 1  
            Qy = zeros(bifur,1);                            %transverse flow at the junction
            H123 = zeros(bifur,1);                          %depth at junction point
            dzdy = zeros(bifur,1);                          %transverse slope at junction
            qsy = zeros(bifur,1);                           %width averaged transverse sediment exchange at nodal point
            Qsy = zeros(bifur,1);                           %transverse sediment exchange at nodal point
          elseif switchndpt == 2
            rat_Qtot(bifur,1)= zeros(bifur,1); 
          end
         %shields number for nodalpoint calculation
          tau_bolla1 = zeros(bifur,1);
          tau_bolla2 = zeros(bifur,1);
          tau_bolla3 = zeros(bifur,1);


        end
    end

    %ndlpt relation Bolla for confluence
    if Nconf>0
        for conf = 1:Nconf
          if switchndpt == 0 || switchndpt == 1 
            Qy_conf = zeros(conf,1);
            H123_conf = zeros(conf,1);
            dzdy_conf = zeros(conf,1);
            qsy_conf = zeros(conf,1);
            Qsy_conf = zeros(conf,1);
          elseif switchndpt == 2
            rat_Qtot(bifur,1)= zeros(conf,1); 
          end 
         %shields number for nodalpoint calculation
         tau_bolla1_conf = zeros(conf,1);
         tau_bolla2_conf = zeros(conf,1);
         tau_bolla3_conf = zeros(conf,1);


        end
    end
    ini_dia_sedtan_simpang_1
end