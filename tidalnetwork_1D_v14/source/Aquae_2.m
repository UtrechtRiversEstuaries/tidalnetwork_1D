for channel=1:Nc
  %% in-channel HD

        A{channel}=B{channel}'.*(-eta{channel}(1,:)'+Z{channel}(:,1));
        P{channel}=B{channel}';%+((2.*-eta{channel}(1,:)')+(2*Z{channel}(:,1)));
        

        
        for i = 1:Nx{channel}-1

        %channel configurations for preissmann scheme
            Btheta{channel}(i) =  ((B{channel}(i+1)+B{channel}(i))/2)  ;  %assumed constant width through time at each location
            A12{channel}(i)   =  (A{channel}(i+1)+A{channel}(i))/2;

        %known variables for continuity equation to calculate next spatial
        %and temporal grid
            C1{channel}(i) = Btheta{channel}(i)./(2*dt);
            D1{channel}(i) = theta/dx;
            E1{channel}(i) = ((Btheta{channel}(i)/(2*dt))*Z{channel}(i+1,1)) + ((Btheta{channel}(i)/(2*dt))*Z{channel}(i,1))...
                                 - (((1-theta)/dx)*(Q{channel}(i+1,1)-Q{channel}(i,1)));

        %known variables for momentum equation to calculate next spatial
        %and temporal grid
            CC{channel}(i) = g*A12{channel}(i)*theta/dx;
            DD{channel}(i) = (1/(2*dt)) + (adv*theta*Q{channel}(i+1,1)/(dx*A{channel}(i+1)))...
                + ((Cd*abs(Q{channel}(i+1))*P{channel}(i+1))/(((A{channel}(i+1))^2)*2));
            EE{channel}(i) = (1/(2*dt)) - (adv*theta*Q{channel}(i,1)/(dx*A{channel}(i)))...
                + ((Cd*abs(Q{channel}(i,1))*P{channel}(i))/(((A{channel}(i))^2)*2));
            FF{channel}(i) = (Q{channel}(i+1,1)/(2*dt)) + (Q{channel}(i,1)/(2*dt))...
                - ((g*A12{channel}(i)/dx)*(1-theta)*(Z{channel}(i+1,1)-Z{channel}(i,1)))...
                -(adv*(1-theta)*(Q{channel}(i+1,1)/(dx*A{channel}(i+1)))*(Q{channel}(i+1,1)))...
                +(adv*(1-theta)*(Q{channel}(i,1)/(dx*A{channel}(i)))*(Q{channel}(i,1)));
            
         %known variables for double sweep thomas algorithm (Liggett &
         %Cunge, 1975)
            MM{channel}(i) = (1+(DD{channel}(i)/EE{channel}(i)))...
                /((CC{channel}(i)/EE{channel}(i))-(C1{channel}(i)/D1{channel}(i)));
            NN{channel}(i) = -((FF{channel}(i)/EE{channel}(i))+(E1{channel}(i)/D1{channel}(i)))...
                /((CC{channel}(i)/EE{channel}(i))-(C1{channel}(i)/D1{channel}(i)));
            MN{channel}(i) = ((CC{channel}(i)/EE{channel}(i)) + (C1{channel}(i)/D1{channel}(i)))...
                /((CC{channel}(i)/EE{channel}(i)) - (C1{channel}(i)/D1{channel}(i)));
            
            
            ML{channel}(i) = ((FF{channel}(i)/DD{channel}(i))-(E1{channel}(i)/D1{channel}(i)))...
                /((CC{channel}(i)/DD{channel}(i))-(C1{channel}(i)/D1{channel}(i)));  
            
            LL{channel}(i) = (-1-(EE{channel}(i)/DD{channel}(i)))...
                /((CC{channel}(i)/DD{channel}(i))-(C1{channel}(i)/D1{channel}(i)));
       
            NL{channel}(i) = ((CC{channel}(i)/DD{channel}(i))+(C1{channel}(i)/D1{channel}(i)))...
                /((CC{channel}(i)/DD{channel}(i))-(C1{channel}(i)/D1{channel}(i)));
            
          %calculate known variable for linear assumption (Q=KZi+Mi+NiZ1)
            if i==1
                K{channel}(i+1) = ((-C1{channel}(i)/D1{channel}(i)) - (CC{channel}(i)/EE{channel}(i)))...
                    / (1+(DD{channel}(i)/EE{channel}(i)));
                M{channel}(i+1) = ((FF{channel}(i)/EE{channel}(i)) + (E1{channel}(i)/D1{channel}(i)))...
                    / (1+(DD{channel}(i)/EE{channel}(i)));
                N{channel}(i+1) = ((CC{channel}(i)/EE{channel}(i)) - (C1{channel}(i)/D1{channel}(i)))...
                    / (1+(DD{channel}(i)/EE{channel}(i)));
            else
                K{channel}(i+1) = - ( C1{channel}(i) + (C1{channel}(i)*MN{channel}(i)) - D1{channel}(i)*K{channel}(i)*MN{channel}(i))...
                    /((C1{channel}(i)*MM{channel}(i))+D1{channel}(i)-(D1{channel}(i)*K{channel}(i)*MM{channel}(i)));
                M{channel}(i+1) = ((-C1{channel}(i)*NN{channel}(i))+(D1{channel}(i)*K{channel}(i)*NN{channel}(i))+(D1{channel}(i)*M{channel}(i))+E1{channel}(i))...
                    /((C1{channel}(i)*MM{channel}(i))+D1{channel}(i)-(D1{channel}(i)*K{channel}(i)*MM{channel}(i)));
                N{channel}(i+1) = D1{channel}(i)*N{channel}(i)...
                    /((C1{channel}(i)*MM{channel}(i))+D1{channel}(i)-(D1{channel}(i)*K{channel}(i)*MM{channel}(i)));
            end
            
        end  
          for i = Nx{channel}-1:-1:1
              
                if i==Nx{channel}-1
                    KK{channel}(i) = (((CC{channel}(i)/DD{channel}(i))+(C1{channel}(i)/D1{channel}(i))) ...
                                         / (1+(EE{channel}(i)/DD{channel}(i))));
                    MK{channel}(i) = (((FF{channel}(i)/DD{channel}(i))-(E1{channel}(i)/D1{channel}(i))) ...
                                         / (1+(EE{channel}(i)/DD{channel}(i))));
                    NK{channel}(i) = (((-CC{channel}(i)/DD{channel}(i))+(C1{channel}(i)/D1{channel}(i))) ...
                                         / (1+(EE{channel}(i)/DD{channel}(i))));

                else

                   KK{channel}(i) = ((C1{channel}(i)*NL{channel}(i))+C1{channel}(i)+(D1{channel}(i)*KK{channel}(i+1)*NL{channel}(i)))...
                       /(D1{channel}(i)-(C1{channel}(i)*LL{channel}(i))-(D1{channel}(i)*KK{channel}(i+1)*LL{channel}(i)));
                   MK{channel}(i) = ((C1{channel}(i)*ML{channel}(i))+(D1{channel}(i)*KK{channel}(i+1)*ML{channel}(i))+(D1{channel}(i)*MK{channel}(i+1))-E1{channel}(i))...
                       /(D1{channel}(i)-(C1{channel}(i)*LL{channel}(i))-(D1{channel}(i)*KK{channel}(i+1)*LL{channel}(i)));
                   NK{channel}(i) = D1{channel}(i)*NK{channel}(i+1)...
                       /(D1{channel}(i)-(C1{channel}(i)*LL{channel}(i))-(D1{channel}(i)*KK{channel}(i+1)*LL{channel}(i)));

                end
            
         end     
    end
    
%% calculate Q and Z at boundaries
      %calculate matrix of independent variables
      for channel = 1:Nc
         if any (channel  == Upbound)
            F_bound(((Nb*(channel-1))+1):(Nb*channel)) = [Q{channel}(1,2);-MK{channel}(1);-M{channel}(Nx{channel});0];  
         elseif any (channel == seabound)
            F_bound(((Nb*(channel-1))+1):(Nb*channel)) = [0;-MK{channel}(1);-M{channel}(Nx{channel});Z{channel}(Nx{channel},2)];
         else
            F_bound(((Nb*(channel-1))+1):(Nb*channel)) = [0;-MK{channel}(1);-M{channel}(Nx{channel});0];
         end
      end

%% junction boundary(es) HD
      %matrix of known variables
%       for channel = 1:Nc
%         M_bound((channel*Nb)-2,(channel*Nb)-3:channel*Nb) = [-1,KK{channel}(1),0,NK{channel}(1)];
%         M_bound((channel*Nb)-1,(channel*Nb)-2:channel*Nb) = [N{channel}(Nx{channel}),-1,K{channel}(Nx{channel})];
%       end
      
      %%WARNING modified matrix structure first if you want to add more
      %%channel
%       for channel = 1:Nc
%           if any (channel == bifurDowbound)
%               M_bound(channel*Nb,channel*Nb:channel*Nb+2) = [1,0,-1];
%               M_bound((channel*Nb)+1,(channel*Nb)-1:((channel*Nb)+Nb)+1) = [1,0,-1,0,0,0,-1];
%           elseif any (channel == seabound)
%               M_bound(channel*Nb,channel*Nb) = 1;
%           
%            %modification is necessary for multiple bifurcation case
%             if channel == seabound(end)
%               M_bound(channel*Nb-3,(channel*Nb)-2) = -1;
%               M_bound(channel*Nb-3,bifurDowbound*Nb) = 1;
%             end
%           end
%       end
        
    for channel = 1:Nc
        if any (channel == Upbound)
            chan_bound(1,(Nb*(channel-1))+1) = 1;
        elseif any (channel == bifurUpbound)
            [temp_x,temp_y] = find(tobif==channel);
            cek_y = find(temp_y==1);
            temp_y(cek_y)=[];
            temp_x(cek_y)=[];
            chan_1 = tobif(temp_x,1);
            
            if temp_y == 2
                chan_3 = tobif(temp_x,3);
                % to satisfy Q1 = Q2+Q3
                chan_bound(1,(Nb*(chan_1-1))+3)=1;          %channel 1 in this connection
                chan_bound(1,(Nb*(chan_3-1))+1)=-1;         %channel 3 in this connection
                chan_bound(1,(Nb*(channel-1))+1)=-1;        %channel 2 in this connection
            elseif temp_y == 3
                chan_1 = tobif(temp_x,1);
                %to satisfy Z1 = Z3
                chan_bound(1,(Nb*(chan_1-1))+4)=1;
                chan_bound(1,(Nb*(channel-1))+2)=-1;
            end
        elseif any(channel == confUpbound)
            [temp_x,temp_y] = find(tocon==channel);
            cek_y = find(temp_y~=3);
            temp_y(cek_y)=[];
            temp_x(cek_y)=[];
            chan_1 = tocon(temp_x,1);
            chan_2 = tocon(temp_x,2);
            % to satisfy Q1 = Q2+Q3
            chan_bound(1,(Nb*(chan_1-1))+3)=1;          %channel 1 in this connection
            chan_bound(1,(Nb*(chan_2-1))+3)=1;          %channel 2 in this connection
            chan_bound(1,(Nb*(channel-1))+1)=-1;        %channel 3 in this connection
            
        end
        
        chan_bound(2,((Nb*(channel-1))+1):(Nb*channel))=[-1,KK{channel}(1),0,NK{channel}(1)];
        chan_bound(3,((Nb*(channel-1))+1):(Nb*channel))=[0,N{channel}(Nx{channel}),-1,K{channel}(Nx{channel})];
        chan_bound(4,((Nb*(channel-1))+1):(Nb*channel))=[0,0,0,1];
        
        if any(channel == bifurDowbound)
            [temp_x,temp_y] = find(tobif==channel);
            cek_y = find(temp_y~=1);
            temp_y(cek_y)=[];
            temp_x(cek_y)=[];
            chan_2 = tobif(temp_x,2);
            %satistify Z1=Z2
            chan_bound(4,(Nb*(chan_2-1))+2) = -1;
        elseif any(channel == confDowbound)
            [temp_x,temp_y] = find(tocon==channel);
            cek_y = find(temp_y==3);
            temp_y(cek_y)=[];
            temp_x(cek_y)=[];
            chan_3 = tocon(temp_x,3);
            %satisfy Z1&2 = Z3
            chan_bound(4,(Nb*(chan_3-1))+2) = -1;
        end
       M_bound(((Nb*(channel-1))+1):(Nb*channel),:) = chan_bound;
       chan_bound = zeros(Nb,Nb*Nc);
    end
      ZQ_bound = M_bound\F_bound;
      
      for channel =  1:Nc
            Q{channel}(1,2) = ZQ_bound((channel*Nb)-3);
            Z{channel}(1,2) = ZQ_bound((channel*Nb)-2);
            Q{channel}(Nx{channel},2) = ZQ_bound((channel*Nb)-1);
            Z{channel}(Nx{channel},2) = ZQ_bound((channel*Nb));
      end
      
      for channel = Nc:-1:1
          for i = Nx{channel}-1:-1:2

              Z{channel}(i,2) = (MM{channel}(i)*Q{channel}(i+1,2))+(MN{channel}(i)*Z{channel}(i+1,2))+NN{channel}(i);
              Q{channel}(i,2) = (K{channel}(i)*Z{channel}(i,2))+M{channel}(i)+(N{channel}(i)*Z{channel}(1,2));
              
          end
      end