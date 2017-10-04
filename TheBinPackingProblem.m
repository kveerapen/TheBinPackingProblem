% Load test data %

%C = 100;           The test data should be plugged into the Command Window. 
%items = 100;
%lambda = 5;
%mu = 500;
%sig = 1;
 

disp('This is the 1-D Bin Packing Heuristic Algorithm');
disp('1 - Uniform Distribution');
disp('2 - Normal Distribution');
disp('3 - Poisson Distribution');
chos=input('Please choose which distribution you would like to test: ');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNIFORM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DISTRIBUTION
    
if chos==1,
    W = randi([0,C],1,items)                                                 %Uniform Random Number Generator
    n = length(W);                                                           % n is the number of items in the list
    Bins = C*ones(1,n);                                                      % Bins is represented as a 1 by n matrix with all ones, since no item has packed it yet.
    i = 1:C;                                                                 % i is the index from 1 to the capacity of the bin
    bin_number = 1;                                                          % we start with the first bin
    A = zeros(n,n);                                                          % A is the answer matrix that initially starts off with all zeros
    a = n;                                                              
    flags = zeros(1,n);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('1 - Ordered List(Decreasing)');
    disp('2 - Non-Ordered List');
    chos1=input('Please choose which arrangement you would like to test: ');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chos1==1, 
        
        Wn = sort(W,'descend')                                                % Orders the list of items, W, into descending order 
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('0 - Both Algorithms');
        disp('1 - First Fit Method');
        disp('2 - Best Fit Method');
        chos2=input('Now please choose which method you want to Bin Pack: ');
        if isempty(chos2)
            chos2=0;
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if chos2==0 || chos2==1,
            tic                                                               % This function calculates the time taken from the start 'tic' to the end of the algorithm 'toc'
            aa = A; Bins1 = Bins;                                             % The variables are assigned different variables so that it doesnt clash with the other algorithms
           
            
            if Wn(1)> Bins1(bin_number)                                       % Because the first item would be the highest value in the list of objects since its in dec order, 
                                                                              % if the largest item in the list doesnt fit into the bin, the packing cannot be carried out.
                display('The packing is impossible')
            else
                k=1;
   
                while k<=n                                                    % The While loop gets executed by using k iterations from 1 to n
        
                    if  Wn(k)<=Bins1(bin_number)                                                                 % First of all we start by selecting the first term in the list of objects, to see whether it will fit into the first bin
                        aa(bin_number,k)=k;                                   % If the item fits, we print out the first answer matrix A with the first term in the first row which is bin 1, and so on.
                        Bins1(bin_number)=Bins1(bin_number) - Wn(k);          % We now determine the residual capacity of the bin after the item has been packed into it, so we know the remaining space left in the bins.
                        k=k+1;                                                % We then move onto the next item in the list and carry out the same procedure again until all the items are packed
                    else
                        bin_number=bin_number+1;                              % If the item doesnt fit, it will open up a new bin and pack it into that
                    end;
        %display(Wn(k:end));                                                  % This will display the remaining weights in the list everytime an object has been packed.
        
                end
            end
            
   Answer1 = aa;                                                               % Final Answer Matrix
   display(Answer1)
   formatSpec = 'Uniform Distribution/First Fit Decreasing/Number of Bins = %d.';
   Conclusion = sprintf(formatSpec,bin_number);                                % Conclusion, prints out the number of bins required to pack all items
   display(Conclusion)
   
   EmptySpace = (1 - (sum(W)/(bin_number*C)))*100                              % This is the empty space percentage calculation described in the thesis
   
           
   toc
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if chos2==0 || chos2==2,
            
            tic
            ab = A; Bins2 = Bins; bin_number=0;
             if Wn(1)> Bins2(1)                                                 % Because the first item would be the highest value in the list of objects since its in dec order,  
                 display('The packing is impossible')                           % if it doesnt fit into the bin, the packing cannot be carried out.display('The packing is impossible')
             else
                while a~=0
                    bin_number=bin_number + 1;    
                    for j=1:n                                                   % The While loop gets executed by using k iterations from 1 to n
                        if  flags(j)==0 && Wn(j)<=Bins2(bin_number)             % First of all we start by selecting the first term in the list of objects, to see whether it will fit into the first bin
                            ab(bin_number,j)=j;                                 % If the item fits, we print out the first answer matrix A with the first term in the first row which is bin 1, and so on.
                            a = a - 1;
                            Bins2(bin_number) = Bins2(bin_number) - Wn(j);
                            %Wn(j)= Wn(j) + 10*C;   
                            flags(j)=1;                                         % We now determine the residual capacity of the bin after the item has been packed into it, so we know the remaining space left in the bins.
                        end                                                     % We then move onto the next item in the list and carry out the same procedure again until all the items are packed
                    end                                                         % If the item doesnt fit, it will open up a new bin and pack it into that1
                end
            end
        Answer2 = ab;
        display(Answer2)
        formatSpec = 'Uniform Distribution/Best Fit Decreasing/Number of Bins = %d.';
        Conclusion = sprintf(formatSpec,bin_number);
        display(Conclusion)
        EmptySpace = (1 - (sum(W)/(bin_number*C)))*100
        toc
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chos1==2,
    n = length(W);
    Bins = C*ones(1,n);
    i = 1:C;
    bin_number = 1; 
    A = zeros(n,n);
    a = n;
    flags = zeros(1,n);
    Wn = sort(W,'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('0 - Both Algorithms');
        disp('1 - First Fit Method');
        disp('2 - Best Fit Method');
        chos2=input('Now please choose which method you want to Bin Pack: ');
        if isempty(chos2)
            chos2=0;
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if chos2==0 || chos2==1,
            
            tic
            ac = A; Bins3 = Bins; 
            if Wn(1)> Bins3(bin_number)                      
                                                            
                display('The packing is impossible')
            else
                k=1;
   
                while k<=n                                  
        
                    if  W(k)<=Bins3(bin_number)            
                        ac(bin_number,k)=k;                
                        Bins3(bin_number)=Bins3(bin_number) - W(k); 
                        k=k+1;                            
                    else
                        bin_number=bin_number+1;          
                    end;
            %display(Wn(k:end));                
        
                end
            end
            
        Answer3 = ac;
        display(Answer3)
        formatSpec = 'Uniform Distribution/First Fit/Number of Bins = %d.';
        Conclusion = sprintf(formatSpec,bin_number);
        display(Conclusion) 
        EmptySpace = (1 - (sum(W)/(bin_number*C)))*100
        toc
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if chos2==0 || chos2==2,
            tic
            ad = A; Bins4 = Bins; bin_number=0;
            
             if Wn(1)> Bins4(1)                                                         % Because the first item would be the highest value in the list of objects since its in dec order, 
                  
                 display('The packing is impossible')                                   % if it doesnt fit into the bin, the packing cannot be carried out.display('The packing is impossible')
            
             else
   
                while a~=0
                    bin_number=bin_number + 1;    
                    for j=1:n                                                           % The While loop gets executed by using k iterations from 1 to n
                        if  flags(j)==0 && W(j)<=Bins4(bin_number)                      % First of all we start by selecting the first term in the list of objects, to see whether it will fit into the first bin
                            ad(bin_number,j)=j;                                         % If the item fits, we print out the first answer matrix A with the first term in the first row which is bin 1, and so on.
                            a = a - 1;
                            Bins4(bin_number) = Bins4(bin_number) - W(j);
                            %Wn(j)= Wn(j) + 10*C;   
                            flags(j)=1;                                                 % We now determine the residual capacity of the bin after the item has been packed into it, so we know the remaining space left in the bins.
                        end                                                             % We then move onto the next item in the list and carry out the same procedure again until all the items are packed
                    end                                                                 % If the item doesnt fit, it will open up a new bin and pack it into that1
   
                end
   
            
             end
             
        Answer4 = ad;
        display(Answer4)
        formatSpec = 'Uniform Distribution/Best Fit/Number of Bins = %d.';
        Conclusion = sprintf(formatSpec,bin_number);
        display(Conclusion)
        EmptySpace = (1 - (sum(W)/(bin_number*C)))*100
        toc
        end
    end
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NORMAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DISTRIBUTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if chos==2,
        
    W = normrnd(mu,sig,[1 items])
    n = length(W);
    Bins = C*ones(1,n);
    i = 1:C;
    bin_number = 1; 
    A = zeros(n,n);
    a = n;
    flags = zeros(1,n);
    
    disp('1 - Ordered List(Decreasing)');
    disp('2 - Non-Ordered List');
    chos1=input('Please choose which arrangement you would like to test: ');
   
    
    if chos1==1, 
        
        Wn = sort(W,'descend')
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('0 - Both Algorithms');
        disp('1 - First Fit Method');
        disp('2 - Best Fit Method');
        chos2=input('Now please choose which method you want to Bin Pack: ');
        if isempty(chos2)
            chos2=0;
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if chos2==0 || chos2==1,
            tic
            ae = A; Bins5 = Bins; 
            if Wn(1)> Bins5(bin_number)                                         % Because the first item would be the highest value in the list of objects since its in dec order, 
                                                                                % if it doesnt fit into the bin, the packing cannot be carried out.
                display('The packing is impossible')
            else
                k=1;
   
                while k<=n                                                      % The While loop gets executed by using k iterations from 1 to n
        
                    if  Wn(k)<=Bins5(bin_number)                                % First of all we start by selecting the first term in the list of objects, to see whether it will fit into the first bin
                        ae(bin_number,k)=k;                                     % If the item fits, we print out the first answer matrix A with the first term in the first row which is bin 1, and so on.
                        Bins5(bin_number)=Bins5(bin_number) - Wn(k);            % We now determine the residual capacity of the bin after the item has been packed into it, so we know the remaining space left in the bins.
                        k=k+1;                                                  % We then move onto the next item in the list and carry out the same procedure again until all the items are packed
                    else
                        bin_number=bin_number+1;                                % If the item doesnt fit, it will open up a new bin and pack it into that
                    end;
        %display(Wn(k:end));                                                    % This will display the remaining weights in the list everytime an object has been packed.
        
                end
            end
            
   Answer5 = ae;
   display(Answer5)
   formatSpec = 'Normal Distribution/First Fit Decreasing/Number of Bins = %d.';
   Conclusion = sprintf(formatSpec,bin_number);
   display(Conclusion) 
   EmptySpace = (1 - (sum(W)/(bin_number*C)))*100     
   toc
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if chos2==0 || chos2==2,
            tic
            af = A; Bins6 = Bins; bin_number=0;
            
             if Wn(1)> Bins6(1)                                                         % Because the first item would be the highest value in the list of objects since its in dec order, 
                  
                 display('The packing is impossible')                                   % if it doesnt fit into the bin, the packing cannot be carried out.display('The packing is impossible')
            
             else
   
                while a~=0
                    bin_number=bin_number + 1;    
                    for j=1:n                                                           % The While loop gets executed by using k iterations from 1 to n
                        if  flags(j)==0 && Wn(j)<=Bins6(bin_number)                     % First of all we start by selecting the first term in the list of objects, to see whether it will fit into the first bin
                            af(bin_number,j)=j;                                         % If the item fits, we print out the first answer matrix A with the first term in the first row which is bin 1, and so on.
                            a = a - 1;
                            Bins6(bin_number) = Bins6(bin_number) - Wn(j);
                            %Wn(j)= Wn(j) + 10*C;   
                            flags(j)=1;                                                 % We now determine the residual capacity of the bin after the item has been packed into it, so we know the remaining space left in the bins.
                        end                                                             % We then move onto the next item in the list and carry out the same procedure again until all the items are packed
                    end                                                                 % If the item doesnt fit, it will open up a new bin and pack it into that1
   
                end
   
            
             end
             
        Answer6 = af;
        display(Answer6)
        formatSpec = 'Normal Distribution/Best Fit Decreasing/Number of Bins = %d.';
        Conclusion = sprintf(formatSpec,bin_number);
        display(Conclusion)
        EmptySpace = (1 - (sum(W)/(bin_number*C)))*100
        toc
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chos1==2,
    n = length(W);
    Bins = C*ones(1,n);
    i = 1:C;
    bin_number = 1; 
    A = zeros(n,n);
    a = n;
    flags = zeros(1,n);
    Wn = sort(W,'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('0 - Both Algorithms');
        disp('1 - First Fit Method');
        disp('2 - Best Fit Method');
        chos2=input('Now please choose which method you want to Bin Pack: ');
        if isempty(chos2)
            chos2=0;
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if chos2==0 || chos2==1,
            tic
            ag = A; Bins7 = Bins; 
            if Wn(1)> Bins7(bin_number)                                     % Because the first item would be the highest value in the list of objects since its in dec order, 
                                                                            % if it doesnt fit into the bin, the packing cannot be carried out.
                display('The packing is impossible')
            else
                k=1;
   
                while k<=n                                                   % The While loop gets executed by using k iterations from 1 to n
        
                    if  W(k)<=Bins7(bin_number)                                 % First of all we start by selecting the first term in the list of objects, to see whether it will fit into the first bin
                        ag(bin_number,k)=k;                                  % If the item fits, we print out the first answer matrix A with the first term in the first row which is bin 1, and so on.
                        Bins7(bin_number)=Bins7(bin_number) - W(k);             % We now determine the residual capacity of the bin after the item has been packed into it, so we know the remaining space left in the bins.
                        k=k+1;                                              % We then move onto the next item in the list and carry out the same procedure again until all the items are packed
                    else
                        bin_number=bin_number+1;                            % If the item doesnt fit, it will open up a new bin and pack it into that
                    end;
        %display(Wn(k:end));                                                % This will display the remaining weights in the list everytime an object has been packed.
        
                end
            end
            
        Answer7 = ag;
        display(Answer7)
        formatSpec = 'Normal Distribution/First Fit/Number of Bins = %d.';
        Conclusion = sprintf(formatSpec,bin_number);
        display(Conclusion) 
        EmptySpace = (1 - (sum(W)/(bin_number*C)))*100
        toc
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if chos2==0 || chos2==2,
            tic
            ah = A; Bins8 = Bins; bin_number=0;
            
             if Wn(1)> Bins8(1)                                                 % Because the first item would be the highest value in the list of objects since its in dec order, 
                  
                 display('The packing is impossible')                           % if it doesnt fit into the bin, the packing cannot be carried out.display('The packing is impossible')
            
             else
   
                while a~=0
                    bin_number=bin_number + 1;    
                    for j=1:n                                                    % The While loop gets executed by using k iterations from 1 to n
                        if  flags(j)==0 && W(j)<=Bins8(bin_number)                            % First of all we start by selecting the first term in the list of objects, to see whether it will fit into the first bin
                            ah(bin_number,j)=j;                                 % If the item fits, we print out the first answer matrix A with the first term in the first row which is bin 1, and so on.
                            a = a - 1;
                            Bins8(bin_number) = Bins8(bin_number) - W(j);
                            %Wn(j)= Wn(j) + 10*C;   
                            flags(j)=1;                                          % We now determine the residual capacity of the bin after the item has been packed into it, so we know the remaining space left in the bins.
                        end                                                    % We then move onto the next item in the list and carry out the same procedure again until all the items are packed
                    end                                                         % If the item doesnt fit, it will open up a new bin and pack it into that1
   
                end
   
            
             end
             
        Answer8 = ah;
        display(Answer8)
        formatSpec = 'Normal Distribution/Best Fit/Number of Bins = %d.';
        Conclusion = sprintf(formatSpec,bin_number);
        display(Conclusion)
        EmptySpace = (1 - (sum(W)/(bin_number*C)))*100
        toc
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POISSON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DISTRIBUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if chos==3,

        
    W = poissrnd(lambda,1,items)
    n = length(W);
    Bins = C*ones(1,n);
    i = 1:C;
    bin_number = 1; 
    A = zeros(n,n);
    a = n;
    flags = zeros(1,n);
    
    disp('1 - Ordered List(Decreasing)');
    disp('2 - Non-Ordered List');
    chos1=input('Please choose which arrangement you would like to test: ');
   
    
    if chos1==1, 
        
        Wn = sort(W,'descend')
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('0 - Both Algorithms');
        disp('1 - First Fit Method');
        disp('2 - Best Fit Method');
        chos2=input('Now please choose which method you want to Bin Pack: ');
        if isempty(chos2)
            chos2=0;
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if chos2==0 || chos2==1,
            tic
            ai = A; Bins9 = Bins; 
            if Wn(1)> Bins9(bin_number)                                             % Because the first item would be the highest value in the list of objects since its in dec order, 
                                                                                    % if it doesnt fit into the bin, the packing cannot be carried out.
                display('The packing is impossible')
            else
                k=1;
   
                while k<=n                                                              % The While loop gets executed by using k iterations from 1 to n
        
                    if  Wn(k)<=Bins9(bin_number)                                    % First of all we start by selecting the first term in the list of objects, to see whether it will fit into the first bin
                        ai(bin_number,k)=k;                                         % If the item fits, we print out the first answer matrix A with the first term in the first row which is bin 1, and so on.
                        Bins9(bin_number)=Bins9(bin_number) - Wn(k);                % We now determine the residual capacity of the bin after the item has been packed into it, so we know the remaining space left in the bins.
                        k=k+1;                                                       % We then move onto the next item in the list and carry out the same procedure again until all the items are packed
                    else
                        bin_number=bin_number+1;                                    % If the item doesnt fit, it will open up a new bin and pack it into that
                    end;
                   %display(Wn(k:end));                                             % This will display the remaining weights in the list everytime an object has been packed.
        
                end
            end
            
   Answer9 = ai;
   display(Answer9)
   formatSpec = 'Poisson Distribution/First Fit Decreasing/Number of Bins = %d.';
   Conclusion = sprintf(formatSpec,bin_number);
   display(Conclusion) 
   EmptySpace = (1 - (sum(W)/(bin_number*C)))*100     
   toc
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if chos2==0 || chos2==2,
            tic
            aj = A; Bins10 = Bins; bin_number=0;
            
             if Wn(1)> Bins10(1)                                                            % Because the first item would be the highest value in the list of objects since its in dec order, 
                  
                 display('The packing is impossible')                                       % if it doesnt fit into the bin, the packing cannot be carried out.display('The packing is impossible')
            
             else
   
                while a~=0
                    bin_number=bin_number + 1;    
                    for j=1:n                                                                   % The While loop gets executed by using k iterations from 1 to n
                        if  flags(j)==0 && Wn(j)<=Bins10(bin_number)                            % First of all we start by selecting the first term in the list of objects, to see whether it will fit into the first bin
                            aj(bin_number,j)=j;                                             % If the item fits, we print out the first answer matrix A with the first term in the first row which is bin 1, and so on.
                            a = a - 1;
                            Bins10(bin_number) = Bins10(bin_number) - Wn(j);
                            %Wn(j)= Wn(j) + 10*C;   
                            flags(j)=1;                                                     % We now determine the residual capacity of the bin after the item has been packed into it, so we know the remaining space left in the bins.
                        end
                    end                                                                   % We then move onto the next item in the list and carry out the same procedure again until all the items are packed
                end                                                                        % If the item doesnt fit, it will open up a new bin and pack it into that1
             end
        Answer10 = aj;
        display(Answer10)
        formatSpec = 'Poisson Distribution/Best Fit Decreasing/Number of Bins = %d.';
        Conclusion = sprintf(formatSpec,bin_number);
        display(Conclusion)
        EmptySpace = (1 - (sum(W)/(bin_number*C)))*100
        toc
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chos1==2,
    n = length(W);
    Bins = C*ones(1,n);
    i = 1:C;
    bin_number = 1; 
    A = zeros(n,n);
    a = n;
    flags = zeros(1,n);
    Wn = sort(W,'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('0 - Both Algorithms');
        disp('1 - First Fit Method');
        disp('2 - Best Fit Method');
        chos2=input('Now please choose which method you want to Bin Pack: ');
        if isempty(chos2)
            chos2=0;
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if chos2==0 || chos2==1,
            tic
            ak = A; Bins11 = Bins; 
            if Wn(1)> Bins11(bin_number)                     % Because the first item would be the highest value in the list of objects since its in dec order, 
                                                            % if it doesnt fit into the bin, the packing cannot be carried out.
                display('The packing is impossible')
            else
                k=1;
   
                while k<=n                                                      % The While loop gets executed by using k iterations from 1 to n
        
                    if  W(k)<=Bins11(bin_number)                                % First of all we start by selecting the first term in the list of objects, to see whether it will fit into the first bin
                        ak(bin_number,k)=k;                                     % If the item fits, we print out the first answer matrix A with the first term in the first row which is bin 1, and so on.
                        Bins11(bin_number)=Bins11(bin_number) - W(k);           % We now determine the residual capacity of the bin after the item has been packed into it, so we know the remaining space left in the bins.
                        k=k+1;                                                  % We then move onto the next item in the list and carry out the same procedure again until all the items are packed
                    else
                        bin_number=bin_number+1;                                % If the item doesnt fit, it will open up a new bin and pack it into that
                    end;
        %display(Wn(k:end));                                                    % This will display the remaining weights in the list everytime an object has been packed.
        
                end
            end
            
        Answer11 = ak;
        display(Answer11)
        formatSpec = 'Poisson Distribution/First Fit/Number of Bins = %d.';
        Conclusion = sprintf(formatSpec,bin_number);
        display(Conclusion) 
        EmptySpace = (1 - (sum(W)/(bin_number*C)))*100
        toc
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if chos2==0 || chos2==2,
            tic
            al = A; Bins12 = Bins; bin_number=0;
            
             if Wn(1)> Bins12(1)                                                % Because the first item would be the highest value in the list of objects since its in dec order, 
                    
                 display('The packing is impossible')                                   % if it doesnt fit into the bin, the packing cannot be carried out.display('The packing is impossible')
            
             else
   
                while a~=0
                    bin_number=bin_number + 1;    
                    for j=1:n                                                       % The While loop gets executed by using k iterations from 1 to n
                        if  flags(j)==0 && W(j)<=Bins12(bin_number)                            % First of all we start by selecting the first term in the list of objects, to see whether it will fit into the first bin
                            al(bin_number,j)=j;                                     % If the item fits, we print out the first answer matrix A with the first term in the first row which is bin 1, and so on.
                            a = a - 1;
                            Bins12(bin_number) = Bins12(bin_number) - W(j);
                            %Wn(j)= Wn(j) + 10*C;   
                            flags(j)=1;                                              % We now determine the residual capacity of the bin after the item has been packed into it, so we know the remaining space left in the bins.
                        end                                                          % We then move onto the next item in the list and carry out the same procedure again until all the items are packed
                    end                                                              % If the item doesnt fit, it will open up a new bin and pack it into that1
   
                end
   
            
             end
             
        Answer12 = al;
        display(Answer12)
        formatSpec = 'Poisson Distribution/Best Fit/Number of Bins = %d.';
        Conclusion = sprintf(formatSpec,bin_number);
        display(Conclusion)
        EmptySpace = (1 - (sum(W)/(bin_number*C)))*100
        toc
        end
    end
end

    
%end
