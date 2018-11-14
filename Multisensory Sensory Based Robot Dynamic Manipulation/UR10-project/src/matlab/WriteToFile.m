function  WriteToFile( fid, input_matrix, name ,varargin )
%WRITETOFILE Summary of this function goes here
%   Detailed explanation goes here

sz = size(input_matrix);

lineLength = 3000;
if nargin == 3
    if length(sz) > 2
        for iii = 1:sz(3)
            fprintf(fid,strcat(name,'(:,:,',num2str(iii),')= ')) ;
            fprintf(fid,'[');
            for ii = 1:sz(1)
                counter = 0;
                for i = 1:sz(2)
                    
                    charVec = char(input_matrix(ii,i,iii));
                    counter = counter +length(charVec)+1;
                    if counter < 2000
                        if i == sz(2)
                            fprintf(fid,'%s',input_matrix(ii,i,iii));
                        else
                            fprintf(fid,'%s,',input_matrix(ii,i,iii));
                        end
                    else
                        counter = 0;
                        fprintf(fid,' ...');
                        fprintf(fid,'\n');
                        if i == sz(2)
                            fprintf(fid,'%s',input_matrix(ii,i,iii));
                        else
                            fprintf(fid,'%s,',input_matrix(ii,i,iii));
                        end
                    end
                end
                
                if ii ~= sz(1)
                    fprintf(fid,';');
                    fprintf(fid,'\n');
                else
                end
                
            end
            fprintf(fid,'];');
            
            fprintf(fid,'\n');
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        
    else
        
        
        fprintf(fid,strcat(name,'= ')) ;
        fprintf(fid,'[');
        for ii = 1:sz(1)
            counter = 0;
            for i = 1:sz(2)
                bigTerm = false;
                charVec = char(input_matrix(ii,i));
                if length(charVec) > lineLength
                    bigTerm = true;
                end
                counter = counter +length(charVec)+1;
                if counter < lineLength
                    if i == sz(2)
                        fprintf(fid,'%s',input_matrix(ii,i));
                    else
                        fprintf(fid,'%s,',input_matrix(ii,i));
                    end
                else
                    counter = 0;
                    if bigTerm == false
%                         fprintf(fid,' ...');
                        fprintf(fid,'\n');
                        if i == sz(2)
                            fprintf(fid,'%s',input_matrix(ii,i));
                        else
                            fprintf(fid,'%s,',input_matrix(ii,i));
                        end
                    else
                        while bigTerm
                            
                            for ind = lineLength:-1:1
                                if strcmp(charVec(ind),'+') || strcmp(charVec(ind),'-')
                                    break
                                end
                                assert(ind~=1,'Big Term witout operators');
                            end
                            a = charVec(1:ind-1);
                            b = charVec(ind:end);

                            fprintf(fid,'%s',a);
%                             fprintf(fid,' ...');
                            fprintf(fid,'\n');
%                             fprintf(fid,'%s',b);

                            
                             
                            if length(b) < lineLength
                                bigTerm = false;
                                if i == sz(2)
                                    fprintf(fid,'%s',b);
                                else
                                    fprintf(fid,'%s,',b);
                                end
                            else
                                charVec = b;
                            end
                                
                        end
                    end
                end

            end
            
            if ii ~= sz(1)
                fprintf(fid,';');
                fprintf(fid,'\n');
            else
            end
            
        end
        fprintf(fid,'];');
        
        fprintf(fid,'\n');
    end
    
elseif (nargin == 4 && strcmp(varargin{1},'noSpace'))
    
    if length(sz) > 2
        for iii = 1:sz(3)
            fprintf(fid,strcat(name,'(:,:,',num2str(iii),')= ')) ;
            fprintf(fid,' ');
            for ii = 1:sz(1)
                counter = 0;
                for i = 1:sz(2)
                    
                    charVec = char(input_matrix(ii,i,iii));
                    counter = counter +length(charVec)+1;
                    if counter < 2000
                        if i == sz(2)
                            fprintf(fid,'%s',input_matrix(ii,i,iii));
                        else
                            fprintf(fid,'%s,',input_matrix(ii,i,iii));
                        end
                    else
                        counter = 0;
%                         fprintf(fid,' ...');
                        fprintf(fid,'\n');
                        if i == sz(2)
                            fprintf(fid,'%s',input_matrix(ii,i,iii));
                        else
                            fprintf(fid,'%s,',input_matrix(ii,i,iii));
                        end
                    end
                end
                
                if ii ~= sz(1)
                    fprintf(fid,';');
                    fprintf(fid,'\n');
                else
                end
                
            end
            fprintf(fid,';');
            
            fprintf(fid,'\n');
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        
    else
        
        
        fprintf(fid,strcat(name,'= ')) ;
        fprintf(fid,' ');
        for ii = 1:sz(1)
            counter = 0;
            for i = 1:sz(2)
                bigTerm = false;
                charVec = ccode(input_matrix(ii,i));
                if length(charVec) > lineLength
                    bigTerm = true;
                end
                counter = counter +length(charVec)+1;
                if counter < lineLength
                    if i == sz(2)
                        fprintf(fid,'%s',ccode(input_matrix(ii,i)));
                    else
                        fprintf(fid,'%s,',ccode(input_matrix(ii,i)));
                    end
                else
                    counter = 0;
                    if bigTerm == false
                        fprintf(fid,' ...');
                        fprintf(fid,'\n');
                        if i == sz(2)
                            fprintf(fid,'%s',ccode(input_matrix(ii,i)));
                        else
                            fprintf(fid,'%s,',ccode(input_matrix(ii,i)));
                        end
                    else
                        while bigTerm
                            
                            for ind = lineLength:-1:1
                                if strcmp(charVec(ind),'+') || strcmp(charVec(ind),'-')
                                    break
                                end
                                assert(ind~=1,'Big Term witout operators');
                            end
                            a = charVec(1:ind-1);
                            b = charVec(ind:end);

                            fprintf(fid,'%s',a);
%                             fprintf(fid,' ...');
                            fprintf(fid,'\n');
%                             fprintf(fid,'%s',b);

                            
                             
                            if length(b) < lineLength
                                bigTerm = false;
                                if i == sz(2)
                                    fprintf(fid,'%s',b);
                                else
                                    fprintf(fid,'%s,',b);
                                end
                            else
                                charVec = b;
                            end
                                
                        end
                    end
                end

            end
            
            if ii ~= sz(1)
                fprintf(fid,';');
                fprintf(fid,'\n');
            else
            end
            
        end
        fprintf(fid,';');
        
        fprintf(fid,'\n');
    end
    
elseif strcmp(varargin{1},'verbose')
    
    if length(sz) > 2
        for iii = 1:sz(3)
            for ii = 1:sz(1)
                for i = 1:sz(2)
                    fprintf(fid,strcat(name,'(',num2str(ii),',',num2str(i),',',num2str(iii),')= '));
                    fprintf(fid,'%s',input_matrix(ii,i,iii));
                    fprintf(fid,';');
                    fprintf(fid,'\n');
                end
            end
            fprintf(fid,'\n');
            fprintf(fid,'\n');
        end
        
        
    else
        
        
        for ii = 1:sz(1)
            for i = 1:sz(2)
                fprintf(fid,strcat(name,'(',num2str(ii-1),',',num2str(i-1),')= '));
                fprintf(fid,'%s',ccode(input_matrix(ii,i)));
                fprintf(fid,';');
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        
    end

end

    
if nargin == 4 && strcmp(varargin{1},'noSpace')
    
else
    fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf(fid,'\n');
    fprintf(fid,'\n');

end
end

