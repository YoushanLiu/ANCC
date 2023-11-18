function [resp] = readpaz(filename)

resp.zeros = [];
resp.poles = [];
resp.CONSTANT = [];

fid = fopen(filename, 'r');
    while(~feof(fid))
       line = fgetl(fid);
       line = strtrim(line);
       S = regexp(line, '^[^\*]+.', 'match');
       if (~isempty(S))
           S = regexp(line, '[ \t]*ZEROS[ \t]+\d', 'match');
           if (~isempty(S))
               S = regexp(line, '[ \t]', 'split');
               nz = str2double(S{2});
               resp.zeros = zeros(nz, 1);
               for k = 1:1:nz
                   line = fgetl(fid);
                   line = strtrim(line);
                   S = regexp(line, '[ \t]*[+-]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?[ \t]+[+-]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match');
                   if (~isempty(S))
                       S = regexp(line, '[ \t]', 'split');
                       resp.zeros(k) = str2double(S{1}) + 1i*str2double(S{2});
                   else
                       break;
                   end
               end
           end
           S = regexp(line, '[ \t]*POLES[ \t]+\d', 'match');
           if (~isempty(S))
               S = regexp(line, '[ \t]', 'split');
               np = str2double(S{2});
               resp.poles = zeros(np, 1);
               for k = 1:1:np
                   line = fgetl(fid);
                   line = strtrim(line);
                   S = regexp(line, '[ \t]*[+-]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?[ \t]+[+-]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match');
                   if (~isempty(S))
                       S = regexp(line, '[ \t]', 'split');
                       resp.poles(k) = str2double(S{1}) + 1i*str2double(S{2});
                   else
                       break;
                   end
               end
           end
           S = regexp(line, '[ \t]*CONSTANT[ \t]+[+-]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match');
           if (~isempty(S))
               S = regexp(line, '[ \t]', 'split');
               resp.CONSTANT = str2double(S{2});
           end
       end
    end
fclose(fid);


% if (~exist(filename, 'file'))
%     return
% end
% %indx_slash = strfind(filename, '/');
% %indx_backslash = strfind(filename, '\');
% %if (~isempty(indx_slash))
% %    idx = indx_slash(end);
% %else
% %    idx = indx_backslash(end);
% %end
% %indx_dot = strfind(filename(idx+1:end), '.');
% %Instrument_Code = filename(idx+indx_dot(3)+2:idx+indx_dot(3)+2);
% 
% 
% fid = fopen(filename, 'r');
%    line = fgetl(fid);
%    line = strtrim(line);
%    line = regexp(line, '[\t ]*', 'split');
%    nzeros = round(str2double(line{2}));
% %    if ('H' == Instrument_Code)
% %        resp.zeros = zeros(nzeros+1, 1);
% %    else
%        resp.zeros = zeros(nzeros, 1);
% %    end
%    for j = 1:1:nzeros
%        line = fgetl(fid);
%        line = strtrim(line);
%        line = regexp(line, '[\t ]*', 'split');
%        resp.zeros(j) = str2double(line{1}) + 1i*str2double(line{2});
%    end
% %    if ('H' == Instrument_Code)
% %        resp.zeros(nzeros+1) = 0.0 + 0.0j;
% %    end
%    line = fgetl(fid);
%    line = strtrim(line);
%    line = regexp(line, ' ', 'split');
%    npoles = round(str2double(line{2}));
%    resp.poles = zeros(npoles, 1);
%    for j = 1:1:npoles
%        line = fgetl(fid);
%        line = strtrim(line);
%        line = regexp(line, '[\t ]*', 'split');
%        resp.poles(j) = str2double(line{1}) + 1i*str2double(line{2});
%    end
%    line = fgetl(fid);
%    line = strtrim(line);
%    line = regexp(line, '[\t ]*', 'split');
%    resp.CONSTANT = str2double(line{2});
% fclose(fid);
