function clearData( days )
%CLEARDATA This searces the user directorys and tutorial folders and
%removes the previos data runs that have not been edited in the last 'days'

%Get path name
if nargin<1
    days=0;
end
fclose all;
p = mfilename('fullpath');
p = p(1:length(p)-length(mfilename));
p = p(1:end -8);
removed = 0;
bytes = 0;
folder = 'User';
for fold = [1,2]
    listing = ls(strcat(p,folder,filesep));
    szl = size(listing);
    today = now();
    for i=3:szl(1)
        newpath = strcat(p,folder,filesep,listing(i,:));
        if  isdir(newpath)
            pathtoclear = strcat(newpath,filesep,'Logs',filesep);
            files = dir(pathtoclear);
            for j = 3:length(files)
                if today - datenum(files(j).date) > days || days == 0
                    filetogo = strcat(pathtoclear,filesep,files(j).name);
                    delete(filetogo)
                    removed = removed +1;
                    bytes = bytes + files(j).bytes;
                end
            end
            pathtoclear = strcat(newpath,filesep,'Data',filesep);
            files = dir(pathtoclear);
            for j = 3:length(files)
                if today - datenum(files(j).date) > days || days == 0
                    filetogo = strcat(pathtoclear,filesep,files(j).name);
                    delete(filetogo)
                    removed = removed +1;
                    bytes = bytes + files(j).bytes;
                end
            end

        end
    end
    folder = 'Tutorial';
end
fprintf(1,'Files Removed: %d\n',removed);
fprintf(1,'Bytes Recovered: %d\n',bytes);
end
