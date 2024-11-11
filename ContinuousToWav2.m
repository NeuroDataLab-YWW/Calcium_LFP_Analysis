function ContinuousToWav(Fs,Ref)
switch nargin
    case 0
        Fs = 30000;
        Ref = 0;
    case 1
        Ref = 0;
end
Rdata=0;
if Ref ~= 0
    [filenameR,pathnameR,indexR] = uigetfile('*.continuous','Select One Ref Data');
    if indexR 
    strR = [pathnameR filenameR];
    display('waiting for ref reading...');
    Rdata = readAllContineuousData(strR);
    display('ref reading completed!')
    Rdata = Rdata/1000;
    path = [pathnameR '*.continuous'];
    [filename,pathname,index] = uigetfile(path,'Select Multiple Data Files','MultiSelect','on');
    else
        display('please select the ref file');
        return;
    end
    
else
    [filename,pathname,index] = uigetfile('*.continuous','Select Multiple Data Files','MultiSelect','on');
end
if ~index
    display('please select the data file');
    return; 
end
if iscell(filename) 
        for i = 1:length(filename)
            str = [pathname filename{i}];
            display(filename{i});
            display('reading...');
            samples = readAllContineuousData(str);
            samples = samples/1000-Rdata;
            str2 = [str(1:end-11),'.wav'];
            audiowrite(str2,samples,Fs);
            display('read completed!');
        end
else
        str = [pathname filename];
        display(filename);
        display('reading...');
        samples = readAllContineuousData(str);
        samples = samples/1000-Rdata;
        str2 = [str(1:end-11),'.wav'];
        audiowrite(str2,samples,Fs);
        display('read completed!');
end    
end