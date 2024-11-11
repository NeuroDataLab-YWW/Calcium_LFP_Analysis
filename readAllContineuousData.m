function samples = readAllContineuousData(filename)

fid = fopen(filename);

samples=[];
% timestamp = fread(fid, 1, 'int64',0,'l');
% N = fread(fid, 1, 'uint16',0,'l');
% recordingNumber = fread(fid, 1, 'uint16', 0, 'l');
fseek(fid,1036,'cof');
samples = [samples;fread(fid, 1024, 'int16',0,'b')*0.195];
while ~feof(fid)
    fseek(fid,22,'cof');
    samples = [samples;fread(fid, 1024, 'int16',0,'b')*0.195];    
end
samples(end-4:end) = [];
fclose(fid);
