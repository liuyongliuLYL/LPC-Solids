function A = readBinary(filename)

fid = fopen(filename, 'rb');

% 读取二维 double 数组
A = fread(fid,'double');

% 关闭文件
fclose(fid);

end