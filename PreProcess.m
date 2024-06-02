
clear
fn = 1; 
% [filename, pathname] = uigetfile('*.*', 'Select Image File');
pathname = 'D:\Edu\MSc\MSU\Sci\Darkfield\Tail-less\11uM tailless\11uM tailless\20201125 ch4\';
filename = 'ch4_11uM_tailless.tif';

filePath = matlab.desktop.editor.getActiveFilename;
splf = split(filePath, "\");
splf(max(size(splf))) = [];
fn = split(filename,".");

log_file_name = string(datestr(now, 30)) + "_" + fn(1) + "_Log.txt";
fid = fopen(fullfile(join(splf, "\"), log_file_name), 'a'); %Create log file
disp('File is opening..')
Imtif = tiffread([pathname, filename]);
imtif_num = length(Imtif);
n = 5; % parameter that regulates intensity
fprintf(fid, '%s: %s\n', datestr(now, 31), "File is opened: " + pathname + filename);

%%
% convert struct to 3d array
Idouble = zeros(size(Imtif(1).data,1), size(Imtif(1).data,2), imtif_num);
for ss=1:imtif_num
    Idouble(:,:,ss) = im2double(Imtif(ss).data);
end

sz = size(Idouble);
set(0,'DefaultFigureVisible','on')

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.2, 0.7, 0.7]);
sv = sliceViewer(Idouble);
sv.DisplayRange = [mean2(Idouble(:,:,1)) - 7 * std2(Idouble(:,:,1)), mean2(Idouble(:,:,1)) + 7 * std2(Idouble(:,:,1))];
ls = round(sz(1) / 2 - sz(1) * 0.1);
ws = round(sz(2) / 2 - sz(2) * 0.1);
rec = drawrectangle('Position',[ls, ws, round(sz(1) * 0.1), round(sz(2) * 0.1)], 'FaceAlpha', 0);
disp('Please configure ROI for MT to fit it in, press enter and choose slices of interest')
pause  
posrec = round(rec.Position);
z1 = input('Please enter the number of the first frame ');
z2 = input('Please enter the number of the last frame ');
close(gcf)
fprintf(fid, '%s: %s%d,%d,%d,%d%s\n', datestr(now, 31), ...
    "User selected ROI with coordinates: [", ls, ws, round(sz(1) * 0.1), round(sz(2) * 0.1), "]");
fprintf(fid, '%s: %s\n', datestr(now, 31), ...
    "User selected frames from " + string(z1) + " to " + string(z2));
Idouble = Idouble(posrec(2):posrec(2)+posrec(4)-1, posrec(1):posrec(1)+posrec(3)-1, z1:z2);
imtot = size(Idouble,3);
frames = 1:imtot;
%% denoising
fprintf(fid, '%s: %s\n', datestr(now, 31), "Denoising is started");
disp('Denosing is in progress. Please wait...')
Id = zeros(size(Idouble));
for ii=z1:z2
    Item = Imtif(ii).data;
    Id(:,:,ii-z1+1) = double(Item(posrec(2):posrec(2)+posrec(4)-1, posrec(1):posrec(1)+posrec(3)-1));
end
Imagad = bm4d(Id,'Gauss',0,'mp',0,0);
fprintf(fid, '%s: %s\n', datestr(now, 31), "Denoising is completed");
% sliceViewer(Imagad)
%% contrast enhancment
Id = im2double(uint16(Imagad));
Imagecr = zeros(size(Id));
avg_tot = mean(Id(:));
sigma_tot = 0;
for ii=1:imtot
sigma_tot = sigma_tot + std2(Id(:,:,ii));
end
sigma_tot = sigma_tot / imtot;
avg_prev = avg_tot;
sigma_prev = sigma_tot;
for ii=1:imtot
    Itemp = Id(:,:,ii);
    avg = mean(Itemp, 'all');
    sigma = std2(Itemp);
    if sigma > 0.2*sigma_tot
        avg = avg_prev;
        sigma = sigma_prev;
    end
    avg_prev = avg;
    sigma_prev = sigma;
    Itemp=imadjust(Itemp, [max([0, avg-n*sigma]), min([1, avg+n*sigma])]);
    Imagecr(:,:,ii) = Itemp;
end
disp('Video has been prepared for MT tracking')
fprintf(fid, '%s: %s\n', datestr(now, 31), "Intensity is corrected, video has been prepared for MT tracking");
fclose(fid);
% figure
% sliceViewer(Imagecr)
