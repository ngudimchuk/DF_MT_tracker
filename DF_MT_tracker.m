%%%%%%%%%%% Settings and Options %%%%%%%%%%
calibration_on = 0; % if calibration mode is on (1) then only frames with MT seeds (full MT) should be selected for processing 
                 % parameter 'tube' should be adjusted during the clibration
                 % so that number of PFs in MT seed corresponds to 14 PF
pix = 72e-9; % pixel size in meters
fixed = 1; % (0) user defines MT region each frame, (1) region from 1st frame is used for each subsequent frame
show = 1;  % display (1) or don't display (0) tracking results of tip position for each frame
pad = 4; % defines padding of initial box in pixels
linfit = 40; % number of points behind MT tip used to to find MT axis (linear fit)
len = 40; % length of line along which new tip position (back and forward) will be searched (using linear fit, < linfit)
level_dist = 4; % number of sections within which PF number will be estimated along MT 
rtip = 0.5; % section size to estimate PF number, um
tube = 38; % intensity of full MT used as a reference to calculate PF number at MT tip and shaft (set after calibration)
lpf = 4; % minimum number of PFs that can be detected
dynam_len = 100; % number of coordinates near the tip that are recalculated at each iteration (another part of backbone is memorized)

%%
fid = fopen(fullfile(join(splf, "\"), log_file_name), 'a'); %Open log file again
if calibration_on
    fprintf(fid, '\n%s: %s\n', datestr(now, 31), "Calibration is started");
else
    fprintf(fid, '\n%s: %s\n', datestr(now, 31), "MT tracking is started");
end

fprintf(fid, '%s\n', "Parameters used: ");
fprintf(fid, '%s%d\n', "pix = ", pix);
fprintf(fid, '%s%d\n', "fixed = ", fixed);
fprintf(fid, '%s%d\n', "show = ", show);
fprintf(fid, '%s%d\n', "pad = ", pad);
fprintf(fid, '%s%d\n', "len = ", len);
fprintf(fid, '%s%d\n', "linfit = ", linfit);
fprintf(fid, '%s%d\n', "rtip = ", rtip);
fprintf(fid, '%s%d\n', "level_dist = ", level_dist);
fprintf(fid, '%s%d\n', "tube = ", tube);
fprintf(fid, '%s%d\n', "lpf = ", lpf);
fprintf(fid, '%s%d\n', "dynam_len = ", dynam_len);

acc = 0.1; % spatial sampling limit of algorithm in pixels
bb = ceil(pad);
if bb <= 1; bb = 2; end % insures that box isn't one pixel wide, resulting in error of Gaussian fit

tippos = cell(1,imtot);
backpos = cell(1,imtot);
cumstd = zeros(imtot,1);
mtleum = zeros(imtot,1);
pf_tip_num = zeros(imtot,level_dist);
tip_extens = zeros(imtot,1);
Intensity.data = cell(length(frames),1);
Intensity.coord = cell(length(frames),1);
Shape.xx = cell(length(frames),1);
Shape.yy = cell(length(frames),1);

if calibration_on
    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.2, 0.7, 0.7]);
    sv = sliceViewer(Idouble);
    imfin = input('Please enter the number of the last frame ');
    fprintf(fid, '%s%d\n', "Number of last frame for calibration = ", imfin);
    close(gcf)
else
    imfin = imtot;
end

ii = 1; previi = ii;
xi = []; yi = [];
xmt = []; ymt = [];
Tran = 0;
inc = 0;
roic = 0;
prev_ymt = [];
prev_xmt = [];
backup_xmt = [];
backup_ymt = [];
stdis(1:2) = [5 5];
[sz_x, sz_y, sz_z] = size(Imagad);

for bini = 1:14*4
    pflin(bini) = (bini/14)^2;
end

while ii <= imfin
    %% Define MT Region 
    I = Imagad(:,:,ii); 
    Inorm = Idouble(:,:,ii);
    Idat = Imtif(ii+z1-1).data;
    Idat = im2double(Idat(posrec(2):posrec(2)+posrec(4)-1,posrec(1):posrec(1)+posrec(3)-1));
    
    [Imin,Imax,xi,yi,chkpts,roic] = clickreader(Imagecr(:,:,ii),fixed,ii,xi,yi,tippos,roic,Inorm);
    if ii == 1
        while length(xi) ~= 2
            close all
            disp('Number of user defined points does not equal to 2 points, please redo selection')
            [Imin,Imax,xi,yi,chkpts,roic] = clickreader(Imagecr(:,:,ii),fixed,ii,xi,yi,tippos,roic,Inorm);
        end
        fprintf(fid, 'User selected MT ends with the following coordinates:\n');
        fprintf(fid, '%s%d%s%d%s%d%s%d\n', "x1 = ", xi(1), ", y1 = ", yi(1), ", x2 = ", xi(2), ", y2 = ", yi(2));
    end
    if chkpts
        previi = ii;
    end
    
    
    %% Defines region in which MT is contained
    xrange = sort([xi(1) xi(2)])+[-bb bb];
    yrange = sort([yi(1) yi(2)])+[-bb bb];
    
    if xrange(1) <= 0
        xrange(1) = 1;
    elseif yrange(1) <= 0
        yrange(1) = 1;
    end
    [rI,cI] = size(I);
    if xrange(2) > cI
        xrange(2) = cI;
    elseif yrange(2) > rI
        yrange(2) = rI;
    end
    
    if ii == 1
        %Sub image is defined by MT region and used to calculate center of
        %mass
        Ireg = I(yrange(1):yrange(2), xrange(1):xrange(2));
        rwgt = sum(Ireg(:));
        [rIr,cIr] = size(Ireg);
        
        xmass = 0; ymass = 0;
        for ct=1:cIr
            xmass = xmass + ct*sum(Ireg(:,ct));
        end
        for rt=1:rIr
            ymass = ymass + rt*sum(Ireg(rt,:));
        end
        xcm = xmass/rwgt;
        ycm = ymass/rwgt;
        
        regangle = atand(ycm/xcm);
        
        %MT region is rotated or not based on center of mass
        if regangle>45 || regangle<-45
            thk=1*abs(1/cosd(90-regangle)); % thickness of the MT "looks" bigger
            % as the MT rotates
            T=1; % transposed condition
        else
            thk=1*abs(1/cosd(regangle));
            T=0; % untransposed condition
        end
        lthk = 2 * thk;
        rthk = round(3 * thk);
        mtthk = 2;
        Tran = T;
    end
    
    %% Find MT Spline (backbone) based on Gaussian Fit of MT Region Intensities %%
    compos = [];
    tries = 0;
    flag = 0; % flag that is activated if MT's tip cannot be find
    flag2 = 0; % flag that is activated in case of abrupt change in MT/background intensity (if =1 then skip the frame)
    backup_xmt = prev_xmt; % backups will be used in case if some impurities are detected and this frame is skipped
    backup_ymt = prev_ymt;
    back_yi = yi;
    back_xi = xi;
    if T == 0 %region was not transposed
        while isempty(compos)
            tries = tries+1;
            if tries > 1
                if ind_tip > 4*(tries+1)
                    yi=back_yi;
                    xi=back_xi;
                    yi=round([yi(1);last_coord_y(ind_tip-5*tries+1)]);
                    xi=round([xi(1);last_coord_x(ind_tip-5*tries+1)]);
                else
                    xmt=backup_xmt;
                    ymt=backup_ymt;
                    flag=1;
                    break
                end
            end
            inc = (xi(2)-xi(1))/abs(xi(2)-xi(1)); %defines direction of scan through MT region
            scansz = abs(xi(2)-xi(1)+1);
            xmt = zeros(1,scansz);
            ymt = zeros(1,scansz);
            xmt_dir = cell(2,1);
            ymt_dir = cell(2,1);
            scanfl = 0;
            ishift = 0;
            if scansz > dynam_len && ii > 1
                xi(1) = xi(2) - inc * dynam_len;
            end
            
            for dir = [1,2] % scanning intensity forward and backward
                if scansz > dynam_len && ~isempty(prev_ymt)
                    gauss = cell(1, dynam_len+1);
                    if dir == 1
                        scan = scansz - dynam_len;
                        init_scan = scan - 1;
                        ymt_dir{dir}(1:scansz-dynam_len) = prev_ymt(1:scansz-dynam_len);
                        xmt_dir{dir}(1:scansz-dynam_len) = prev_xmt(1:scansz-dynam_len);
                    else
                        scan = 1;
                        init_scan = scan;
                        ymt_dir{dir}((dynam_len+1):scansz) = flip(prev_ymt(1:scansz-dynam_len));
                        xmt_dir{dir}((dynam_len+1):scansz) = flip(prev_xmt(1:scansz-dynam_len));
                    end
                else
                    gauss = cell(1, scansz);
                    scan = 1;
                    init_scan = scan;
                end
                inc_dir=[inc; -inc];
                
                for ic = xi(dir):inc_dir(dir):xi(dir+(1*dir==1)-(1*dir==2))
                    xmt_dir{dir}(scan) = ic;
                    if scan == 1
                        xdata = (yi(dir)-pad:yi(dir)+pad)';
                        ydata = I(yi(dir)-pad:yi(dir)+pad,ic);
                        ydata = ydata-min(ydata)+1;
                        lowb = [0 yi(dir)-pad lthk min(ydata)];
                        uppb = [max(ydata)+200 yi(dir)+pad rthk max(ydata)];
                        appval = [max(ydata),yi(dir), lthk, min(ydata)];
                        gauss{scan} = gauss_func(xdata, ydata, appval, lowb, uppb);
                        dg = gauss{scan}(2)-yi(dir);
                        prev_st = yi(dir);
                    else
                        if dir==1
                            if ii == 1 || scan > length(prev_ymt)
                                stval = gauss{scan-1}(2);
                            else
                                stval = prev_ymt(scan);
                            end
                        else
                            if ii == 1 || gauss{scan-1}(2) < prev_ymt(end) || scan > length(prev_ymt)  %|| scan <= dynam_len 
                                stval = gauss{scan-1}(2);
                            else 
                                ishift = ishift + 1;
                                stval = flip(prev_ymt);
                                stval = stval(ishift);
                            end
                        end
                        rstval = max([pad+1, round(stval)]);
                        if rstval > sz_x + pad + 1
                            rstval = sz_x - pad - 1;
                        end
                        xdata = (rstval-pad:rstval+pad)';
                        ydata = I(rstval-pad:rstval+pad,ic);
                        ydata = ydata-min(ydata)+1;
                        appval = [max(ydata),rstval,lthk,min(ydata)];
                        if ~isempty(gauss{scan-init_scan})
                            appval = gauss{scan-init_scan};
                        end
                        lowb = [0 rstval-pad lthk min(ydata)];
                        uppb = [max(ydata)+200 rstval+pad rthk max(ydata)];
                        gauss{scan} = gauss_func(xdata,ydata,appval,lowb,uppb);
                        dg = gauss{scan}(2) - stval;
                        prev_st = stval;
                    end
                    
                    iac = pad;
                    while abs(dg) > stdis(dir) % if the new position is more than 3 pix far from the last one, then constrain the region for gaussian fitting
                        iac = iac / 1.5;
                        xmid = find(xdata > prev_st, 1) - 1;
                        leftb = max([1, ceil(xmid-iac)]);
                        rigb = min([length(xdata), ceil(xmid+iac)]);
                        lowb1 = [0 xdata(leftb) lthk min(ydata)];
                        uppb1 = [max(ydata)+200 xdata(rigb) rthk max(ydata)];
                        gauss{scan} = gauss_func(xdata(leftb:rigb), ydata(leftb:rigb), appval, lowb1, uppb1);
                        dg = gauss{scan}(2) - prev_st;
                    end
                    
                    ymt_dir{dir}(scan) = gauss{scan}(2);
                    scan = scan + 1;
                end
            end
            ymt_dir{2} = flip(ymt_dir{2});
            non_dynam = max([1, scansz - dynam_len]);
            compos = non_dynam + find(abs(ymt_dir{1}(non_dynam:end)-ymt_dir{2}(non_dynam:end)) < 2, 1, 'last') - 1; % find where two lines intercept
        end
        if flag == 0
            xmt = xmt_dir{1};
            ymt = smooth([ymt_dir{1}(1:compos-1), ymt_dir{2}(compos:end)], 0.3, 'loess')';
            if ii > 1
                ml = min([length(ymt), length(prev_ymt)]);
                if mean(ymt(1:ml)) - mean(prev_ymt(1:ml))>5
                    flag2 = 1;
                end
            end
            prev_xmt = xmt;
            prev_ymt = ymt;
        end
        
    else  %region was transposed
        while isempty(compos)
            tries = tries+1;
            if tries > 1
                if ind_tip > 4*(tries+1)
                    yi=back_yi;
                    xi=back_xi;
                    yi=round([yi(1); last_coord_y(ind_tip-4*tries+1)]);
                    xi=round([xi(1); last_coord_x(ind_tip-4*tries+1)]);
                else
                    xmt = prev_xmt;
                    ymt = prev_ymt;
                    flag = 1;
                    break
                end
            end
            inc = (yi(2)-yi(1)) / abs(yi(2)-yi(1)); % defines direction of scan through MT region
            scansz = abs(yi(2) - yi(1))+1;
            xmt_dir = cell(2,1);
            ymt_dir = cell(2,1);
            scanfl = 0;
            ishift = 0;
            if scansz > dynam_len && ii > 1
                yi(1) = yi(2) - inc * dynam_len;
            end
            
            for dir = [1,2] % scanning intensity forward and backward
                if scansz > dynam_len && ii > 1
                    gauss = cell(1, dynam_len+1);
                    
                    if dir == 1
                        scan = scansz - dynam_len;
                        init_scan = scan - 1;
                        ymt_dir{dir}(1:scansz-dynam_len) = prev_ymt(1:scansz-dynam_len);
                        xmt_dir{dir}(1:scansz-dynam_len) = prev_xmt(1:scansz-dynam_len);
                    else
                        scan = 1;
                        init_scan = scan;
                        ymt_dir{dir}((dynam_len+1):scansz) = flip(prev_ymt(1:scansz-dynam_len));
                        xmt_dir{dir}((dynam_len+1):scansz) = flip(prev_xmt(1:scansz-dynam_len));
                    end
                else
                    gauss = cell(1, scansz);
                    scan = 1;
                    init_scan = scan;
                end
                inc_dir=[inc; -inc];
                
                
                for ir = yi(dir):inc_dir(dir):yi(dir+(1*dir==1)-(1*dir==2))
                    ymt_dir{dir}(scan) = ir;
                    if scan == 1
                        xdata = (xi(dir)-pad:xi(dir)+pad);
                        ydata = I(ir,xi(dir)-pad:xi(dir)+pad);
                        ydata = ydata-min(ydata)+1;
                        lowb = [0 xi(dir)-pad lthk min(ydata)];
                        uppb = [max(ydata)+200 xi(dir)+pad rthk max(ydata)];
                        appval = [max(ydata),xi(dir),lthk, min(ydata)];
                        gauss{scan} = gauss_func(xdata,ydata,appval,lowb,uppb);
                        dg = gauss{scan}(2) - xi(dir);
                        prev_st = xi(dir);
                    else
                        if dir==1
                            if ii == 1 || scan > length(prev_xmt)
                                stval = gauss{scan-1}(2);
                            else
                                stval = prev_xmt(scan);
                            end
                        else
                            if ii == 1 || gauss{scan-1}(2) < prev_xmt(end-5) || scan > length(prev_xmt)  %|| scan <= dynam_len 
                                stval = gauss{scan-1}(2);
                            else 
                                ishift = ishift + 1;
                                stval = flip(prev_xmt);
                                stval = stval(ishift);
                            end
                        end
                        rstval = max([pad+1, round(stval)]);
                        if rstval > sz_y + pad + 1
                            rstval = sz_y - pad - 1;
                        end
                        xdata = (rstval-pad:rstval+pad);
                        ydata = I(ir, rstval-pad:rstval+pad);
                        ydata = ydata-min(ydata)+1;
                        appval = [max(ydata),rstval,lthk, min(ydata)];
                        if ~isempty(gauss{scan-init_scan})
                            appval = gauss{scan-init_scan};
                        end
                        lowb = [0 rstval-pad lthk min(ydata)];
                        uppb = [max(ydata)+200 rstval+pad rthk max(ydata)];
                        gauss{scan} = gauss_func(xdata,ydata,appval,lowb,uppb);
                        dg = gauss{scan}(2)-stval;
                        prev_st = stval;
                    end
                    
                    iac = pad;
                    while abs(dg)> stdis(dir) % if the new position is more than stdis pix far from the last one, then constrain the region for gaussian fitting
                        iac = iac/1.5;
                        xmid = find(xdata>prev_st,1)-1;
                        leftb = max([1,ceil(xmid-iac)]);
                        rigb = min([length(xdata),ceil(xmid+iac)]);
                        lowb1 = [0 xdata(leftb) lthk min(ydata)];
                        uppb1 = [max(ydata)+200 xdata(rigb) rthk max(ydata)];
                        gauss{scan} = gauss_func(xdata(leftb:rigb),ydata(leftb:rigb),appval,lowb1,uppb1);
                        dg = gauss{scan}(2)-prev_st;
                    end
                    xmt_dir{dir}(scan) = gauss{scan}(2);
                    scan = scan+1;
                end
            end
            
            xmt_dir{2} = flip(xmt_dir{2});
            non_dynam = max([1, scansz - dynam_len]);
            compos = non_dynam + find(abs(xmt_dir{1}(non_dynam:end)-xmt_dir{2}(non_dynam:end)) < 2, 1, 'last') - 1;% find where two lines intercect
        end
        if flag == 0
            ymt = ymt_dir{1};
            xmt = smooth([xmt_dir{1}(1:compos-1) ,xmt_dir{2}(compos:end)], 0.3, 'rloess')';
            if ii > 1
                ml = min([length(xmt), length(prev_xmt)]);
                if mean(xmt(1:ml)) - mean(prev_xmt(1:ml)) > 5
                    flag2 = 1;
                end
            end
            prev_xmt = xmt;
            prev_ymt = ymt;
        end
    end
    
    backpos{ii} = [xmt; ymt];
    yi = back_yi;
    xi = back_xi;
    
    %% Fitting MT Spline (backbone) end to Estimate MT axis (x''-y'')
    xline=cell(1);
    yline=cell(1);
    
    if T == 0 %region was not transposed
        pnts = min([length(xmt)-1, linfit]);
        pmt = polyfit(xmt(end-pnts:end), ymt(end-pnts:end),1); %linear fit of MT axis
        if inc > 0
            xline = ((round(max(xmt)-len)):inc:(round(max(xmt)))+len);
        else
            xline = ((round(min(xmt)+len)):inc:(round(min(xmt))-len));
        end
        
        yline = round(pmt(1).*xline+pmt(2));
    else
        pnts = min([length(xmt)-1, linfit]);
        pmt = polyfit(ymt(end-pnts:end), xmt(end-pnts:end),1); %linear fit of MT axis
        
        if inc > 0
            yline = ((round(max(ymt)-len)):inc:(round(max(ymt)+len)));
        else
            yline = ((round(min(ymt)+len)):inc:(round(min(ymt)-len)));
        end
        
        xline = round(pmt(1).*yline + pmt(2));
    end
    
    %get rid of points along MT axis that lie outside the image
    
    [Idr,Idc] = size(I);
    incol = find(xline>(ceil(thk)) & xline<(Idc-(ceil(thk))));
    inrow = find(yline>(ceil(thk)) & yline<(Idr-(ceil(thk))));
    if isequal(incol,inrow)
        xline = xline(incol);
        yline = yline(incol);
    elseif numel(incol) < numel(inrow)
        xline = xline(incol);
        yline = yline(incol);
    elseif numel(inrow) < numel(incol)
        xline = xline(inrow);
        yline = yline(inrow);
    end

    
    %% Find Brightness Profile along MT axis 
    [~, ad_line] = min(abs(xline - xmt(end)) + abs(yline - ymt(end)));
    xline = round([xmt, xline(ad_line+1:end)]);
    yline = round([ymt, yline(ad_line+1:end)]);
    br = zeros(1, size(xline, 2));
    br_noise = br;
    br_tip = br;
    br_tip_noise = br;
    if T == 0 && inc > 0
        bound = min([30, yline-1, size(Idat,1)-yline-1]);
        for i = 1:length(xline)
            br(i) = mean(Idat((yline(i)-mtthk):(yline(i)+mtthk), xline(i)));
            br_noise(i) = mean(Idat([yline(i)-bound:yline(i)-mtthk,...
                yline(i)+mtthk:yline(i)+bound], xline(i)));
            br_tip(i) = mean(I((yline(i)-mtthk):(yline(i)+mtthk),xline(i)));
            br_tip_noise(i) = mean(I([yline(i)-bound:yline(i)-mtthk,...
                yline(i)+mtthk:yline(i)+bound], xline(i)));
        end
        xx = (xmt-xmt(1)).*pix.*1e9/1000;
        yy = (ymt-ymt(1)).*pix.*1e9/1000;
    elseif T == 0 && inc<0
        bound = min([30, yline-1,size(Idat,1)-yline-1]);
        for i = 1:length(xline)
            br(i) = mean(Idat((yline(i)-mtthk):(yline(i)+mtthk), xline(i)));
            br_noise(i) = mean(Idat([yline(i)-bound:yline(i)-mtthk,...
                yline(i)+mtthk:yline(i)+bound], xline(i)));
            br_tip(i) = mean(I((yline(i)-mtthk):(yline(i)+mtthk), xline(i)));
            br_tip_noise(i) = mean(I([yline(i)-bound:yline(i)-mtthk,...
                yline(i)+mtthk:yline(i)+bound], xline(i)));
        end
        xx = (-(xmt-xmt(1))).*pix.*1e9/1000;
        yy = (-flip(ymt-ymt(1))).*pix.*1e9/1000;
    elseif T == 1 && inc>0
        bound = min([30, xline-1,size(Idat,2)-xline-1]);
        for i = 1:length(xline)
            br(i) = mean(Idat(yline(i),(xline(i)-mtthk):(xline(i)+mtthk)));
            br_noise(i) = mean(Idat(yline(i),[xline(i)-bound:xline(i)-mtthk,...
                xline(i)+mtthk:xline(i)+bound]));
            br_tip(i) = mean(I(yline(i),(xline(i)-mtthk):(xline(i)+mtthk)));
            br_tip_noise(i) = mean(I(yline(i),[xline(i)-bound:xline(i)-mtthk,...
                xline(i)+mtthk:xline(i)+bound]));
        end
        xx = (xmt-xmt(1)).*pix.*1e9/1000;
        yy = -(ymt-ymt(1)).*pix.*1e9/1000;
    elseif T == 1 && inc<0
        bound = min([30, xline-1,size(Idat,2)-xline-1]);
        for i = 1:length(xline)
            br(i) = mean(Idat(yline(i),(xline(i)-mtthk):(xline(i)+mtthk)));
            br_noise(i) = mean(Idat(yline(i),[xline(i)-bound:xline(i)-mtthk,...
                xline(i)+mtthk:xline(i)+bound]));
            br_tip(i) = mean(I(yline(i),(xline(i)-mtthk):(xline(i)+mtthk)));
            br_tip_noise(i) = mean(I(yline(i),[xline(i)-bound:xline(i)-mtthk,...
                xline(i)+mtthk:xline(i)+bound]));
        end
        xx = -(ymt-ymt(1)).*pix.*1e9/1000;
        yy = -(xmt-xmt(1)).*pix.*1e9/1000;
    end
    
    if ii == 1
        m_br_tip_prev = std(br_tip);
        m_br_tip_noise_prev = std(br_tip_noise);
        nonskip = 1;
    else
        m_br_tip_now = std(br_tip);
        m_br_tip_noise_now  = std(br_tip_noise);
        if (m_br_tip_now > m_br_tip_prev * 1.3 || m_br_tip_now < m_br_tip_prev * 0.7) || ...
            (m_br_tip_noise_now > m_br_tip_noise_prev * 1.3 || m_br_tip_noise_now < m_br_tip_noise_prev * 0.7)
            nonskip = 0;
            prev_xmt = backup_xmt;
            prev_ymt = backup_ymt;
        else
            nonskip = 1;
        end
        m_br_tip_prev = m_br_tip_now;
        m_br_tip_noise_prev = m_br_tip_noise_now;
    end
    
    if nonskip
    
        % determine distance along the MT axis at which brightness values were measured
        distmt = sqrt((xline-xline(1)).^2+(yline-yline(1)).^2); %fliplr()

        % background subtraction
        br_noise(br_noise > median(br_noise)+std(br_noise)) = median(br_noise);
        br_bg = (br - br_noise) ./ (max(br) - min(br_noise));
        
        if ii == 1
            intq = prctile(br_tip_noise, 65) - prctile(br_tip_noise, 25);
        end

        if isnan(median(br_tip_noise))
            fi = imgaussfilt(I, 30);
            br_tip_noise = median(fi(:));
        else
            br_tip_noise(br_tip_noise>(prctile(br_tip_noise, 25)+1.5*intq)...
                | br_tip_noise<(prctile(br_tip_noise, 25)-1.5*intq)) = prctile(br_tip_noise, 50);
            sm_br_tip = smooth(br_tip_noise, 0.5, 'rloess')';
            br_tip_noise = sm_br_tip;
        end
        
%         br_bg_tip = (br_tip - br_tip_noise) / prctile(br_tip - br_tip_noise, 80);
        br_bg_tip = (br_tip - br_tip_noise);
        % PF numbers determination
        [steps,~,~] = l1pwc(br_bg_tip, tube/10, 0); %min([20,length(xmt)])
        smint = (steps - max([0, min([tube*(4/14)^2, min(steps)])])) / tube;
        pfam = discretize(abs(smint), pflin);
        pfam(smint <= (1/14)^2) = 0;

        % tip position determination
%         cumpf = [0; cumsum(diff(pfam))];
        stip = max([5, length(xmt) - dynam_len, find(smooth(pfam,10) <= lpf,1)]);
        etip = stip + find(smooth(pfam(stip:end), 10) <= lpf-1, 1) - 1;
        ind_tip = stip + find(pfam(stip:etip) >= max([lpf, min(pfam(stip:etip))]), 1, 'last') - 1; 
        if isempty(ind_tip) || ind_tip < 3
            if stip >= etip
                ind_tip = stip;
            else
                ind_tip = length(pfam);
            end
        end
        full_mt_end = find(pfam(1:ind_tip) >= min([12, max(pfam(1:ind_tip))]), 1, 'last');
        tip_pntx = distmt(ind_tip);
        tippos{ii}(1,1) = xline(ind_tip);
        tippos{ii}(1,2) = yline(ind_tip);

        distum = distmt.*pix.*1e9/1000;
        pfam(isnan(pfam)) = 0;
        tip_pf = full_mt_end + find(pfam(full_mt_end:end) <= lpf, 1) - 1;

        alpf = lpf + 1;
        while isempty(tip_pf)
            alpf = alpf + 1;
            tip_pf = full_mt_end + find(pfam(full_mt_end:end) <= alpf, 1) - 1;
        end

        for ran = 1:level_dist
            sind = find(distum > distum(tip_pf) - ran*rtip, 1);
            eind = find(distum > distum(tip_pf) - (ran-1)*rtip, 1);
            pf_tip_num(ii,ran) = mean(pfam(sind:eind), 'omitnan');
        end

        tip_extens(ii) = distum(tip_pf) - distum(full_mt_end);
        mtleum(ii,1) = distum(tip_pf); % MT length in um
        last_coord_x = xline;
        last_coord_y = yline;
            %% Plot MT Spline (backbone) & axes on original image
            if show == 1 
                imshowpair(Imagecr(:,:,ii), Inorm, 'montage')
                hold on
                plot(xmt, ymt,'Color',[0.92,0.92,0.21],'LineWidth',1)
                scatter(tippos{ii}(1,1), tippos{ii}(1,2), 50, 'Marker','o', 'MarkerFaceColor',...
                    'None', 'MarkerEdgeColor',[0.00, 0.85, 1.00], 'LineWidth',1)
                scatter(tippos{ii}(1,1)+size(Idat,2), tippos{ii}(1,2), 50, 'Marker','o', 'MarkerFaceColor',...
                    'None', 'MarkerEdgeColor',[0.00, 0.85, 1.00], 'LineWidth',1)
                annotation('textbox', 'String', sprintf('Frame # %g',ii), 'position',[0.5,0.02,0.1,0.1], 'LineStyle','none');
                pause(0.001)
                hold off
                delete(findall(gcf, 'type', 'annotation'))
            end

            if flag2 == 0
                lnd = mslowess(distmt'.*pix.*1e9/1000, br_bg', 'Order', 2);
                Intensity.data{ii,1} = br_bg;
                Intensity.coord{ii,1} = distmt.*pix.*1e9/1000;

                Shape.xx{ii,1} = xx;
                Shape.yy{ii,1} = yy;
            else
                tippos{ii}(1,1) = tippos{ii-1}(1,1);
                tippos{ii}(1,2) = tippos{ii-1}(1,2);
                pf_tip_num(ii,:) = NaN;
                tip_extens(ii) = NaN;
            end
    end
    if nonskip == 0
        fprintf(fid, '%s\n', "frame # " + string(ii) + " is skipped");
        tippos{ii}(1,1) = tippos{ii-1}(1,1);
        tippos{ii}(1,2) = tippos{ii-1}(1,2);
        pf_tip_num(ii,:) = NaN;
        tip_extens(ii) = NaN;
    end
    ii = ii+1;
end

if calibration_on
    
    fprintf(fid, '%s: %s\n\n', datestr(now, 31), "Calibration is completed, recommended paramter tube = " + string(round((nanmean(pf_tip_num(1:imfin,level_dist)) / 14) ^ 2 * tube,2)));
else
    fprintf(fid, '%s: %s\n\n', datestr(now, 31), "MT tracking is completed");
end
fclose(fid);
%%
figure
hold on
cm = colormap(hot(size(pf_tip_num(1:imfin,:), 2) * 2));
for dd=1:size(pf_tip_num(1:imfin,:), 2)
scatter(frames(1:imfin), pf_tip_num(1:imfin,dd), 5, 'MarkerFaceAlpha',0.4, 'MarkerEdgeAlpha',0.6, 'MarkerFaceColor',cm(dd,:), 'MarkerEdgeColor',cm(dd,:))
plot(frames(1:imfin), smooth(frames(1:imfin), pf_tip_num(1:imfin,dd), 0.2, 'rlowess'), 'LineWidth',1.5, 'Color',cm(dd,:))
end
xlabel('Frame #')
ylabel('Mean PF number')
xlim([frames(1), frames(imfin)])
set(gca,'FontSize',13)
set(gcf,'Position',[408.2  320.2  685.6  252])
%%
tippos_m = cell2mat(tippos);
tippos_x = tippos_m(1:2:end);
tippos_y = tippos_m(2:2:end);


mtlength = sqrt((tippos_x - xi(1)).^2 + (tippos_y - yi(1)).^2) * pix * 1e9 / 1000;
%%
figure
hold on
plot(frames(1:imfin), mtlength(1:imfin),'.','Color','b','MarkerSize',9)
xlabel('Frame #')
smtl = smooth(frames(1:imfin), mtlength(1:imfin), 5/length(mtlength(1:imfin)),'rlowess');
plot(frames(1:imfin), smtl(1:imfin), 'b')
ylabel('MT length, \mum')
xlim([frames(1), frames(imfin)])
set(gca,'FontSize',13)
set(gcf,'Position',[408.2  320.2  685.6  252])
yyaxis right
plot(frames(1:imfin), tip_extens(1:imfin), '.','Color',[0.99,0.44,0.23])
ylabel('Tip extension, \mum') 
plot(frames(1:imfin), smooth(tip_extens(1:imfin), 15/length(mtlength(1:imfin)),'rlowess'), 'Color', [0.99,0.44,0.23])
ylim([0 1.0*max(tip_extens(1:imfin))])

%% save results in .csv-files
if calibration_on == 0
    fp = join(splf, "\");
    % PF number
    writematrix(frames(1:imfin), fp{1} + "\" + fn(1) + "_" + datestr(now, 30)+ "_" + "PF_number.csv", 'Delimiter', ';')
    for dd=1:size(pf_tip_num(1:imfin,:), 2)
        writematrix(pf_tip_num(1:imfin,dd)', fp{1} + "\" + fn(1) + "_" + datestr(now, 30)+ "_" + "PF_number.csv", 'Delimiter', ';','WriteMode', 'append')
    end
    % Tracked tip
    writematrix(frames(1:imfin), fp{1} + "\" + fn(1) + "_" + datestr(now, 30)+ "_" + "MT_length.csv", 'Delimiter', ';')
    writematrix(mtlength(1:imfin), fp{1} + "\" + fn(1) + "_" + datestr(now, 30)+ "_" + "MT_length.csv", 'Delimiter', ';','WriteMode', 'append')
    % Tip extension
    writematrix(frames(1:imfin), fp{1} + "\" + fn(1) + "_" + datestr(now, 30)+ "_" + "Tip_extension.csv", 'Delimiter', ';')
    writematrix(tip_extens(1:imfin)', fp{1} + "\" + fn(1) + "_" + datestr(now, 30)+ "_" + "Tip_extension.csv", 'Delimiter', ';','WriteMode', 'append')
end
%% clickreader
function [Imin, Imax, xi, yi, chkpts, roic] = clickreader(I, fixed, ii, xi, yi, tippos, roic, Iad)
Irange = max(I(:)) - min(I(:));
Imax = min(I(:)) + Irange;
Imin = min(I(:));
chkpts = 0;
switch (fixed)
    case 1
        if ii == 1 
            figure
            imshowpair(I, Iad, 'montage')
            disp('Define MT analysis region (two clicks)...')
            [xi, yi] = getpts(gcf); chkpts=1;
            xi = round(xi);
            yi = round(yi);
        elseif mod(ii-1, 1)==0
            xi = [xi(1); round(tippos{ii-1}(1,1))];
            yi = [yi(1); round(tippos{ii-1}(1,2))];
        end
    case 0
        if mod(ii-1,1) == 0
            figure('Units','centimeters','Position',[4 0 25.4 19.1])
            pos_mt = [0 0.07 0.3 0.9];
            axes('Position', pos_mt);
            set(gcf, 'Color',[1 1 1])
            imshow(I,[Imin Imax],'InitialMagnification','fit'); disp('Define MT analysis region (two clicks)...')
            [xi,yi] = getpts(gcf); chkpts=1;
            xi = round(xi);
            yi = round(yi);
        end
end
end

%% gauss_func 
function fit_res = gauss_func(xdata, ydata, x, lb, ub)
options = optimset('Display', 'off');
fit_res = lsqcurvefit(@fun_fit, x, xdata, ydata, lb, ub, options);
    function F = fun_fit(x, coordinates)
        F = x(1) .* exp(-(coordinates - x(2)).^2 / (2 * x(3)^2)) + x(4);
    end
end


