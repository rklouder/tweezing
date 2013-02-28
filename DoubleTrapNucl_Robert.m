function DoubleTrapNuclFXupdated(filename, settingsfile, savedir);
%
% Analyzes double trap tweezer data files
% Updated for Jasons Machine
% Stephan Grill, Sept 2004, Dec 3005
%
% USAGE:
%       DoubleTrapJ(filename, settingsfile, savedir)
%
%       filename is the name of the dataset to load, without extension
%       settingsfile is the name of the file where analysis variables are
%           passed.
%           Run DTall to create a sample settingsfile
%       savedir is the directory to save data summary files

close all; %closes all figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIND FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first find all the files

if ~exist('filename')
    disp('ERROR: No datafile specified.');
    return;
end;
% check if flie has an extension, remove it if it does
if isempty(findstr(filename,'.'))
    basefilename = filename;
else
    basefilename = filename(1:findstr(filename,'.')-1);
end;



disp(' ');
disp(['=>Analyzing dataset ' basefilename]);

%Check if a saving directory was specified and create it if it does not exist
if ~exist('savedir')
    disp('  No saving directory specified, using "./Summary".');
    savedir='Summary';
end;
if exist(savedir) ~= 7 %check if the directory exists
    mkdir(savedir);
    disp(['  Saving directory "./' savedir '" created.']);
end;

% all the files that will be used
f.txtfilename   = strcat(basefilename,'.txt');
f.datfilename   = strcat(basefilename,'.2d');
f.c1filename    = strcat(basefilename,'.2c1');
f.c2filename    = strcat(basefilename,'.2c2');

f.pdffilename   = strcat(savedir,'/',basefilename,'_s1.pdf');
f.pdf2filename   = strcat(savedir,'/',basefilename,'_s2.pdf');
f.fvpfilename   = strcat(savedir,'/',basefilename,'_fvp.mat');


%This is a hack since DTloadbinCH and DTloadbinJ use different formats for
%the filename
% HCH 2007-02-25
[pathstr, basefilename1, ext, versn] = fileparts(filename);
% HCH end
f.matfilename    = strcat(savedir,'/',basefilename1,'.mat');
f.FXfilename    = strcat(savedir,'/',basefilename1,'_FX.mat');
f.fvpfilename   = strcat(savedir,'/',basefilename1,'_fvp.mat');
%this is were the results are stored
an.file = basefilename;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD SETTINGSFILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load settings

%check if settings file was specified
if ~exist('settingsfile')
    disp('  No settings-file specified, using "settings.txt".');
    settingsfile='settings.txt';
end;

% load them
s = DTloadset(settingsfile,0,0);

%set arrays to zero
clear d; clear av; clear cal; clear event; clear mc;
cal.comment = ''; %just because not all files have this
cal.driftcorr = 1; %just because not all files have this



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the shit



if exist(f.fvpfilename)==0
    if exist(f.matfilename )==0
        % load binary data
        [d av cal] = DTloadbinCH(f,s);
    else
        disp(['->Loading "' f.matfilename  '"...']);
        load(f.matfilename );

        if ~exist('av')
            disp('ERROR. Bad file.');
            return;
        end;

        %just because not all files have these fields...
        if ~isfield(cal, 'comment'), cal.comment = ''; end;
        if ~isfield(cal, 'driftcorr'), cal.driftcorr = 1; end;

        if exist('d')
            disp(['  ',num2str(length(d.T1.X)),' lines of raw and averaged data read in.']);
        else
            if exist('av.d.T1.X')
                disp(['  ',num2str(length(av.d.T1.X)),' lines of averaged data read in.']);
            end;
        end;
    end;
    %only look at Y direction, remove X directions
    cal2 = cal;
    av2 = av;
    clear cal.rec; clear av;

    cal.rec.rate        = cal2.rec.rate;
    cal.rec.T1.offset   = cal2.rec.T1.Y.offset;
    cal.rec.T1.dr       = cal2.rec.T1.Y.dr;
    cal.rec.T1.k        = cal2.rec.T1.Y.k;

    cal.rec.T1.X.offset = cal2.rec.T1.X.offset;
    cal.rec.T1.X.dr     = cal2.rec.T1.X.dr;
    cal.rec.T1.X.k      = cal2.rec.T1.X.k;

    cal.rec.T2.offset   = cal2.rec.T2.Y.offset;
    cal.rec.T2.dr       = cal2.rec.T2.Y.dr;
    cal.rec.T2.k        = cal2.rec.T2.Y.k;

    cal.rec.T2.X.offset = cal2.rec.T2.X.offset;
    cal.rec.T2.X.dr     = cal2.rec.T2.X.dr;
    cal.rec.T2.X.k      = cal2.rec.T2.X.k;


    av.d.T1             = av2.d.T1.Y;
    av.d.T2             = av2.d.T2.Y;
    av.d.T2pos          = av2.d.T2.Ypos;
    av.d.time           = av2.d.time;

    % hang on to X data to control for flow
    av.d.T2x            = av2.d.T2.X;
    av.d.T1x            = av2.d.T1.X;

    clear cal2; clear av2;
    %only have y directions left

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALC FORCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Force...

    av.d.force.T1 = -(av.d.T1 - cal.rec.T1.offset) / cal.rec.T1.dr * cal.rec.T1.k;
    av.d.force.T2 = (av.d.T2 - cal.rec.T2.offset) / cal.rec.T2.dr * cal.rec.T2.k;

    av.d.force.T1x = (av.d.T1x - cal.rec.T1.X.offset) / cal.rec.T1.X.dr * cal.rec.T1.X.k;
    av.d.force.T2x = (av.d.T2x - cal.rec.T2.X.offset) / cal.rec.T2.X.dr * cal.rec.T2.X.k;

    av.d.avXforce = (av.d.force.T1x + av.d.force.T2x) / 2;
    av.d.avforce = (av.d.force.T1 + av.d.force.T2) / 2;

    %added trap displacements
    av.d.addedtd = (av.d.T2 - cal.rec.T2.offset) / cal.rec.T2.dr - (av.d.T1 - cal.rec.T1.offset) / cal.rec.T1.dr;

    av.fx = av.d;
    av.fx.Extension = av.fx.T2pos - av.fx.addedtd;
    av.fx.pos.T2 = av.fx.T2pos;
    av.fx = DTsplitFX(av.fx);

    save  (f.fvpfilename, 'av')

else
    disp(['->Loading "' f.fvpfilename  '"...']);
    load(f.fvpfilename);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loop through each pull %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist(f.FXfilename)==0
    nr=0;
    nf=0;
    sgtime=s.sgtime;
    for i=1: length(av.fx.split),

        handles.fig1 = figure(1); clf;
        hold on;

        if(av.fx.split(i).direction > 0),
            nr=nr+1;
            plot(av.fx.split(i).ext,av.fx.split(i).force,'-k');
            [rip(nr)]=[DTripLBupdated(i, av)];
        elseif(av.fx.split(i).direction < 0);

            nf=nf+1;
            plot(av.fx.split(i).ext,av.fx.split(i).force,'-r');
            [fold(nf)]=[DTripLBupdated(i, av)];
        end
    end

    save  (f.FXfilename, 'rip', 'fold')


else
    disp(['->Loading "' f.FXfilename  '"...']);
    load(f.FXfilename);
end

%%%%%%%%%%%%%%%%%% Plot the summary figure %%%%%%%%%%%%%%%%%%%%%%%
figure(5); clf; hold on;
i1=0;
i2=0;

subplot(2,1,1) %subplot with all unfold and refold sizes and forces
title([' WT Nucleosome ' basefilename1],'Interpreter','None');
axis off
dtext = 0.11; 
pytext = 1;
rytext = 1;
pfxtext = 0.3;
rxtext=0.54;
rfxtext = 0.9;
hold on

text(0.2, pytext, 'unfold')
pytext = pytext - dtext;
text(0, pytext, 'size upper(nm)')
%text(0.18, pytext, 'size lower(nm)')
text(pfxtext, pytext, 'force(pN)')
pytext = pytext - dtext;


text(0.7, rytext, 'refold')
rytext = rytext - dtext;
text(rxtext, rytext, 'size upper(nm)')
%text(rxtext + 0.18, rytext, 'size lower(nm)')
text(rfxtext, rytext, 'force(pN)')
rytext = rytext - dtext;


for i=1:length(av.fx.split),
    if(av.fx.split(i).direction > 0),
        i1=i1+1;
        if ~isempty(rip(i1).sizeMin)%pim changed
                hold on
                text(0, pytext, num2str(round(rip(i1).sizeMax*10)/10))%pim changed
                %text(0.18, pytext, num2str(round(rip(i1).sizeMin*10)/10))%pim changed
                text(pfxtext, pytext, num2str(round(rip(i1).forcerMin*10)/10))
                pytext = pytext - dtext;
        end

            elseif(av.fx.split(i).direction < 0),
                i2=i2+1;
                if ~isempty(fold(i2).sizeMin)                
                     hold on
                     text(rxtext, rytext, num2str(round(fold(i2).sizeMax*10)/10))
                     %text(rxtext + 0.18, rytext, num2str(round(fold(i2).sizeMin*10)/10))
                     text(rfxtext, rytext, num2str(round(fold(i2).forcerMax*10)/10))
                     
                      rytext = rytext - dtext;
                end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2)
i1=0;
i2=0;

offsetL=1100;
spacing=50;
for i=1:length(av.fx.split),
    if(av.fx.split(i).direction > 0),%for going up, the part that commented out will plot the lower rip size
        i1=i1+1;
        plot(av.fx.split(i).ext+(i-1)*spacing+offsetL,av.fx.split(i).force,'-k');
        hold on
        if ~isempty(rip(i1).sizeMin)%pim changed
            %plot ( rip(i1).extFoldMin+(i-1)*spacing+offsetL,rip(i1).forcerMin, 'm.')%plot the rip at forcerMin
            %plot ( rip(i1).extUnfoldMin+(i-1)*spacing+offsetL, rip(i1).forcerMin, 'm.')%plot the rip at forcerMin
            plot ( rip(i1).extFoldMax+(i-1)*spacing+offsetL,rip(i1).forcerMax, 'g.')
            plot ( rip(i1).extUnfoldMax+(i-1)*spacing+offsetL, rip(i1).forcerMax, 'g.')
            %for m=1:length(rip(i1).extUnfoldMin)
            %plot ( [rip(i1).extUnfoldMin(m)+(i-1)*spacing+offsetL rip(i1).extFoldMin(m)+(i-1)*spacing+offsetL],...
            %    [rip(i1).forcerMin(m) rip(i1).forcerMin(m)], 'Color', 'm')
            %end
            for m=1:length(rip(i1).extUnfoldMax)
            plot ( [rip(i1).extUnfoldMax(m)+(i-1)*spacing+offsetL rip(i1).extFoldMax(m)+(i-1)*spacing+offsetL],...
                [rip(i1).forcerMax(m) rip(i1).forcerMax(m)], 'Color', 'g')
            end
        end
        
    elseif(av.fx.split(i).direction < 0),%for going down
        i2=i2+1;
        plot(av.fx.split(i).ext+(i-2)*spacing+offsetL,av.fx.split(i).force,'-r');
                hold on
        if ~isempty(fold(i2).sizeMin)%pim changed
            %plot ( fold(i2).extFoldMin+(i-2)*spacing+offsetL,fold(i2).forcerMin, 'c.')
            %plot ( fold(i2).extUnfoldMin+(i-2)*spacing+offsetL, fold(i2).forcerMin, 'c.')
            plot ( fold(i2).extFoldMax+(i-2)*spacing+offsetL,fold(i2).forcerMax, 'b.')
            plot ( fold(i2).extUnfoldMax+(i-2)*spacing+offsetL, fold(i2).forcerMax, 'b.')
            %for m=1:length(fold(i2).extUnfoldMin)
            %plot ( [fold(i2).extUnfoldMin(m)+(i-2)*spacing+offsetL fold(i2).extFoldMin(m)+(i-2)*spacing+offsetL],...
            %    [fold(i2).forcerMin(m) fold(i2).forcerMin(m)], 'Color', 'c')
            %end
            for m=1:length(fold(i2).extUnfoldMax)
            plot ( [fold(i2).extUnfoldMax(m)+(i-2)*spacing+offsetL fold(i2).extFoldMax(m)+(i-2)*spacing+offsetL],...
                [fold(i2).forcerMax(m) fold(i2).forcerMax(m)], 'Color', 'b')
            end
        end
    end
end

xlabel('extension (nm)');
ylabel('force (pN)');
ylim([0 30]); %force limits

f.pdfFXfilename   = strcat(savedir,filesep,basefilename1,'_FX.pdf'); % LB change
saveas(5, f.pdfFXfilename); %LB change


disp('=>Done.');
