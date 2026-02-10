% reading ansys surface data
clear all
directory = "D:\TGS380_CFD_transient\TGS380_028";
fileIdentifier = "pressure";
dt = 1.0345E-05; % each iterations
colNoOfPressure = 4;

read_write_flag = input("read the raw data or load already saved data? 1 for read or 0 for load: ");

fileInfo = dir(fullfile(directory,[fileIdentifier+"*.csv"]));
nFiles = length(fileInfo);
fileName = fullfile(fileInfo(1).folder,fileInfo(1).name);
info = cfx_getBlockInfo(fileName);
nNodes = info.node.dataEndLine-info.node.dataStartLine+1;
time = zeros(nFiles,1);
iter = zeros(nFiles,1);
pressure = zeros(nFiles,nNodes);
S = cfx_readWithReadmatrix(fileName,info);
if(read_write_flag==1)
    for iFile = 1:nFiles
        fileName = fullfile(fileInfo(iFile).folder,fileInfo(iFile).name);
        fprintf("reading = %s, %d of %d\n",fileInfo(iFile).name,iFile,nFiles);
        P = cfx_readWithReadmatrix_onlyPressure(fileName, info, colNoOfPressure);
        a = split(fileInfo(iFile).name,"_");
        itert = a{2}(1:end-4);
        iter(iFile) = str2num(itert);
        time(iFile) = iter(iFile)*dt;
        pressure(iFile,:) = P.nodes.pressure_kg_m_1_s_2;
    end
    save(fullfile(directory,fileIdentifier),"pressure","time","S")
    for iFile = 1:nFiles
        p = pressure(iFile,:)';
        data2write = [S.nodes.x_m,S.nodes.y_m,S.nodes.z_m,p];
        fname = sprintf("m_%s_%07d.csv",fileIdentifier,iter(iFile));
        disp("writing "+fname)
        fullname = fullfile(directory,fname);
        writematrix(data2write,fullname)
    end

else
    %%
    load(fullfile(directory,fileIdentifier),"pressure","time","S")
end
TR = triangulation(S.faces,[S.nodes.x_m,S.nodes.y_m,S.nodes.z_m]);
mean_P = mean(mean(pressure));
avg_P_surf = mean(pressure);
if(read_write_flag==1)
    p = avg_P_surf';
    data2write = [S.nodes.x_m,S.nodes.y_m,S.nodes.z_m,p];
    fname = sprintf("m_%s_avg.csv",fileIdentifier);
    disp("writing "+fname)
    fullname = fullfile(directory,fname);
    writematrix(data2write,fullname)
end


%% 
figure; ax = axes; hold(ax,'on');
trisurf(TR,avg_P_surf,"edgeColor","none");
axis equal; view(-90,180);
camlight; lighting gouraud; 
colormap("jet");
colorbar;
title('Avg Pressure');

%%
[F,Fx,Fy,Fz] = pressureForcesTriSurface(S, pressure, ...
    'PressureLocation','vertex', ...
    'ChunkSize', 50, ...
    'ForceSign', +1);
plot(time,Fx,time,Fy,time,Fz); legend({"Fx","F_{axial}","Fy"})
xlabel("time [s]")
ylabel("Force")
figure;
surfAvgP = mean(pressure,2);
plot(time,surfAvgP,time,pressure(:,1)',time,pressure(:,100)',time,pressure(:,500)',time,pressure(:,1000)',time,pressure(:,5000)',time,pressure(:,10000)'); hold on;
%%
time = time-time(1);
dt_solid = 1e-4;
time_new = 0:dt_solid:3e-3;
[Fx_new, Fy_new, Fz_new, time_new] = interpForcesToNewTime(time, Fx, Fy, Fz, time_new);
figure;
plot(time,Fx,time,Fy,time,Fz); legend({"Fx","F_{axial}","Fy"})
hold on
plot(time_new,Fx_new,"-*",time_new,Fy_new,"-*",time_new,Fz_new,"-*"); legend({"Fx","F_{axial}","Fy"})
xlabel("time [s]")
ylabel("Force")


%%
if(0)
    % Two-view, side-by-side animation (same TR / same pressure), written to one GIF

    %
    gifFile = makePressureGif(TR, pressure, time, mean_P, directory, ...
        'GifName','animation1.gif', ...
        'StartIndex',1, ...
        'EndIndex',nFiles, ...
        'View',[0 180], ...
        'UseCameraSpin',true);

    gifFile = makePressureGif(TR, pressure, time, mean_P, directory, ...
        'GifName','animation3.gif', ...
        'StartIndex',1, ...
        'EndIndex',nFiles, ...
        'View',[0 -90], ...
        'UseCameraSpin',false);

    gifFile = makePressureGifTwoViews(TR, pressure, time, mean_P, directory, ...
        'GifName','animation11.gif', ...
        'StartIndex',1, ...
        'EndIndex',nFiles, ...
        'ViewLeft',[0 -90], ...
        'ViewRight',[-180 90]);


end
%% Frequency
Prandom1 = pressure(:,1);
Prandom2 = pressure(:,100);
Prandom3 = pressure(:,1000);
Prandom4 = pressure(:,10000);
Prandom5 = pressure(:,20000);
Prandom6 = pressure(:,2000);
Prandom7 = pressure(:,200);
Prandom8 = pressure(:,20);
tData = timetable(Fx,Fy,Fz,Prandom1,Prandom2,Prandom3,Prandom4,Prandom5,Prandom6,Prandom7,Prandom8,'RowTimes',seconds(time));
tDataAll = timetable(pressure,'RowTimes',seconds(time));





%% separating the faces
[faceLabel, faceGroups, boundaryEdges] = splitSurfacesByAngle(TR, 20);
fprintf('Found %d surface patches.\n', numel(faceGroups));
F = TR.ConnectivityList;
P = TR.Points;

%% original list of patches

figure; ax = axes; hold(ax,'on');
h = patch('Faces',F,'Vertices',P, ...
          'FaceVertexCData', faceLabel, ... % color per vertex is expected
          'FaceColor','flat','EdgeColor','none');
axis equal; view(3); % camlight; lighting gouraud; colorbar;
title('Surface patches (by dihedral angle)');

%% merge some faces together
innner_core_patches = 33; % user
ring_core_patches = [34:37,38,39:54]; % user
plate_patches = [1:32]; %user
mergeGroups = { ring_core_patches,plate_patches };
faceLabel2  = mergePatchGroups(faceLabel, mergeGroups, 'Renumber', true);

figure; ax = axes; hold(ax,'on');
h = patch('Faces',F,'Vertices',P, ...
          'FaceVertexCData', faceLabel2, ... % color per vertex is expected
          'FaceColor','flat','EdgeColor','none');
axis equal; view(3);
camlight; lighting gouraud; 
colorbar;
title('Surface patches (by dihedral angle)');

patchIDs =  1;
F = TR.ConnectivityList;
P = TR.Points;
faceKeep = ismember(faceLabel2, patchIDs);
Fsub = F(faceKeep,:);
figure;
patch('Faces', Fsub, 'Vertices', P, ...
      'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
axis equal; view(3); camlight headlight; lighting gouraud;
title(sprintf('Patches: %s', mat2str(patchIDs)));

patchIDs =  2;
F = TR.ConnectivityList;
P = TR.Points;
faceKeep = ismember(faceLabel2, patchIDs);
Fsub = F(faceKeep,:);
figure;
patch('Faces', Fsub, 'Vertices', P, ...
      'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
axis equal; view(3); camlight headlight; lighting gouraud;
title(sprintf('Patches: %s', mat2str(patchIDs)));


patchIDs =  3;
F = TR.ConnectivityList;
P = TR.Points;
faceKeep = ismember(faceLabel2, patchIDs);
Fsub = F(faceKeep,:);
figure;
patch('Faces', Fsub, 'Vertices', P, ...
      'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
axis equal; view(3); camlight headlight; lighting gouraud;
title(sprintf('Patches: %s', mat2str(patchIDs)));
%% hide specific faces
% patchesToHide = [[1:10],[11:30]];
% TRvis = hidePatches(TR, faceLabel, patchesToHide);
% % Replot (you can still color by faceLabel(keep) if you want)
% keep = ~ismember(faceLabel, patchesToHide);
% figure;
% plotSurfacePatches(TRvis, faceLabel(keep), 'LabelPatches', false);
% title('Hidden patches removed');
%% show specific faces

patchIDs =  [10]; % user
F = TR.ConnectivityList;
P = TR.Points;
faceKeep = ismember(faceLabel, patchIDs);
Fsub = F(faceKeep,:);
figure;
patch('Faces', Fsub, 'Vertices', P, ...
      'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
axis equal; view(3); camlight headlight; lighting gouraud;
title(sprintf('Patches: %s', mat2str(patchIDs)));

%%

[centerMean, centerArea] = patchCenterXYZ(TR, faceLabel, 10);

fprintf('Mean-of-centroids center:  (%.6g, %.6g, %.6g)\n', centerMean);
fprintf('Area-weighted center:      (%.6g, %.6g, %.6g)\n', centerArea);

% Optional: show it in the plot
hold on;
plot3(centerArea(1), centerArea(2), centerArea(3), 'ko', 'MarkerFaceColor','k');
text(centerArea(1), centerArea(2), centerArea(3), '  patch 3 center', 'Color','k');

%%

patch_with_issue = 10; % user
% get the coordinate from the previous section
[faceLabel3, newIDs,splitInfo] = splitPatchByRadiusXZ(TR, faceLabel, patch_with_issue, 0.0140, ...
    'CenterXZ', [0.0873354 0.0931166], ...
    'Method', 'centroid');  % or 'vertexFraction'

insidePatchID  = newIDs.InsideID;
outsidePatchID = newIDs.OutsideID;

disp([insidePatchID outsidePatchID])

% Plot the updated patches
figure; ax = axes; hold(ax,'on');
h = patch('Faces',F,'Vertices',P, ...
          'FaceVertexCData', faceLabel3, ... % color per vertex is expected
          'FaceColor','flat','EdgeColor','none');
axis equal; view(3); % camlight; lighting gouraud; colorbar;
title('Surface patches (by dihedral angle)');


%%
patchIDs =  [55]; % new patch % user
F = TR.ConnectivityList;
P = TR.Points;
faceKeep = ismember(faceLabel3, patchIDs);
Fsub = F(faceKeep,:);
figure;
patch('Faces', Fsub, 'Vertices', P, ...
      'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
axis equal; view(3); camlight headlight; lighting gouraud;
title(sprintf('Patches: %s', mat2str(patchIDs)));


%%
innner_core_patches = 33; %user
ring_core_patches = [34:37,38,39:55]; %user
plate_patches = [1:32]; %user
mergeGroups = { ring_core_patches,plate_patches};
faceLabel_merged  = mergePatchGroups(faceLabel3, mergeGroups, 'Renumber', true);

figure; ax = axes; hold(ax,'on');
h = patch('Faces',F,'Vertices',P, ...
          'FaceVertexCData', faceLabel_merged, ... % color per vertex is expected
          'FaceColor','flat','EdgeColor','none');
axis equal; view(3);
camlight; lighting gouraud; 
colorbar;
title('Surface patches (by dihedral angle)');

% Plot all patches (unique labels) in separate figures
plotPatchesSeparateFigures(TR, faceLabel_merged);


%% calculating forces for each patches

% out = pressureForcesByPatches(S, faceLabel_merged, pressure, ...
%     'PressureLocation','vertex', ...
%     'ForceSign', +1);

O   = [ 0.08734  -0.23179   0.09312 ];
Zp  = [ 0.08734   1.00000   0.09312 ];
XZp = [ 0.00000  -0.23179   0.09312 ];

% ptest = ones(1, size(pressure,2));

out = pressureForcesByPatches(S, faceLabel_merged, pressure, ...
    'PressureLocation','vertex', ...
    'ForceSign', +1, ...
    'ComputeTorque', true, ...
    'TorqueReference', O);

[R, O, ex, ey, ez] = coordFrameFromAxisPoints(O, Zp, XZp);
plotTRwithFrameVectors(TR, O, ex, ey, ez, 'FaceAlpha', 0.25, 'ShowEdges', false);
out = convertForcesnTorqueToLocalCoor(out,R);

Fx_sum = sum(out.F_aboutX, 2);   % Nt x 1
Fy_sum = sum(out.F_aboutY, 2);
Fz_sum = sum(out.F_aboutZ, 2);

maxErr = max(vecnorm([Fx_sum - Fx;Fy_sum - Fy;Fz_sum - Fz], 2, 2));
fprintf('Max ||Ftot - sum(patches)|| over time = %.6g\n', maxErr);

figure; hold on; grid on;
plot(time, Fx_sum, 'LineWidth', 1.5);
plot(time, Fy_sum, 'LineWidth', 1.5);
plot(time, Fz_sum, 'LineWidth', 1.5);
xlabel('Time, t');
ylabel('Force');
title('Sum of forces across all patches vs time - transformed');
legend('Sum Fx','Sum Fy','Sum Fz','Location','best');

%%
% out = pressureForcesByPatches(...);  % computed earlier
plotPatchForcesOverTime(time, out);  % Fx, Fy, Fz, |F| for all patches

% Only Fx and Fz, and only patches 2,5,9:
% plotPatchForcesOverTime(t, out, 'Components',{'Fx','Fz'}, 'PatchIDs',[2 5 9]);
% Separate figures per component:
% plotPatchForcesOverTime(t, out, 'SeparateFigures',true);
% plotPatchForcesOverTime(t, out, 'Components',{'Fx','Fz'});

%% 
Tx_sum = sum(out.T_aboutX, 2);   % Nt x 1
Ty_sum = sum(out.T_aboutY, 2);
Tz_sum = sum(out.T_aboutZ, 2);

figure; hold on; grid on;
plot(time, Tx_sum, 'LineWidth', 1.5);
plot(time, Ty_sum, 'LineWidth', 1.5);
plot(time, Tz_sum, 'LineWidth', 1.5);
xlabel('Time, t');
ylabel('Torque');
title('Sum of forces across all patches vs time - Local');
legend('Sum Tx','Sum Ty','Sum Tz')
%%

% out = pressureForcesByPatches(...);   % already computed
% t is your Nt x 1 time vector

patchNames = ["plate", "inner", "ring"];  % must match numel(out.PatchIDs)
writePatchForcesCSV(fullfile(fileInfo(1).folder,"all_forces.csv"), time, out, patchNames);

%% Avg surf Pressure
