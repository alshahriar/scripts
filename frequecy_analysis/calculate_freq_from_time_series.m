clear all
[file_to_process, test_dir] = uigetfile({'*.xlsx'},'Select Files to process','MultiSelect','on');
file_names = dir(fullfile(test_dir,file_to_process));
nFiles = length(file_names);

%% from user
row_no_header = 1;
data_start_from_line = 2;
row_no_of_units = 0;
trim_data_flag = 0;
row_no_of_description = 0;
samplingRate = 1/(abs(49.987534-49.862877)*1e-3);
samplingRate = 1/1e-4;


for iFile=1:nFiles
    file2read = fullfile(file_names(iFile).folder,file_names(iFile).name);

    sheets = sheetnames(file2read);
    for i=1:length(sheets)
        disp(string(i)+": "+string(sheets(i)))
    end
    sheetNo = input("which sheet to read: ");
    sheetRead  = sheets(sheetNo);
    opts = detectImportOptions(file2read,"Sheet",sheetRead);
    opts.VariableNamesRange = ['A',char(string(row_no_header))];
    opts.DataRange = ['A',char(string(data_start_from_line))];
    %opts.Delimiter = ",";
    opts.VariableNamingRule="preserve";
    %
    if(row_no_of_units~=0)
        opts.VariableUnitsLine = row_no_of_units;
    end    
    if(row_no_of_description~=0)
        opts.VariableDescriptionsLine = row_no_of_description;
    end
    %opts.DataLines = data_start_from_line;
    opts.Sheet = char(sheetRead);
    raw_data = readtimetable(file2read,opts,"SampleRate",samplingRate);
    %raw_data = readtable(file2read,opts);
    raw_data = rmmissing(raw_data);
    if(~isempty(raw_data))
        % raw_data.timeStamp = string(raw_data.timeStamp);
        raw_data_height = height(raw_data);
        raw_data_width = width(raw_data);
        %raw_data = add_units_from_csv(raw_data,file2read,row_no_of_units);
        if(trim_data_flag==1)
            var_selected = search_select_var_to_trim(raw_data);
            for i= 1:width(var_selected)
                figure(i)
                plot(raw_data.(var_selected{1,i}))
                title(var_selected{1,i})
                ax = gca;
                figure_axis_min_max_limit(ax)
            end

            trimS = input("starting row : ");
            trimE = input("ending row - last value ="+raw_data_height+" : ");

            for i= 1:width(var_selected)
                clf(figure(i))
            end

            trimmed_data = raw_data;
            trimmed_data([1:trimS, trimE:raw_data_height],:) = [];
        else
            trimmed_data = raw_data;
        end

        %data_stat = summarizeTableStats(trimmed_data);
        %data_stat_T = rows2vars(data_stat);
        %data_stat_T.OriginalVariableNames = [];
        %data_stat_T.Properties.RowNames = data_stat.Properties.VariableNames;
        %data_stat_T.Properties.VariableNames = data_stat.Properties.RowNames;
        

        %new_row = data_stat("mean",:);
        %new_row.Properties.RowNames = {char("File_"+iFile)};
        % if(iFile==1)
        %     data_mean = data_stat([],:);
        %     data_mean = [data_mean;new_row];
        % else
        %     data_mean = [data_mean;new_row];
        % end

        %sum_data = summary(trimmed_data);
        % writetable(trimmed_data,fullfile(file_names(iFile).folder,"processed_"+file_names(iFile).name),"WriteRowNames",true,"WriteVariableNames",true)

        % sheet_title = simplify_the_string(file_names(iFile).name);
        % sheet_title = string(sheet_title(1:25))+"_"+iFile;
        % writetable(data_stat,fullfile(file_names(iFile).folder,"sum_processed.xlsx"),"Sheet",sheet_title,"WriteRowNames",true)

        % if(0)
        %     writetable(cooling_duty,fullfile(file_names(iFile).folder,"cooling_duty_"+file_names(iFile).name),"WriteRowNames",true)
        %     cooling_duty = groupsummary(trimmed_data,"[CC]Cooling Mode (0=None; 1=Inverter; 2=Motor; 3=Both)");
        % end
    else
        disp("The following file is empty")
        disp(file_names(iFile).name)
    end
end
% write reading summary
% writetable(data_mean,fullfile(file_names(iFile).folder,"sum_processed.xlsx"),"Sheet","summary_mean","WriteRowNames",true)
% write_the_summary_of_data_files(file_names)

%% sensor connections
% if(0)
%     % EF230234
%     sensor_info = readtable("EF230234.csv");
% elseif(0)
%       %EF230180 is missing    
%     sensor_info = readtable("EF240436.csv");
% elseif(1)
%     sensor_info = readtable("EF240488.csv");
% end
% 
% x_axis_col = "Compressor 1 suction saturation temp"; % SST
% % x_axis_col = "Compressor 1 discharge temp (average of 4x sensors)"; % SDT
% for i=1:height(sensor_info)
%     if(string(sensor_info.channel{i})~="-")
%         figure(i)
%         plot(data_mean.(x_axis_col),data_mean.(string(sensor_info.channel{i})),"*")
%         xlabel(x_axis_col)
%         ylabel(string(sensor_info.location{i}));
%         folder2write = fullfile(file_names(iFile).folder,simplify_the_string(x_axis_col));
%         if(~exist(folder2write,"dir")) 
%             mkdir(folder2write);
%         end
%         saveas(gcf,fullfile(folder2write,[simplify_the_string(string(sensor_info.channel{i}))+".png"]))
%     end
% end


%% find the sensor column
% colIdxs = [];
% sensor = table();
% for i=1:height(sensor_info)
%     if(string(sensor_info.channel{i})~="-")
%         colIdx = find(strcmp(data_mean.Properties.VariableNames, string(sensor_info.channel{i})));
%         colIdxs = [colIdxs colIdx];
%         newCol = data_mean(:,colIdx);
%         newCol.Properties.VariableNames = string(sensor_info.location{i});
%         sensor = [sensor newCol];
%     else
%         sensor.(string(sensor_info.location{i})) = cell(height(sensor),1);
%     end
% end
% %% find other important column
% col_index = [91, 94, 95, 113, 228, 229, 235, 237, 245]; %EF234
% importantT = data_mean(:,col_index);
% sensor = [sensor importantT];
% writetable(sensor,fullfile(file_names(iFile).folder,"sum_processed.xlsx"),"Sheet","sensor_mean","WriteRowNames",true)
%% correlation coefficient
% if(0)
%     [R,PValue] = corrplot(sensor);
%     Rn = tril(table2array(R),-1);
%     %%
%     Rheatmap = standardizeMissing(Rn,0);
%     map = [1 0 0;
%         0.9 0.3 0.3;
%         0.9 0.6 0.6;
%         1 1 1;
%         0.6 0.6 0.9;
%         0.3 0.3 0.9;
%         0 0 1];
%     h = heatmap(Rheatmap,Colormap=map,ColorLimits=[-1 1]);
%     h.Title = "Correlation Coefficients";
%     h.YLabel = "First Variable";
%     h.XLabel = "Second Variable";
% end