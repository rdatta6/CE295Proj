clear, clc
temp_data = readtable('San Juan Temperature Data .csv'); 
temp_data = table2timetable(temp_data); 


[Group, Hour] = findgroups(temp_data.DATE.Hour);
AvgValue = splitapply(@mean, temp_data(:, {'TEMP'}), Group);
tResult = table(Hour,AvgValue)

tResult = table2array(tResult)
tResult(1, 2) = 77

tResult = (tResult(:, 2) - 32).*(5/9)