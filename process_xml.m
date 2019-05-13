function [avgFR, frametime] = process_xml(raw)

sep = strfind(raw,'<Frame relativeTime');
frameNum = numel(sep);
for i = 1:frameNum
    curr_idx = strfind(raw(sep(i):sep(i)+100),'"');
    toGet = sep(i)+curr_idx(1);
    if i == 1 %because first entry for "relativetime" is 0
        frametime(i) = str2num(raw(1,toGet));
    else
        frametime(i) = str2num(raw(1,toGet:toGet+8));
    end
end

exactPeriod = (frametime(end) - frametime(1))/frameNum;
avgFR = 1/exactPeriod;