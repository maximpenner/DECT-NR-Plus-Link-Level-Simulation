function [] = clf_all()

    % https://de.mathworks.com/matlabcentral/answers/478668-how-to-clear-not-close-all-the-opened-figures

    FigList = findall(groot, 'Type', 'figure');

    for iFig = 1:numel(FigList)
        try
            clf(FigList(iFig));
        catch
            % Nothing to do
        end
    end

end
