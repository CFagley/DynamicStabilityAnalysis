function varargout = Figure(varargin)
% Call by:
% Figure Example_Name
fig_name = varargin{1};
if isempty(findobj('name',fig_name))
    if nargout >0
        varargout{1} = figure('name',fig_name);
    else
        figure('name',fig_name)
    end
else
figure(findobj('name',fig_name))
hold all
end
hold all
