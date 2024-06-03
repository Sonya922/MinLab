function convertWNtoIndex(Wavenumber, index_val)

% index_val = [920,2200,3300,3400,3800]';
index_array = "Index_" + index_val;
for num = 1: size(index_array,1)
% [index_name,index_return] = deal([],[]);
index_name = index_array(num);
index_return = find(abs(Wavenumber-index_val(num))<1);
if isempty(index_return)
    index_return = find(abs(Wavenumber-(index_val(num)-1))<1);
end
index_return = index_return(1);
assignin('base',index_name, index_return)
end

end