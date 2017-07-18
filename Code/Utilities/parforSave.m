function parforSave( fName, accy, reg, predicted, true, model )

% A helper function that enables saving data from inside a parfor loop.
%
% Copyright, Henryk Blasinski 2017.

    save(fName,'accy','reg','predicted','true','model');

end

