function res_bool =CheckPowerofTwo(in)
if in > 0
    if (bitand(in,(in-1)) == 0)
        res_bool=true;
    else
        res_bool=false;
    end
else
    error('The input should be greater than 0');
end
end