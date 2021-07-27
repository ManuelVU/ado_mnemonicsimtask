function F = setFigure(F, position, dual)
%SETFIGURE sets a figure in a way that copes with dual monitor changes
%    setFigure(F, position, dual)
%    F is figure handle
%    position is normalized [left bottom width height'
%    dual is 'upDown' or 'leftRight' depending on dual monitor arrangement
tmp = get(0, 'MonitorPositions');
if tmp(3) == 1920 && tmp(4) == 1080
    dual = 'single';
elseif tmp(3) == 1920
    dual = 'upDown';
elseif tmp(4) == 1080
    dual = 'leftRight';
else
    warning('not sure about monitors');
end
switch dual
    case 'single'
        vec(1) = (tmp(3)-tmp(1))*position(1);
        vec(2) = (tmp(4)-tmp(2))*position(2);
        vec(3) = (tmp(3)-tmp(1))*position(3);
        vec(4) = (tmp(4)-tmp(2))*position(4);
    case 'upDown'
        vec(1) = (tmp(3)-tmp(1))*position(1);
        vec(2) = (tmp(4)-tmp(2))*position(2)/2;
        vec(3) = (tmp(3)-tmp(1))*position(3);
        vec(4) = (tmp(4)-tmp(2))*position(4)/2;
    case 'leftRight'
        vec(1) = (tmp(3)-tmp(1))*position(1)/2;
        vec(2) = (tmp(4)-tmp(2))*position(2);
        vec(3) = (tmp(3)-tmp(1))*position(3)/2;
        vec(4) = (tmp(4)-tmp(2))*position(4);
end
set(F, ...
    'color'             , 'w'                , ...
    'units'             , 'pixels'           , ...
    'position'          ,  vec               , ...
    'papersize'         , [10 6.75]          , ...
    'paperpositionmode' , 'auto'             );

