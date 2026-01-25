function [] = annotation(x, y, str)    
    dim = [x, y, 0.1, 0.1];
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', LineStyle='none', Interpreter='none', FontName='FixedWidth');
end