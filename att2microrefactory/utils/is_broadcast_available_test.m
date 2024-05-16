function brdcst = is_broadcast_available_test()
    % Explicitly declare and initialize the global variable
    global is_broadcast_available_test_var;
    
    % Initialize the variable to false by default
    is_broadcast_available_test_var = true;
    
    % Check if the variable has been previously set to true
    if ~isempty(is_broadcast_available_test_var)
        % Broadcast is available only from release 2016b
        vers   = version('-release');
        versn  = round(str2double(vers(1:4)));
        versl  = vers(5);
        if(versn>=2016 && versl=='b')
            is_broadcast_available_test_var = true;
        end
    end
    brdcst = is_broadcast_available_test_var;
end