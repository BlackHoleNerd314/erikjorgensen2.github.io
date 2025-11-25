function [] = Seconds2Gregorian(s0)



    s = floor(s0);

    % --- Constants ---
    SECS_PER_DAY = 86400; % (60 * 60 * 24)

    eon0 = 3652425*SECS_PER_DAY;
    total_seconds = mod(s,eon0); 
    epoch = floor(s/eon0);
    
    
    % --- 1. Separate Days and Time ---
    
    % Get the total number of full days
    total_days = floor(total_seconds / SECS_PER_DAY);
    
    % Get the remaining seconds within the current day
    seconds_today = rem(total_seconds, SECS_PER_DAY);

    % Calculate day of the week

    Weekday = mod(total_days, 7) + 1;
    
    % --- 2. Calculate Hour, Minute, Second ---
    
    H = floor(seconds_today / 3600);
    seconds_today = rem(seconds_today, 3600);
    
    Min = floor(seconds_today / 60);
    S = rem(seconds_today, 60);
    
    % --- 3. Calculate Year ---
    % We loop through each year, subtracting the number of days in that
    % year from our total, until we find the correct year.
    
    Y = 0; % Start at Year 0
    days_remaining = total_days;
    
    while true
        % Check if the *current* year (Y) is a leap year
        is_leap = (mod(Y, 4) == 0 && mod(Y, 100) ~= 0) || (mod(Y, 400) == 0);
        
        if is_leap
            days_in_this_year = 366;
        else
            days_in_this_year = 365;
        end
        
        % If we have more days remaining than are in this year,
        % subtract them and move to the next year.
        if days_remaining >= days_in_this_year
            days_remaining = days_remaining - days_in_this_year;
            Y = Y + 1;
        else
            % We found the correct year.
            break;
        end
    end
    
    % --- 4. Calculate Month and Day ---
    % 'days_remaining' is now the 0-indexed day of the year (0-364 or 0-365)
    
    % Define the days in each month
    days_in_months_normal = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    days_in_months_leap =   [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    
    % Check if our *final* year (Y) is a leap year to get month lengths
    is_leap = (mod(Y, 4) == 0 && mod(Y, 100) ~= 0) || (mod(Y, 400) == 0);
    
    if is_leap
        month_lengths = days_in_months_leap;
    else
        month_lengths = days_in_months_normal;
    end
    
    % Loop through the months to find the correct one
    M = 1; % Start at Month 1 (January)
    
    % 'D' will be the 1-indexed day of the year
    D_of_year = days_remaining + 1; 
    
    while D_of_year > month_lengths(M)
        D_of_year = D_of_year - month_lengths(M);
        M = M + 1;
    end
    
    % The remaining 'D_of_year' is our final day of the month
    D = D_of_year;

  % print to screen

  disp([epoch,Y,M,D,Weekday,H,Min,S]);

end
    
