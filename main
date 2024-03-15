%Input of the pn junction diode characteristics  data. Name/change variables as you see fit.
% Conversion from milliamps to amps for diodes F through J
FCurrent = millianpsF / 1000;
GCurrent = millianpsG / 1000;
HCurrent = millianpsH / 1000;
ICurrent = millianpsI / 1000;
JCurrent = millianpsJ / 1000;

% Assigning voltage variables for diodes F through J
FVoltage = VoltsF;
GVoltage = VoltsG;
HVoltage = VoltsH;
IVoltage = VoltsI/7;
JVoltage = VoltsJ/7;

% Filter out values smaller than 0.0001 amps for each diode and calculate the logarithm of the filtered current
Fvalid_indices = FCurrent > 0.0001;
Fcurrent_filtered = FCurrent(Fvalid_indices);
Fvoltage_filtered = FVoltage(Fvalid_indices);
Flog_current_filtered = log10(Fcurrent_filtered);

Gvalid_indices = GCurrent > 0.0001;
Gcurrent_filtered = GCurrent(Gvalid_indices);
Gvoltage_filtered = GVoltage(Gvalid_indices);
Glog_current_filtered = log10(Gcurrent_filtered);

Hvalid_indices = HCurrent > 0.0001;
Hcurrent_filtered = HCurrent(Hvalid_indices);
Hvoltage_filtered = HVoltage(Hvalid_indices);
Hlog_current_filtered = log10(Hcurrent_filtered);

Ivalid_indices = ICurrent > 0.0001;
Icurrent_filtered = ICurrent(Ivalid_indices);
Ivoltage_filtered = IVoltage(Ivalid_indices);
Ilog_current_filtered = log10(Icurrent_filtered);

Jvalid_indices = JCurrent > 0.0001;
Jcurrent_filtered = JCurrent(Jvalid_indices);
Jvoltage_filtered = JVoltage(Jvalid_indices);
Jlog_current_filtered = log10(Jcurrent_filtered);

% Constants
e = 1.602176634e-19; % Charge of an electron in Coulombs
k = 1.380649e-23; % Boltzmann's constant in J/K
T = 273.15 + 25; % Temperature in Kelvin, assuming 25 degrees Celsius

% Dataset initialization for diodes F through J
datasets = {'F', 'G', 'H', 'I', 'J'};
log_currents_filtered = {Flog_current_filtered, Glog_current_filtered, Hlog_current_filtered, Ilog_current_filtered, Jlog_current_filtered};
voltages_filtered = {Fvoltage_filtered, Gvoltage_filtered, Hvoltage_filtered, Ivoltage_filtered, Jvoltage_filtered};

% Prepare for plotting
figure;
hold on;

% Iterate through each dataset for diodes F through J
for datasetIndex = 1:length(datasets)
    x = voltages_filtered{datasetIndex};
    y = log_currents_filtered{datasetIndex};
    N = length(x);
    
    % Preallocate the slope array with NaN values
    slopes = NaN(N, 1);
    
    % Calculate slopes using finite differences
    for i = 2:N-1
        slopes(i) = (y(i+1) - y(i-1)) / (x(i+1) - x(i-1));
    end
    slopes(1) = (y(2) - y(1)) / (x(2) - x(1)); % First point
    slopes(N) = (y(N) - y(N-1)) / (x(N) - x(N-1)); % Last point

    %Next chunks of code can be copied/deleted depending on the amount of pn junction characteristics are being inputted. 
    %Adapt the log-current value (-2.49 in first instance) depending on the current density value for the pn junction diode in question. 
    if datasetIndex == 1
    % Find the index closest to log-current of -2.49
    [~, closestIndex] = min(abs(y + 2.49));
    
    % Slope and voltage at the closest index
    voltageAtMinus1 = x(closestIndex);
    slopeAtMinus1 = slopes(closestIndex);
    
    % Calculate ideality factor
    ideality_factor = (e*log10(exp(1)))/ (slopeAtMinus1 * k * T);
    
    % Plotting the filtered data
    plot(x, y, 'LineWidth', 2, 'DisplayName', sprintf('%s Data', datasets{datasetIndex}));
    
    % Plot the linear trend line
    xRange = linspace(min(x), max(x)+20, 100);
    yLine = slopeAtMinus1 * (xRange - voltageAtMinus1) + (-2.49);
    plot(xRange, yLine, '--', 'LineWidth', 1.5, 'DisplayName', sprintf('%s Trend Line', datasets{datasetIndex}));
    ylim([-4 -1])

    % Display calculated values
    fprintf('For plot %s at log-current = -2.49:\n', datasets{datasetIndex});
    fprintf('Voltage: %f V\n', voltageAtMinus1);
    fprintf('Slope: %f\n', slopeAtMinus1);
    fprintf('Ideality Factor: %f\n\n', ideality_factor);
    end

    if datasetIndex == 2
        % Find the index closest to log-current of -1.64
        [~, closestIndex] = min(abs(y + 1.64));
        
        % Slope and voltage at the closest index
        voltageAtMinus2 = x(closestIndex);
        slopeAtMinus2 = slopes(closestIndex);
        
        % Calculate ideality factor
        ideality_factor = (e*log10(exp(1)))/ (slopeAtMinus2 * k * T);
        
        % Plotting the filtered data
        plot(x, y, 'LineWidth', 2, 'DisplayName', sprintf('%s Data', datasets{datasetIndex}));
        
        % Plot the linear trend line
        xRange = linspace(min(x), max(x)+20, 100);
        yLine = slopeAtMinus2 * (xRange - voltageAtMinus2) + (-1.64);
        plot(xRange, yLine, '--', 'LineWidth', 1.5, 'DisplayName', sprintf('%s Trend Line', datasets{datasetIndex}));
        ylim([-4 -1])

        % Display calculated values
        fprintf('For plot %s at log-current = -1.64:\n', datasets{datasetIndex});
        fprintf('Voltage: %f V\n', voltageAtMinus2);
        fprintf('Slope: %f\n', slopeAtMinus2);
        fprintf('Ideality Factor: %f\n\n', ideality_factor);
    end

    if datasetIndex == 3
        % Find the index closest to log-current of -1.49
        [~, closestIndex] = min(abs(y + 1.49));
        
        % Slope and voltage at the closest index
        voltageAtMinus3 = x(closestIndex);
        slopeAtMinus3 = slopes(closestIndex);
        
        % Calculate ideality factor
        ideality_factor = (e*log10(exp(1)))/ (slopeAtMinus3 * k * T);
        
        % Plotting the filtered data
        plot(x, y, 'LineWidth', 2, 'DisplayName', sprintf('%s Data', datasets{datasetIndex}));
        
        % Plot the linear trend line
        xRange = linspace(min(x), max(x)+20, 100);
        yLine = slopeAtMinus3 * (xRange - voltageAtMinus3) + (-1.49);
        plot(xRange, yLine, '--', 'LineWidth', 1.5, 'DisplayName', sprintf('%s Trend Line', datasets{datasetIndex}));
        ylim([-4 -1])

        % Display calculated values
        fprintf('For plot %s at log-current = -1.49:\n', datasets{datasetIndex});
        fprintf('Voltage: %f V\n', voltageAtMinus3);
        fprintf('Slope: %f\n', slopeAtMinus3);
        fprintf('Ideality Factor: %f\n\n', ideality_factor);
    end
    if datasetIndex == 4
        % Find the index closest to log-current of -2.64
        [~, closestIndex] = min(abs(y + 2.64));
        
        % Slope and voltage at the closest index
        voltageAtMinus4 = x(closestIndex);
        slopeAtMinus4 = slopes(closestIndex);
        
        % Calculate ideality factor
        ideality_factor = (e*log10(exp(1)))/ (slopeAtMinus4 * k * T);
        
        % Plotting the filtered data
        plot(x, y, 'LineWidth', 2, 'DisplayName', sprintf('%s Data', datasets{datasetIndex}));
        
        % Plot the linear trend line
        xRange = linspace(min(x), max(x)+20, 100);
        yLine = slopeAtMinus4 * (xRange - voltageAtMinus4) + (-2.64);
        plot(xRange, yLine, '--', 'LineWidth', 1.5, 'DisplayName', sprintf('%s Trend Line', datasets{datasetIndex}));
        ylim([-4 -1])

        % Display calculated values
        fprintf('For plot %s at log-current = -2.64:\n', datasets{datasetIndex});
        fprintf('Voltage: %f V\n', voltageAtMinus4);
        fprintf('Slope: %f\n', slopeAtMinus4);
        fprintf('Ideality Factor: %f\n\n', ideality_factor);
    end

    if datasetIndex == 5
        % Find the index closest to log-current of -2.64
        [~, closestIndex] = min(abs(y + 2.64));
        
        % Slope and voltage at the closest index
        voltageAtMinus5 = x(closestIndex);
        slopeAtMinus5 = slopes(closestIndex);
        
        % Calculate ideality factor
        ideality_factor = (e*log10(exp(1)))/ (slopeAtMinus5 * k * T);
        
        % Plotting the filtered data
        plot(x, y, 'LineWidth', 2, 'DisplayName', sprintf('%s Data', datasets{datasetIndex}));
        
        % Plot the linear trend line
        xRange = linspace(min(x), max(x)+20, 100);
        yLine = slopeAtMinus5 * (xRange - voltageAtMinus5) + (-2.64);
        plot(xRange, yLine, '--', 'LineWidth', 1.5, 'DisplayName', sprintf('%s Trend Line', datasets{datasetIndex}));
        ylim([-4 -1])

        % Display calculated values
        fprintf('For plot %s at log-current = -2.64:\n', datasets{datasetIndex});
        fprintf('Voltage: %f V\n', voltageAtMinus5);
        fprintf('Slope: %f\n', slopeAtMinus5);
        fprintf('Ideality Factor: %f\n\n', ideality_factor);
    end
end

hold off;
xlabel('Voltage (V)');
ylabel('Log_{10}(Current (A))');
title('Diode Characteristics: Log_{10}(Current) vs Voltage');
legend('show');
grid on;
