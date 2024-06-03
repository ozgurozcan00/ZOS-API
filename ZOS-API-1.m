%Öncelikle zemax'ı MatLab üzerinden interaktif bir şekilde kontrol etmek
%için aşağıdaki kod bloğu çalıştırılır
%Interactive Extension tıklanıp çalıştırıldıktan sonra kod çalıştırılır
%Status: Connected olmalıdır

if ~exist('instance', 'var')
    instance = 0;
else
    try
        instance = int32(instance);
    catch
        instance = 0;
        warning('Invalid parameter {instance}');
    end


% Initialize the OpticStudio connection
TheApplication = InitConnection(instance);
if isempty(TheApplication)
    % failed to initialize a connection
    TheApplication = 'Failed to connect to OpticStudio';
else
    import ZOSAPI.*;
end
end

function app = InitConnection(instance)

import System.Reflection.*;

% Find the installed version of OpticStudio.

% This method assumes the helper dll is in the .m file directory.
% p = mfilename('fullpath');
% [path] = fileparts(p);
% p = strcat(path, '\', 'ZOSAPI_NetHelper.dll' );
% NET.addAssembly(p);

% This uses a hard-coded path to OpticStudio
NET.addAssembly('C:\Program Files\Zemax OpticStudio\ZOS-API\Libraries\ZOSAPI_NetHelper.dll');

success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize();
% Note -- uncomment the following line to use a custom initialization path
% success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize('C:\Program Files\OpticStudio\');
if success == 1
    LogMessage(strcat('Found OpticStudio at: ', char(ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory())));
else
    app = [];
    return;
end

% Now load the ZOS-API assemblies
NET.addAssembly(AssemblyName('ZOSAPI_Interfaces'));
NET.addAssembly(AssemblyName('ZOSAPI'));

% Create the initial connection class
TheConnection = ZOSAPI.ZOSAPI_Connection();

% Attempt to create a Standalone connection

% NOTE - if this fails with a message like 'Unable to load one or more of
% the requested types', it is usually caused by try to connect to a 32-bit
% version of OpticStudio from a 64-bit version of MATLAB (or vice-versa).
% This is an issue with how MATLAB interfaces with .NET, and the only
% current workaround is to use 32- or 64-bit versions of both applications.
app = TheConnection.ConnectAsExtension(instance);
if isempty(app)
   HandleError('Failed to connect to OpticStudio!');
end
if ~app.IsValidLicenseForAPI
	app.CloseApplication();
    HandleError('License check failed!');
    app = [];
end

end

function LogMessage(msg)
disp(msg);
end

function HandleError(error)
ME = MXException(error);
throw(ME);
end

%İnteraktif bir şekilde zemax'a uygulama yazmak üzere aşağıdaki
%Metodlar ve değişkenler çağırılır

%Matlab ortamında Zemax değişkenleri için bir sistem oluşturulur
TheSystem = TheApplication.PrimarySystem;
TheSystem.MakeSequential();


%SystemData metodu optik sistem için çağırılır
TheSystemData = TheSystem.SystemData;

%Açıklık tipi ve boyutu için aşağıdaki metodlar çağırılır
TheSystemData.Aperture.ApertureType=ZOSAPI.SystemData.ZemaxApertureType.EntrancePupilDiameter;
TheSystemData.Aperture.ApertureValue=10;

%Field Of View,dalga boyu ve materyal ataması için aşağıdaki metodlar ile ayarlanır
sysField = TheSystem.SystemData.Fields;
NewField_2=sysField.AddField(0,3.0,0);
 %Add 3 wavelengths: F,d,C
slPreset = TheSystemData.Wavelengths.SelectWavelengthPreset(ZOSAPI.SystemData.WavelengthPreset.FdC_Visible);
TheSystem.SystemData.MaterialCatalogs.AddCatalog("SCHOTT");


%Lens yaratmak üzere aşağıdaki yüzey metodları eklenir

TheLDE = TheSystem.LDE;
TheLDE.InsertNewSurfaceAt(1);
TheLDE.InsertNewSurfaceAt(1);
TheLDE.InsertNewSurfaceAt(1);

Surf_1 = TheLDE.GetSurfaceAt(1);
Surf_2 = TheLDE.GetSurfaceAt(2);
Surf_3 = TheLDE.GetSurfaceAt(3);
Surf_4 = TheLDE.GetSurfaceAt(4);
Surf_1.IsStop = true;

Surf_1.Thickness=5;
Surf_2.Thickness=5;
Surf_2.Radius=100;
Surf_2.Material="N-BK7";
Surf_3.Thickness=3;
Surf_3.Radius=-30;
Surf_3.Material="F2";
Surf_4.Radius=-80;

QuickFocus = TheSystem.Tools.OpenQuickFocus();
QuickFocus.RunAndWaitForCompletion();
QuickFocus.Close();

%Optimiziasyon
TheLDE = TheSystem.LDE;
tools = TheSystem.Tools;
tools.RemoveAllVariables();
Surface_Last = TheLDE.GetSurfaceAt(TheLDE.NumberOfSurfaces - 2);
Solver = Surface_Last.RadiusCell.CreateSolveType(ZOSAPI.Editors.SolveType.FNumber);
Solver.S_FNumber_.FNumber = 3.1415;
Surface_Last.RadiusCell.SetSolveData(Solver);

tools.SetAllRadiiVariable();            
Surface2 = TheLDE.GetSurfaceAt(2);
Surface3 = TheLDE.GetSurfaceAt(3);

Surface2.ThicknessCell.MakeSolveVariable();
Surface3.ThicknessCell.MakeSolveVariable();

TheMFE = TheSystem.MFE;
OptWizard = TheMFE.SEQOptimizationWizard;

OptWizard.Data = 1;
OptWizard.OverallWeight=1;
OptWizard.Ring=4;
OptWizard.IsGlassUsed= true;
OptWizard.GlassMin = 3.0;
OptWizard.GlassMax = 10.0;
OptWizard.GlassEdge = 3.0;
OptWizard.IsAirUsed = true;
OptWizard.AirMin = 0.5;
OptWizard.AirMax = 100.0;
OptWizard.AirEdge = 0.5;
OptWizard.Apply();

LocalOpt = TheSystem.Tools.OpenLocalOptimization();
if ~isempty(LocalOpt)
    LocalOpt.Algorithm = ZOSAPI.Tools.Optimization.OptimizationAlgorithm.DampedLeastSquares;
    LocalOpt.Cycles = ZOSAPI.Tools.Optimization.OptimizationCycles.Automatic;
    LocalOpt.NumberOfCores = 8;
    fprintf('Local Optimization...\n');
    fprintf('Initial Merit Function %6.3f\n', LocalOpt.InitialMeritFunction);
    LocalOpt.RunAndWaitForCompletion();
    fprintf('Final Merit Function %6.3f\n', LocalOpt.CurrentMeritFunction);
    LocalOpt.Close();
end

QuickFocus = TheSystem.Tools.OpenQuickFocus();
QuickFocus.RunAndWaitForCompletion();
QuickFocus.Close();

HammerOptimTimeInSeconds = 60;
HammerOpt = TheSystem.Tools.OpenHammerOptimization();

if ~isempty(HammerOpt)
    HammerOpt.Algorithm = ZOSAPI.Tools.Optimization.OptimizationAlgorithm.DampedLeastSquares;
    HammerOpt.NumberOfCores = 8;
    fprintf('Hammer Optimization for %i seconds...\n', HammerOptimTimeInSeconds);
    fprintf('Initial Merit Function %6.3f\n',HammerOpt.InitialMeritFunction);
    HammerOpt.RunAndWaitWithTimeout(HammerOptimTimeInSeconds);
    fprintf('Final Merit Function %6.3f\n', HammerOpt.CurrentMeritFunction);
    HammerOpt.Cancel();
    HammerOpt.WaitForCompletion();
    HammerOpt.Close();
end

QuickFocus = TheSystem.Tools.OpenQuickFocus();
QuickFocus.RunAndWaitForCompletion();
QuickFocus.Close();

%Analiz
TheRayFan = TheSystem.Analyses.New_Analysis(ZOSAPI.Analysis.AnalysisIDM.RayFan);
SD = TheSystem.Analyses.New_Analysis(ZOSAPI.Analysis.AnalysisIDM.SeidelDiagram);
FSD = TheSystem.Analyses.New_Analysis(ZOSAPI.Analysis.AnalysisIDM.FocalShiftDiagram);
Results = TheRayFan.GetResults();

newSpot = TheSystem.Analyses.New_StandardSpot();
fprintf('Spot has analysis specific settings? %i\n', newSpot.HasAnalysisSpecificSettings);
spotSet = newSpot.GetSettings();
spotSet.RayDensity = 15;
newSpot.ApplyAndWaitForCompletion();
analysis = TheSystem.Analyses.New_Analysis(ZOSAPI.Analysis.AnalysisIDM.ShadedModel);

newWin = TheSystem.Analyses.New_FftMtf(); 
newWin_Settings.MaximumFrequency = 50;
newWin_Settings.SampleSize = ZOSAPI.Analysis.SampleSizes.S_256x256;
newWin.ApplyAndWaitForCompletion();
newWin_Results = newWin.GetResults();

figure('position', [50, 150, 900, 600])
hold on;
grid on;

dataSeries = newWin_Results.DataSeries;
cc=lines(double(newWin_Results.NumberOfDataSeries));

for gridN=1:newWin_Results.NumberOfDataSeries
    data = dataSeries(gridN);
    y = data.YData.Data.double;
    x = data.XData.Data.double;
    plot(x,y(:,1),'-','color',cc(gridN,:));
    plot(x,y(:,2),':','color',cc(gridN,:));
end

title("FFT MTF");
xlabel('Spatial Frequency in cycles per mm');
ylabel('Modulus of the OTF');
legend('0^\circ tangential', '0^\circ sagittal', '14^\circ tangential', ...
    '14^\circ sagittal', '20^\circ tangential', '20^\circ sagittal');

