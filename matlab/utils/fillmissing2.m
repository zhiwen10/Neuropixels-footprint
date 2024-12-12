function [F,TF] = fillmissing2(A,method,varargin)
%FILLMISSING2   Fill missing entries in two dimensions
%    First argument must be a numeric matrix.
%    F = FILLMISSING2(A,INTERP) fills NaN entries in A using the
%    interpolation method specified by INTERP, which must be:
%        "nearest"  - Nearest (in terms of index) non-missing entry.
%        "linear"   - Linear interpolation of non-missing entries.
%        "natural"  - Natural neighbor interpolation of non-missing entries.
%        "cubic"    - Cubic interpolation of non-missing entries.
%        "v4"       - Biharmonic spline interpolation of non-missing entries.
%    Filling using interpolation methods is supported only for double and
%    single type data.
%
%    F = FILLMISSING2(A,MOV,K) fills NaN entries using a centered, 2D
%    moving window formed from the neighboring non-missing entries. K
%    specifies the window length in both the X and Y directions and must
%    be a positive scalar. MOV specifies the moving window method, which
%    must be:
%      "movmean"   - Moving mean of neighboring non-missing entries.
%      "movmedian" - Moving median of neighboring non-missing entries.
%
%    F = FILLMISSING2(A,MOV,{KX,KY}) uses a moving window defined by KX,
%    the window length in the X (rows) direction and KY, the window
%    length in the Y (columns) direction. KX and KY may either be positive
%    scalar values, or two-element vectors with the form: [NB NF], where NB
%    indicates the number of elements before the missing value that are
%    included in the moving window, and NF is the number ahead of the 
%    missing value that are included.
%
%    Optional arguments:
%
%   F = FILLMISSING2(A,METHOD,...,MissingLocations=M) specifies the
%   missing data locations according to the logical array M. Elements of M
%   that are true indicate missing data in the corresponding element of A.
%   M must be the same size as A.
%
%   F = FILLMISSING2(A,METHOD,...,SamplePoints={X,Y}) also specifies the
%   sample points X and Y used by the fill method. X contains the sample 
%   points in the X-direction (each sample point corresponds to a row) 
%   and Y contains the sample points in the Y-direction (each sample point
%   corresponds to a column). X and Y must be floating-point, duration, or 
%   datetime vectors. Each sample points vector must be sorted and cannot 
%   contain duplicate values. By default, FILLMISSING2 uses data sampled
%   uniformly at points X = [1 2 3 ... ], Y = [1 2 3 ... ].
%
%    [F,TF] = FILLMISSING2(A,...) also returns a logical array TF
%    indicating which entries of A were filled to make F. TF is the same
%    size as F.
%
%   Examples:
%
%    % Fill missing entries using 2D linear interpolation:
%      A = [ 17   NaN     1     8    15
%           NaN     5     7    14    16
%             4     6    13    20   NaN
%            10    12    19   NaN     3
%            11    18   NaN     2     9]
%      B = fillmissing2(A,"linear")
%
%    % Fill missing entries using a 2D moving mean, with a square, 2x2
%    % window centered on each query point:
%      A = [ 17    24   NaN     8    15
%            23     5     7    14    16
%           NaN     6    13    20    22
%            10    12    19    21   NaN
%            11    18    25   NaN     9]
%      B = fillmissing2(A,"movmean",2)
%  
%    % Fill missing entries using a 2D moving median, with a rectangular
%    % window that includes two points to the left of the query point, one
%    % point to the right, three points above the query point and no points
%    % below:
%      A = magic(6);
%      A([1 10 12 23 30]) = NaN
%      B = fillmissing2(A,"movmedian",{[2,1],[3,0]})
%
%    % Fill missing entries using 2D cubic interpolation with sample points:
%      x = 2:2:20;
%      y = 0:0.5:2;
%      A = sin(x)'*y;
%      B = fillmissing2(A,"cubic",SamplePoints={x,y})
%
%    % Fill missing entries with a moving median using a 3x3 window, with 
%    % missing locations supplied:
%      A = randi(10,[5,5],'int32')
%      ma = A > 8
%      B = fillmissing2(A,"movmedian",3,MissingLocations=ma)
%
%   See also ISMISSING, STANDARDIZEMISSING, RMMISSING, ISNAN, ISNAT
%            ISOUTLIER, FILLMISSING, FILLOUTLIERS, RMOUTLIERS, SMOOTHDATA

%   Copyright 2022-2023 The MathWorks, Inc.

arguments
    A
    method {mustBeTextScalar}
end

arguments (Repeating)
    varargin
end
% Validate data
if ~(isnumeric(A) && ismatrix(A))
    error(message("MATLAB:fillmissing2:InvalidDataShapeOrType"));
end

% Validate method
validMethods = ["movmean","movmedian","nearest","linear","natural", ...
    "cubic","v4"];
indMethod = matlab.internal.math.checkInputName(method,validMethods);
if sum(indMethod) ~= 1
    % Also catch ambiguities for fillmissing2(A,"n")
    error(message("MATLAB:fillmissing2:MethodInvalid"));
end
method = validMethods{indMethod};

if ~(indMethod(1) || indMethod(2)) % method is not "movmean" or "movmedian"
    [F,TF] = fillmissing2Interp(A,method,varargin{:});
else
    [F,TF] = fillmissing2Mov(A,method,varargin{:});
end
if issparse(A)
    F = sparse(F);
    TF = sparse(TF);
end
end

%--------------------------------------------------------------------------

function [F,TF] = fillmissing2Interp(A,method,options)
arguments
    A
    method
    options.MissingLocations {mustBeA(options.MissingLocations,"logical"), ...
        validateMissingLocations(A,options.MissingLocations)}
    options.SamplePoints {mustBeA(options.SamplePoints,"cell"), ...
        matlab.internal.math.validate2DSamplePoints(A,options.SamplePoints)}
end

if ~isfloat(A)
    error(message("MATLAB:fillmissing2:InterpOnlyFloats"));
end

if isfield(options,"MissingLocations")
    ma = options.MissingLocations;
    allMissingLocations = ma | ismissing(A);
else
    ma = ismissing(A);
    allMissingLocations = ma;
end

spMustBeDouble = true;
sp = generateOrConvertSP(options,size(A),spMustBeDouble);

F = A;
TF = false(size(A));

% Quick return when there are no missing values or empties
if ~any(ma,"all") || isempty(A)
    return
end

[X,Y] = find(~allMissingLocations);
if isempty(X)
    % Quick return when all values are missing
    return
end
X = sp{1}(X);
Y = sp{2}(Y);
[xq,yq] = find(ma);
xq = sp{1}(xq);
yq = sp{2}(yq);

fillVals = griddata(X,Y,double(A(~allMissingLocations)),xq,yq,method);
if ~isempty(fillVals)
    F(ma) = fillVals;
    TF(ma) = ~ismissing(fillVals);
end
end

%--------------------------------------------------------------------------

function [F,TF] = fillmissing2Mov(A,method,window,options)
arguments
    A
    method
    window {matlab.internal.math.validate2DWindow}
    options.MissingLocations (:,:) {mustBeA(options.MissingLocations,"logical"), ...
        validateMissingLocations(A,options.MissingLocations)}
    options.SamplePoints {mustBeA(options.SamplePoints,"cell"), ...
        matlab.internal.math.validate2DSamplePoints(A,options.SamplePoints,window)}
end

if isfield(options,"MissingLocations")
    ma = options.MissingLocations;
else
    ma = ismissing(A);
end
if ~isfield(options,"SamplePoints") && isduration(window)
    error(message("MATLAB:gridded2DData:SPForTimeWindow"));
end

spMustBeDouble = false;
sp = generateOrConvertSP(options,size(A),spMustBeDouble);

if ~iscell(window)
    window = {window,window};
end

% The type of the sample points will dictate the type of the window
[window{1}, inclusiveUpperBoundX] = matlab.internal.math.convertAndSplitWindow(window{1},sp{1});
[window{2}, inclusiveUpperBoundY] = matlab.internal.math.convertAndSplitWindow(window{2},sp{2});

if strcmp(method,"movmean")
    [F,TF] = matlab.internal.math.fillmovmean2(full(A),window{1},window{2},inclusiveUpperBoundX,inclusiveUpperBoundY,sp{1},sp{2},full(ma));
else
    [F,TF] = matlab.internal.math.fillmovmedian2(full(A),window{1},window{2},inclusiveUpperBoundX,inclusiveUpperBoundY,sp{1},sp{2},full(ma));
end
end

%--------------------------------------------------------------------------

function validateMissingLocations(A,ma)
if ~isequal(size(ma),size(A))
    error(message("MATLAB:fillmissing2:MissingLocationsInvalid"));
end
end

%--------------------------------------------------------------------------

function sp = generateOrConvertSP(options,sizeA,spMustBeDouble)
if isfield(options,"SamplePoints")
    sp = options.SamplePoints;
    sp{1} = matlab.internal.math.convertSamplePoints(sp{1},spMustBeDouble);
    sp{2} = matlab.internal.math.convertSamplePoints(sp{2},spMustBeDouble);
else
    sp = {1:sizeA(1),1:sizeA(2)};
end
end
