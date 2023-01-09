function [ ApodizationSpectrum ] = OCTFileGetApodizationSpectrum( handle )
% OCTFILEGETINTENSITY  Get the ApodizationSpectrum data from an .oct file.
%   data = OCTFILEGETINTENSITY( handle, dataName ) Get the ApodizationSpectrum data from an .oct file.
%
    ApodizationSpectrum = OCTFileGetRealData( handle, 'data\ApodizationSpectrum.data' );
end