function blkStruct = slblocks
  
  % Basic definition of the library
  Browser.Library = 'fod';
  Browser.Name    = 'FOMCON Toolbox';

  % Backwards compatibility
  blkStruct.Name = 'FOMCON Toolbox';
  blkStruct.OpenFcn = 'fod';
  blkStruct.MaskDisplay = '';
  
  % Output the structure
  blkStruct.Browser = Browser;