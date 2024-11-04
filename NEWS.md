# dynatopGIS 0.3.0.1011

- Output models in the format of v0.3 of `dynatop`
- Adjust workflow to take catchment map outline as initial input
- Adjust convert_channel to allow for merging of waterbodies and network simplification


# dynatopGIS 0.2.4

- Switch (with suboptimal code) to `terra` package over `raster` to fix
  handling of crs in meta data
- Removal of rgdal from suggests

# dynatopGIS 0.2.3

- Calls to rgdal / rgeos replaced by calls to raster

# dynatopGIS 0.2.2

- fixed bug in computation of channel area in the model
- fixed bug in reloading meta data when these was a classification
- added run-off attenuation parameters for hillslope HRU in create_model method

# dynatopGIS 0.2.1

- Tidying up of vignettes and package imports for CRAN

# dynatopGIS 0.2.0.9101

## Breaking changes

- Previous classify method broken into two parts
  - `classify` to cut a layer into classes
  - `combine_classes` to combine classes using cantor pairing and burns
  
- Output of model altered to reflect changes in `dynatop` v0.2.0.9101

# dynatopGIS 0.2.0.9030

## Breaking changes

- Incorporated rainfall and pet weight into the model output

# dynatopGIS 0.2.0.9030

## Breaking changes

- New format for flow linkages in a separate table

## other changes

- Added ability to compute weights for gridded inputs

# dynatopGIS 0.2.0.9001

## Breaking changes

- Altered model structure output to match change to surface zone celerity from
  time constant in dynatop v0.2.0.9001

# dynatopGIS 0.2.0.9000

## Breaking Changes

- Altered analysis of catchment to provide data for the bottom cross section
  of the HSU. This effects:
    - Derivation of flow directions and fractions for the HSU
	- Estimation of the width is now based on contour length

## Other changes

- More efficient implementation of property computation using few passes of the
  data
  
# dynatop 0.11

## Other Changes

- Data sets for vignette split off into dynatopData package
- Under the hood changes to package build 

# dynatopGIS 0.1

## Breaking Changes

- Code base reformulated in an Object orientated form using the `R6`
  package. Except for change below algorithms as for v0.0.4
- Output format of model changed to reflect `dynatop` v0.1

## New Features

- Reformulation allows for data and code to saved in a single object allowing full
  reproducibility
- Additional plotting functions
- Model creation algorithm changes to:
    - take sequence for most downslope band (rather then upslope)
	- Adjusts sequence rather then merging or dropping areas (still
  experimental) to ensure connectivity

# dynatopGIS 0.0.4

## Breaking Changes

- Uses a single multilayer raster file rather then a folder of raster files
- New algorithms will not produce identical results to v0.0.3	

## New Features

- New sink fill algorithm which should use fewer passes of the DEM
	- Fills from the sink with lowest possible value first.
	- Drops non-connected cells on the edge of the catchment
- Rewritten computation of HSU properties based on Quinn et al. 1995 and
  implemented in two passes of the DEM

# dynatopGIS 0.0.3

## Bug fixes

- Various minor bug

# dynatopGIS 0.0.2

## New features
- Added ability to continue a run of sink_fill by starting with previous
  filled dem 
  
## Bug fixes
- Model generation - missing average gradient property added and class of data frame columns altered
- Linking Rcpp function - alter NAMESPACE so should be included

## Other changes
- Vignette know includes running of model and burning in a class

# dynatopGIS 0.0.1

## Other changes
- Removed vignette data from package to reduce size. Now on github
- Minor spelling and format changes to vignette

# dynatopGIS 0.0.0.9000

## Context
This package is the result of an almost complete rewrite of the dynatopmodel package
formally on CRAN and the associated development code (not in the public
domain).

This package contains the code for extracting the GIS summaries required for
dynamic TOPMODEL and helpers for inserting default parameters. The package [dynatop]{https://waternumbers.github.io/dynatop/} contains the
code for performing simulations.

## New Features
- Change in approach to writing all intermediate files to
      - produce a record of analysis
	  - allow for future developments with larger size GIS data
- Porting of sink_fill and atb calculations from original 'topmodel' CRAN package
  to Rcpp
- Improved method of determining river network intersection with DEM from shape files
- Changes to sink filling and atb calculations to recognise burnt in river network
  and correctly split flow by gradient
- Alterations to output mew style dynamic TOPMODEL object
- Vignettes added documenting function usage

## Regressions of note
- No processing of distance histogram for river routing
